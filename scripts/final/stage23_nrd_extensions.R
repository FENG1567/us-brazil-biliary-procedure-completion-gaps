#!/usr/bin/env Rscript

# NRD extensions: hospital-year variation, associated national burden,
# and COVID-era resilience. All estimates are observational associations.

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(lme4)
  library(parallel)
  library(splines)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("usage: stage23_nrd_extensions.R PROJECT OUTDIR [NBOOT] [CORES]")
project <- normalizePath(args[[1L]], mustWork = TRUE)
outdir <- args[[2L]]
nboot <- if (length(args) >= 3L) as.integer(args[[3L]]) else 500L
cores <- if (length(args) >= 4L) as.integer(args[[4L]]) else 8L
if (!is.finite(nboot) || nboot < 100L) stop("NBOOT must be at least 100")
if (!is.finite(cores) || cores < 1L || cores > 8L) stop("CORES must be 1--8")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
setDTthreads(cores)
set.seed(20260923L)

bin <- function(x) as.integer(tolower(as.character(x)) %in% c("1", "true", "t", "yes"))
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
fref <- function(x) {
  z <- as.character(x)
  z[is.na(z) | z == "" | z %in% c("-99", "-9", "-8", "-6")] <- "unknown"
  tt <- sort(table(z), decreasing = TRUE)
  factor(z, levels = c(names(tt)[1L], sort(setdiff(names(tt), names(tt)[1L]))))
}
mode_value <- function(x) {
  z <- as.character(x); z <- z[!is.na(z) & z != ""]
  if (!length(z)) return(NA_character_)
  names(sort(table(z), decreasing = TRUE))[1L]
}
effn <- function(w) {
  w <- w[is.finite(w) & w > 0]
  if (!length(w)) return(NA_real_)
  sum(w)^2 / sum(w^2)
}

# The hospital analyses must use the same source risk set as the primary
# analysis: this day-3 landmark parquet with events reconstructed from all
# five linked-admission files; do not substitute the earlier day-3
# eligibility flag.
input <- file.path(project, "stage12_revision", "results", "nrd_day3_landmark_v3.0.parquet")
linked_glob <- file.path(project, "derived", "analysis_ready_v1.0", "stage9", "nrd_*_index_patient_admissions_v1.0.parquet")
stage20_dir <- file.path(project, "outputs", "stage20_cgh_final_20260916", "production_v1")
stage20_boot <- file.path(stage20_dir, "stage20_nrd_hospital_year_bootstrap_replicates.csv")
stage20_log <- file.path(stage20_dir, "stage20_nrd_run.log")
if (!file.exists(input) || !file.exists(stage20_boot) || !file.exists(stage20_log)) stop("required primary-analysis inputs are missing")
linked_files <- Sys.glob(linked_glob)
if (length(linked_files) != 5L) stop("expected five linked-admission parquet files, found ", length(linked_files))

index_cols <- c(
  "cohort", "source_row", "index_id", "NRD_VisitLink", "NRD_DaysToEvent_num", "index_discharge_day",
  "target_primary", "urgent_admission", "observed_prior_cholecystectomy", "definitive_completion_primary",
  "treatment_day", "treatment_timing_unknown", "index_died", "LOS_num", "DISCWT_num", "YEAR", "NRD_STRATUM", "hospital_cluster_id",
  "HOSP_NRD", "AGE_num", "FEMALE_num", "AWEEKEND_num", "ELECTIVE_num", "HCUP_ED_num",
  "PAY1", "ZIPINC_QRTL", "PL_NCHS", "HOSP_BEDSIZE", "HOSP_UR_TEACH", "H_CONTRL",
  "HOSP_URCAT4", "DMONTH", "prior_90d_disease_volume",
  "biliary_admissions_d3_93", "biliary_cost_d3_93", "index_cost_2022"
)
raw <- as.data.table(read_parquet(input, col_select = index_cols))
idx <- raw[, .(
  row_id = .I, index_id_key = as.character(index_id), cohort, visit = as.character(NRD_VisitLink),
  year = as.integer(num(YEAR)), idx_source = as.character(source_row), idx_admit = num(NRD_DaysToEvent_num),
  idx_discharge = num(index_discharge_day)
)]
adm <- rbindlist(lapply(linked_files, function(p) as.data.table(read_parquet(p))), fill = TRUE)
need_adm <- c("source_row", "NRD_VisitLink", "NRD_DaysToEvent", "ELECTIVE", "SAMEDAYEVENT", "REHABTRANSFER", "DIED", "YEAR", "related_dx_any")
if (length(setdiff(need_adm, names(adm)))) stop("linked-admission columns missing: ", paste(setdiff(need_adm, names(adm)), collapse = ", "))
adm[, `:=`(
  visit = as.character(NRD_VisitLink), year = as.integer(num(YEAR)), adm_source = as.character(source_row),
  adm_day = num(NRD_DaysToEvent), adm_elective = bin(ELECTIVE), adm_sameday = num(SAMEDAYEVENT) > 0,
  adm_rehab = bin(REHABTRANSFER), adm_death = bin(DIED), adm_biliary = bin(related_dx_any)
)]
candidate <- merge(idx, adm[, .(visit, year, adm_source, adm_day, adm_elective, adm_sameday, adm_rehab, adm_death, adm_biliary)],
  by = c("visit", "year"), allow.cartesian = TRUE)
candidate <- candidate[adm_source != idx_source & is.finite(adm_day) & adm_day > idx_discharge]
candidate[, day_from_admission := adm_day - idx_admit]
full_followup <- candidate[, .(
  biliary_day = suppressWarnings(min(day_from_admission[adm_biliary == 1L & !adm_sameday & adm_rehab == 0L], na.rm = TRUE)),
  followup_death_day = suppressWarnings(min(day_from_admission[adm_death == 1L], na.rm = TRUE))
), by = .(row_id, index_id_key, cohort)]
full_followup[!is.finite(biliary_day), biliary_day := NA_real_]
full_followup[!is.finite(followup_death_day), followup_death_day := NA_real_]
nonb_followup <- candidate[adm_biliary == 0L & adm_elective == 0L & !adm_sameday & adm_rehab == 0L,
  .(nonbiliary_day = suppressWarnings(min(day_from_admission))), by = .(row_id, index_id_key, cohort)]
nonb_followup[!is.finite(nonbiliary_day), nonbiliary_day := NA_real_]
raw[, row_id := .I]
raw[full_followup, `:=`(biliary_day = i.biliary_day, followup_death_day = i.followup_death_day), on = .(row_id, index_id = index_id_key, cohort)]
raw[nonb_followup, nonbiliary_day := i.nonbiliary_day, on = .(row_id, index_id = index_id_key, cohort)]
raw[, index_death_day := fifelse(bin(index_died) == 1L, num(LOS_num), NA_real_)]
raw[, death_day := pmin(index_death_day, followup_death_day, na.rm = TRUE)]
raw[!is.finite(death_day), death_day := NA_real_]
raw[, `:=`(
  target_primary_b = bin(target_primary), urgent_b = bin(urgent_admission), prior_chole = bin(observed_prior_cholecystectomy),
  definitive = bin(definitive_completion_primary), trt_day = num(treatment_day), timing_unknown = bin(treatment_timing_unknown),
  base_weight = num(DISCWT_num),
  age = num(AGE_num), female = bin(FEMALE_num), weekend = bin(AWEEKEND_num),
  elective = bin(ELECTIVE_num), ed = bin(HCUP_ED_num),
  payer = fref(PAY1), income = fref(ZIPINC_QRTL), rurality = fref(PL_NCHS),
  bedsize = fref(HOSP_BEDSIZE), teaching = fref(HOSP_UR_TEACH), ownership = fref(H_CONTRL),
  urbanicity = fref(HOSP_URCAT4), year_f = fref(YEAR), month_f = fref(DMONTH),
  prior_volume = num(prior_90d_disease_volume),
  year = as.integer(num(YEAR)), cluster = as.character(hospital_cluster_id),
  boot_stratum = interaction(YEAR, NRD_STRATUM, drop = TRUE),
  biliary_admissions = num(biliary_admissions_d3_93), biliary_cost = num(biliary_cost_d3_93)
)]
base_ok <- raw$target_primary_b == 1L & raw$urgent_b == 1L &
  (raw$cohort != "acute_cholecystitis" | raw$prior_chole == 0L) & raw$timing_unknown == 0L &
  (!is.finite(raw$death_day) | raw$death_day > 3L) & (!is.finite(raw$biliary_day) | raw$biliary_day > 3L) &
  is.finite(raw$base_weight) & raw$base_weight > 0 & !is.na(raw$cluster) & raw$cluster != "" & is.finite(num(raw$LOS_num)) & num(raw$LOS_num) >= 0
d <- raw[base_ok]
d[, `:=`(
  treatment = as.integer(definitive == 1L & is.finite(trt_day) & trt_day <= 3L),
  biliary90 = as.integer(is.finite(biliary_day) & biliary_day > 3L & biliary_day <= 93L),
  negative90 = as.integer(is.finite(nonbiliary_day) & nonbiliary_day > 3L & nonbiliary_day <= 93L)
)]
if (any(d[, .N, by = cohort]$N < 100L) || any(d[, uniqueN(treatment), by = cohort]$V1 != 2L)) stop("day-3 risk-set support failed")
d[!is.finite(age), age := median(age, na.rm = TRUE)]
d[!is.finite(prior_volume), prior_volume := 0]
d[, model_weight := base_weight / mean(base_weight)]
expected_day3 <- data.table(cohort = c("acute_cholecystitis", "choledocholithiasis"), expected_n = c(233764L, 117293L))
observed_day3 <- d[, .(observed_n = .N, early_n = sum(treatment), biliary90_events = sum(biliary90), nonbiliary90_events = sum(negative90)), by = cohort]
cohort_check <- merge(expected_day3, observed_day3, by = "cohort", all = TRUE)
cohort_check[, `:=`(
  landmark_input = normalizePath(input),
  linked_glob = linked_glob, linked_file_count = length(linked_files), event_reconstruction = "full linked-admission reconstruction",
  matches_primary_day3 = observed_n == expected_n
)]
fwrite(cohort_check, file.path(outdir, "stage23_nrd_cohort_check.csv"))
if (!all(cohort_check$matches_primary_day3)) stop("hospital extension day-3 counts do not match the primary analysis")

# A monotone log term is used for prior disease volume: it remains defined in
# sparse strata where an internally estimated spline has no distinct knots.
base_terms <- c(
  "ns(age,3)", "female", "weekend", "elective", "ed", "payer", "income", "rurality",
  "bedsize", "teaching", "ownership", "urbanicity", "year_f", "month_f", "log1p(prior_volume)"
)

fit_hospital_year <- function(x, outcome, cohort_name) {
  x <- copy(x)
  x[, y := as.integer(get(outcome))]
  terms <- base_terms
  if (cohort_name == "choledocholithiasis") terms <- c(terms, "prior_chole")
  form <- as.formula(paste("y ~", paste(terms, collapse = " + "), "+ (1|cluster)"))
  fit <- glmer(
    form, data = x, family = binomial(), weights = model_weight, nAGQ = 0,
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5), calc.derivs = TRUE)
  )
  x[, `:=`(
    p_fixed = as.numeric(predict(fit, type = "response", re.form = NA)),
    p_cond = as.numeric(predict(fit, type = "response"))
  )]
  national <- weighted.mean(x$y, x$base_weight)
  random <- as.data.table(ranef(fit, condVar = TRUE)$cluster, keep.rownames = "cluster")
  setnames(random, "(Intercept)", "random_intercept")
  postvar <- attr(ranef(fit, condVar = TRUE)$cluster, "postVar")
  random[, posterior_variance := as.numeric(postvar[1, 1, ])]
  var_re <- as.numeric(VarCorr(fit)$cluster[1, 1])
  hospital <- x[, .(
    sample_n = .N,
    weighted_n = sum(base_weight),
    events = sum(y),
    weighted_events = sum(base_weight * y),
    observed_rate = weighted.mean(y, base_weight),
    expected = sum(base_weight * p_fixed),
    predicted = sum(base_weight * p_cond),
    year = unique(year)[1L],
    hospital_id = mode_value(HOSP_NRD),
    bedsize = mode_value(bedsize),
    teaching = mode_value(teaching),
    ownership = mode_value(ownership),
    urbanicity = mode_value(urbanicity)
  ), by = cluster]
  hospital <- merge(hospital, random, by = "cluster", all.x = TRUE)
  hospital[, `:=`(
    risk_standardized_rate = pmin(1, pmax(0, national * predicted / pmax(expected, 1e-8))),
    reliability = pmin(1, pmax(0, 1 - posterior_variance / pmax(var_re, 1e-12))),
    cohort = cohort_name,
    outcome = outcome,
    eligible_n20 = sample_n >= 20,
    eligible_stable = sample_n >= 20 & events >= 5 & (sample_n - events) >= 5
  )]
  opt <- fit@optinfo
  gradient <- if (!is.null(opt$derivs$gradient)) max(abs(opt$derivs$gradient)) else NA_real_
  diagnostics <- data.table(
    cohort = cohort_name, outcome = outcome, n = nrow(x), events = sum(x$y),
    hospital_years = uniqueN(x$cluster), random_intercept_variance = var_re,
    median_odds_ratio = exp(qnorm(.75) * sqrt(2 * var_re)), singular = isSingular(fit, tol = 1e-4),
    max_abs_gradient = gradient, national_weighted_rate = national,
    hospital_years_n20 = sum(hospital$eligible_n20), hospital_years_stable = sum(hospital$eligible_stable),
    rsr_p10 = quantile(hospital[eligible_n20 == TRUE]$risk_standardized_rate, .10, na.rm = TRUE),
    rsr_p50 = quantile(hospital[eligible_n20 == TRUE]$risk_standardized_rate, .50, na.rm = TRUE),
    rsr_p90 = quantile(hospital[eligible_n20 == TRUE]$risk_standardized_rate, .90, na.rm = TRUE),
    interpretation = "risk-standardized hospital-year distribution; not a quality ranking"
  )
  list(fit = fit, hospital = hospital, diagnostics = diagnostics)
}

hospital_results <- list()
for (cohort_name in c("acute_cholecystitis", "choledocholithiasis")) {
  x <- d[cohort == cohort_name]
  for (outcome in c("treatment", "biliary90")) {
    key <- paste(cohort_name, outcome, sep = "__")
    hospital_results[[key]] <- fit_hospital_year(x, outcome, cohort_name)
    cat(sprintf("HOSPITAL_MODEL_DONE cohort=%s outcome=%s\n", cohort_name, outcome)); flush.console()
  }
}
hospital_rates <- rbindlist(lapply(hospital_results, `[[`, "hospital"), fill = TRUE)
hospital_diagnostics <- rbindlist(lapply(hospital_results, `[[`, "diagnostics"), fill = TRUE)
fwrite(hospital_rates, file.path(outdir, "stage23_nrd_hospital_year_risk_standardized_rates.csv"))
fwrite(hospital_diagnostics, file.path(outdir, "stage23_nrd_hospital_year_model_diagnostics.csv"))
saveRDS(lapply(hospital_results, `[[`, "fit"), file.path(outdir, "stage23_nrd_hospital_year_models.rds"), compress = FALSE)

completion_hosp <- hospital_rates[cohort == "acute_cholecystitis" & outcome == "treatment" & eligible_n20 == TRUE]
readmission_hosp <- hospital_rates[cohort == "acute_cholecystitis" & outcome == "biliary90" & eligible_n20 == TRUE,
  .(cluster, readmission_rsr = risk_standardized_rate, readmission_reliability = reliability)]
paired <- merge(
  completion_hosp[, .(cluster, year, weighted_n, completion_rsr = risk_standardized_rate, completion_reliability = reliability)],
  readmission_hosp, by = "cluster", all = FALSE
)
correlation <- data.table(
  cohort = "acute_cholecystitis",
  hospital_years = nrow(paired),
  spearman_completion_readmission = cor(paired$completion_rsr, paired$readmission_rsr, method = "spearman"),
  pearson_completion_readmission = cor(paired$completion_rsr, paired$readmission_rsr, method = "pearson"),
  interpretation = "ecological hospital-year association; not a patient-level or causal effect"
)
fwrite(correlation, file.path(outdir, "stage23_nrd_hospital_year_completion_readmission_correlation.csv"))

# Associated national burden under peer-benchmark gap closure. The completion
# model is treated as fixed; hospital-year resampling and the primary-analysis
# target-trial bootstrap propagate sampling and association uncertainty.
stage20_reps <- fread(stage20_boot)[
  cohort == "acute_cholecystitis" & outcome == "biliary_readmission" & day == 90 & finite == TRUE
]
if (nrow(stage20_reps) < 1800L) stop("insufficient primary-analysis bootstrap replicates")
mean_cost_per_admission <- d[
  cohort == "acute_cholecystitis" & is.finite(biliary_cost) & biliary_cost >= 0 & biliary_admissions > 0,
  sum(base_weight * biliary_cost) / sum(base_weight * biliary_admissions)
]
closures <- c(.25, .50, .75, 1.00)
peer_p75 <- quantile(completion_hosp$risk_standardized_rate, .75, na.rm = TRUE)
point_shortfall <- completion_hosp[, sum(weighted_n * pmax(peer_p75 - risk_standardized_rate, 0))]
point_rd <- stage20_reps[, median(difference_early_minus_not_early)]
point_rows <- rbindlist(lapply(closures, function(closure) data.table(
  cohort = "acute_cholecystitis", scenario = paste0(round(100 * closure), "_percent_peer_gap_closure"),
  closure_fraction = closure, peer_benchmark = peer_p75,
  associated_newly_completed = closure * point_shortfall,
  associated_biliary_readmissions = closure * point_shortfall * abs(point_rd),
  associated_inpatient_cost_2022_usd = closure * point_shortfall * abs(point_rd) * mean_cost_per_admission
)))

burden_reps <- rbindlist(lapply(seq_len(max(1000L, nboot)), function(replicate_id) {
  sampled <- completion_hosp[, .SD[sample.int(.N, .N, replace = TRUE)], by = year]
  benchmark <- quantile(sampled$risk_standardized_rate, .75, na.rm = TRUE)
  shortfall <- sampled[, sum(weighted_n * pmax(benchmark - risk_standardized_rate, 0))]
  rd <- sample(stage20_reps$difference_early_minus_not_early, 1L)
  rbindlist(lapply(closures, function(closure) data.table(
    replicate = replicate_id, closure_fraction = closure, peer_benchmark = benchmark,
    associated_newly_completed = closure * shortfall,
    associated_biliary_readmissions = closure * shortfall * abs(rd),
    associated_inpatient_cost_2022_usd = closure * shortfall * abs(rd) * mean_cost_per_admission
  )))
}))
burden_summary <- burden_reps[, .(
  newly_completed_lcl = quantile(associated_newly_completed, .025),
  newly_completed_ucl = quantile(associated_newly_completed, .975),
  readmissions_lcl = quantile(associated_biliary_readmissions, .025),
  readmissions_ucl = quantile(associated_biliary_readmissions, .975),
  cost_lcl = quantile(associated_inpatient_cost_2022_usd, .025),
  cost_ucl = quantile(associated_inpatient_cost_2022_usd, .975)
), by = closure_fraction]
burden <- merge(point_rows, burden_summary, by = "closure_fraction", all.x = TRUE)
burden[, `:=`(
  bootstrap_replicates = max(1000L, nboot),
  negative_control_guard = "residual_selection_detected",
  interpretation = "potentially addressable associated burden under a peer-benchmark scenario; not preventable causal counts or guaranteed savings"
)]
fwrite(burden, file.path(outdir, "stage23_nrd_peer_benchmark_associated_burden.csv"))
fwrite(burden_reps, file.path(outdir, "stage23_nrd_peer_benchmark_burden_replicates.csv"))

# COVID-era resilience: re-estimate the overlap population separately in the
# pre-pandemic, acute-pandemic, and recovery periods.
acute <- copy(d[cohort == "acute_cholecystitis"])
acute[, era := fifelse(year <= 2019L, "2018-2019_pre", fifelse(year == 2020L, "2020_acute", "2021-2022_recovery"))]
# Calendar year is intentionally omitted within era so the 2020-only stratum
# has the same estimand and a full-rank design as the two multi-year strata.
era_terms <- setdiff(base_terms, "year_f")
ps_formula <- as.formula(paste("treatment ~", paste(era_terms, collapse = " + ")))
fit_era <- function(x) {
  fit <- suppressWarnings(glm(ps_formula, data = x, family = binomial(), weights = base_weight / mean(base_weight)))
  if (!isTRUE(fit$converged)) stop("era PS failed")
  x[, ps := pmin(pmax(as.numeric(predict(fit, type = "response")), .001), .999)]
  x[, ow := base_weight * ifelse(treatment == 1L, 1 - ps, ps)]
  x
}
era_point <- rbindlist(lapply(unique(acute$era), function(era_name) {
  x <- fit_era(copy(acute[era == era_name]))
  data.table(
    era = era_name, n = nrow(x), weighted_population = sum(x$base_weight),
    weighted_completion = weighted.mean(x$treatment, x$base_weight),
    overlap_risk_early = weighted.mean(x$biliary90[x$treatment == 1L], x$ow[x$treatment == 1L]),
    overlap_risk_not_early = weighted.mean(x$biliary90[x$treatment == 0L], x$ow[x$treatment == 0L]),
    risk_difference = weighted.mean(x$biliary90[x$treatment == 1L], x$ow[x$treatment == 1L]) -
      weighted.mean(x$biliary90[x$treatment == 0L], x$ow[x$treatment == 0L]),
    ess_early = effn(x$ow[x$treatment == 1L]), ess_not_early = effn(x$ow[x$treatment == 0L]),
    ps_outside_005_095 = mean(x$ps < .05 | x$ps > .95)
  )
}))

bootstrap_era_one <- function(x, seed) {
  # mclapply owns the eight-process bootstrap budget; prevent each child from
  # inheriting the parent data.table thread setting and oversubscribing it.
  setDTthreads(1L)
  set.seed(seed)
  clusters <- unique(x[, .(boot_stratum, cluster)])
  drawn <- clusters[, .(draw = sample(cluster, .N, replace = TRUE)), by = boot_stratum]
  counts <- drawn[, .(mult = .N), by = .(boot_stratum, draw)]
  multiplier <- counts$mult[match(paste(x$boot_stratum, x$cluster), paste(counts$boot_stratum, counts$draw))]
  multiplier[is.na(multiplier)] <- 0
  w <- x$base_weight * multiplier
  ok <- w > 0 & is.finite(w)
  # glm evaluates a weights expression in its data frame.  Materialize the
  # cluster-bootstrap weight to keep its evaluation and the resampled cohort
  # aligned in both serial and forked runs.
  boot_data <- copy(x[ok])
  boot_data[, boot_weight := w[ok]]
  fit <- suppressWarnings(glm(ps_formula, data = boot_data, family = binomial(),
    weights = boot_weight / mean(boot_weight)))
  # Cluster resamples can make a factor column aliased.  This is handled by
  # glm's QR prediction machinery and is not a failed replicate; require a
  # converged fit and finite standardized predictions instead.
  if (!isTRUE(fit$converged)) return(c(early = NA, not_early = NA, rd = NA, completion = NA))
  ps <- pmin(pmax(as.numeric(predict(fit, newdata = x, type = "response")), .001), .999)
  if (any(!is.finite(ps))) return(c(early = NA, not_early = NA, rd = NA, completion = NA))
  ow <- w * ifelse(x$treatment == 1L, 1 - ps, ps)
  early <- weighted.mean(x$biliary90[x$treatment == 1L], ow[x$treatment == 1L])
  not_early <- weighted.mean(x$biliary90[x$treatment == 0L], ow[x$treatment == 0L])
  c(early = early, not_early = not_early, rd = early - not_early, completion = weighted.mean(x$treatment[ok], w[ok]))
}

era_boot <- list()
for (era_name in unique(acute$era)) {
  x <- copy(acute[era == era_name])
  seeds <- 202609230L + seq_len(nboot) + match(era_name, sort(unique(acute$era))) * 100000L
  if (.Platform$OS.type == "unix" && cores > 1L) {
    estimates <- mclapply(seeds, function(seed) bootstrap_era_one(x, seed), mc.cores = cores, mc.preschedule = TRUE)
  } else {
    estimates <- lapply(seeds, function(seed) bootstrap_era_one(x, seed))
  }
  z <- as.data.table(do.call(rbind, estimates))
  z[, `:=`(era = era_name, replicate = .I)]
  era_boot[[era_name]] <- z
  cat(sprintf("COVID_BOOT_DONE era=%s success=%d/%d\n", era_name, sum(is.finite(z$rd)), nboot)); flush.console()
}
era_boot <- rbindlist(era_boot)
era_ci <- era_boot[, .(
  bootstrap_success = sum(is.finite(rd)),
  risk_early_lcl = quantile(early[is.finite(early)], .025), risk_early_ucl = quantile(early[is.finite(early)], .975),
  risk_not_early_lcl = quantile(not_early[is.finite(not_early)], .025), risk_not_early_ucl = quantile(not_early[is.finite(not_early)], .975),
  rd_lcl = quantile(rd[is.finite(rd)], .025), rd_ucl = quantile(rd[is.finite(rd)], .975),
  completion_lcl = quantile(completion[is.finite(completion)], .025), completion_ucl = quantile(completion[is.finite(completion)], .975)
), by = era]
era_final <- merge(era_point, era_ci, by = "era", all.x = TRUE)
era_final[, interpretation := "era-stratified observational association; no pandemic causal effect"]
fwrite(era_final, file.path(outdir, "stage23_nrd_covid_resilience.csv"))
fwrite(era_boot, file.path(outdir, "stage23_nrd_covid_resilience_bootstrap_replicates.csv"))

gate <- data.table(
  check = c(
    "primary_day3_cohort_match", "four_hospital_models", "hospital_models_finite", "hospital_model_gradient",
    "burden_replicates", "covid_bootstrap_success", "cost_positive"
  ),
  pass = c(
    all(cohort_check$matches_primary_day3),
    nrow(hospital_diagnostics) == 4L,
    all(is.finite(hospital_diagnostics$median_odds_ratio)),
    all(hospital_diagnostics$max_abs_gradient < .05 | is.na(hospital_diagnostics$max_abs_gradient)),
    uniqueN(burden_reps$replicate) >= 1000L,
    all(era_final$bootstrap_success >= .90 * nboot),
    is.finite(mean_cost_per_admission) && mean_cost_per_admission > 0
  )
)
fwrite(gate, file.path(outdir, "stage23_nrd_extensions_gate.csv"))
if (!all(gate$pass)) stop("NRD extension gate failed")
cat("NRD_EXTENSIONS_PASS\n")
