#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(parallel)
  library(splines)
  library(survey)
})

options(survey.lonely.psu = "adjust")
args <- commandArgs(trailingOnly = TRUE)
project <- normalizePath(if (length(args) >= 1) args[[1]] else ".", mustWork = TRUE)
nboot <- if (length(args) >= 2) as.integer(args[[2]]) else 500L
cores <- if (length(args) >= 3) as.integer(args[[3]]) else 8L
cores <- max(1L, min(8L, cores))
if (nboot < 100L) stop("nboot must be at least 100")
setDTthreads(cores)

input <- file.path(project, "derived", "analysis_ready_v1.0", "stage11", "nrd_day3_landmark_v2.0.parquet")
outdir <- file.path(project, "results", "stage11")
qcdir <- file.path(project, "qc", "stage11")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(qcdir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(input))

binary <- function(x) {
  z <- tolower(as.character(x))
  as.numeric(z %in% c("1", "true", "t", "yes"))
}

fref <- function(x) {
  z <- as.character(x)
  z[is.na(z) | z == "" | z %in% c("-9", "-8", "-6", "-99")] <- "unknown"
  tab <- sort(table(z), decreasing = TRUE)
  factor(z, levels = c(names(tab)[1], sort(setdiff(unique(z), names(tab)[1]))))
}

weighted_var <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]
  if (!length(x)) return(NA_real_)
  m <- sum(w * x) / sum(w)
  sum(w * (x - m)^2) / sum(w)
}

effective_n <- function(w) {
  w <- w[is.finite(w) & w > 0]
  if (!length(w)) return(NA_real_)
  sum(w)^2 / sum(w^2)
}

d <- as.data.table(read_parquet(input))
d[, `:=`(
  age = as.numeric(AGE_num), female = binary(FEMALE_num), weekend = binary(AWEEKEND_num),
  elective = binary(ELECTIVE_num), ed = binary(HCUP_ED_num),
  cholangitis = binary(dx_any_cholangitis), obstruction = binary(dx_any_biliary_obstruction),
  pancreatitis = binary(dx_any_biliary_acute_pancreatitis),
  payer = fref(PAY1), income = fref(ZIPINC_QRTL), rurality = fref(PL_NCHS),
  bedsize = fref(HOSP_BEDSIZE), teaching = fref(HOSP_UR_TEACH),
  ownership = fref(H_CONTRL), urbanicity = fref(HOSP_URCAT4),
  chole_state = fref(known_chole_state), year_f = fref(YEAR),
  prior_volume = as.numeric(prior_90d_disease_volume), base_weight = as.numeric(DISCWT_num),
  treatment = binary(early_completion_day3),
  biliary = binary(biliary_readmission_d3_93), nonbiliary = binary(nonbiliary_readmission_d3_93),
  death = binary(observed_inpatient_death_d3_93), composite = binary(biliary_or_death_d3_93),
  total_cost = as.numeric(total_inpatient_cost_d3_93),
  boot_stratum = interaction(YEAR, NRD_STRATUM, drop = TRUE),
  cluster = factor(hospital_cluster_id)
)]
d[!is.finite(age), age := median(age, na.rm = TRUE)]
d[!is.finite(prior_volume), prior_volume := 0]
d[!is.finite(base_weight) | base_weight <= 0, base_weight := NA_real_]

base_terms <- c(
  "ns(age,3)", "female", "weekend", "elective", "ed", "payer", "income", "rurality",
  "cholangitis", "obstruction", "pancreatitis", "bedsize", "teaching", "ownership",
  "urbanicity", "year_f", "log1p(prior_volume)"
)

terms_for <- function(cohort_name) {
  if (cohort_name == "choledocholithiasis") c(base_terms, "chole_state") else base_terms
}

make_formula <- function(lhs, cohort_name, treatment_name = NULL) {
  rhs <- terms_for(cohort_name)
  if (!is.null(treatment_name)) rhs <- c(treatment_name, rhs)
  as.formula(paste(lhs, "~", paste(rhs, collapse = " + ")))
}

safe_glm <- function(formula, data, family, weights) {
  data <- as.data.frame(data)
  data$.analysis_weights <- as.numeric(weights)
  suppressWarnings(glm(formula, data = data, family = family, weights = .analysis_weights,
                       control = glm.control(maxit = 50)))
}

balance_table <- function(x, cohort_name, mm, w_before, w_after) {
  calc_smd <- function(v, w) {
    a <- x$treatment
    w1 <- w * (a == 1); w0 <- w * (a == 0)
    m1 <- sum(w1 * v) / sum(w1); m0 <- sum(w0 * v) / sum(w0)
    sdpool <- sqrt((weighted_var(v[a == 1], w[a == 1]) + weighted_var(v[a == 0], w[a == 0])) / 2)
    if (!is.finite(sdpool) || sdpool == 0) return(0)
    (m1 - m0) / sdpool
  }
  data.table(
    cohort = cohort_name,
    variable = colnames(mm),
    smd_before = apply(mm, 2, calc_smd, w = w_before),
    smd_after = apply(mm, 2, calc_smd, w = w_after)
  )
}

prepare_analysis <- function(x, cohort_name) {
  ps_formula <- make_formula("treatment", cohort_name)
  ps_fit <- safe_glm(ps_formula, x, binomial(), x$base_weight / mean(x$base_weight, na.rm = TRUE))
  if (!isTRUE(ps_fit$converged)) stop("Propensity model did not converge")
  x[, ps := pmin(pmax(as.numeric(predict(ps_fit, type = "response")), 0.001), 0.999)]
  if (any(!is.finite(x$ps))) stop("Non-finite propensity score")
  x[, overlap_weight := base_weight * fifelse(treatment == 1, 1 - ps, ps)]
  x[, analysis_weight := overlap_weight / mean(overlap_weight, na.rm = TRUE)]
  mm <- model.matrix(ps_formula, data = x)[, -1, drop = FALSE]
  list(data = x, ps_fit = ps_fit,
       balance = balance_table(x, cohort_name, mm, x$base_weight, x$overlap_weight))
}

estimate_prepared <- function(prepared, cohort_name, outcome, estimate_atu = FALSE, estimate_rr = TRUE) {
  x <- copy(prepared$data)
  des <- svydesign(ids = ~hospital_cluster_id, strata = ~boot_stratum,
                   weights = ~analysis_weight, data = x, nest = TRUE)
  rd_fit <- svyglm(as.formula(paste(outcome, "~ treatment")), design = des, family = gaussian())
  rd <- unname(coef(rd_fit)["treatment"]); rdse <- sqrt(vcov(rd_fit)["treatment", "treatment"])
  lrr <- NA_real_; lrrse <- NA_real_
  if (estimate_rr) {
    rr_fit <- svyglm(as.formula(paste(outcome, "~ treatment")), design = des,
                     family = quasipoisson(link = "log"))
    lrr <- unname(coef(rr_fit)["treatment"]); lrrse <- sqrt(vcov(rr_fit)["treatment", "treatment"])
  }
  y <- x[[outcome]]
  r1 <- weighted.mean(y[x$treatment == 1], x$overlap_weight[x$treatment == 1])
  r0 <- weighted.mean(y[x$treatment == 0], x$overlap_weight[x$treatment == 0])
  support <- x$ps >= 0.05 & x$ps <= 0.95
  atu <- x$treatment == 0 & support
  out_fit <- NULL; atu_r0 <- NA_real_; atu_r1 <- NA_real_
  if (estimate_atu) {
    out_formula <- make_formula(outcome, cohort_name, "treatment")
    out_fit <- safe_glm(out_formula, x, binomial(), x$base_weight / mean(x$base_weight, na.rm = TRUE))
    if (!isTRUE(out_fit$converged)) stop("Outcome model did not converge")
    x1 <- copy(x); x1[, treatment := 1]
    x0 <- copy(x); x0[, treatment := 0]
    m1 <- pmin(pmax(as.numeric(predict(out_fit, newdata = x1, type = "response")), 1e-6), 1 - 1e-6)
    m0 <- pmin(pmax(as.numeric(predict(out_fit, newdata = x0, type = "response")), 1e-6), 1 - 1e-6)
    atu_r0 <- weighted.mean(m0[atu], x$base_weight[atu])
    atu_r1 <- weighted.mean(m1[atu], x$base_weight[atu])
  }

  estimate <- data.table(
    cohort = cohort_name, outcome = outcome, estimand = "day3_overlap_population",
    n = nrow(x), events = sum(y), early_completed = sum(x$treatment == 1),
    not_early_completed = sum(x$treatment == 0),
    risk_not_early = r0, risk_early = r1, risk_difference = rd,
    rd_lcl = rd - 1.96 * rdse, rd_ucl = rd + 1.96 * rdse,
    risk_ratio = ifelse(estimate_rr, exp(lrr), r1 / r0),
    rr_lcl = ifelse(estimate_rr, exp(lrr - 1.96 * lrrse), NA_real_),
    rr_ucl = ifelse(estimate_rr, exp(lrr + 1.96 * lrrse), NA_real_),
    ps_min = min(x$ps), ps_p01 = quantile(x$ps, .01), ps_p99 = quantile(x$ps, .99), ps_max = max(x$ps),
    outside_005_095 = mean(!support), ess_early = effective_n(x$overlap_weight[x$treatment == 1]),
    ess_not_early = effective_n(x$overlap_weight[x$treatment == 0]),
    atu_support_n = sum(atu), atu_support_weighted_n = sum(x$base_weight[atu]),
    atu_model_risk_not_early = atu_r0, atu_model_risk_if_early = atu_r1,
    atu_model_rd = atu_r1 - atu_r0,
    interpretation = "fixed_day3_landmark_association_not_randomized_effect"
  )
  list(estimate = estimate, data = x, balance = prepared$balance,
       ps_fit = prepared$ps_fit, outcome_fit = out_fit)
}

fit_primary <- function(x, cohort_name, outcome, estimate_atu = TRUE, estimate_rr = TRUE) {
  prepared <- prepare_analysis(x, cohort_name)
  estimate_prepared(prepared, cohort_name, outcome, estimate_atu = estimate_atu, estimate_rr = estimate_rr)
}

primary_results <- list(); balances <- list(); fitted_primary <- list()
for (cc in c("acute_cholecystitis", "choledocholithiasis")) {
  cat(sprintf("PRIMARY_START %s\n", cc)); flush.console()
  x <- copy(d[cohort == cc & binary(primary_day3_eligible) == 1 & is.finite(base_weight)])
  prepared <- prepare_analysis(copy(x), cc)
  fit <- estimate_prepared(prepared, cc, "biliary", estimate_atu = TRUE)
  primary_results[[paste(cc, "biliary")]] <- fit$estimate
  balances[[cc]] <- fit$balance
  fitted_primary[[cc]] <- fit
  for (outcome in c("composite", "death", "nonbiliary")) {
    primary_results[[paste(cc, outcome)]] <-
      estimate_prepared(prepared, cc, outcome, estimate_atu = FALSE)$estimate
  }
  cat(sprintf("PRIMARY_DONE %s\n", cc)); flush.console()
}
primary <- rbindlist(primary_results, fill = TRUE)
balance <- rbindlist(balances, fill = TRUE)
fwrite(primary, file.path(outdir, "nrd_day3_primary_estimates.csv"))
fwrite(balance, file.path(qcdir, "nrd_day3_balance.csv"))

# Baseline characteristics before and after overlap weighting.
baseline_vars <- c("age", "prior_volume", "female", "weekend", "elective", "ed",
                   "cholangitis", "obstruction", "pancreatitis")
baseline_factors <- c("payer", "income", "rurality", "bedsize", "teaching",
                      "ownership", "urbanicity", "year_f")
baseline_rows <- list()
for (cc in names(fitted_primary)) {
  x <- fitted_primary[[cc]]$data
  for (v in baseline_vars) {
    for (weighted in c(FALSE, TRUE)) {
      w <- if (weighted) x$overlap_weight else x$base_weight
      for (a in 0:1) {
        keep <- x$treatment == a
        baseline_rows[[length(baseline_rows) + 1L]] <- data.table(
          cohort = cc, variable = v, treatment = a,
          weighting = ifelse(weighted, "overlap", "survey"),
          mean = weighted.mean(x[[v]][keep], w[keep]), n = sum(keep), weighted_n = sum(w[keep])
        )
      }
    }
  }
  if (cc == "choledocholithiasis") baseline_factors <- unique(c(baseline_factors, "chole_state"))
  for (v in baseline_factors) {
    for (level in levels(x[[v]])) {
      indicator <- as.numeric(x[[v]] == level)
      for (weighted in c(FALSE, TRUE)) {
        w <- if (weighted) x$overlap_weight else x$base_weight
        for (a in 0:1) {
          keep <- x$treatment == a
          baseline_rows[[length(baseline_rows) + 1L]] <- data.table(
            cohort = cc, variable = paste(v, level, sep = ":"), treatment = a,
            weighting = ifelse(weighted, "overlap", "survey"),
            mean = weighted.mean(indicator[keep], w[keep]), n = sum(keep), weighted_n = sum(w[keep])
          )
        }
      }
    }
  }
}
fwrite(rbindlist(baseline_rows), file.path(outdir, "nrd_day3_baseline_by_treatment.csv"))

# Outcome subtype analyses use the identical primary cohort and propensity weights.
subtype_map <- c(
  recurrent_acute_d3_93 = "Acute cholecystitis",
  recurrent_cbd_d3_93 = "Choledocholithiasis",
  recurrent_cholangitis_d3_93 = "Cholangitis",
  recurrent_pancreatitis_d3_93 = "Biliary pancreatitis",
  repeat_ercp_d3_93 = "Repeat therapeutic ERCP",
  urgent_surgery_d3_93 = "Urgent surgery"
)
subtype_rows <- list()
for (cc in names(fitted_primary)) {
  x <- fitted_primary[[cc]]$data
  des <- svydesign(ids = ~hospital_cluster_id, strata = ~boot_stratum,
                   weights = ~analysis_weight, data = x, nest = TRUE)
  for (v in names(subtype_map)) {
    x[, subtype_y := binary(get(v))]
    des <- update(des, subtype_y = x$subtype_y)
    fit <- svyglm(subtype_y ~ treatment, design = des, family = gaussian())
    rd <- unname(coef(fit)["treatment"]); se <- sqrt(vcov(fit)["treatment", "treatment"])
    r1 <- weighted.mean(x$subtype_y[x$treatment == 1], x$overlap_weight[x$treatment == 1])
    r0 <- weighted.mean(x$subtype_y[x$treatment == 0], x$overlap_weight[x$treatment == 0])
    subtype_rows[[length(subtype_rows) + 1L]] <- data.table(
      cohort = cc, outcome = subtype_map[[v]], events = sum(x$subtype_y),
      risk_not_early = r0, risk_early = r1, risk_difference = rd,
      rd_lcl = rd - 1.96 * se, rd_ucl = rd + 1.96 * se
    )
  }
}
fwrite(rbindlist(subtype_rows), file.path(outdir, "nrd_day3_subtype_estimates.csv"))

# Prespecified sensitivity definitions. Every variant retains a common, fixed landmark.
variant_specs <- list(
  broad_eligibility_day3 = list(day = 3, elig = "broad", exposure = "primary"),
  age18_79_day3 = list(day = 3, elig = "primary_age79", exposure = "primary"),
  exclude_transfer_day3 = list(day = 3, elig = "primary_transfer_clean", exposure = "primary"),
  principal_all_admissions_day3 = list(day = 3, elig = "principal_all", exposure = "primary"),
  primary_day2 = list(day = 2, elig = "dynamic", exposure = "primary"),
  primary_day5 = list(day = 5, elig = "dynamic", exposure = "primary"),
  primary_day7 = list(day = 7, elig = "dynamic", exposure = "primary"),
  strict_procedure_day3 = list(day = 3, elig = "primary", exposure = "strict"),
  broad_procedure_day3 = list(day = 3, elig = "primary", exposure = "broad")
)

make_variant_data <- function(dd, cc, spec) {
  L <- spec$day
  z <- copy(dd[cohort == cc & is.finite(base_weight)])
  base_elig <- binary(z$target_primary) == 1 & binary(z$urgent_admission) == 1 &
    binary(z$no_definite_prior_chole_for_acute) == 1 & binary(z$treatment_timing_unknown) == 0 &
    (is.na(z$death_day_from_admission) | z$death_day_from_admission > L) &
    (is.na(z$biliary_day_from_admission) | z$biliary_day_from_admission > L)
  if (spec$elig %in% c("primary", "primary_age79", "primary_transfer_clean", "dynamic")) {
    base_elig <- base_elig & binary(z$no_coded_contraindication_proxy) == 1
  }
  if (spec$elig == "primary_age79") base_elig <- base_elig & z$age <= 79
  if (spec$elig == "primary_transfer_clean") base_elig <- base_elig & binary(z$index_transfer_or_continuation_proxy) == 0
  if (spec$elig == "principal_all") {
    base_elig <- binary(z$target_primary) == 1 & binary(z$no_definite_prior_chole_for_acute) == 1 &
      binary(z$treatment_timing_unknown) == 0 & binary(z$no_coded_contraindication_proxy) == 1 &
      (is.na(z$death_day_from_admission) | z$death_day_from_admission > L) &
      (is.na(z$biliary_day_from_admission) | z$biliary_day_from_admission > L)
  }
  z <- z[base_elig]
  dayvar <- if (spec$exposure == "strict") z$strict_treatment_day else if (spec$exposure == "broad") z$broad_treatment_day else z$treatment_day
  z[, treatment := as.numeric(!is.na(dayvar) & dayvar <= L)]
  z[, biliary := as.numeric(!is.na(biliary_day_from_admission) &
                              biliary_day_from_admission > L &
                              biliary_day_from_admission <= L + 90)]
  z
}

sensitivity_rows <- list()
for (cc in c("acute_cholecystitis", "choledocholithiasis")) {
  for (nm in names(variant_specs)) {
    spec <- variant_specs[[nm]]
    # Total-only is meaningful for acute; broad therapeutic ERCP is meaningful for CBD.
    if (nm == "strict_procedure_day3" && cc == "choledocholithiasis") next
    if (nm == "broad_procedure_day3" && cc == "acute_cholecystitis") next
    x <- make_variant_data(d, cc, spec)
    if (nrow(x) < 1000 || length(unique(x$treatment)) < 2) next
    cat(sprintf("SENSITIVITY_START %s %s n=%d\n", cc, nm, nrow(x))); flush.console()
    fit <- fit_primary(x, cc, "biliary", estimate_atu = FALSE, estimate_rr = FALSE)$estimate
    fit[, variant := nm]
    sensitivity_rows[[paste(cc, nm)]] <- fit
    cat(sprintf("SENSITIVITY_DONE %s %s\n", cc, nm)); flush.console()
  }
}
sensitivity <- rbindlist(sensitivity_rows, fill = TRUE)
fwrite(sensitivity, file.path(outdir, "nrd_day3_sensitivity_estimates.csv"))

# Treatment and outcome frequencies by coded contraindication proxy. These are
# descriptive diagnostics, not strata in which treatment candidacy is assumed exchangeable.
proxy_rows <- list()
for (cc in c("acute_cholecystitis", "choledocholithiasis")) {
  x <- copy(d[cohort == cc & binary(broad_day3_eligible) == 1 & is.finite(base_weight)])
  x[, proxy_stratum := fifelse(binary(dx_any_cirrhosis_proxy) == 1 & binary(dx_any_coagulation_disorder_proxy) == 1, "both",
                        fifelse(binary(dx_any_cirrhosis_proxy) == 1, "cirrhosis_only",
                        fifelse(binary(dx_any_coagulation_disorder_proxy) == 1, "coagulopathy_only", "neither")))]
  proxy_rows[[cc]] <- x[, .(
    n = .N, weighted_n = sum(base_weight),
    early_completion = weighted.mean(treatment, base_weight),
    biliary_event = weighted.mean(biliary, base_weight),
    observed_inpatient_death = weighted.mean(death, base_weight),
    biliary_or_death = weighted.mean(composite, base_weight)
  ), by = .(cohort, proxy_stratum)]
}
fwrite(rbindlist(proxy_rows), file.path(outdir, "nrd_day3_contraindication_proxy_strata.csv"))

# Clinically anchored single-confounder sensitivity grid. It is a bias analysis,
# not a claim that unmeasured confounding has been removed.
qba <- rbindlist(lapply(c("acute_cholecystitis", "choledocholithiasis"), function(cc) {
  z <- primary[cohort == cc & outcome == "biliary"]
  grid <- CJ(prevalence_early = c(.05, .10, .20),
             prevalence_not_early = c(.10, .20, .40, .60),
             confounder_outcome_rr = c(2, 4, 8, 12))
  grid <- grid[prevalence_not_early >= prevalence_early]
  grid[, observed_rr := z$risk_ratio]
  grid[, imbalance_factor := (1 + prevalence_early * (confounder_outcome_rr - 1)) /
                              (1 + prevalence_not_early * (confounder_outcome_rr - 1))]
  grid[, corrected_rr := observed_rr / imbalance_factor]
  grid[, corrected_rd := z$risk_not_early * corrected_rr - z$risk_not_early]
  grid[, `:=`(cohort = cc,
              clinical_anchor = "unmeasured treatment contraindication or frailty more prevalent among patients without early completion",
              interpretation = "sensitivity_grid_not_identified_effect")]
  grid
}))
fwrite(qba, file.path(outdir, "nrd_day3_quantitative_bias_grid.csv"))

# Directional bounds for incompletely observed out-of-hospital death. The grid
# adds hypothetical death among patients without an observed biliary/death event.
death_bounds <- rbindlist(lapply(c("acute_cholecystitis", "choledocholithiasis"), function(cc) {
  z <- primary[cohort == cc & outcome == "composite"]
  g <- CJ(unobserved_death_early = c(0, .005, .01, .02, .04),
          unobserved_death_not_early = c(0, .005, .01, .02, .04))
  g[, adjusted_composite_early := z$risk_early + (1 - z$risk_early) * unobserved_death_early]
  g[, adjusted_composite_not_early := z$risk_not_early + (1 - z$risk_not_early) * unobserved_death_not_early]
  g[, adjusted_composite_rd := adjusted_composite_early - adjusted_composite_not_early]
  g[, `:=`(cohort = cc,
            interpretation = "directional scenario for unobserved out-of-hospital death; not an identified competing-risk estimate")]
  g
}))
fwrite(death_bounds, file.path(outdir, "nrd_day3_out_of_hospital_death_scenarios.csv"))

# Point burden and net inpatient facility-cost estimates in the same untreated,
# common-support population used by the outcome model.
burden_point <- list()
for (cc in names(fitted_primary)) {
  fit <- fitted_primary[[cc]]
  x <- fit$data
  support <- x$ps >= 0.05 & x$ps <= 0.95
  atu <- x$treatment == 0 & support
  cost_ok <- is.finite(x$index_cost_2022) & x$index_cost_2022 > 0 &
    is.finite(x$total_cost) & x$total_cost > 0 &
    x$priced_biliary_admissions_d3_93 == x$biliary_admissions_d3_93
  cost_formula <- make_formula("total_cost", cc, "treatment")
  cost_fit <- safe_glm(cost_formula, x[cost_ok], Gamma(link = "log"),
                       x$base_weight[cost_ok] / mean(x$base_weight[cost_ok]))
  x1 <- copy(x); x1[, treatment := 1]
  x0 <- copy(x); x0[, treatment := 0]
  c1 <- as.numeric(predict(cost_fit, newdata = x1, type = "response"))
  c0 <- as.numeric(predict(cost_fit, newdata = x0, type = "response"))
  cost_rd_atu <- weighted.mean(c1[atu] - c0[atu], x$base_weight[atu])
  z <- primary[cohort == cc & outcome == "biliary"]
  for (closure in c(.25, .50, .75)) {
    new_complete <- z$atu_support_weighted_n * closure
    burden_point[[paste(cc, closure)]] <- data.table(
      cohort = cc, scenario = paste0(round(100 * closure), "_percent_of_supported_not_early_group"),
      weighted_supported_not_early = z$atu_support_weighted_n,
      modeled_additional_early_completions = new_complete,
      associated_biliary_readmissions = -new_complete * z$atu_model_rd,
      net_inpatient_facility_cost_change_2022usd = new_complete * cost_rd_atu,
      atu_model_rd = z$atu_model_rd, atu_total_cost_rd = cost_rd_atu,
      interpretation = "model_based_associated_scenario_in_same_eligible_common_support_population_not_causal_savings"
    )
  }
}
fwrite(rbindlist(burden_point), file.path(outdir, "nrd_day3_aligned_burden_point.csv"))

# Design-consistent stratified hospital-year bootstrap with propensity and outcome
# models refitted in every replicate.
sample_multiplier <- function(x) {
  clusters <- unique(x[, .(boot_stratum, cluster)])
  draws <- clusters[, .(sampled = sample(cluster, .N, replace = TRUE)), by = boot_stratum]
  counts <- draws[, .N, by = sampled]
  ans <- counts$N[match(x$cluster, counts$sampled)]
  ans[is.na(ans)] <- 0
  ans
}

resolve_beta <- function(fit, design, treatment_required = FALSE) {
  beta <- fit$coefficients
  if (!isTRUE(fit$converged) || any(is.infinite(beta)) || all(is.na(beta))) return(NULL)
  if (treatment_required && (is.na(beta["treatment"]) || !"treatment" %in% names(beta))) return(NULL)
  beta[is.na(beta)] <- 0
  reconstructed <- fit$family$linkinv(drop(design %*% beta))
  if (any(!is.finite(reconstructed)) || max(abs(reconstructed - fit$fitted.values), na.rm = TRUE) > 1e-6) return(NULL)
  beta
}

boot_metric_names <- c(
  as.vector(rbind(
    paste0(c("biliary", "composite", "death", "nonbiliary"), ".risk_early"),
    paste0(c("biliary", "composite", "death", "nonbiliary"), ".risk_not_early")
  )),
  "biliary.rd", "composite.rd", "death.rd", "nonbiliary.rd",
  "atu.rd", "atu.weighted_n", "atu.cost_rd",
  "burden25", "burden50", "burden75",
  "netcost25", "netcost50", "netcost75",
  "weighted_population", "weighted_early", "weighted_not_early"
)

failed_boot <- function() setNames(rep(NA_real_, length(boot_metric_names)), boot_metric_names)

fit_boot <- function(x, X, Xout, multiplier) {
  w <- x$base_weight * multiplier
  active <- is.finite(w) & w > 0
  if (sum(active) < ncol(Xout) + 20 || length(unique(x$treatment[active])) < 2) return(failed_boot())
  xa <- x[active]; wa <- w[active]
  Xa <- X[active, , drop = FALSE]; Xoa <- Xout[active, , drop = FALSE]
  a <- xa$treatment
  psfit <- glm.fit(Xa, a, family = binomial(), weights = wa)
  psb <- resolve_beta(psfit, Xa)
  if (is.null(psb)) return(failed_boot())
  ps <- pmin(pmax(plogis(drop(Xa %*% psb)), .001), .999)
  ow <- wa * ifelse(a == 1, 1 - ps, ps)
  metric <- c()
  for (yname in c("biliary", "composite", "death", "nonbiliary")) {
    y <- xa[[yname]]
    metric <- c(metric,
      setNames(weighted.mean(y[a == 1], ow[a == 1]), paste0(yname, ".risk_early")),
      setNames(weighted.mean(y[a == 0], ow[a == 0]), paste0(yname, ".risk_not_early")))
  }
  y <- xa$biliary
  outfit <- glm.fit(Xoa, y, family = binomial(), weights = wa)
  outb <- resolve_beta(outfit, Xoa, treatment_required = TRUE)
  if (is.null(outb)) return(failed_boot())
  trtcol <- which(colnames(Xoa) == "treatment")
  X1 <- Xoa; X0 <- Xoa; X1[, trtcol] <- 1; X0[, trtcol] <- 0
  m1 <- pmin(pmax(plogis(drop(X1 %*% outb)), 1e-6), 1 - 1e-6)
  m0 <- pmin(pmax(plogis(drop(X0 %*% outb)), 1e-6), 1 - 1e-6)
  support <- ps >= .05 & ps <= .95
  atu <- a == 0 & support
  if (sum(atu) < 100) return(failed_boot())
  atu_rd <- weighted.mean(m1[atu] - m0[atu], wa[atu])
  atu_n <- sum(wa[atu])

  cost_ok <- is.finite(xa$index_cost_2022) & xa$index_cost_2022 > 0 &
    is.finite(xa$total_cost) & xa$total_cost > 0 &
    xa$priced_biliary_admissions_d3_93 == xa$biliary_admissions_d3_93
  if (sum(cost_ok) < ncol(Xoa) + 20) return(failed_boot())
  costfit <- glm.fit(Xoa[cost_ok, , drop = FALSE], xa$total_cost[cost_ok],
                     family = Gamma(link = "log"), weights = wa[cost_ok])
  costb <- resolve_beta(costfit, Xoa[cost_ok, , drop = FALSE], treatment_required = TRUE)
  if (is.null(costb)) return(failed_boot())
  c1 <- exp(drop(X1 %*% costb)); c0 <- exp(drop(X0 %*% costb))
  cost_rd <- weighted.mean(c1[atu] - c0[atu], wa[atu])
  metric <- c(metric,
    biliary.rd = unname(metric["biliary.risk_early"] - metric["biliary.risk_not_early"]),
    composite.rd = unname(metric["composite.risk_early"] - metric["composite.risk_not_early"]),
    death.rd = unname(metric["death.risk_early"] - metric["death.risk_not_early"]),
    nonbiliary.rd = unname(metric["nonbiliary.risk_early"] - metric["nonbiliary.risk_not_early"]),
    atu.rd = atu_rd, atu.weighted_n = atu_n, atu.cost_rd = cost_rd,
    burden25 = -atu_n * .25 * atu_rd, burden50 = -atu_n * .50 * atu_rd,
    burden75 = -atu_n * .75 * atu_rd,
    netcost25 = atu_n * .25 * cost_rd, netcost50 = atu_n * .50 * cost_rd,
    netcost75 = atu_n * .75 * cost_rd,
    weighted_population = sum(wa), weighted_early = sum(wa * a), weighted_not_early = sum(wa * (1 - a))
  )
  stopifnot(identical(names(metric), boot_metric_names))
  metric
}

cl <- makePSOCKcluster(cores, outfile = "")
boot_rows <- list()
tryCatch({
  invisible(clusterEvalQ(cl, {
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
    suppressPackageStartupMessages(library(data.table)); setDTthreads(1L); NULL
  }))
  clusterExport(cl, c("sample_multiplier", "fit_boot", "resolve_beta",
                      "boot_metric_names", "failed_boot"), envir = environment())
  for (cc in c("acute_cholecystitis", "choledocholithiasis")) {
    cat(sprintf("BOOTSTRAP_START %s nboot=%d\n", cc, nboot)); flush.console()
    x <- fitted_primary[[cc]]$data
    ps_formula <- make_formula("treatment", cc)
    out_formula <- make_formula("biliary", cc, "treatment")
    X <- model.matrix(ps_formula, data = x)
    Xout <- model.matrix(out_formula, data = x)
    clusterExport(cl, c("x", "X", "Xout"), envir = environment())
    reps <- parLapply(cl, seq_len(nboot), function(b) {
      set.seed(20260914 + b + ifelse(x$cohort[1] == "acute_cholecystitis", 0L, 100000L))
      fit_boot(x, X, Xout, sample_multiplier(x))
    })
    tab <- as.data.table(do.call(rbind, reps))
    tab[, `:=`(cohort = cc, replicate = .I)]
    boot_rows[[cc]] <- tab
    cat(sprintf("BOOTSTRAP_DONE %s\n", cc)); flush.console()
  }
}, finally = try(stopCluster(cl), silent = TRUE))

boot <- rbindlist(boot_rows, fill = TRUE)
fwrite(boot, file.path(outdir, "nrd_day3_design_bootstrap_replicates.csv"))
metric_cols <- setdiff(names(boot), c("cohort", "replicate"))
boot_summary <- boot[, rbindlist(lapply(metric_cols, function(v) {
  z <- get(v); z <- z[is.finite(z)]
  data.table(metric = v, replicates = length(z), mean = mean(z), se = sd(z),
             lcl = quantile(z, .025), ucl = quantile(z, .975))
})), by = cohort]
fwrite(boot_summary, file.path(outdir, "nrd_day3_design_bootstrap_summary.csv"))

max_smd <- balance[, max(abs(smd_after), na.rm = TRUE), by = cohort]
gate <- data.table(
  status = ifelse(min(boot_summary$replicates) >= .90 * nboot && max(max_smd$V1) < .10, "PASS", "REVIEW"),
  requested_bootstrap = nboot,
  minimum_successful_bootstrap = min(boot_summary$replicates),
  maximum_primary_smd = max(max_smd$V1),
  primary_population = "fixed day-3 landmark, principal urgent, no coded cirrhosis/coagulopathy, known PRDAY",
  adjustment_timing = "baseline or admission-recorded only; no LOS, discharge disposition, APR severity, APR mortality, or post-treatment ICU variables",
  variance = "survey linearization primary plus hospital-year resampling within NRD year-stratum with propensity and outcome models refit",
  interpretation = "association; residual clinical candidacy confounding remains"
)
fwrite(gate, file.path(qcdir, "nrd_day3_analysis_gate.csv"))
cat(sprintf("NRD_%s nboot=%d min_success=%d max_smd=%.5f\n",
            gate$status, nboot, gate$minimum_successful_bootstrap, gate$maximum_primary_smd))
if (gate$status != "PASS") quit(status = 2)
