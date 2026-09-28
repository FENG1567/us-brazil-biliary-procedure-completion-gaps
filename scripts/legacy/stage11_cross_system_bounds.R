#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(arrow); library(data.table); library(parallel); library(splines)})

args <- commandArgs(trailingOnly = TRUE)
project <- normalizePath(if (length(args) >= 1) args[[1]] else ".", mustWork = TRUE)
nboot <- if (length(args) >= 2) as.integer(args[[2]]) else 500L
cores <- max(1L, min(8L, if (length(args) >= 3) as.integer(args[[3]]) else 8L))
setDTthreads(cores)
outdir <- file.path(project, "results", "stage11"); qcdir <- file.path(project, "qc", "stage11")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE); dir.create(qcdir, recursive = TRUE, showWarnings = FALSE)

binary <- function(x) as.numeric(tolower(as.character(x)) %in% c("1", "true", "t", "yes"))
nrd_path <- file.path(project, "derived", "analysis_ready_v1.0", "nrd_index_cohorts_v1.0.parquet")
sus_path <- file.path(project, "derived", "analysis_ready_v1.0", "datasus_index_cohorts_cnes_v1.0.parquet")

nrd <- as.data.table(read_parquet(nrd_path))
nrd <- nrd[target_primary == TRUE & as.integer(YEAR) %in% 2021:2022 &
             !dx_any_biliary_or_pancreatic_malignancy & !dx_any_cirrhosis_proxy &
             !dx_any_coagulation_disorder_proxy &
             (cohort != "acute_cholecystitis" | known_chole_state != "observed_prior_cholecystectomy")]
nrd[, `:=`(
  system = "United States", year = as.integer(YEAR), age = as.numeric(AGE_num),
  female = binary(FEMALE_num), cholangitis = binary(dx_any_cholangitis),
  obstruction = binary(dx_any_biliary_obstruction), pancreatitis = binary(dx_any_biliary_acute_pancreatitis),
  outcome_lower = fifelse(cohort == "acute_cholecystitis", binary(pr_any_cholecystectomy_total),
                          binary(pr_any_ercp_clearance_proxy | pr_any_surgical_cbd_clearance)),
  outcome_upper = fifelse(cohort == "acute_cholecystitis",
                          binary(pr_any_cholecystectomy_total | pr_any_cholecystectomy_partial),
                          binary(pr_any_ercp_therapeutic | pr_any_surgical_cbd_clearance)),
  raw_weight = as.numeric(DISCWT_num), cluster = paste0("US:", hospital_cluster_id),
  boot_stratum = paste0("US:", YEAR, ":", NRD_STRATUM)
)]

sus <- as.data.table(read_parquet(sus_path))
sus <- sus[cross_system_primary_eligible == TRUE & file_year %in% 2021:2022 & !malignancy_any &
             !cirrhosis_any & !coagulopathy_any &
             (cohort != "acute_cholecystitis" | !prior_chole_proxy_any)]
sus[, `:=`(
  system = "Brazil", year = as.integer(file_year), age = as.numeric(age_years),
  female = as.numeric(as.character(SEXO) == "3"), cholangitis = binary(cholangitis_any),
  obstruction = binary(biliary_obstruction_any), pancreatitis = binary(biliary_pancreatitis_any),
  outcome_lower = binary(definitive_completion_narrow),
  outcome_upper = fifelse(cohort == "acute_cholecystitis",
                          binary(rd_cholecystectomy | op_cholecystectomy),
                          binary(rd_therapeutic_ercp | op_therapeutic_ercp |
                                 rd_clearance_explicit | op_clearance_explicit)),
  raw_weight = 1, cluster = paste0("BR:", file_year, ":", CNES),
  boot_stratum = paste0("BR:", file_year, ":", file_uf)
)]

keep <- c("cohort", "system", "year", "age", "female", "cholangitis", "obstruction", "pancreatitis",
          "outcome_lower", "outcome_upper", "raw_weight", "cluster", "boot_stratum")
all_data <- rbindlist(list(nrd[, ..keep], sus[, ..keep]), use.names = TRUE)
age_replacement <- median(all_data$age, na.rm = TRUE)
female_replacement <- median(all_data$female, na.rm = TRUE)
if (!is.finite(age_replacement) || !is.finite(female_replacement)) stop("No finite cross-system imputation value")
set(all_data, which(!is.finite(all_data$age)), "age", age_replacement)
set(all_data, which(!is.finite(all_data$female)), "female", female_replacement)
all_data[, `:=`(system_us = as.numeric(system == "United States"), year_2022 = as.numeric(year == 2022))]

sample_multiplier <- function(d) {
  cl <- unique(d[, .(boot_stratum, cluster)])
  dr <- cl[, .(sampled = sample(cluster, .N, replace = TRUE)), by = boot_stratum]
  ct <- dr[, .N, by = sampled]
  z <- ct$N[match(d$cluster, ct$sampled)]; z[is.na(z)] <- 0; z
}

fit_bounds <- function(d, X, mult = rep(1, nrow(d)), details = FALSE) {
  wr <- d$raw_weight * mult
  active <- is.finite(wr) & wr > 0
  if (sum(active) < ncol(X) + 20 || length(unique(d$system_us[active])) < 2) return(NULL)
  we <- wr
  for (s in 0:1) {
    ii <- active & d$system_us == s
    if (!any(ii) || sum(wr[ii]) <= 0) return(NULL)
    we[ii] <- wr[ii] / sum(wr[ii])
  }
  fit <- suppressWarnings(glm.fit(X[active, , drop = FALSE], d$system_us[active], family = binomial(), weights = we[active]))
  b <- fit$coefficients
  if (!isTRUE(fit$converged) || any(is.infinite(b)) || all(is.na(b))) return(NULL)
  b[is.na(b)] <- 0
  ps <- rep(NA_real_, nrow(d)); ps[active] <- pmin(pmax(plogis(drop(X[active, , drop = FALSE] %*% b)), .001), .999)
  support <- active & ps >= .05 & ps <= .95
  ow <- rep(0, nrow(d))
  ow[support] <- we[support] * ifelse(d$system_us[support] == 1, 1 - ps[support], ps[support])
  if (any(!is.finite(ow)) || sum(ow) <= 0) return(NULL)
  vals <- c()
  for (defn in c("lower", "upper")) {
    y <- d[[paste0("outcome_", defn)]]
    for (s in c("us", "br")) {
      ii <- support & d$system_us == ifelse(s == "us", 1, 0)
      vals[paste(s, defn, sep = "_")] <- weighted.mean(y[ii], ow[ii])
    }
  }
  vals["gap_lower_bound"] <- vals["us_lower"] - vals["br_upper"]
  vals["gap_upper_bound"] <- vals["us_upper"] - vals["br_lower"]
  vals["weighted_support_us"] <- weighted.mean(support[d$system_us == 1], wr[d$system_us == 1])
  vals["weighted_support_br"] <- weighted.mean(support[d$system_us == 0], wr[d$system_us == 0])
  vals["record_support_us"] <- weighted.mean(support[d$system_us == 1], mult[d$system_us == 1])
  vals["record_support_br"] <- weighted.mean(support[d$system_us == 0], mult[d$system_us == 0])
  if (!details) return(vals)
  list(values = vals, ps = ps, support = support, overlap_weight = ow, equal_weight = we)
}

wvar <- function(x, w) { m <- weighted.mean(x, w); sum(w * (x - m)^2) / sum(w) }
smd <- function(x, g, w) {
  m1 <- weighted.mean(x[g == 1], w[g == 1]); m0 <- weighted.mean(x[g == 0], w[g == 0])
  sp <- sqrt((wvar(x[g == 1], w[g == 1]) + wvar(x[g == 0], w[g == 0])) / 2)
  if (!is.finite(sp) || sp == 0) 0 else (m1 - m0) / sp
}

point_rows <- list(); boot_rows <- list(); balance_rows <- list()
cl <- makePSOCKcluster(cores, outfile = "")
tryCatch({
  invisible(clusterEvalQ(cl, {Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1"); library(data.table); setDTthreads(1L); NULL}))
  clusterExport(cl, c("sample_multiplier", "fit_bounds"), envir = environment())
  for (cc in c("acute_cholecystitis", "choledocholithiasis")) {
    d <- copy(all_data[cohort == cc])
    X <- model.matrix(~ ns(age, 3) + female + cholangitis + obstruction + pancreatitis + year_2022, d)
    pt <- fit_bounds(d, X, details = TRUE)
    if (!is.list(pt)) stop("Cross-system point model failed")
    nms <- names(pt$values)
    clusterExport(cl, c("d", "X", "nms"), envir = environment())
    reps <- parLapply(cl, seq_len(nboot), function(b) {
      set.seed(20260914 + b + ifelse(d$cohort[1] == "acute_cholecystitis", 0L, 100000L))
      z <- fit_bounds(d, X, sample_multiplier(d), details = FALSE)
      if (is.null(z)) setNames(rep(NA_real_, length(nms)), nms) else z
    })
    bt <- as.data.table(do.call(rbind, reps)); bt[, `:=`(cohort = cc, replicate = .I)]
    boot_rows[[cc]] <- bt
    for (v in nms) {
      z <- bt[[v]]; z <- z[is.finite(z)]
      point_rows[[paste(cc, v)]] <- data.table(
        cohort = cc, estimand = v, estimate = pt$values[[v]], lcl = quantile(z, .025),
        ucl = quantile(z, .975), bootstrap_success = length(z),
        interpretation = "2021-2022 common-support coding-bound comparison; not a causal health-system effect"
      )
    }
    mm <- X[, -1, drop = FALSE]
    balance_rows[[cc]] <- data.table(
      cohort = cc, variable = colnames(mm),
      smd_before = apply(mm, 2, smd, g = d$system_us, w = pt$equal_weight),
      smd_after = apply(mm, 2, smd, g = d$system_us, w = pt$overlap_weight)
    )
  }
}, finally = try(stopCluster(cl), silent = TRUE))

points <- rbindlist(point_rows); boots <- rbindlist(boot_rows); balance <- rbindlist(balance_rows)
fwrite(points, file.path(outdir, "cross_system_2021_2022_coding_bounds.csv"))
fwrite(boots, file.path(outdir, "cross_system_2021_2022_bootstrap.csv"))
fwrite(balance, file.path(qcdir, "cross_system_2021_2022_balance.csv"))
gate <- data.table(
  status = ifelse(min(points$bootstrap_success) >= .90 * nboot && max(abs(balance$smd_after)) < .10, "PASS", "REVIEW"),
  requested_bootstrap = nboot, minimum_success = min(points$bootstrap_success),
  maximum_smd_after = max(abs(balance$smd_after)),
  period = "2021-2022 in both systems",
  outcome_mapping = "lower and upper procedure-code definitions; only the resulting interval is interpreted",
  target = "measured covariate overlap population after 0.05-0.95 system-propensity restriction",
  interpretation = "bounded measurement comparison, not causal decomposition"
)
fwrite(gate, file.path(qcdir, "cross_system_2021_2022_gate.csv"))
cat(sprintf("CROSS_SYSTEM_%s min_success=%d max_smd=%.4f\n", gate$status, gate$minimum_success, gate$maximum_smd_after))
if (gate$status != "PASS") quit(status = 2)
