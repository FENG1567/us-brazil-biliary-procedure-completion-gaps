#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(parallel)
  library(splines)
})

args <- commandArgs(trailingOnly = TRUE)
project <- normalizePath(if (length(args) >= 1) args[[1]] else ".", mustWork = TRUE)
nboot <- if (length(args) >= 2) as.integer(args[[2]]) else 500L
cores <- if (length(args) >= 3) as.integer(args[[3]]) else 8L
if (!is.finite(nboot) || nboot < 100) stop("nboot must be at least 100")
cores <- max(1L, min(8L, cores))
set.seed(20260912)
outdir <- file.path(project, "results", "stage26_cross_system")
qcdir <- file.path(project, "qc", "stage26_cross_system")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(qcdir, recursive = TRUE, showWarnings = FALSE)

binary <- function(x) as.numeric(tolower(as.character(x)) %in% c("1", "true", "t", "yes"))

nrd_path <- file.path(project, "derived", "analysis_ready_v1.0", "nrd_index_cohorts_v1.0.parquet")
sus_path <- file.path(project, "derived", "analysis_ready_v1.0", "datasus_index_cohorts_cnes_v1.0.parquet")
stopifnot(file.exists(nrd_path), file.exists(sus_path))

nrd <- as.data.table(read_parquet(nrd_path, col_select = c(
  "cohort", "target_primary", "definitive_completion_primary", "YEAR", "AGE_num", "FEMALE_num",
  "dx_any_cholangitis", "dx_any_biliary_obstruction", "dx_any_biliary_acute_pancreatitis",
  "dx_any_cirrhosis_proxy", "dx_any_coagulation_disorder_proxy", "DISCWT_num", "NRD_STRATUM",
  "hospital_cluster_id", "dx_any_biliary_or_pancreatic_malignancy"
)))
nrd <- nrd[target_primary == TRUE & !dx_any_biliary_or_pancreatic_malignancy]
nrd[, `:=`(
  system = "United States", outcome = binary(definitive_completion_primary), year = as.integer(YEAR),
  age = as.numeric(AGE_num), female = binary(FEMALE_num), cholangitis = binary(dx_any_cholangitis),
  obstruction = binary(dx_any_biliary_obstruction), pancreatitis = binary(dx_any_biliary_acute_pancreatitis),
  cirrhosis = binary(dx_any_cirrhosis_proxy), coagulopathy = binary(dx_any_coagulation_disorder_proxy),
  raw_weight = as.numeric(DISCWT_num),
  cluster = paste0("US:", hospital_cluster_id), boot_stratum = paste0("US:", YEAR, ":", NRD_STRATUM)
)]

sus <- as.data.table(read_parquet(sus_path, col_select = c(
  "cohort", "cross_system_primary_eligible", "definitive_completion_primary", "file_year", "age_years",
  "SEXO", "cholangitis_any", "biliary_obstruction_any", "biliary_pancreatitis_any", "cirrhosis_any",
  "coagulopathy_any", "CNES", "file_uf", "malignancy_any"
)))
sus <- sus[cross_system_primary_eligible == TRUE & !malignancy_any]
sus[, `:=`(
  system = "Brazil", outcome = binary(definitive_completion_primary), year = as.integer(file_year),
  age = as.numeric(age_years), female = as.numeric(as.character(SEXO) == "3"),
  cholangitis = binary(cholangitis_any), obstruction = binary(biliary_obstruction_any),
  pancreatitis = binary(biliary_pancreatitis_any), cirrhosis = binary(cirrhosis_any),
  coagulopathy = binary(coagulopathy_any), raw_weight = 1,
  cluster = paste0("BR:", file_year, ":", CNES), boot_stratum = paste0("BR:", file_year, ":", file_uf)
)]

common_names <- c("cohort", "system", "outcome", "year", "age", "female", "cholangitis", "obstruction",
                  "pancreatitis", "cirrhosis", "coagulopathy", "raw_weight", "cluster", "boot_stratum")
all_data <- rbindlist(list(nrd[, ..common_names], sus[, ..common_names]), use.names = TRUE)
all_data[!is.finite(age), age := median(age, na.rm = TRUE)]
all_data[, year_f := factor(year)]
all_data[, system_us := as.numeric(system == "United States")]

weighted_var <- function(x, w) {
  m <- sum(w * x) / sum(w)
  sum(w * (x - m)^2) / sum(w)
}

smd <- function(x, group, w) {
  w1 <- w * (group == 1); w0 <- w * (group == 0)
  m1 <- sum(w1 * x) / sum(w1); m0 <- sum(w0 * x) / sum(w0)
  den <- sqrt((weighted_var(x[group == 1], w[group == 1]) + weighted_var(x[group == 0], w[group == 0])) / 2)
  if (!is.finite(den) || den == 0) 0 else (m1 - m0) / den
}

sample_cluster_multiplier <- function(d) {
  clusters <- unique(d[, .(boot_stratum, cluster)])
  draws <- clusters[, .(sampled = sample(cluster, .N, replace = TRUE)), by = boot_stratum]
  counts <- draws[, .N, by = sampled]
  ans <- counts$N[match(d$cluster, counts$sampled)]
  ans[is.na(ans)] <- 0
  ans
}

fit_components <- function(d, X, multiplier = rep(1, nrow(d)), return_details = FALSE) {
  weight_raw <- d$raw_weight * multiplier
  active <- is.finite(weight_raw) & weight_raw > 0
  if (sum(active) < ncol(X) || any(table(d$system[active]) == 0)) {
    return(rep(NA_real_, 8))
  }
  weight_equal <- weight_raw
  for (ss in unique(d$system)) {
    ii <- d$system == ss
    weight_equal[ii] <- weight_raw[ii] / sum(weight_raw[ii]) * nrow(d) / 2
  }

  resolve_aliased <- function(fit, fit_design, prediction_design = fit_design) {
    beta <- fit$coefficients
    if (!isTRUE(fit$converged) || isTRUE(fit$boundary) ||
        any(is.infinite(beta)) || all(is.na(beta))) return(NULL)
    aliased <- which(is.na(beta))
    if (length(aliased)) {
      estimable <- which(!is.na(beta))
      if (!length(estimable)) return(NULL)
      estimable_qr <- qr(fit_design[, estimable, drop = FALSE])
      for (jj in aliased) {
        mapping <- qr.coef(estimable_qr, fit_design[, jj])
        if (any(!is.finite(mapping))) return(NULL)
        fit_residual <- fit_design[, jj] - drop(fit_design[, estimable, drop = FALSE] %*% mapping)
        prediction_residual <- prediction_design[, jj] -
          drop(prediction_design[, estimable, drop = FALSE] %*% mapping)
        tolerance <- 1e-8 * max(1, max(abs(fit_design[, jj]), na.rm = TRUE))
        # An aliased parameter may be zeroed only if the same linear dependency
        # holds on the population where predictions are standardized.  This
        # prevents arbitrary extrapolation to covariate patterns absent from a
        # system-specific bootstrap sample.
        if (max(abs(fit_residual), na.rm = TRUE) > tolerance ||
            max(abs(prediction_residual), na.rm = TRUE) > tolerance) return(NULL)
      }
      beta[aliased] <- 0
    }
    reconstructed <- plogis(drop(fit_design %*% beta))
    if (any(!is.finite(reconstructed)) ||
        max(abs(reconstructed - fit$fitted.values), na.rm = TRUE) > 1e-8) {
      return(NULL)
    }
    predicted <- plogis(drop(prediction_design %*% beta))
    if (any(!is.finite(predicted))) return(NULL)
    list(beta = beta, predicted = predicted)
  }

  X_active <- X[active, , drop = FALSE]
  ps_fit <- glm.fit(X_active, d$system_us[active], family = binomial(), weights = weight_equal[active])
  ps_resolved <- resolve_aliased(ps_fit, X_active)
  if (is.null(ps_resolved)) return(rep(NA_real_, 8))
  ps <- rep(NA_real_, nrow(d))
  ps[active] <- pmin(pmax(ps_resolved$predicted, 0.001), 0.999)
  support <- active & ps >= 0.05 & ps <= 0.95
  target_weight <- weight_equal * ifelse(d$system_us == 1, 1 - ps, ps) * support
  if (!any(support) || !all(is.finite(target_weight[support])) || sum(target_weight[support]) <= 0) {
    return(rep(NA_real_, 8))
  }
  us <- active & d$system_us == 1; br <- active & d$system_us == 0
  if (length(unique(d$outcome[us])) != 2L || length(unique(d$outcome[br])) != 2L) {
    return(rep(NA_real_, 8))
  }
  us_fit <- glm.fit(X[us, , drop = FALSE], d$outcome[us], family = binomial(), weights = weight_raw[us])
  br_fit <- glm.fit(X[br, , drop = FALSE], d$outcome[br], family = binomial(), weights = weight_raw[br])
  us_resolved <- resolve_aliased(us_fit, X[us, , drop = FALSE], X[support, , drop = FALSE])
  br_resolved <- resolve_aliased(br_fit, X[br, , drop = FALSE], X[support, , drop = FALSE])
  if (is.null(us_resolved) || is.null(br_resolved)) return(rep(NA_real_, 8))
  pred_us <- us_resolved$predicted
  pred_br <- br_resolved$predicted
  raw_us <- weighted.mean(d$outcome[us], weight_raw[us])
  raw_br <- weighted.mean(d$outcome[br], weight_raw[br])
  std_us <- weighted.mean(pred_us, target_weight[support])
  std_br <- weighted.mean(pred_br, target_weight[support])
  raw_gap <- raw_us - raw_br
  residual_gap <- std_us - std_br
  composition <- raw_gap - residual_gap
  # Record-level overlap in a cluster bootstrap replicate must use the
  # resampling multiplicity as its empirical frequency.  Counting unsampled
  # zero-weight records in the denominator would center the bootstrap
  # distribution far below the point estimate.
  record_overlap_fraction <- weighted.mean(support[active], multiplier[active])
  values <- c(raw_us, raw_br, raw_gap, std_us, std_br, residual_gap, composition, record_overlap_fraction)
  if (!return_details) return(values)
  list(values = values, ps = ps, support = support, target_weight = target_weight, X = X,
       weight_equal = weight_equal, weight_raw = weight_raw)
}

point_rows <- list(); balance_rows <- list(); bootstrap_rows <- list(); support_rows <- list()
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
setDTthreads(1L)
for (cc in c("acute_cholecystitis", "choledocholithiasis")) {
    d <- copy(all_data[cohort == cc])
    # Use the common 2018-2022 calendar window. The upstream cohort's
    # cross_system_primary_eligible flag already excludes years in which the
    # disease-specific therapeutic-ERCP phenotype is not observable.
    if (cc == "choledocholithiasis") d <- d[year %in% 2018:2022]
    d[, year_f := droplevels(factor(year))]
    formula <- ~ ns(age, 3) + female + cholangitis + obstruction + pancreatitis + cirrhosis + coagulopathy + year_f
    X <- model.matrix(formula, data = d)
    point <- fit_components(d, X, return_details = TRUE)
    if (!is.list(point)) stop("Cross-system point model failed before bootstrap")
  vals <- point$values
  names(vals) <- c("raw_us", "raw_brazil", "raw_gap", "standardized_us", "standardized_brazil",
                   "residual_system_gap", "measured_composition_component", "record_overlap_fraction")

    # Sequential bootstrap is deliberate here. The record-level design matrix
    # is large and PSOCK serialization can fail with SIGPIPE on this dataset.
    boot <- lapply(seq_len(nboot), function(b) {
      set.seed(20260912 + b + ifelse(cc == "acute_cholecystitis", 0, 100000))
      mult <- sample_cluster_multiplier(d)
      fit_components(d, X, mult, return_details = FALSE)
    })
  boot <- do.call(rbind, boot)
  colnames(boot) <- names(vals)
  ok <- complete.cases(boot)
  if (sum(ok) < 0.9 * nboot) stop("More than 10 percent of cross-system bootstrap replicates failed")
  boot <- boot[ok, , drop = FALSE]
  ci <- apply(boot, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  point_table <- data.table(
    database = "NRD_vs_SIH-SUS", cohort = cc, estimand = names(vals), estimate = as.numeric(vals),
    lcl = ci[1, names(vals)], ucl = ci[2, names(vals)], bootstrap_replicates = nrow(boot),
    interpretation = "common_overlap_standardization_not_pure_causal_system_effect"
  )
  point_rows[[cc]] <- point_table
  bootstrap_rows[[cc]] <- data.table(cohort = cc, replicate = seq_len(nrow(boot)), as.data.table(boot))

  mm <- point$X[, -1, drop = FALSE]
  balance <- data.table(
    database = "NRD_vs_SIH-SUS", cohort = cc, variable = colnames(mm),
    smd_before = apply(mm, 2, smd, group = d$system_us, w = point$weight_equal),
    smd_after = apply(mm, 2, smd, group = d$system_us, w = point$target_weight)
  )
  balance_rows[[cc]] <- balance
    # Materialize row-aligned point-model diagnostics before grouping by
    # system. Referencing point$support or point$ps directly inside a grouped
    # data.table expression would retain full-cohort length and mismatch each
    # system-specific weight vector.
    d[, `:=`(point_support = point$support, point_ps = point$ps)]
    support_rows[[cc]] <- d[, .(
      n = .N,
      weighted_support_fraction = weighted.mean(point_support, raw_weight),
      record_support_fraction = mean(point_support),
      ps_p01 = quantile(point_ps, 0.01), ps_p50 = quantile(point_ps, 0.50), ps_p99 = quantile(point_ps, 0.99)
    ), by = system][, `:=`(database = "NRD_vs_SIH-SUS", cohort = cc)]
}

results <- rbindlist(point_rows, fill = TRUE)
balance <- rbindlist(balance_rows, fill = TRUE)
bootstraps <- rbindlist(bootstrap_rows, fill = TRUE)
support <- rbindlist(support_rows, fill = TRUE)
fwrite(results, file.path(outdir, "cross_system_completion_decomposition.csv"))
fwrite(bootstraps, file.path(outdir, "cross_system_decomposition_bootstrap.csv"))
fwrite(balance, file.path(qcdir, "cross_system_overlap_balance.csv"))
fwrite(support, file.path(qcdir, "cross_system_overlap_support.csv"))

status <- if (all(abs(balance$smd_after) < 0.10) && all(support$weighted_support_fraction >= 0.80)) "PASS" else "REVIEW_BOUNDED_OVERLAP"
gate <- data.table(status = status, bootstrap_replicates = nboot,
                   max_abs_smd_after = max(abs(balance$smd_after)),
                   min_weighted_support = min(support$weighted_support_fraction),
                   note = "A REVIEW result retains bounded overlap-population estimates but prohibits a single headline decomposed percentage.")
fwrite(gate, file.path(qcdir, "cross_system_decomposition_gate.csv"))
cat(sprintf("CROSS_SYSTEM_%s bootstrap=%d max_smd=%.4f min_support=%.4f\n",
            status, nboot, gate$max_abs_smd_after, gate$min_weighted_support))
