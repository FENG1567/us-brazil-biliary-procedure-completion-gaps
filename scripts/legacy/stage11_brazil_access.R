#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(splines)
})

args <- commandArgs(trailingOnly = TRUE)
project <- normalizePath(if (length(args)) args[[1]] else ".", mustWork = TRUE)
setDTthreads(max(1L, min(8L, if (length(args) >= 2) as.integer(args[[2]]) else 8L)))
input <- file.path(project, "derived", "analysis_ready_v1.0", "stage9",
                   "datasus_care_flow_network_geography_v1.0.parquet")
outdir <- file.path(project, "results", "stage11")
qcdir <- file.path(project, "qc", "stage11")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(qcdir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(input))

binary <- function(x) as.numeric(tolower(as.character(x)) %in% c("1", "true", "t", "yes"))
fmiss <- function(x) {
  z <- as.character(x); z[is.na(z) | z == ""] <- "missing"
  ref <- names(which.max(table(z)))
  factor(z, levels = c(ref, sort(setdiff(unique(z), ref))))
}

# HC1 cluster-robust sandwich covariance implemented locally.  For two-way
# clustering we use the Cameron-Gelbach-Miller inclusion-exclusion form:
# origin + destination - origin-by-destination intersection.  Negative
# numerical eigenvalues are truncated at zero, matching fix=TRUE behavior.
cluster_meat <- function(score, cluster) {
  g <- as.integer(factor(cluster))
  s <- rowsum(score, g, reorder = FALSE)
  G <- nrow(s); n <- nrow(score); k <- ncol(score)
  if (G <= 1) stop("Insufficient clusters")
  (G / (G - 1)) * ((n - 1) / max(1, n - k)) * crossprod(s)
}

psd_fix <- function(v) {
  v <- (v + t(v)) / 2
  ee <- eigen(v, symmetric = TRUE)
  ee$values[ee$values < 0] <- 0
  out <- ee$vectors %*% (ee$values * t(ee$vectors))
  dimnames(out) <- dimnames(v)
  out
}

glm_cluster_vcovs <- function(fit, origin, destination) {
  b <- coef(fit); keep <- is.finite(b)
  X <- model.matrix(fit)[, keep, drop = FALSE]
  if (length(origin) != nrow(X)) {
    omitted <- if (is.null(fit$na.action)) integer() else as.integer(fit$na.action)
    used <- setdiff(seq_along(origin), omitted)
    origin <- origin[used]; destination <- destination[used]
  }
  if (length(origin) != nrow(X) || length(destination) != nrow(X)) {
    stop("Cluster vectors do not align with the fitted model frame")
  }
  working_w <- as.numeric(fit$weights)
  bread <- qr.solve(crossprod(X, X * working_w), diag(ncol(X)), tol = 1e-10)
  score <- X * as.numeric(fit$y - fit$fitted.values)
  m_orig <- cluster_meat(score, origin)
  m_dest <- cluster_meat(score, destination)
  m_both <- cluster_meat(score, interaction(origin, destination, drop = TRUE))
  inner <- list(
    two_way_origin_destination = m_orig + m_dest - m_both,
    destination_hospital = m_dest,
    residence_municipality = m_orig
  )
  lapply(inner, function(meat) {
    small <- psd_fix(bread %*% meat %*% bread)
    full <- matrix(0, nrow = length(b), ncol = length(b),
                   dimnames = list(names(b), names(b)))
    full[keep, keep] <- small
    full
  })
}

standardize <- function(fit, vc, x, variable, low, high) {
  x0 <- copy(x); x1 <- copy(x)
  if (is.factor(x[[variable]])) {
    lev <- levels(x[[variable]])
    x0[[variable]] <- factor(low, levels = lev)
    x1[[variable]] <- factor(high, levels = lev)
  } else {
    x0[[variable]] <- low; x1[[variable]] <- high
  }
  tt <- delete.response(terms(fit)); bn <- names(coef(fit))
  X0 <- model.matrix(tt, x0)[, bn, drop = FALSE]
  X1 <- model.matrix(tt, x1)[, bn, drop = FALSE]
  b <- coef(fit); b[is.na(b)] <- 0
  p0 <- plogis(drop(X0 %*% b)); p1 <- plogis(drop(X1 %*% b))
  w <- rep(1 / nrow(x), nrow(x))
  g <- colSums((X1 * (p1 * (1 - p1)) - X0 * (p0 * (1 - p0))) * w)
  se <- sqrt(max(0, drop(t(g) %*% vc[bn, bn, drop = FALSE] %*% g)))
  rd <- mean(p1 - p0)
  data.table(variable = variable, low = as.character(low), high = as.character(high),
             risk_low = mean(p0), risk_high = mean(p1), risk_difference = rd,
             rd_se = se, rd_lcl = rd - 1.96 * se, rd_ucl = rd + 1.96 * se)
}

d <- as.data.table(read_parquet(input))
d <- d[cnes_match_type != "unmatched" & cross_system_primary_eligible == TRUE & prior_window_complete == TRUE]
d[, `:=`(
  y = binary(definitive_completion_primary), age = as.numeric(age_years),
  female = as.numeric(as.character(SEXO) == "3"), race = fmiss(RACA_COR),
  admission_type = fmiss(CAR_INT), cholangitis = binary(cholangitis_any),
  obstruction = binary(biliary_obstruction_any), pancreatitis = binary(biliary_pancreatitis_any),
  cirrhosis = binary(cirrhosis_any), coagulopathy = binary(coagulopathy_any),
  structural_readiness = as.numeric(structural_readiness_count),
  local_prior_capability = binary(local_prior_capability),
  destination_prior_operations = as.numeric(destination_prior_operations),
  log_destination_prior_operations = log1p(as.numeric(destination_prior_operations)),
  flow = factor(care_flow_category, levels = c("same_municipality", "within_state_external", "interstate")),
  log_distance = log1p(as.numeric(centroid_distance_km)),
  distance_missing = as.numeric(!is.finite(as.numeric(centroid_distance_km))),
  outward_dependence = as.numeric(origin_outward_dependence),
  origin_hhi = as.numeric(origin_destination_hhi),
  hhi_missing = as.numeric(!is.finite(as.numeric(origin_destination_hhi))),
  log_catchment = log1p(as.numeric(destination_catchment_breadth)),
  destination_external_share = as.numeric(destination_external_share),
  year_f = fmiss(file_year), state = fmiss(file_uf),
  origin_cluster = interaction(file_year, origin_muni, drop = TRUE),
  destination_cluster = interaction(file_year, destination_hospital, drop = TRUE)
)]
imputation_vars <- c("age", "structural_readiness", "destination_prior_operations",
                     "log_destination_prior_operations", "log_distance", "outward_dependence",
                     "origin_hhi", "log_catchment", "destination_external_share")
pre_imputation_missingness <- rbindlist(lapply(imputation_vars, function(v) {
  d[, .(variable = v, missing_before_imputation = sum(!is.finite(get(v))), records = .N), by = cohort]
}))
for (v in imputation_vars) {
  replacement <- median(d[[v]], na.rm = TRUE)
  if (!is.finite(replacement)) stop(sprintf("No finite imputation value for %s", v))
  set(d, which(!is.finite(d[[v]])), v, replacement)
}
d[!is.finite(distance_missing), distance_missing := 1]
d[!is.finite(hhi_missing), hhi_missing := 1]
d[!is.finite(female), female := median(female, na.rm = TRUE)]
d <- d[!is.na(flow)]

analysis_vars <- c(
  "y", "age", "female", "race", "admission_type", "cholangitis", "obstruction",
  "pancreatitis", "cirrhosis", "coagulopathy", "structural_readiness",
  "local_prior_capability", "log_destination_prior_operations", "flow",
  "log_distance", "distance_missing", "outward_dependence", "origin_hhi",
  "hhi_missing", "log_catchment", "destination_external_share", "year_f", "state",
  "origin_cluster", "destination_cluster"
)
d[, analysis_complete := complete.cases(.SD), .SDcols = analysis_vars]
missingness <- d[, .(
  records = .N,
  completion = mean(y),
  local_flow = mean(flow == "same_municipality"),
  interstate_flow = mean(flow == "interstate")
), by = .(cohort, analysis_complete)]
missingness <- merge(
  CJ(cohort = unique(d$cohort), analysis_complete = c(FALSE, TRUE), unique = TRUE),
  missingness, by = c("cohort", "analysis_complete"), all.x = TRUE
)
missingness[is.na(records), `:=`(records = 0L, completion = 0, local_flow = 0, interstate_flow = 0)]
missing_variable_counts <- rbindlist(lapply(analysis_vars, function(v) {
  d[, .(variable = v, missing_after_handling = sum(is.na(get(v))), records = .N), by = cohort]
}))
missing_variable_counts <- merge(missing_variable_counts, pre_imputation_missingness,
                                 by = c("cohort", "variable", "records"), all.x = TRUE)
missing_variable_counts[is.na(missing_before_imputation), missing_before_imputation := 0L]
d <- d[analysis_complete == TRUE]

patient <- "ns(age,3) + female + race + admission_type + cholangitis + obstruction + pancreatitis + cirrhosis + coagulopathy + year_f + state"
forms <- list(
  B0_patient = as.formula(paste("y ~", patient)),
  B1_structure = as.formula(paste("y ~ structural_readiness +", patient)),
  B2_prior_performance = as.formula(paste("y ~ structural_readiness + local_prior_capability + log_destination_prior_operations +", patient)),
  B3_care_flow = as.formula(paste("y ~ structural_readiness + local_prior_capability + log_destination_prior_operations + flow + log_distance + distance_missing + outward_dependence + origin_hhi + hhi_missing + log_catchment + destination_external_share +", patient))
)

coef_rows <- list(); contrast_rows <- list(); model_rows <- list(); flow_counts <- list()
for (cc in c("acute_cholecystitis", "choledocholithiasis")) {
  x <- copy(d[cohort == cc])
  stopifnot(nrow(x) > 1000, uniqueN(x$origin_cluster) > 20, uniqueN(x$destination_cluster) > 20)
  flow_counts[[cc]] <- x[, .(n = .N, completion = mean(y)), by = .(cohort, care_flow_category)]
  for (mn in names(forms)) {
    fit <- suppressWarnings(glm(forms[[mn]], data = x, family = binomial(), control = glm.control(maxit = 75)))
    stopifnot(isTRUE(fit$converged))
    stopifnot(nobs(fit) == nrow(x))
    vcs <- glm_cluster_vcovs(fit, x$origin_cluster, x$destination_cluster)
    vc_dest <- vcs$destination_hospital
    vc_orig <- vcs$residence_municipality
    vc_two <- vcs$two_way_origin_destination
    for (vc_name in c("two_way_origin_destination", "destination_hospital", "residence_municipality")) {
      vc <- switch(vc_name, two_way_origin_destination = vc_two,
                   destination_hospital = vc_dest, residence_municipality = vc_orig)
      se <- sqrt(diag(vc)); b <- coef(fit)
      coef_rows[[paste(cc, mn, vc_name)]] <- data.table(
        cohort = cc, model = mn, inference = vc_name, term = names(b), estimate = b,
        se = se, odds_ratio = exp(b), or_lcl = exp(b - 1.96 * se),
        or_ucl = exp(b + 1.96 * se), p_value = 2 * pnorm(-abs(b / se))
      )
    }
    model_rows[[paste(cc, mn)]] <- data.table(
      cohort = cc, model = mn, n = nrow(x), completed = sum(x$y),
      origins = uniqueN(x$origin_cluster), destinations = uniqueN(x$destination_cluster),
      parameters = length(coef(fit)), aic = AIC(fit), converged = fit$converged
    )
    if (mn == "B3_care_flow") {
      q_struct <- quantile(x$structural_readiness, c(.25, .75), na.rm = TRUE)
      if (isTRUE(all.equal(unname(q_struct[[1]]), unname(q_struct[[2]])))) {
        q_struct <- quantile(x$structural_readiness, c(.10, .90), na.rm = TRUE)
      }
      if (isTRUE(all.equal(unname(q_struct[[1]]), unname(q_struct[[2]])))) {
        q_struct <- quantile(x$structural_readiness, c(.01, .99), na.rm = TRUE)
      }
      if (isTRUE(all.equal(unname(q_struct[[1]]), unname(q_struct[[2]])))) {
        stop(sprintf("Structural readiness has no usable contrast in %s", cc))
      }
      q_ops <- quantile(x$destination_prior_operations, c(.25, .75), na.rm = TRUE)
      cs <- list(
        standardize(fit, vc_two, x, "structural_readiness", q_struct[[1]], q_struct[[2]]),
        standardize(fit, vc_two, x, "local_prior_capability", 0, 1),
        standardize(fit, vc_two, x, "log_destination_prior_operations", log1p(q_ops[[1]]), log1p(q_ops[[2]])),
        standardize(fit, vc_two, x, "flow", "same_municipality", "within_state_external"),
        standardize(fit, vc_two, x, "flow", "same_municipality", "interstate")
      )
      z <- rbindlist(cs, fill = TRUE)
      z[, `:=`(cohort = cc, model = mn, inference = "two_way_origin_destination")]
      contrast_rows[[cc]] <- z
    }
  }
}

coefs <- rbindlist(coef_rows, fill = TRUE)
coefs[, fdr_p_value := NA_real_]
coefs[model == "B3_care_flow" & inference == "two_way_origin_destination" &
        grepl("structural_readiness|local_prior_capability|log_destination_prior_operations|flow|log_distance|outward_dependence|origin_hhi|log_catchment|destination_external_share", term),
      fdr_p_value := p.adjust(p_value, method = "BH")]
contrasts <- rbindlist(contrast_rows, fill = TRUE)
models <- rbindlist(model_rows, fill = TRUE)
flows <- rbindlist(flow_counts, fill = TRUE)

# Reviewer-requested care-flow sensitivities.  The first uses a linear
# probability model after exact within-origin-year demeaning, so care-flow
# coefficients are identified only by admissions from the same residence
# municipality and year that reached different destination categories.  The
# covariance is clustered on both origin-year and destination-year.  The
# second excludes interstate movements, the records most compatible with
# long-distance referral selection, and re-estimates the within-state contrast.
within_origin_lpm <- function(x, cohort_name) {
  x <- copy(x)
  x[, origin_flow_levels := uniqueN(flow), by = origin_cluster]
  x <- x[origin_flow_levels > 1]
  stopifnot(nrow(x) > 1000)
  X <- model.matrix(~ flow + structural_readiness + local_prior_capability +
                      log_destination_prior_operations + log_distance + distance_missing +
                      outward_dependence + origin_hhi + hhi_missing + log_catchment +
                      destination_external_share + ns(age,3) + female + race + admission_type +
                      cholangitis + obstruction + pancreatitis + cirrhosis + coagulopathy + year_f,
                    data = x)
  gid <- as.integer(factor(x$origin_cluster))
  counts <- tabulate(gid)
  Xbar <- rowsum(X, gid, reorder = FALSE) / counts
  Xw <- X - Xbar[gid, , drop = FALSE]
  ybar <- rowsum(matrix(x$y, ncol = 1), gid, reorder = FALSE)[, 1] / counts
  yw <- x$y - ybar[gid]
  keep <- apply(Xw, 2, function(z) is.finite(sd(z)) && sd(z) > 1e-10)
  Xw <- Xw[, keep, drop = FALSE]
  fit <- lm.fit(Xw, yw)
  ok <- is.finite(fit$coefficients)
  Xv <- Xw[, ok, drop = FALSE]
  beta <- fit$coefficients[ok]
  names(beta) <- colnames(Xv)
  resid <- yw - drop(Xv %*% beta)
  score <- Xv * resid
  bread <- qr.solve(crossprod(Xv), diag(ncol(Xv)), tol = 1e-10)
  meat_origin <- cluster_meat(score, x$origin_cluster)
  meat_destination <- cluster_meat(score, x$destination_cluster)
  meat_intersection <- cluster_meat(score, interaction(x$origin_cluster, x$destination_cluster, drop = TRUE))
  vc <- bread %*% (meat_origin + meat_destination - meat_intersection) %*% bread
  wanted <- intersect(c("flowwithin_state_external", "flowinterstate"), names(beta))
  se <- sqrt(pmax(0, diag(vc)[match(wanted, names(beta))]))
  data.table(
    cohort = cohort_name, sensitivity = "within_origin_year_lpm",
    contrast = wanted, n = nrow(x), origins = uniqueN(x$origin_cluster),
    destinations = uniqueN(x$destination_cluster), estimate = beta[wanted], se = se,
    lcl = beta[wanted] - 1.96 * se, ucl = beta[wanted] + 1.96 * se,
    interpretation = "within residence-municipality-year adjusted risk-difference association; two-way clustered; not a referral effect"
  )
}

sensitivity_rows <- list()
for (cc in c("acute_cholecystitis", "choledocholithiasis")) {
  x <- copy(d[cohort == cc])
  sensitivity_rows[[paste(cc, "within")]] <- within_origin_lpm(x, cc)
  xr <- droplevels(x[flow != "interstate"])
  fit <- suppressWarnings(glm(forms$B3_care_flow, data = xr, family = binomial(), control = glm.control(maxit = 75)))
  if (!isTRUE(fit$converged)) stop("Same-state sensitivity did not converge")
  stopifnot(nobs(fit) == nrow(xr))
  vc <- glm_cluster_vcovs(fit, xr$origin_cluster, xr$destination_cluster)$two_way_origin_destination
  z <- standardize(fit, vc, xr, "flow", "same_municipality", "within_state_external")
  sensitivity_rows[[paste(cc, "same_state")]] <- data.table(
    cohort = cc, sensitivity = "exclude_interstate", contrast = "within_state_external_vs_same_municipality",
    n = nrow(xr), origins = uniqueN(xr$origin_cluster), destinations = uniqueN(xr$destination_cluster),
    estimate = z$risk_difference, se = z$rd_se, lcl = z$rd_lcl, ucl = z$rd_ucl,
    interpretation = "marginal care-flow association after excluding interstate movements; two-way clustered"
  )
}
sensitivity <- rbindlist(sensitivity_rows, fill = TRUE)
fwrite(coefs, file.path(outdir, "brazil_access_block_coefficients.csv"))
fwrite(contrasts, file.path(outdir, "brazil_access_marginal_contrasts.csv"))
fwrite(models, file.path(qcdir, "brazil_access_model_diagnostics.csv"))
fwrite(missingness, file.path(qcdir, "brazil_access_complete_case_diagnostics.csv"))
fwrite(missing_variable_counts, file.path(qcdir, "brazil_access_variable_missingness.csv"))
fwrite(flows, file.path(outdir, "brazil_care_flow_descriptives.csv"))
fwrite(sensitivity, file.path(outdir, "brazil_access_sensitivity.csv"))

primary <- coefs[model == "B3_care_flow" & inference == "two_way_origin_destination"]
gate <- data.table(
  status = ifelse(all(models$converged) && all(models$completed / models$parameters >= 20) &&
                    all(is.finite(contrasts$risk_difference)) && nrow(primary) > 0 &&
                    all(is.finite(sensitivity$estimate)), "PASS", "REVIEW"),
  cohorts = uniqueN(models$cohort), models = nrow(models),
  minimum_events_per_parameter = min(models$completed / models$parameters),
  uncertainty = "Cameron-Gelbach-Miller two-way cluster-robust covariance by residence municipality and destination hospital; each one-way covariance reported as sensitivity",
  interpretation = "residence-to-destination care-flow association; within-origin and interstate-exclusion sensitivities; not a verified interhospital transfer pathway or causal mediation"
)
fwrite(gate, file.path(qcdir, "brazil_access_gate.csv"))
cat(sprintf("BRAZIL_ACCESS_%s contrasts=%d models=%d\n", gate$status, nrow(contrasts), nrow(models)))
if (gate$status != "PASS") quit(status = 2)
