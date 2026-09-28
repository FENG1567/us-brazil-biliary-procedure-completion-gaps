#!/usr/bin/env Rscript

# Brazil extension: road travel time, geographic equity, and
# associational capacity/referral policy scenarios. No causal effects are claimed.

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(splines)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("usage: stage23_brazil_geoequity_policy.R INPUT_PARQUET OUTDIR")
input <- args[[1L]]
outdir <- args[[2L]]
if (!file.exists(input)) stop("missing input parquet: ", input)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
setDTthreads(min(8L, max(1L, getDTthreads())))
set.seed(20260923L)

bin <- function(x) as.integer(tolower(as.character(x)) %in% c("1", "true", "t", "yes"))
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
fref <- function(x, preferred = NULL) {
  z <- as.character(x)
  z[is.na(z) | z == "" | z %in% c("-99", "-9", "-8", "-6")] <- "unknown"
  tab <- sort(table(z), decreasing = TRUE)
  first <- if (!is.null(preferred) && preferred %in% names(tab)) preferred else names(tab)[1L]
  factor(z, levels = c(first, sort(setdiff(names(tab), first))))
}
road_band <- function(minutes, same_municipality) {
  # data.table::fifelse requires its test to be scalar or the full vector
  # length; referral scenarios deliberately pass a scalar FALSE.
  same_municipality <- rep_len(same_municipality, length(minutes))
  factor(
    fifelse(same_municipality, "same_municipality",
      fifelse(minutes <= 60, "external_le_60",
        fifelse(minutes <= 120, "external_61_120",
          fifelse(minutes <= 240, "external_121_240", "external_gt_240")))),
    levels = c("same_municipality", "external_le_60", "external_61_120", "external_121_240", "external_gt_240")
  )
}

cluster_meat <- function(score, cluster) {
  g <- as.integer(factor(cluster))
  G <- max(g)
  if (G <= 1L) stop("insufficient clusters")
  sums <- rowsum(score, g, reorder = FALSE)
  crossprod(sums) * (G / (G - 1))
}
nearest_psd <- function(v) {
  v <- (v + t(v)) / 2
  ee <- eigen(v, symmetric = TRUE)
  ee$values[ee$values < 0] <- 0
  ee$vectors %*% (ee$values * t(ee$vectors))
}

fit_cluster_logit <- function(data, formula, weight, origin_cluster, destination_cluster) {
  frame <- model.frame(formula, data = data, na.action = na.fail)
  y <- model.response(frame)
  full <- model.matrix(formula, frame)
  w <- weight / mean(weight)
  qr_fit <- qr(sqrt(w) * full, tol = 1e-10)
  keep <- sort(qr_fit$pivot[seq_len(qr_fit$rank)])
  x <- full[, keep, drop = FALSE]
  fit <- suppressWarnings(glm.fit(x, y, family = binomial(), weights = w, control = glm.control(maxit = 100L)))
  if (!isTRUE(fit$converged) || any(!is.finite(fit$coefficients))) stop("cluster logit did not converge")
  mu <- pmin(pmax(fit$fitted.values, 1e-8), 1 - 1e-8)
  info <- crossprod(x, x * as.numeric(w * mu * (1 - mu)))
  bread <- solve(info)
  score <- x * as.numeric(w * (y - mu))
  meat_origin <- cluster_meat(score, origin_cluster)
  meat_destination <- cluster_meat(score, destination_cluster)
  meat_intersection <- cluster_meat(score, interaction(origin_cluster, destination_cluster, drop = TRUE))
  vcov_two_way <- nearest_psd(bread %*% (meat_origin + meat_destination - meat_intersection) %*% bread)
  list(
    coefficients = fit$coefficients, vcov = vcov_two_way, formula = formula,
    kept_columns = colnames(x), converged = fit$converged, rank = ncol(x),
    full_columns = ncol(full), fitted = mu, n = nrow(data), events = sum(y),
    origin_clusters = uniqueN(origin_cluster), destination_clusters = uniqueN(destination_cluster)
  )
}

design_for <- function(model, newdata) {
  matrix_full <- model.matrix(model$formula, model.frame(model$formula, data = newdata, na.action = na.fail))
  absent <- setdiff(model$kept_columns, colnames(matrix_full))
  if (length(absent)) stop("new-data design is missing columns: ", paste(absent, collapse = ", "))
  matrix_full[, model$kept_columns, drop = FALSE]
}
standardized_prediction <- function(model, newdata, weights) {
  x <- design_for(model, newdata)
  eta <- drop(x %*% model$coefficients)
  p <- plogis(eta)
  risk <- weighted.mean(p, weights)
  gradient <- colSums(x * as.numeric(weights * p * (1 - p))) / sum(weights)
  list(risk = risk, gradient = gradient, probability = p, design = x)
}
contrast_result <- function(label, low, high, vcov, cohort, analysis, note) {
  difference <- high$risk - low$risk
  gradient <- high$gradient - low$gradient
  se <- sqrt(max(0, drop(t(gradient) %*% vcov %*% gradient)))
  data.table(
    cohort = cohort, analysis = analysis, contrast = label,
    risk_low = low$risk, risk_high = high$risk, risk_difference = difference,
    rd_se = se, rd_lcl = difference - 1.96 * se, rd_ucl = difference + 1.96 * se,
    interpretation = note
  )
}
scenario_result <- function(label, baseline, scenario_risk, scenario_gradient, vcov, cohort, eligible_n, allocation_fraction, note) {
  difference <- scenario_risk - baseline$risk
  gradient <- scenario_gradient - baseline$gradient
  se <- sqrt(max(0, drop(t(gradient) %*% vcov %*% gradient)))
  data.table(
    cohort = cohort, scenario = label, baseline_risk = baseline$risk,
    scenario_risk = scenario_risk, risk_difference = difference,
    rd_se = se, rd_lcl = difference - 1.96 * se, rd_ucl = difference + 1.96 * se,
    eligible_admissions = eligible_n, mean_capacity_allocation_fraction = allocation_fraction,
    interpretation = note
  )
}

raw <- as.data.table(read_parquet(input))
raw[, `:=`(
  eligible = bin(cross_system_primary_eligible) == 1L & bin(prior_window_complete) == 1L,
  outcome = bin(definitive_completion_primary),
  age10 = (num(age_years) - 65) / 10,
  female = as.integer(as.character(SEXO) %in% c("3", "F", "2")),
  cholangitis = bin(cholangitis_any), obstruction = bin(biliary_obstruction_any),
  pancreatitis = bin(biliary_pancreatitis_any), cirrhosis = bin(cirrhosis_any),
  coagulopathy = bin(coagulopathy_any),
  structural_readiness = num(structural_readiness_count),
  local_capability = bin(local_prior_capability),
  destination_capability = bin(destination_prior_capability),
  log_destination_operations = log1p(pmax(num(destination_prior_operations), 0)),
  outward_dependence = num(origin_outward_dependence), origin_hhi = num(origin_destination_hhi),
  log_catchment = log1p(pmax(num(destination_catchment_breadth), 0)),
  destination_external_share = num(destination_external_share),
  road_minutes = num(road_duration_minutes), road_km = num(road_distance_km),
  route_ok = as.character(route_status) %in% c("ok", "same_municipality"),
  route_flag = bin(route_plausibility_flag),
  same_municipality = as.character(origin_muni) == as.character(destination_muni),
  gdp_q = factor(as.integer(num(origin_gdp_pc_quartile)), levels = 1:4, labels = paste0("Q", 1:4)),
  population_q = factor(as.integer(num(origin_population_quartile)), levels = 1:4, labels = paste0("Q", 1:4)),
  macroregion = fref(origin_macroregion, "Southeast"),
  year_f = fref(file_year), month_f = fref(file_month), race = fref(RACA_COR),
  admission_type = fref(CAR_INT), complexity = fref(COMPLEX), specialty = fref(ESPEC),
  flow = fref(care_flow_category, "same_municipality"),
  origin_cluster = interaction(file_year, origin_muni, drop = TRUE),
  destination_cluster = interaction(file_year, destination_hospital, drop = TRUE),
  raw_weight = 1.0,
  year = as.integer(num(file_year)), origin = as.character(origin_muni),
  destination = as.character(destination_muni), destination_id = as.character(destination_hospital)
)]
raw[!is.finite(age10), age10 := median(age10, na.rm = TRUE)]
raw[!is.finite(structural_readiness), structural_readiness := 0]
raw[, hhi_missing := as.integer(!is.finite(origin_hhi))]
raw[, origin_hhi := fifelse(is.finite(origin_hhi), origin_hhi, median(origin_hhi[is.finite(origin_hhi)], na.rm = TRUE)), by = .(cohort, year)]
raw[!is.finite(outward_dependence), outward_dependence := 0]
raw[!is.finite(destination_external_share), destination_external_share := 0]
raw[, travel_band := road_band(raw$road_minutes, raw$same_municipality)]

flow_check <- raw[, .(
  records = .N,
  eligible_records = sum(eligible),
  invalid_route = sum(eligible & !route_ok),
  route_plausibility_excluded = sum(eligible & route_ok & route_flag == 1L),
  missing_equity = sum(eligible & (is.na(gdp_q) | is.na(population_q) | is.na(macroregion))),
  analytic_records = sum(eligible & route_ok & route_flag == 0L & !is.na(gdp_q) & !is.na(population_q))
), by = cohort]
fwrite(flow_check, file.path(outdir, "stage23_brazil_geoequity_flow_check.csv"))

d <- raw[eligible & route_ok & route_flag == 0L & !is.na(gdp_q) & !is.na(population_q)]
d <- droplevels(d)
if (nrow(d) < 10000L) stop("insufficient geoequity analytic records")

patient_terms <- c(
  "age10", "I(age10^2)", "I(age10^3)", "female", "cholangitis", "obstruction", "pancreatitis", "cirrhosis", "coagulopathy",
  "year_f", "month_f", "race", "admission_type", "complexity", "specialty"
)
network_terms <- c(
  "travel_band", "gdp_q", "population_q", "macroregion", "local_capability", "destination_capability",
  "log_destination_operations", "I(log_destination_operations^2)", "structural_readiness", "I(structural_readiness^2)",
  "flow", "outward_dependence",
  "origin_hhi", "hhi_missing", "log_catchment", "destination_external_share"
)
formula_main <- as.formula(paste("outcome ~", paste(c(network_terms, patient_terms), collapse = " + ")))

models <- list(); coefficients <- list(); contrasts <- list(); scenarios <- list(); diagnostics <- list()
for (cohort_name in c("acute_cholecystitis", "choledocholithiasis")) {
  x <- droplevels(copy(d[cohort == cohort_name]))
  if (nrow(x) < 5000L) stop("insufficient cohort records: ", cohort_name)
  model <- fit_cluster_logit(x, formula_main, x$raw_weight, x$origin_cluster, x$destination_cluster)
  models[[cohort_name]] <- model
  se <- sqrt(diag(model$vcov))
  coefficients[[cohort_name]] <- data.table(
    cohort = cohort_name, term = names(model$coefficients), estimate = model$coefficients, se = se,
    lcl = model$coefficients - 1.96 * se, ucl = model$coefficients + 1.96 * se,
    odds_ratio = exp(model$coefficients), or_lcl = exp(model$coefficients - 1.96 * se),
    or_ucl = exp(model$coefficients + 1.96 * se), inference = "two_way_origin_destination_cluster_robust"
  )
  diagnostics[[cohort_name]] <- data.table(
    cohort = cohort_name, n = model$n, events = model$events, origins = model$origin_clusters,
    destinations = model$destination_clusters, rank = model$rank, full_columns = model$full_columns,
    outcome_rate = mean(x$outcome), route_median_minutes = median(x[same_municipality == FALSE]$road_minutes),
    route_p90_minutes = quantile(x[same_municipality == FALSE]$road_minutes, .90),
    road_centroid_spearman = cor(x[same_municipality == FALSE]$road_minutes, x[same_municipality == FALSE]$centroid_distance_km,
      method = "spearman", use = "complete.obs"),
    interpretation = "associational model with two-way clustered inference"
  )

  # Common-population marginal contrasts.
  compare_factor <- function(variable, low_level, high_level, label) {
    low_data <- copy(x); high_data <- copy(x)
    low_data[[variable]] <- factor(low_level, levels = levels(x[[variable]]))
    high_data[[variable]] <- factor(high_level, levels = levels(x[[variable]]))
    low <- standardized_prediction(model, low_data, x$raw_weight)
    high <- standardized_prediction(model, high_data, x$raw_weight)
    contrast_result(label, low, high, model$vcov, cohort_name, "geographic_equity",
      "model-standardized association in the measured-covariate population; not a causal effect")
  }
  contrasts[[paste0(cohort_name, "_gdp")]] <- compare_factor("gdp_q", "Q1", "Q4", "GDP_per_capita_Q4_vs_Q1")
  contrasts[[paste0(cohort_name, "_pop")]] <- compare_factor("population_q", "Q1", "Q4", "population_Q4_vs_Q1")
  contrasts[[paste0(cohort_name, "_travel")]] <- compare_factor("travel_band", "external_le_60", "external_gt_240", "road_time_gt240_vs_le60_minutes")
  for (region in setdiff(levels(x$macroregion), "Southeast")) {
    contrasts[[paste(cohort_name, "region", region, sep = "_")]] <- compare_factor(
      "macroregion", "Southeast", region, paste0(region, "_vs_Southeast")
    )
  }
  low_local <- copy(x); high_local <- copy(x)
  low_local[, local_capability := 0L]; high_local[, local_capability := 1L]
  contrasts[[paste0(cohort_name, "_local")]] <- contrast_result(
    "local_prior_capability_yes_vs_no",
    standardized_prediction(model, low_local, x$raw_weight),
    standardized_prediction(model, high_local, x$raw_weight),
    model$vcov, cohort_name, "network_capability",
    "model-standardized association; capability status is not randomized"
  )

  # Origin-year alternative: the fastest observed capable, non-interstate
  # destination for that origin and year.  Only observed network options enter.
  alternatives <- x[
    destination_capability == 1L & !same_municipality & flow != "interstate" & is.finite(road_minutes),
    .SD[order(road_minutes, -log_destination_operations)][1L],
    by = .(origin, year)
  ][, .(
    origin, year,
    alt_destination_id = destination_id,
    alt_destination = destination,
    alt_road_minutes = road_minutes,
    alt_structural_readiness = structural_readiness,
    alt_log_destination_operations = log_destination_operations,
    alt_log_catchment = log_catchment,
    alt_destination_external_share = destination_external_share
  )]
  capacity <- x[, .(observed_destination_volume = .N), by = .(year, destination_id)]
  alternatives <- merge(alternatives, capacity, by.x = c("year", "alt_destination_id"), by.y = c("year", "destination_id"), all.x = TRUE)
  alternatives[, additional_capacity := pmax(1L, floor(.20 * observed_destination_volume))]
  x <- merge(x, alternatives, by = c("origin", "year"), all.x = TRUE, sort = FALSE)
  x[, referral_eligible := local_capability == 0L & destination_capability == 0L & !is.na(alt_destination_id) &
      is.finite(alt_road_minutes) & (same_municipality | alt_road_minutes <= road_minutes)]
  demand <- x[referral_eligible == TRUE, .(demand = .N, additional_capacity = max(additional_capacity)),
    by = .(year, alt_destination_id)]
  demand[, allocation_fraction := pmin(1, additional_capacity / pmax(demand, 1))]
  x[demand, allocation_fraction := i.allocation_fraction, on = .(year, alt_destination_id)]
  x[is.na(allocation_fraction), allocation_fraction := 0]

  baseline <- standardized_prediction(model, x, x$raw_weight)
  scenario_a_data <- copy(x)
  scenario_a_data[local_capability == 0L, local_capability := 1L]
  scenario_a <- standardized_prediction(model, scenario_a_data, x$raw_weight)
  scenarios[[paste0(cohort_name, "_A")]] <- scenario_result(
    "A_local_capability_expansion", baseline, scenario_a$risk, scenario_a$gradient,
    model$vcov, cohort_name, sum(x$local_capability == 0L), 1,
    "associational standardization setting observed local capability to present; not a causal forecast"
  )

  referral_data <- copy(x)
  referral_data[referral_eligible == TRUE, `:=`(
    destination_capability = 1L,
    structural_readiness = alt_structural_readiness,
    log_destination_operations = alt_log_destination_operations,
    log_catchment = alt_log_catchment,
    destination_external_share = alt_destination_external_share,
    road_minutes = alt_road_minutes,
    same_municipality = FALSE,
    travel_band = road_band(alt_road_minutes, FALSE),
    flow = factor("within_state_external", levels = levels(flow))
  )]
  referral_full <- standardized_prediction(model, referral_data, x$raw_weight)
  allocation <- x$allocation_fraction
  scenario_b_probability <- baseline$probability + allocation * (referral_full$probability - baseline$probability)
  scenario_b_risk <- weighted.mean(scenario_b_probability, x$raw_weight)
  scenario_b_gradient <- colSums(
    baseline$design * as.numeric(x$raw_weight * baseline$probability * (1 - baseline$probability)) +
      allocation * (
        referral_full$design * as.numeric(x$raw_weight * referral_full$probability * (1 - referral_full$probability)) -
          baseline$design * as.numeric(x$raw_weight * baseline$probability * (1 - baseline$probability))
      )
  ) / sum(x$raw_weight)
  scenarios[[paste0(cohort_name, "_B")]] <- scenario_result(
    "B_capacity_constrained_referral_reconfiguration", baseline, scenario_b_risk, scenario_b_gradient,
    model$vcov, cohort_name, sum(x$referral_eligible), mean(allocation[x$referral_eligible]),
    "associational scenario using the fastest observed capable in-state destination and a 20% observed-volume capacity cap"
  )

  combined_data <- copy(referral_data)
  combined_data[local_capability == 0L, local_capability := 1L]
  combined_full <- standardized_prediction(model, combined_data, x$raw_weight)
  local_base <- scenario_a
  scenario_c_probability <- local_base$probability + allocation * (combined_full$probability - local_base$probability)
  scenario_c_risk <- weighted.mean(scenario_c_probability, x$raw_weight)
  local_derivative <- local_base$design * as.numeric(x$raw_weight * local_base$probability * (1 - local_base$probability))
  combined_derivative <- combined_full$design * as.numeric(x$raw_weight * combined_full$probability * (1 - combined_full$probability))
  scenario_c_gradient <- colSums(local_derivative + allocation * (combined_derivative - local_derivative)) / sum(x$raw_weight)
  scenarios[[paste0(cohort_name, "_C")]] <- scenario_result(
    "C_combined_capability_and_referral", baseline, scenario_c_risk, scenario_c_gradient,
    model$vcov, cohort_name, sum(x$local_capability == 0L | x$referral_eligible),
    mean(allocation[x$referral_eligible]),
    "combined associational scenario; not an intervention effect or guaranteed policy yield"
  )

  scenario_feasibility <- x[, .(
    cohort = cohort_name,
    analytic_admissions = .N,
    origins = uniqueN(origin),
    observed_destinations = uniqueN(destination_id),
    origins_with_observed_capable_alternative = uniqueN(origin[!is.na(alt_destination_id)]),
    referral_eligible_admissions = sum(referral_eligible),
    capacity_allocated_admission_equivalents = sum(allocation),
    median_allocation_fraction = median(allocation[referral_eligible]),
    policy_guard = "observed-network common-support and 20-percent destination-volume cap"
  )]
  fwrite(scenario_feasibility, file.path(outdir, paste0("stage23_", cohort_name, "_policy_feasibility.csv")))
  cat(sprintf("BRAZIL_MODEL_DONE cohort=%s n=%d\n", cohort_name, nrow(x))); flush.console()
}

coefficient_table <- rbindlist(coefficients, fill = TRUE)
contrast_table <- rbindlist(contrasts, fill = TRUE)
scenario_table <- rbindlist(scenarios, fill = TRUE)
diagnostic_table <- rbindlist(diagnostics, fill = TRUE)
fwrite(coefficient_table, file.path(outdir, "stage23_brazil_geoequity_model_coefficients.csv"))
fwrite(contrast_table, file.path(outdir, "stage23_brazil_geoequity_marginal_contrasts.csv"))
fwrite(scenario_table, file.path(outdir, "stage23_brazil_policy_scenarios.csv"))
fwrite(diagnostic_table, file.path(outdir, "stage23_brazil_geoequity_model_diagnostics.csv"))
saveRDS(models, file.path(outdir, "stage23_brazil_geoequity_models.rds"), compress = FALSE)

gate <- data.table(
  check = c(
    "two_models", "adequate_clusters", "full_rank_reduction_recorded", "finite_contrasts",
    "finite_scenarios", "route_exclusion_below_10pct", "policy_language_guard"
  ),
  pass = c(
    length(models) == 2L,
    all(diagnostic_table$origins > 20 & diagnostic_table$destinations > 20),
    all(diagnostic_table$rank <= diagnostic_table$full_columns),
    all(is.finite(contrast_table$risk_difference) & is.finite(contrast_table$rd_se)),
    all(is.finite(scenario_table$risk_difference) & is.finite(scenario_table$rd_se)),
    all(flow_check$route_plausibility_excluded / pmax(flow_check$eligible_records, 1) < .10),
    all(grepl("associational|not a causal", scenario_table$interpretation))
  )
)
fwrite(gate, file.path(outdir, "stage23_brazil_geoequity_policy_gate.csv"))
if (!all(gate$pass)) stop("Brazil geoequity/policy gate failed")
cat("BRAZIL_GEOEQUITY_POLICY_PASS\n")
