# Statistical analysis plan: extension analyses

Finalized on 2026-09-17 before inspecting the extension model estimates.

## Objective

Raise the manuscript's evidentiary ceiling without converting observational associations into causal claims. The extension analyses integrate: (1) the day-3 NRD target-trial framework and bias-control analyses; (2) a measurement-aware US-Brazil decomposition; (3) Brazil referral-network, road-travel-time, and geographic-equity analyses; (4) hospital-year risk-standardized variation and nationally weighted attributable burden; and (5) explicitly associational policy scenarios. The previously explored hospital-capacity-change quasi-experiment is excluded from all extension analyses and publication files.

## Analysis populations and prespecified estimands

### NRD day-3 emulation

- Source: 2018-2022 NRD day-3 landmark dataset used for the primary analysis.
- Primary cohort: adults with acute cholecystitis, April-September index admissions, at least 90 days of observed look-back, alive and free of the study outcome at day 3, no observed prior cholecystectomy, no malignant biliary obstruction, and no coded cirrhosis/coagulopathy proxy in the conservative eligibility definition.
- Exploratory cohort: choledocholithiasis using the prespecified gallbladder-status-aware completion phenotype.
- Exposure strategy: definitive completion by day 3 versus not completed by day 3.
- Primary outcome: observed biliary readmission from day 4 through day 93 after index admission.
- Primary contrast: overlap-population standardized absolute risk difference at 90 days. Hospital-year cluster bootstrap remains the uncertainty procedure.
- Negative-control outcome: nonbiliary nonelective readmission over the same window.
- Sensitivity analyses: alternative day-2/day-4 landmarks; quantitative unmeasured-confounding grid; COVID-era strata; multiple recurrent biliary admissions using Andersen-Gill with robust clustering and PWP total-time for the first three events.
- Interpretation: observational target-trial emulation. No randomized-treatment or causal-effect language.

### Cross-system completion analysis

- Common calendar period: 2021-2022.
- Systems: NRD (United States) and SIH-SUS/CNES (Brazil), modeled separately and aligned only through a common measured-covariate support population.
- Primary output: system-specific standardized completion proportions plus the residual measured-covariate-aligned contrast.
- Measurement guard: report restrictive and broader procedure-code definitions as a coding-definition range. Do not pool patient records or describe the residual contrast as a causal system effect.
- Decomposition: raw difference = measured-composition component + residual system contrast, with cluster bootstrap intervals.

### Brazil referral network, road access, and geographic equity

- Source: SIH-SUS/CNES care-flow data with IBGE 2022 municipality centroids; official IBGE municipality population and GDP; OpenStreetMap routing through the public OSRM driving profile with a persistent route cache.
- Unit: eligible index admission; municipality-pair road metrics are attached deterministically.
- Road metrics: routed driving time and distance between origin and destination municipality centroids. Same-municipality travel time is coded 0 and analyzed separately from routed external referrals.
- Routing quality gates: request status, snapped-coordinate displacement, route plausibility versus great-circle distance, missing-route fraction, deterministic cache, source URL/template, and retrieval timestamp.
- Equity dimensions: origin-municipality GDP per capita quartile, population-size quartile, macroregion, local prior procedural capability, and road-travel-time categories.
- Outcome: inpatient definitive completion under the prespecified disease-specific phenotype.
- Main estimands: marginal standardized completion probabilities and absolute risk differences across travel-time and socioeconomic strata, with destination-hospital clustering and sensitivity to origin clustering.
- Interpretation: geographic and network associations, not causal referral effects.

### Hospital-year variation and national burden

- Hospital-year cells are the reporting units; patient-level models retain prespecified case-mix adjustment and partial pooling.
- Report distributions of risk-standardized completion and observed readmission, median odds ratios or variance-partition summaries where estimable, and reliability/volume screens. Do not publish hospital league tables.
- National burden uses NRD discharge weights and propagates patient sampling plus hospital-year bootstrap uncertainty.
- Report numbers of incomplete definitive-treatment opportunities, observed biliary readmissions, and inpatient costs associated with scenarios that close 25%, 50%, 75%, or 100% of the modeled completion gap.
- Interpretation: potentially addressable associated burden and scenario-based estimates, not preventable causal counts or guaranteed savings.

### Policy scenarios

- Scenario A: improve local procedural capability among origin municipalities lacking prior capability.
- Scenario B: redirect eligible external referrals to an observed in-state destination with higher prior procedural capability and no longer modeled road time, subject to observed destination-volume capacity constraints.
- Scenario C: combine A and B.
- Standardize outcome predictions under each scenario while holding patient covariates fixed.
- Propagate model and cluster-resampling uncertainty. Include explicit feasibility counts and common-support restrictions.
- Interpretation: associational policy scenarios only; no counterfactual causal claim.

## Multiplicity and reporting hierarchy

The NRD acute-cholecystitis day-3 90-day biliary-readmission risk difference is the sole primary clinical contrast. Cross-system, geographic-equity, hospital-variation, and policy-scenario analyses are prespecified secondary or policy analyses. Choledocholithiasis, recurrent events, COVID strata, negative control, and quantitative-bias analyses are supportive or sensitivity analyses. Exact confidence intervals and standardized absolute risks take priority over p-values; no threshold-based claim of discovery will be made.

## Missing data and model diagnostics

- Preserve explicit missing categories for administrative categorical covariates as in the primary analysis.
- For municipality socioeconomic data, do not impute absent official values into the primary analysis; provide a complete-case sensitivity and missingness profile.
- For road routing failures, retain a missing-route indicator, report failure reasons, and run a great-circle-distance sensitivity.
- Require covariate balance after overlap weighting below an absolute standardized mean difference of 0.10, effective sample size reporting, positivity/common-support diagnostics, bootstrap success of at least 90%, and model convergence/finite-estimate checks.

## Publication integration

The main manuscript will contain exactly eight numbered display items. New analyses will be consolidated into multi-panel figures and dense but readable tables; recurrent-event, COVID, quantitative-bias, routing-QA, and extended model diagnostics will be placed in the Supplement. All source data and analysis scripts will be archived with the publication materials.
