# Aggregate Source Data

This directory contains non-identifying aggregate data used to prepare the displayed figures and selected summary tables. The files contain model estimates, cumulative-incidence summaries, quantitative-bias-analysis grids, distribution summaries, geographic contrasts, and policy-scenario outputs.

The directory does not contain patient-level NRD or Brazilian admission records. Fine-grained NRD hospital-year files with exact small sample or event counts are intentionally excluded from the public release. Approved NRD users can regenerate those intermediate files by running `scripts/final/stage23_nrd_extensions.R` within their secure environment.

The included aggregate files must not be interpreted as hospital rankings, causal effects, verified patient-level referral trajectories, or substitutes for the controlled source data. Variable definitions and interpretation limits are documented in the manuscript, figure contract, and statistical analysis plans.
