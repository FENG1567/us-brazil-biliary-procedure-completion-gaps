# NRD cohort-harmonization note

## Issue identified during submission integration

The original NRD hospital-extension analysis read the earlier `nrd_day3_landmark_v2.0.parquet` and used its conservative `primary_day3_eligible` flag. That source incorporated coded cirrhosis/coagulopathy proxy exclusions. In contrast, the primary analysis used the broad day-3 risk set reconstructed from the `nrd_day3_landmark_v3.0.parquet` plus all five linked-admission files. The two approaches were therefore not the same analytic cohort, while the extension burden component drew the risk-difference distribution from the primary analysis.

This was identified during submission integration, before manuscript result integration. The correction was not selected according to the direction or magnitude of any extension result.

## Correction

The final NRD extension treats the broad primary-analysis cohort as authoritative. It uses:

- Landmark input: `nrd_day3_landmark_v3.0.parquet`;
- all five `nrd_*_index_patient_admissions_v1.0.parquet` linked-admission files;
- the primary analysis's day-3 broad eligibility, recorded early-completion definition, and full linked-admission reconstruction of biliary, nonbiliary, and inpatient-death events.

The code has hard gates requiring exact day-3 cohort counts of 233,764 acute-cholecystitis and 117,293 choledocholithiasis admissions. The hospital-year models, hospital-level correlation, peer-P75 shortfall, associated-burden resampling, and COVID sensitivity all use this harmonized risk set. The burden bootstrap continues to propagate the primary-analysis day-3 risk-difference bootstrap distribution, now matched to the same cohort definition.

## Consequence for prior extension outputs

The prior NRD hospital-extension outputs based on 204,723 acute-cholecystitis and 100,701 choledocholithiasis admissions are superseded and are retained only in an archive. They must not be used in the manuscript, figures, tables, supplement, or interpretation. Brazil extension outputs are unaffected.
