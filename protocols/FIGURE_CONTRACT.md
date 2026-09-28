# Figure contract

Backend: R only (`ggplot2`, `patchwork`, `svglite`, `cairo_pdf`, `ragg`).
Target: Clinical Gastroenterology and Hepatology, 183 mm width, vector PDF as the submitted format.

## Figure 1

Core conclusion: the revised design uses separate disease cohorts, an explicit fixed day-3 landmark, and noninterchangeable roles for the three databases.
Archetype: schematic-led composite.
Panel map: database capability matrix; cohort attrition; landmark timeline; estimand/DAG schematic.
Reviewer risk: readers must see the common time zero, the day-4 to day-93 window, exclusions, later treatment handling, and nonpooling of patient-level records.

## Figure 2

Core conclusion: early disease-specific completion is associated with lower subsequent biliary admission, while mortality, composite, and nonbiliary results reveal the scope and remaining bias.
Archetype: quantitative grid with the primary risk difference as the hero panel.
Panel map: adjusted risks; risk differences for four outcomes; balance; landmark/definition sensitivities.
Reviewer risk: all intervals and target populations must be explicit; the nonbiliary outcome is a falsification signal, not a certification test.

## Figure 3

Core conclusion: recurrence subtypes and model-based burden scenarios are clinically material but bounded to the same eligible untreated common-support population.
Archetype: quantitative grid.
Panel map: recurrence subtypes; 25/50/75% scenario events; net 93-day inpatient facility cost; unmeasured-confounding/death scenarios.
Reviewer risk: cohorts cannot be added; scenario counts are associations, not preventable events; costs include the index stay and eligible follow-up biliary admissions.

## Figure 4

Core conclusion: Brazilian structural resources, prior performance, and care-flow patterns make distinct contributions, and cross-system comparisons remain intervals defined by common support and coding uncertainty.
Archetype: asymmetric mixed-modality quantitative figure.
Panel map: care-flow composition; structural/prior-performance/care-flow marginal contrasts with two-way clustered intervals; 2021-2022 cross-system coding bounds; aligned MIMIC-IV ancillary day-3 estimates.
Reviewer risk: care flow is not verified referral; Brazil uncertainty must reflect origin and destination; MIMIC is not validation; cross-system intervals are not causal system effects.

Statistics: 95% survey-linearization or design/cluster-bootstrap intervals for NRD; two-way origin-destination cluster-robust intervals for Brazil; record-bootstrap intervals for MIMIC; 2021-2022 stratified cluster-bootstrap intervals for cross-system bounds.
Source data: every plotted panel is written to a dedicated CSV under the figure source-data directory.
Image integrity: no photographic or microscopy content; all plotted marks are vector-native and derived from aggregate outputs.
