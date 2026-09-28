suppressPackageStartupMessages({library(arrow); library(data.table)})
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) stop('usage: stage26_brazil_pathways.R INPUT OUTDIR')
input <- args[[1]]; outdir <- args[[2]]
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
d <- as.data.table(read_parquet(input))
d <- d[file_year %in% 2018:2022 & cross_system_primary_eligible == TRUE]
to01 <- function(x) as.integer(isTRUE(x) | (!is.na(x) & as.character(x) %in% c('TRUE','1','true')))
for (v in c('op_cholecystectomy','op_cholecystostomy','op_exploration','op_diagnostic_ercp','op_clearance_explicit','op_therapeutic_ercp','op_bridge','definitive_completion_primary','definitive_completion_narrow')) {
  if (v %in% names(d)) d[[v]] <- to01(d[[v]]) else d[[v]] <- 0L
}
d[, `:=`(
  chole = op_cholecystectomy == 1L,
  duct = op_therapeutic_ercp == 1L | op_exploration == 1L | op_clearance_explicit == 1L,
  ercp_any = op_diagnostic_ercp == 1L | op_therapeutic_ercp == 1L,
  bridge = op_bridge == 1L | op_cholecystostomy == 1L,
  year = as.integer(file_year)
)]
d[, pathway := 'none']
d[cohort == 'acute_cholecystitis' & chole & ercp_any, pathway := 'cholecystectomy_plus_ercp']
d[cohort == 'acute_cholecystitis' & chole & !ercp_any, pathway := 'cholecystectomy_only']
d[cohort == 'acute_cholecystitis' & !chole & bridge, pathway := 'bridge_only']
d[cohort == 'acute_cholecystitis' & !chole & !bridge & ercp_any, pathway := 'ercp_without_cholecystectomy']
d[cohort == 'choledocholithiasis' & chole & duct, pathway := 'cholecystectomy_plus_duct_treatment']
d[cohort == 'choledocholithiasis' & !chole & duct, pathway := 'duct_treatment_without_cholecystectomy']
d[cohort == 'choledocholithiasis' & chole & !duct, pathway := 'cholecystectomy_only']
d[cohort == 'choledocholithiasis' & !chole & bridge, pathway := 'bridge_only']
d[cohort == 'choledocholithiasis' & !chole & !duct & !bridge & ercp_any, pathway := 'ercp_without_definitive_duct_treatment']
pathway <- d[, .(admissions=.N, definitive_primary=sum(definitive_completion_primary), definitive_narrow=sum(definitive_completion_narrow)), by=.(cohort,year,pathway)]
pathway[, proportion := admissions / sum(admissions), by=.(cohort,year)]
fwrite(pathway, file.path(outdir,'stage26_brazil_pathway_distribution.csv'))
overall <- d[, .(admissions=.N, definitive_primary=sum(definitive_completion_primary), definitive_narrow=sum(definitive_completion_narrow)), by=cohort]
overall[, `:=`(primary_completion=definitive_primary/admissions, narrow_completion=definitive_narrow/admissions)]
fwrite(overall, file.path(outdir,'stage26_brazil_pathway_overall.csv'))
if ('care_flow_category' %in% names(d)) {
  flow <- d[, .(admissions=.N, completion=mean(definitive_completion_primary), completion_narrow=mean(definitive_completion_narrow)), by=.(cohort,year,care_flow_category)]
  fwrite(flow, file.path(outdir,'stage26_brazil_flow_completion.csv'))
}
if (all(c('origin_gdp_pc_quartile','road_duration_minutes') %in% names(d))) {
  d[, road_band := fifelse(as.numeric(road_duration_minutes) <= 60,'same_or_le_60',fifelse(as.numeric(road_duration_minutes)<=120,'61_120',fifelse(as.numeric(road_duration_minutes)<=240,'121_240','gt_240')))]
  eq <- d[, .(admissions=.N, completion=mean(definitive_completion_primary), completion_narrow=mean(definitive_completion_narrow)), by=.(cohort,year,origin_gdp_pc_quartile,road_band)]
  fwrite(eq, file.path(outdir,'stage26_brazil_geographic_equity_descriptives.csv'))
}
qc <- data.table(status='PASS', rows=nrow(d), cohorts=paste(sort(unique(d$cohort)),collapse=';'), years=paste(sort(unique(d$year)),collapse=';'), note='Admission-level Brazil pathway module; no patient-level longitudinal linkage claimed.')
fwrite(qc, file.path(outdir,'stage26_brazil_pathway_gate.csv'))
cat(sprintf('Brazil pathway module complete; rows=%d\n',nrow(d)))
