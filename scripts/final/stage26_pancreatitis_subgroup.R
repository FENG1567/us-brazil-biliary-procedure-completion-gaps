#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(arrow); library(data.table); library(survey); library(splines)})
options(survey.lonely.psu="adjust")
args <- commandArgs(trailingOnly=TRUE)
project <- normalizePath(if(length(args)) args[[1]] else ".", mustWork=TRUE)
outdir <- file.path(project, "results", "stage26_pancreatitis_subgroup")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
input <- file.path(project, "derived", "analysis_ready_v1.0", "stage11", "nrd_day3_landmark_v2.0.parquet")
stopifnot(file.exists(input))
binary <- function(x) {
  raw <- tolower(trimws(as.character(x)))
  z <- suppressWarnings(as.numeric(raw))
  z[raw %in% c("true","t","yes","y")] <- 1
  z[raw %in% c("false","f","no","n")] <- 0
  z[is.na(z)] <- 0
  as.numeric(z > 0)
}
fref <- function(x) { x <- as.character(x); x[is.na(x) | x=="" | x %in% c("-9","-8","-6","-99")] <- "unknown"; factor(x) }
d <- as.data.table(read_parquet(input))
# The day-3 primary risk set is the broad day-3 eligibility set:
# principal target diagnosis, urgent admission, no observed prior cholecystectomy
# for the acute-cholecystitis cohort, known treatment timing, survival and no
# biliary event through day 3. The primary analysis intentionally did not apply the
# coded cirrhosis/coagulopathy proxy to this risk set.
d <- d[binary(broad_day3_eligible)==1 & cohort %in% c("acute_cholecystitis","choledocholithiasis")]
d[, `:=`(
  treatment=binary(early_completion_day3), biliary90=binary(biliary_readmission_d3_93),
  pancreatitis=binary(dx_any_biliary_acute_pancreatitis), age=as.numeric(AGE_num),
  female=binary(FEMALE_num), weekend=binary(AWEEKEND_num), ed=binary(HCUP_ED_num),
  cholangitis=binary(dx_any_cholangitis), obstruction=binary(dx_any_biliary_obstruction),
  cirrhosis=binary(dx_any_cirrhosis_proxy), coagulopathy=binary(dx_any_coagulation_disorder_proxy),
  prior_volume=as.numeric(prior_90d_disease_volume), weight=as.numeric(DISCWT_num),
  payer=fref(PAY1), income=fref(ZIPINC_QRTL), rurality=fref(PL_NCHS),
  severity=fref(APRDRG_Severity), mortality=fref(APRDRG_Risk_Mortality),
  bedsize=fref(HOSP_BEDSIZE), teaching=fref(HOSP_UR_TEACH), control=fref(H_CONTRL),
  chole_state=fref(known_chole_state), year_f=fref(YEAR)
)]
d[!is.finite(age), age := median(age, na.rm=TRUE)]
d[!is.finite(prior_volume), prior_volume := 0]
d[!is.finite(weight) | weight <= 0, weight := median(weight[is.finite(weight) & weight > 0])]
weighted_rate <- function(x,w) sum(x*w)/sum(w)
rows <- list(); interactions <- list(); diag <- list()
for (cc in c("acute_cholecystitis","choledocholithiasis")) {
  x <- copy(d[cohort==cc])
  form <- treatment ~ ns(age,3)+female+weekend+ed+payer+income+rurality+cholangitis+obstruction+pancreatitis+cirrhosis+coagulopathy+severity+mortality+bedsize+teaching+control+chole_state+year_f+ns(log1p(prior_volume),3)
  if (cc=="choledocholithiasis") form <- treatment ~ ns(age,3)+female+weekend+ed+payer+income+rurality+cholangitis+obstruction+pancreatitis+cirrhosis+coagulopathy+severity+mortality+bedsize+teaching+control+chole_state+year_f+log1p(prior_volume)
  fit <- glm(form, data=x, family=binomial(), weights=weight/mean(weight))
  x[, ps := pmin(pmax(as.numeric(predict(fit,type="response")),0.001),0.999)]
  x[, ow := weight * fifelse(treatment==1,1-ps,ps)]
  x[, analysis_w := ow/mean(ow)]
  for (pp in 0:1) for (tt in 0:1) {
    y <- x[pancreatitis==pp & treatment==tt]
    rows[[length(rows)+1L]] <- data.table(cohort=cc, pancreatitis=pp, treatment=tt, n=nrow(y), weighted_n=sum(y$weight), completion_rate=weighted_rate(y$treatment,y$weight), biliary90_rate=weighted_rate(y$biliary90,y$weight), overlap_biliary90_rate=weighted_rate(y$biliary90,y$ow))
  }
  for (pp in 0:1) {
    y <- x[pancreatitis==pp]
    r1 <- weighted_rate(y$biliary90[y$treatment==1], y$ow[y$treatment==1]); r0 <- weighted_rate(y$biliary90[y$treatment==0], y$ow[y$treatment==0])
    rows[[length(rows)+1L]] <- data.table(cohort=cc, pancreatitis=pp, treatment=9, n=nrow(y), weighted_n=sum(y$weight), completion_rate=weighted_rate(y$treatment,y$weight), biliary90_rate=NA_real_, overlap_biliary90_rate=r1-r0)
  }
  subgroup_effect0 <- weighted_rate(x$biliary90[x$pancreatitis==0 & x$treatment==1], x$ow[x$pancreatitis==0 & x$treatment==1]) - weighted_rate(x$biliary90[x$pancreatitis==0 & x$treatment==0], x$ow[x$pancreatitis==0 & x$treatment==0])
  subgroup_effect1 <- weighted_rate(x$biliary90[x$pancreatitis==1 & x$treatment==1], x$ow[x$pancreatitis==1 & x$treatment==1]) - weighted_rate(x$biliary90[x$pancreatitis==1 & x$treatment==0], x$ow[x$pancreatitis==1 & x$treatment==0])
  des <- svydesign(ids=~hospital_cluster_id, strata=~NRD_STRATUM, weights=~analysis_w, data=x, nest=TRUE)
  intfit <- tryCatch(svyglm(biliary90 ~ treatment*pancreatitis, design=des, family=quasibinomial()), error=function(e) NULL)
  co <- if(is.null(intfit)) NULL else coef(summary(intfit)); term <- "treatment:pancreatitis"
  term_row <- if(!is.null(co) && term %in% rownames(co)) co[term,,drop=FALSE] else NULL
  interactions[[cc]] <- data.table(cohort=cc, interaction_term=term, interaction_log_odds=if(is.null(term_row)) NA_real_ else unname(term_row[1,"Estimate"]), interaction_se=if(is.null(term_row)) NA_real_ else unname(term_row[1,"Std. Error"]), interaction_p=if(is.null(term_row)) NA_real_ else unname(term_row[1,"Pr(>|t|)"]), effect_no_pancreatitis=subgroup_effect0, effect_with_pancreatitis=subgroup_effect1, difference_in_effects=subgroup_effect1-subgroup_effect0, n=nrow(x), pancreatitis_n=sum(x$pancreatitis), pancreatitis_weighted_n=sum(x$weight[x$pancreatitis==1]), interaction_status=if(is.null(term_row)) "not_estimable_in_survey_model" else "estimated")
  diag[[cc]] <- data.table(cohort=cc, overlap_ps_min=min(x$ps), overlap_ps_p01=quantile(x$ps,.01), overlap_ps_p99=quantile(x$ps,.99), max_weight=max(x$ow), n=nrow(x))
}
fwrite(rbindlist(rows, fill=TRUE), file.path(outdir,"pancreatitis_subgroup_rates.csv"))
fwrite(rbindlist(interactions, fill=TRUE), file.path(outdir,"pancreatitis_interaction.csv"))
fwrite(rbindlist(diag, fill=TRUE), file.path(outdir,"pancreatitis_diagnostics.csv"))
