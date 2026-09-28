#!/usr/bin/env Rscript

# Landmark-specific weighted cumulative-incidence analysis.
# Estimates observed associations, not causal treatment effects.

required_env <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
  "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
  "R_DATATABLE_NUM_THREADS", "ARROW_NUM_THREADS", "ARROW_CPU_COUNT",
  "RCPP_PARALLEL_NUM_THREADS")
bad_env <- required_env[!Sys.getenv(required_env, "") %in% c("1", "1L")]
if (length(bad_env)) stop("thread environment must be fixed before R starts: ", paste(bad_env, collapse=", "))
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
  BLIS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1", NUMEXPR_NUM_THREADS="1",
  R_DATATABLE_NUM_THREADS="1", ARROW_NUM_THREADS="1", ARROW_CPU_COUNT="1",
  RCPP_PARALLEL_NUM_THREADS="1")
suppressPackageStartupMessages({library(arrow); library(data.table); library(splines)})
setDTthreads(1L)
if ("set_cpu_count" %in% getNamespaceExports("arrow")) arrow::set_cpu_count(1L)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4L) stop("usage: stage20_nrd_cif_reanalysis.R INPUT LINKED_GLOB OUTDIR NBOOT")
input <- args[[1L]]
linked_glob <- args[[2L]]
outdir <- args[[3L]]
nboot <- as.integer(args[[4L]])
test_mode <- identical(Sys.getenv("STAGE20_TEST_MODE", "0"), "1")
if (!test_mode && !identical(nboot, 2000L)) stop("production mode requires exactly 2000 bootstrap replicates")
if (test_mode && (nboot < 2L || nboot > 50L)) stop("test-mode NBOOT must be 2--50")
if (!file.exists(input)) stop("missing input parquet: ", input)
linked_files <- Sys.glob(linked_glob)
if (length(linked_files) != 5L) stop("expected five linked-admission parquet files, found ", length(linked_files))
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
workers <- as.integer(Sys.getenv("STAGE20_BOOTSTRAP_WORKERS", "8"))
if (!is.finite(workers) || workers < 1L || workers > 8L) stop("bootstrap worker count must be 1--8")
workers <- if (.Platform$OS.type == "unix") min(workers, parallel::detectCores()) else 1L
set.seed(20260920L)

bin <- function(x) as.integer(tolower(as.character(x)) %in% c("1","true","t","yes"))
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
fref <- function(x) {
  z <- as.character(x); z[is.na(z) | z == "" | z %in% c("-99","-9","-8","-6")] <- "unknown"
  tt <- sort(table(z), decreasing=TRUE); factor(z, levels=c(names(tt)[1L], sort(setdiff(names(tt), names(tt)[1L]))))
}
require_columns <- function(x, cols) {
  miss <- setdiff(cols, names(x)); if (length(miss)) stop("missing required input columns: ", paste(miss, collapse=", "))
}
effn <- function(w) { w <- w[is.finite(w) & w > 0]; if (!length(w)) NA_real_ else sum(w)^2/sum(w^2) }
wvar <- function(x,w) {
  ok <- is.finite(x)&is.finite(w)&w>0; if (!any(ok)) return(NA_real_)
  x<-x[ok]; w<-w[ok]; m<-sum(w*x)/sum(w); sum(w*(x-m)^2)/sum(w)
}

d <- as.data.table(read_parquet(input))
need <- c("cohort","source_row","index_id","NRD_VisitLink","YEAR","NRD_DaysToEvent_num","index_discharge_day",
  "target_primary","urgent_admission","observed_prior_cholecystectomy","proxy_cholecystectomy_to_index",
  "definitive_completion_primary","treatment_day","treatment_timing_unknown","index_died","biliary_day_from_admission",
  "death_day_from_admission","LOS_num","DISCWT_num","NRD_STRATUM","hospital_cluster_id","AGE_num","FEMALE_num",
  "AWEEKEND_num","ELECTIVE_num","HCUP_ED_num","PAY1","ZIPINC_QRTL","PL_NCHS","HOSP_BEDSIZE",
  "HOSP_UR_TEACH","H_CONTRL","HOSP_URCAT4","DMONTH","prior_90d_disease_volume",
  "pr_any_cholecystectomy_total","pr_any_cholecystectomy_partial","pr_any_ercp_clearance_proxy",
  "pr_any_surgical_cbd_clearance","pr_any_ercp_therapeutic","pr_any_ercp_diagnostic",
  "pr_any_surgical_bile_duct_exploration","prday_min_cholecystectomy_total","prday_min_cholecystectomy_partial",
  "prday_min_ercp_clearance_proxy","prday_min_surgical_cbd_clearance","prday_min_ercp_therapeutic",
  "prday_unknown_cholecystectomy_total","prday_unknown_cholecystectomy_partial",
  "prday_unknown_ercp_clearance_proxy","prday_unknown_surgical_cbd_clearance","YEAR","DMONTH")
require_columns(d, unique(need))
d[, `:=`(
  row_id=.I, idx_source=as.character(source_row), index_id_key=as.character(index_id), visit=as.character(NRD_VisitLink),
  year=as.integer(num(YEAR)), idx_admit=num(NRD_DaysToEvent_num), idx_discharge=num(index_discharge_day),
  target_primary_b=bin(target_primary), urgent_b=bin(urgent_admission), prior_chole=bin(observed_prior_cholecystectomy),
  z9049_proxy=bin(proxy_cholecystectomy_to_index), definitive=bin(definitive_completion_primary),
  trt_day=num(treatment_day), timing_unknown=bin(treatment_timing_unknown), index_death=bin(index_died),
  source_biliary_day=num(biliary_day_from_admission), source_death_day=num(death_day_from_admission),
  los=num(LOS_num), base_weight=num(DISCWT_num),
  cluster=as.character(hospital_cluster_id), boot_stratum=interaction(YEAR,NRD_STRATUM,drop=TRUE),
  age=num(AGE_num), female=bin(FEMALE_num), weekend=bin(AWEEKEND_num), elective=bin(ELECTIVE_num), ed=bin(HCUP_ED_num),
  payer=fref(PAY1), income=fref(ZIPINC_QRTL), rurality=fref(PL_NCHS), bedsize=fref(HOSP_BEDSIZE),
  teaching=fref(HOSP_UR_TEACH), ownership=fref(H_CONTRL), urbanicity=fref(HOSP_URCAT4),
  year_f=fref(YEAR), month_f=fref(DMONTH), prior_volume=num(prior_90d_disease_volume)
)]
d[!is.finite(age), age := median(age,na.rm=TRUE)]
d[!is.finite(prior_volume), prior_volume := 0]

# Reconstruct biliary, nonbiliary and death events from all linked admissions.
adm <- rbindlist(lapply(linked_files, function(p) as.data.table(read_parquet(p))), fill=TRUE)
require_columns(adm, c("source_row","NRD_VisitLink","NRD_DaysToEvent","ELECTIVE","SAMEDAYEVENT","REHABTRANSFER","DIED","YEAR","related_dx_any"))
adm[, `:=`(visit=as.character(NRD_VisitLink), year=as.integer(num(YEAR)), adm_source=as.character(source_row),
  adm_day=num(NRD_DaysToEvent), adm_elective=bin(ELECTIVE), adm_sameday=num(SAMEDAYEVENT)>0,
  adm_rehab=bin(REHABTRANSFER), adm_death=bin(DIED), adm_biliary=bin(related_dx_any))]
idx_key <- unique(d[,.(row_id,index_id_key,cohort,visit,year,idx_source,idx_admit,idx_discharge)])
candidate <- merge(idx_key, adm[,.(visit,year,adm_source,adm_day,adm_elective,adm_sameday,adm_rehab,adm_death,adm_biliary)],
  by=c("visit","year"), allow.cartesian=TRUE)
candidate <- candidate[adm_source != idx_source & is.finite(adm_day) & adm_day > idx_discharge]
candidate[, day_from_admission := adm_day-idx_admit]
full_followup <- candidate[,.(
  biliary_day_full=suppressWarnings(min(day_from_admission[adm_biliary==1L & !adm_sameday & adm_rehab==0L],na.rm=TRUE)),
  followup_death_day_full=suppressWarnings(min(day_from_admission[adm_death==1L],na.rm=TRUE))
),by=.(row_id,index_id_key,cohort)]
full_followup[!is.finite(biliary_day_full),biliary_day_full:=NA_real_]
full_followup[!is.finite(followup_death_day_full),followup_death_day_full:=NA_real_]
d[full_followup,`:=`(biliary_day=i.biliary_day_full,followup_death_day=i.followup_death_day_full),on=.(row_id,index_id_key,cohort)]
d[,index_death_day:=fifelse(index_death==1L,los,NA_real_)]
d[,death_day:=pmin(index_death_day,followup_death_day,na.rm=TRUE)]
d[!is.finite(death_day),death_day:=NA_real_]
event_reconstruction_check <- d[,.(
  n=.N,
  source_biliary_events_d3_93=sum(is.finite(source_biliary_day)&source_biliary_day>3&source_biliary_day<=93),
  recomputed_biliary_events_d3_93=sum(is.finite(biliary_day)&biliary_day>3&biliary_day<=93),
  biliary_day_mismatch=sum((is.finite(source_biliary_day)!=is.finite(biliary_day)) |
    (is.finite(source_biliary_day)&is.finite(biliary_day)&source_biliary_day!=biliary_day)),
  source_deaths_d3_93=sum(is.finite(source_death_day)&source_death_day>3&source_death_day<=93),
  recomputed_deaths_d3_93=sum(is.finite(death_day)&death_day>3&death_day<=93),
  death_day_mismatch=sum((is.finite(source_death_day)!=is.finite(death_day)) |
    (is.finite(source_death_day)&is.finite(death_day)&source_death_day!=death_day))
),by=cohort]
nonb_candidate <- candidate[adm_biliary==0L & adm_elective==0L & !adm_sameday & adm_rehab==0L]
nonb_by_landmark <- rbindlist(lapply(2:4, function(L) nonb_candidate[day_from_admission>L & day_from_admission<=L+90,
  .(nonbiliary_after_landmark=min(day_from_admission)), by=.(row_id,index_id_key,cohort)][,landmark:=L]))
if (nrow(nonb_by_landmark)) setkey(nonb_by_landmark,row_id,landmark)
rm(adm,candidate,nonb_candidate,full_followup); invisible(gc())
fwrite(event_reconstruction_check,file.path(outdir,"stage20_nrd_event_reconstruction_check.csv"))

ps_formula <- function(cohort) {
  terms <- c("ns(age,3)","female","weekend","elective","ed","payer","income","rurality","bedsize","teaching",
    "ownership","urbanicity","year_f","month_f","log1p(prior_volume)")
  if (cohort=="choledocholithiasis") terms <- c(terms,"prior_chole")
  as.formula(paste("treatment ~",paste(terms,collapse=" + ")))
}
design_matrix <- function(form,x) {
  mm<-model.matrix(delete.response(terms(form)),data=as.data.frame(x),na.action=na.pass)
  storage.mode(mm)<-"double"; mm[!is.finite(mm)]<-0; mm
}
fit_ps <- function(mm,treatment,w) {
  ok<-is.finite(w)&w>0&is.finite(treatment)
  if (sum(ok)<100L || length(unique(treatment[ok]))!=2L) return(list(ok=FALSE,reason="insufficient exposure support"))
  active<-rep(TRUE,ncol(mm)); if(ncol(mm)>1L) active[-1L]<-vapply(2:ncol(mm),function(j) diff(range(mm[ok,j]))>1e-12,logical(1))
  idx<-which(active); q<-qr(sqrt(w[ok])*mm[ok,idx,drop=FALSE],tol=1e-10)
  if(q$rank<1L) return(list(ok=FALSE,reason="rank-zero PS design"))
  keep<-idx[sort(q$pivot[seq_len(q$rank)])]
  fit<-suppressWarnings(glm.fit(mm[ok,keep,drop=FALSE],treatment[ok],family=binomial(),weights=w[ok],
    control=glm.control(epsilon=1e-8,maxit=100L)))
  if(!isTRUE(fit$converged)||any(!is.finite(fit$coefficients))) return(list(ok=FALSE,reason="nonconvergent/nonfinite PS"))
  eta<-drop(mm[,keep,drop=FALSE]%*%fit$coefficients)
  if(any(!is.finite(eta))) return(list(ok=FALSE,reason="nonfinite PS predictor"))
  list(ok=TRUE,ps=pmin(pmax(plogis(eta),.001),.999),rank=q$rank,columns=paste(colnames(mm)[keep],collapse="|"))
}
smd_table <- function(x,mm,w0,w1) {
  rbindlist(lapply(seq_len(ncol(mm)),function(j){
    z<-mm[,j]; a<-x$treatment==1L
    calc<-function(w){m1<-weighted.mean(z[a],w[a]);m0<-weighted.mean(z[!a],w[!a]);den<-sqrt((wvar(z[a],w[a])+wvar(z[!a],w[!a]))/2);ifelse(is.finite(den)&&den>0,(m1-m0)/den,0)}
    data.table(variable=colnames(mm)[j],smd_before=calc(w0),smd_after=calc(w1))
  }))
}

make_landmark <- function(cohort_name,L,timing_mode="primary") {
  x<-copy(d[cohort==cohort_name])
  base_ok<-x$target_primary_b==1L & x$urgent_b==1L & (cohort_name!="acute_cholecystitis" | x$prior_chole==0L) &
    (!is.finite(x$death_day)|x$death_day>L) & (!is.finite(x$biliary_day)|x$biliary_day>L) &
    is.finite(x$base_weight)&x$base_weight>0 & !is.na(x$cluster)&x$cluster!="" & is.finite(x$los)&x$los>=0
  if(timing_mode=="primary") keep<-base_ok & x$timing_unknown==0L else keep<-base_ok
  x<-x[keep]
  if(timing_mode=="primary") x[,treatment:=as.integer(definitive==1L & is.finite(trt_day) & trt_day<=L)]
  if(timing_mode=="early_bound") x[,treatment:=as.integer((definitive==1L & is.finite(trt_day)&trt_day<=L) | (timing_unknown==1L & definitive==1L))]
  if(timing_mode=="late_bound") x[,treatment:=as.integer(definitive==1L & is.finite(trt_day)&trt_day<=L & timing_unknown==0L)]
  if(length(unique(x$treatment))!=2L) stop("exposure support failed: ",cohort_name," L",L," ",timing_mode)
  x[,landmark:=as.integer(L)]
  if(nrow(nonb_by_landmark)) x[nonb_by_landmark[landmark==L], nonbiliary_day := i.nonbiliary_after_landmark, on=.(row_id,landmark)] else x[,nonbiliary_day:=NA_real_]
  x
}

prepare <- function(x,cohort_name,L,timing_mode="primary") {
  form<-ps_formula(cohort_name);mm<-design_matrix(form,x);fit<-fit_ps(mm,x$treatment,x$base_weight)
  if(!fit$ok) stop("point PS failed: ",cohort_name," L",L," ",fit$reason)
  x[,ps:=fit$ps];x[,ow:=base_weight*ifelse(treatment==1L,1-ps,ps)]
  if(any(!is.finite(x$ow)|x$ow<=0)) stop("invalid overlap weights")
  bal<-smd_table(x,mm,x$base_weight,x$ow)
  list(x=x,mm=mm,form=form,rank=fit$rank,columns=fit$columns,balance=bal,cohort=cohort_name,landmark=L,timing_mode=timing_mode)
}

event_vectors <- function(x,event_col,L,horizon=90L) {
  ev<-x[[event_col]]; de<-x$death_day
  valid_ev<-is.finite(ev)&ev>L&ev<=L+horizon
  valid_de<-is.finite(de)&de>L&de<=L+horizon
  death_first<-valid_de & (!valid_ev | de<=ev)
  event_first<-valid_ev & !death_first
  status<-integer(nrow(x));status[event_first]<-1L;status[death_first]<-2L
  t<-rep(as.integer(horizon),nrow(x));t[event_first]<-as.integer(ev[event_first]-L);t[death_first]<-as.integer(de[death_first]-L)
  list(time=t,status=status)
}
weighted_cif <- function(time,status,w,horizon=90L) {
  if(any(!is.finite(w)|w<=0)) stop("CIF requires positive finite weights")
  S<-1;F1<-0;F2<-0;out<-vector("list",horizon+1L)
  out[[1L]]<-data.table(day=0L,survival=S,cif_event=F1,cif_death=F2,risk_weight=sum(w),event_weight=0,death_weight=0)
  for(t in seq_len(horizon)){
    Y<-sum(w[time>=t]);d1<-sum(w[time==t&status==1L]);d2<-sum(w[time==t&status==2L])
    if(Y>0){F1<-F1+S*d1/Y;F2<-F2+S*d2/Y;S<-S*(1-(d1+d2)/Y)}
    out[[t+1L]]<-data.table(day=t,survival=S,cif_event=F1,cif_death=F2,risk_weight=Y,event_weight=d1,death_weight=d2)
  }
  z<-rbindlist(out)
  if(any(!is.finite(z$cif_event)|z$cif_event< -1e-10|z$cif_event+z$cif_death>1+1e-8)) stop("CIF invariant failed")
  z
}
estimate_prepared <- function(p,weight_col="ow") {
  x<-p$x;all<-list();ep<-list()
  for(outcome in c("biliary_day","nonbiliary_day")) for(g in 0:1){
    z<-x[treatment==g];vec<-event_vectors(z,outcome,p$landmark);cur<-weighted_cif(vec$time,vec$status,z[[weight_col]])
    cur[,`:=`(cohort=p$cohort,landmark=p$landmark,outcome=ifelse(outcome=="biliary_day","biliary_readmission","nonbiliary_nonelective_readmission"),
      treatment=ifelse(g==1L,"early_completion","not_early_completion"),weight_type=weight_col)]
    all[[length(all)+1L]]<-cur
    pts<-if(outcome=="biliary_day") c(30L,60L,90L) else 90L
    ep[[length(ep)+1L]]<-cur[day%in%pts,.(cohort,landmark,outcome,treatment,day,cif_event,cif_death,survival)]
  }
  curve<-rbindlist(all);end<-rbindlist(ep)
  wide<-dcast(end,cohort+landmark+outcome+day~treatment,value.var=c("cif_event","cif_death","survival"))
  wide[,difference_early_minus_not_early:=cif_event_early_completion-cif_event_not_early_completion]
  list(curve=curve,endpoint=wide)
}

# Point estimates for independent day-2, day-3 and day-4 risk sets.
prepared<-list();curves<-list();endpoints<-list();balances<-list();diagnostics<-list();flows<-list();admin<-list()
for(cc in c("acute_cholecystitis","choledocholithiasis")) for(L in 2:4){
  x<-make_landmark(cc,L,"primary");p<-prepare(x,cc,L);est<-estimate_prepared(p)
  prepared[[paste(cc,L)]]<-p;curves[[paste(cc,L)]]<-est$curve;endpoints[[paste(cc,L)]]<-est$endpoint
  balances[[paste(cc,L)]]<-copy(p$balance)[,`:=`(cohort=cc,landmark=L)]
  diagnostics[[paste(cc,L)]]<-data.table(cohort=cc,landmark=L,n=nrow(x),n_early=sum(x$treatment==1L),n_not_early=sum(x$treatment==0L),
    design_weight=sum(x$base_weight),overlap_weight=sum(x$ow),ps_min=min(x$ps),ps_p01=quantile(x$ps,.01),ps_p99=quantile(x$ps,.99),ps_max=max(x$ps),
    ess_early=effn(x$ow[x$treatment==1L]),ess_not_early=effn(x$ow[x$treatment==0L]),max_abs_smd_before=max(abs(p$balance$smd_before)),
    max_abs_smd_after=max(abs(p$balance$smd_after)),ps_rank=p$rank,ps_columns=p$columns)
  allcc<-d[cohort==cc];basepre<-allcc[target_primary_b==1L&urgent_b==1L&(cc!="acute_cholecystitis"|prior_chole==0L)&
    (!is.finite(death_day)|death_day>L)&(!is.finite(biliary_day)|biliary_day>L)&is.finite(base_weight)&base_weight>0&!is.na(cluster)&cluster!=""&is.finite(los)&los>=0]
  flows[[paste(cc,L)]]<-data.table(cohort=cc,landmark=L,pre_timing_n=nrow(basepre),timing_unknown_excluded=sum(basepre$timing_unknown==1L),
    primary_n=nrow(x),early_n=sum(x$treatment==1L),not_early_n=sum(x$treatment==0L))
  # NRD_DaysToEvent is a synthetic person-level linkage clock, not a calendar
  # day-of-year. Full within-year follow-up is guaranteed at month resolution
  # by the discharge-month restriction (April--September): even a
  # September 30 discharge plus 90 days remains in the same calendar year.
  discharge_month<-as.integer(num(x$DMONTH));full<-is.finite(discharge_month)&discharge_month>=4L&discharge_month<=9L
  admin[[paste(cc,L)]]<-data.table(cohort=cc,landmark=L,n=nrow(x),full_window_n=sum(full),incomplete_window_n=sum(!full),
    full_window_weight=sum(x$base_weight[full]),incomplete_window_weight=sum(x$base_weight[!full]))
  cat(sprintf("POINT_DONE cohort=%s landmark=%d n=%d\n",cc,L,nrow(x)));flush.console()
}
curves_dt<-rbindlist(curves);endpoints_dt<-rbindlist(endpoints);balance_dt<-rbindlist(balances);diag_dt<-rbindlist(diagnostics)
fwrite(curves_dt,file.path(outdir,"stage20_nrd_landmark_cif_curves.csv"))
fwrite(endpoints_dt,file.path(outdir,"stage20_nrd_landmark_cif_endpoints.csv"))
fwrite(balance_dt,file.path(outdir,"stage20_nrd_propensity_balance.csv"))
fwrite(diag_dt,file.path(outdir,"stage20_nrd_propensity_diagnostics.csv"))
fwrite(rbindlist(flows),file.path(outdir,"stage20_nrd_patient_flow.csv"))
fwrite(rbindlist(admin),file.path(outdir,"stage20_nrd_administrative_window_check.csv"))

# Procedure-day missingness and classification bounds at day 3.
d[, target_proc := ifelse(cohort=="acute_cholecystitis",bin(pr_any_cholecystectomy_total)|bin(pr_any_cholecystectomy_partial),
  bin(pr_any_ercp_clearance_proxy)|bin(pr_any_surgical_cbd_clearance))]
d[, missing_pattern := fifelse(target_proc==0L,"no_target_procedure",fifelse(is.finite(trt_day)&timing_unknown==0L,"known_min_no_unknown_slot",
  fifelse(is.finite(trt_day)&timing_unknown==1L,"known_min_plus_unknown_slot","all_relevant_days_unknown")))]
missingness_check<-d[,.(n=.N,design_weight=sum(base_weight[is.finite(base_weight)&base_weight>0],na.rm=TRUE)),by=.(cohort,year,teaching,missing_pattern)]
fwrite(missingness_check,file.path(outdir,"stage20_nrd_procedure_day_missingness.csv"))
bound_rows<-list();bound_diag<-list()
for(cc in c("acute_cholecystitis","choledocholithiasis")) for(mode in c("early_bound","late_bound")){
  x<-make_landmark(cc,3L,mode);p<-prepare(x,cc,3L,mode);e<-estimate_prepared(p)$endpoint[outcome=="biliary_readmission"&day==90L]
  e[,timing_mode:=mode];bound_rows[[paste(cc,mode)]]<-e
  bound_diag[[paste(cc,mode)]]<-data.table(cohort=cc,timing_mode=mode,n=nrow(x),early_n=sum(x$treatment==1L),not_early_n=sum(x$treatment==0L),
    max_abs_smd_after=max(abs(p$balance$smd_after)),ess_early=effn(x$ow[x$treatment==1L]),ess_not_early=effn(x$ow[x$treatment==0L]))
}
fwrite(rbindlist(bound_rows),file.path(outdir,"stage20_nrd_timing_unknown_bounds.csv"))
fwrite(rbindlist(bound_diag),file.path(outdir,"stage20_nrd_timing_unknown_bounds_diagnostics.csv"))

# Choledocholithiasis operation phenotype before timing-based exclusion. This
# table is descriptive and intentionally retains unknown PRDAY records.
cbd<-make_landmark("choledocholithiasis",3L,"early_bound")
cbd[,exposure_status:=fcase(timing_unknown==1L,"unknown_timing",
  definitive==1L&is.finite(trt_day)&trt_day<=3L,"recorded_by_day3",default="not_recorded_by_day3")]
cbd[,primary_component:=fcase(bin(pr_any_ercp_clearance_proxy)==1L&bin(pr_any_surgical_cbd_clearance)==1L,"both_recorded",
  bin(pr_any_ercp_clearance_proxy)==1L,"ercp_clearance_proxy_only",bin(pr_any_surgical_cbd_clearance)==1L,"surgical_cbd_clearance_only",default="neither_by_record")]
cbd[,gallbladder_phenotype:=fcase(prior_chole==1L,"strictly_observed_prior_cholecystectomy",z9049_proxy==1L,"index_recorded_Z9049_proxy_only",default="neither_observed")]
cbd[,procedure_day_band:=fcase(!is.finite(trt_day),"unknown",trt_day==0,"0",trt_day==1,"1",trt_day==2,"2",trt_day==3,"3",trt_day==4,"4",trt_day>4,">4",default="other")]
cbd[,`:=`(therapeutic_ercp=as.character(bin(pr_any_ercp_therapeutic)),diagnostic_ercp=as.character(bin(pr_any_ercp_diagnostic)),
  surgical_exploration=as.character(bin(pr_any_surgical_bile_duct_exploration)))]
phen_summary<-function(x,column,dimension){
  ans<-x[,{
    v<-event_vectors(.SD,"biliary_day",3L);z<-weighted_cif(v$time,v$status,base_weight)
    list(n=.N,design_weight=sum(base_weight),biliary_cif90=z[day==90L]$cif_event,observed_inpatient_death_cif90=z[day==90L]$cif_death)
  },by=c("exposure_status",column)]
  setnames(ans,column,"category");ans[,dimension:=dimension];setcolorder(ans,c("dimension","category","exposure_status","n","design_weight","biliary_cif90","observed_inpatient_death_cif90"));ans
}
phen<-rbindlist(list(
  phen_summary(cbd,"primary_component","primary recorded component"),
  phen_summary(cbd,"gallbladder_phenotype","gallbladder-status evidence"),
  phen_summary(cbd,"procedure_day_band","recorded procedure day"),
  phen_summary(cbd,"therapeutic_ercp","therapeutic ERCP descriptor"),
  phen_summary(cbd,"diagnostic_ercp","diagnostic ERCP descriptor"),
  phen_summary(cbd,"surgical_exploration","surgical exploration descriptor")
),fill=TRUE)
fwrite(phen,file.path(outdir,"stage20_choledocholithiasis_procedure_phenotype.csv"))

# Hospital-year bootstrap with PS refit for the day-3 primary and broad-cohort negative control.
sample_mult <- function(x){
  psu<-unique(x[,.(boot_stratum,cluster)]);draw<-psu[,.(draw=sample(cluster,.N,replace=TRUE)),by=boot_stratum]
  cnt<-draw[,.(multiplicity=.N),by=.(boot_stratum,draw)];ans<-cnt$multiplicity[match(paste(x$boot_stratum,x$cluster),paste(cnt$boot_stratum,cnt$draw))]
  ans[is.na(ans)]<-0;as.numeric(ans)
}
boot_one <- function(p,seed){
  set.seed(seed);x<-p$x;mult<-sample_mult(x);bw<-x$base_weight*mult;fit<-fit_ps(p$mm,x$treatment,bw)
  if(!fit$ok)return(data.table(finite=FALSE,reason=fit$reason))
  ow<-bw*ifelse(x$treatment==1L,1-fit$ps,fit$ps)
  if(any(!is.finite(ow))||sum(ow[x$treatment==1L])<=0||sum(ow[x$treatment==0L])<=0)return(data.table(finite=FALSE,reason="invalid bootstrap weights"))
  rows<-list()
  for(outcome in c("biliary_day","nonbiliary_day"))for(g in 0:1){
    keep<-x$treatment==g & ow>0
    z<-x[keep];v<-event_vectors(z,outcome,3L);cif<-weighted_cif(v$time,v$status,ow[keep])
    pts<-if(outcome=="biliary_day")c(30L,60L,90L)else 90L
    rows[[length(rows)+1L]]<-cif[day%in%pts,.(outcome=ifelse(outcome=="biliary_day","biliary_readmission","nonbiliary_nonelective_readmission"),
      treatment=g,day,cif_event,cif_death)]
  }
  z<-rbindlist(rows);wide<-dcast(z,outcome+day~treatment,value.var=c("cif_event","cif_death"))
  setnames(wide,c("cif_event_0","cif_event_1","cif_death_0","cif_death_1"),c("risk_not_early","risk_early","death_not_early","death_early"),skip_absent=TRUE)
  wide[,`:=`(difference_early_minus_not_early=risk_early-risk_not_early,finite=TRUE,reason="",ps_rank=fit$rank)]
  wide
}
boot_cohort <- function(p,seedbase){
  one<-function(i){z<-tryCatch(boot_one(p,seedbase+i),error=function(e)data.table(finite=FALSE,reason=conditionMessage(e)))
    if(!"outcome"%in%names(z))z<-CJ(outcome=c("biliary_readmission","nonbiliary_nonelective_readmission"),day=c(30L,60L,90L))[outcome=="nonbiliary_nonelective_readmission"&day!=90L,day:=NA_integer_][!is.na(day)][,`:=`(risk_not_early=NA_real_,risk_early=NA_real_,death_not_early=NA_real_,death_early=NA_real_,difference_early_minus_not_early=NA_real_,finite=FALSE,reason=z$reason[[1]],ps_rank=NA_integer_)]
    z[,replicate:=i];z}
  ans<-if(.Platform$OS.type=="unix")parallel::mclapply(seq_len(nboot),one,mc.cores=workers,mc.preschedule=TRUE,mc.set.seed=FALSE)else lapply(seq_len(nboot),one)
  rbindlist(ans,fill=TRUE)
}
boot_all<-list()
for(cc in c("acute_cholecystitis","choledocholithiasis")){
  p<-prepared[[paste(cc,3)]];seedbase<-if(cc=="acute_cholecystitis")202609200L else 202619200L
  cat(sprintf("BOOT_START cohort=%s nboot=%d workers=%d\n",cc,nboot,workers));flush.console()
  z<-boot_cohort(p,seedbase);z[,cohort:=cc];boot_all[[cc]]<-z
  cat(sprintf("BOOT_DONE cohort=%s success=%d\n",cc,uniqueN(z[finite==TRUE,replicate])));flush.console()
}
boot<-rbindlist(boot_all,fill=TRUE);fwrite(boot,file.path(outdir,"stage20_nrd_hospital_year_bootstrap_replicates.csv"))
boot_ci<-boot[finite==TRUE,.(bootstrap_total=nboot,bootstrap_success=.N,bootstrap_failure=nboot-.N,bootstrap_success_rate=.N/nboot,
  risk_early_lcl=quantile(risk_early,.025),risk_early_ucl=quantile(risk_early,.975),risk_not_early_lcl=quantile(risk_not_early,.025),risk_not_early_ucl=quantile(risk_not_early,.975),
  difference_lcl=quantile(difference_early_minus_not_early,.025),difference_ucl=quantile(difference_early_minus_not_early,.975),
  difference_sd=sd(difference_early_minus_not_early),difference_unique=uniqueN(difference_early_minus_not_early),min_ps_rank=min(ps_rank),max_ps_rank=max(ps_rank)),by=.(cohort,outcome,day)]
point3<-endpoints_dt[landmark==3L];final<-merge(point3,boot_ci,by=c("cohort","outcome","day"),all.x=TRUE)
fwrite(final,file.path(outdir,"stage20_nrd_day3_primary_bootstrap_ci.csv"))

# Hard gates.
expected<-data.table(cohort=c("acute_cholecystitis","choledocholithiasis"),expected_n=c(233764L,117293L))
observed<-diag_dt[landmark==3L,.(cohort,observed_n=n,max_abs_smd_after,ess_early,ess_not_early)]
gate<-merge(expected,observed,by="cohort",all=TRUE)
bootsum<-boot_ci[,.(bootstrap_min_success=min(bootstrap_success),all_variation=all(is.finite(difference_sd)&difference_sd>0&difference_unique>1L),all_ci=all(is.finite(difference_lcl)&is.finite(difference_ucl)&difference_ucl>difference_lcl)),by=cohort]
gate<-merge(gate,bootsum,by="cohort",all=TRUE)
admin3<-rbindlist(admin)[landmark==3L,.(cohort,incomplete_window_n)]
gate<-merge(gate,admin3,by="cohort",all=TRUE)
gate[,pass:=observed_n==expected_n&max_abs_smd_after<=.10&is.finite(ess_early)&is.finite(ess_not_early)&bootstrap_min_success>=ceiling(.95*nboot)&all_variation&all_ci&incomplete_window_n==0L]
fwrite(gate,file.path(outdir,"stage20_nrd_analysis_gate.csv"))
if(nrow(gate)!=2L||!all(gate$pass))stop("analysis gate failed")
