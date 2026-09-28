from pathlib import Path
import numpy as np
import pandas as pd
import duckdb
import math

ROOT = Path(__file__).resolve().parent.parent.parent
PARQUET = ROOT / "work" / "stage12_final_revision" / "results" / "nrd_day3_landmark_v3.0.parquet"
OUT = ROOT / "work" / "stage26_high_ceiling_reanalysis" / "submission_package" / "input_cache" / "stage26_pancreatitis_subgroup"
OUT.mkdir(parents=True, exist_ok=True)

con = duckdb.connect()
cols = [
    "cohort", "broad_day3_eligible", "early_completion_day3", "biliary_readmission_d3_93",
    "dx_any_biliary_acute_pancreatitis", "AGE_num", "FEMALE_num", "AWEEKEND_num", "HCUP_ED_num",
    "PAY1", "ZIPINC_QRTL", "PL_NCHS", "dx_any_cholangitis", "dx_any_biliary_obstruction",
    "dx_any_cirrhosis_proxy", "dx_any_coagulation_disorder_proxy", "APRDRG_Severity",
    "APRDRG_Risk_Mortality", "HOSP_BEDSIZE", "HOSP_UR_TEACH", "H_CONTRL", "known_chole_state",
    "prior_90d_disease_volume", "DISCWT_num", "YEAR", "hospital_cluster_id"
]
q = "select " + ",".join(cols) + " from read_parquet(?) where broad_day3_eligible and cohort in ('acute_cholecystitis','choledocholithiasis')"
d = con.execute(q, [str(PARQUET)]).fetchdf()

def weighted_mean(y, w):
    return float(np.sum(y*w)/np.sum(w)) if len(y) and np.sum(w)>0 else np.nan

def dummy_matrix(x, cohort):
    # Match the covariate intent of the subgroup analysis with stable, admission-recorded covariates.
    z = pd.DataFrame(index=x.index)
    age = pd.to_numeric(x["AGE_num"], errors="coerce").fillna(pd.to_numeric(x["AGE_num"], errors="coerce").median())
    z["age"] = age
    z["age2"] = age**2
    for c in ["FEMALE_num","AWEEKEND_num","HCUP_ED_num","dx_any_cholangitis","dx_any_biliary_obstruction","dx_any_cirrhosis_proxy","dx_any_coagulation_disorder_proxy"]:
        z[c] = x[c].fillna(False).astype(int)
    cats = ["PAY1","ZIPINC_QRTL","PL_NCHS","APRDRG_Severity","APRDRG_Risk_Mortality","HOSP_BEDSIZE","HOSP_UR_TEACH","H_CONTRL","known_chole_state","YEAR"]
    if cohort == "choledocholithiasis":
        cats = cats + []
    z = pd.concat([z, pd.get_dummies(x[cats].fillna("unknown").astype(str), prefix=cats, drop_first=True, dtype=float)], axis=1)
    prior = np.log1p(pd.to_numeric(x["prior_90d_disease_volume"], errors="coerce").fillna(0).clip(lower=0))
    z["log_prior_volume"] = prior
    return z.astype(float)

def fit_ps(X, y, w):
    X = np.column_stack([np.ones(len(X)), X.to_numpy(dtype=float)])
    y = y.astype(float); w = w.astype(float)
    beta = np.zeros(X.shape[1])
    for _ in range(80):
        eta = np.clip(X @ beta, -30, 30)
        p = 1/(1+np.exp(-eta))
        v = np.maximum(p*(1-p), 1e-7)
        H = (X.T * (w*v)) @ X + np.eye(X.shape[1])*1e-8
        g = X.T @ (w*(y-p))
        step = np.linalg.solve(H, g)
        beta += step
        if np.max(np.abs(step)) < 1e-7: break
    ps = 1/(1+np.exp(-np.clip(X@beta,-30,30)))
    return np.clip(ps, .001, .999)

rates=[]; interactions=[]; diagnostics=[]
for cohort, x in d.groupby("cohort", sort=False):
    x=x.copy().reset_index(drop=True)
    t=x["early_completion_day3"].fillna(False).astype(int).to_numpy()
    y=x["biliary_readmission_d3_93"].fillna(False).astype(int).to_numpy()
    pan=x["dx_any_biliary_acute_pancreatitis"].fillna(False).astype(int).to_numpy()
    w=pd.to_numeric(x["DISCWT_num"], errors="coerce").fillna(1).clip(lower=1e-6).to_numpy()
    ps=fit_ps(dummy_matrix(x, cohort),t,w)
    ow=w*np.where(t==1,1-ps,ps)
    for pp in [0,1]:
        for tt in [0,1]:
            m=(pan==pp)&(t==tt)
            rates.append(dict(cohort=cohort,pancreatitis=pp,treatment=tt,n=int(m.sum()),weighted_n=float(w[m].sum()),completion_rate=weighted_mean(t[m],w[m]),biliary90_rate=weighted_mean(y[m],w[m]),overlap_biliary90_rate=weighted_mean(y[m],ow[m])))
        m=(pan==pp)
        m1=m&(t==1); m0=m&(t==0)
        rates.append(dict(cohort=cohort,pancreatitis=pp,treatment=9,n=int(m.sum()),weighted_n=float(w[m].sum()),completion_rate=weighted_mean(t[m],w[m]),biliary90_rate=np.nan,overlap_biliary90_rate=weighted_mean(y[m1],ow[m1])-weighted_mean(y[m0],ow[m0])))
    # weighted logistic interaction: biliary outcome ~ completion * pancreatitis
    Z=np.column_stack([np.ones(len(x)),t,pan,t*pan])
    b=np.zeros(4)
    for _ in range(80):
        p=1/(1+np.exp(-np.clip(Z@b,-30,30))); v=np.maximum(p*(1-p),1e-7)
        H=(Z.T*(ow*v))@Z+np.eye(4)*1e-8; g=Z.T@(ow*(y-p)); step=np.linalg.solve(H,g); b+=step
        if np.max(np.abs(step))<1e-8: break
    bread=np.linalg.inv((Z.T*(ow*np.maximum(p*(1-p),1e-7)))@Z+np.eye(4)*1e-8)
    score=Z*(ow*(y-p))[:,None]
    meat=np.zeros((4,4))
    for _, idx in x.groupby("hospital_cluster_id", dropna=False).groups.items():
        sg=score[np.asarray(list(idx),dtype=int)].sum(axis=0)
        meat += np.outer(sg,sg)
    cov=bread@meat@bread
    se=float(np.sqrt(max(cov[3,3],0))); stat=float(b[3]/se) if se>0 else np.nan
    pval=float(2*(1-0.5*(1+math.erf(abs(stat)/np.sqrt(2))))) if np.isfinite(stat) else np.nan
    e0=weighted_mean(y[(pan==0)&(t==1)],ow[(pan==0)&(t==1)])-weighted_mean(y[(pan==0)&(t==0)],ow[(pan==0)&(t==0)])
    e1=weighted_mean(y[(pan==1)&(t==1)],ow[(pan==1)&(t==1)])-weighted_mean(y[(pan==1)&(t==0)],ow[(pan==1)&(t==0)])
    interactions.append(dict(cohort=cohort,interaction_term="treatment:pancreatitis",interaction_log_odds=float(b[3]),interaction_se=se,interaction_p=pval,effect_no_pancreatitis=e0,effect_with_pancreatitis=e1,difference_in_effects=e1-e0,n=len(x),pancreatitis_n=int(pan.sum()),pancreatitis_weighted_n=float(w[pan==1].sum()),interaction_status="estimated_local_numpy_weighted_logistic"))
    diagnostics.append(dict(cohort=cohort,overlap_ps_min=float(ps.min()),overlap_ps_p01=float(np.quantile(ps,.01)),overlap_ps_p99=float(np.quantile(ps,.99)),max_weight=float(ow.max()),n=len(x)))

pd.DataFrame(rates).to_csv(OUT/"pancreatitis_subgroup_rates.csv",index=False)
pd.DataFrame(interactions).to_csv(OUT/"pancreatitis_interaction.csv",index=False)
pd.DataFrame(diagnostics).to_csv(OUT/"pancreatitis_diagnostics.csv",index=False)
print(pd.DataFrame(interactions).to_string(index=False))
print(pd.DataFrame(rates).to_string(index=False))
