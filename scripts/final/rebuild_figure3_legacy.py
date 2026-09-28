from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent
PKG = ROOT / 'submission_package'
CACHE = PKG / 'input_cache'
FIG = PKG / 'figures'
SRC = FIG / 'source_data'
FIG.mkdir(exist_ok=True); SRC.mkdir(exist_ok=True)
blue, orange, teal, grey = '#2C7FB8', '#D95F0E', '#238B8D', '#6B7280'
plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.titlesize':8,'axes.labelsize':7,'xtick.labelsize':6.2,'ytick.labelsize':6.2,'legend.fontsize':6,'axes.linewidth':0.6})

hy = pd.read_csv(CACHE/'stage23_brazil_nrd_extensions/results/nrd_extensions/stage23_nrd_hospital_year_model_diagnostics.csv')
rates = pd.read_csv(CACHE/'stage23_brazil_nrd_extensions/results/nrd_extensions/stage23_nrd_hospital_year_risk_standardized_rates.csv')
burden = pd.read_csv(CACHE/'stage23_brazil_nrd_extensions/results/nrd_extensions/stage23_nrd_peer_benchmark_associated_burden.csv')

fig, axs = plt.subplots(2,2,figsize=(7.2,4.72))
fig.subplots_adjust(left=.10,right=.98,bottom=.11,top=.93,wspace=.38,hspace=.42)

# Panel a: US hospital-year distributions by cohort. This avoids presenting
# incompatible US and Brazil completion measures as one directly comparable
# endpoint.
ax=axs[0,0]; z=rates[(rates.outcome=='treatment') & (rates.eligible_n20.astype(str).str.lower().isin(['true','1']))].copy()
for i, cohort in enumerate(['acute_cholecystitis','choledocholithiasis']):
    zz=z[z.cohort==cohort]
    vals=zz.risk_standardized_rate.dropna().to_numpy()
    col=blue if cohort=='acute_cholecystitis' else teal
    ax.violinplot(vals,positions=[i],showmeans=False,showmedians=False,showextrema=False,
                  widths=.55,quantiles=[[.10,.90]],
                  bw_method=.25,points=80)
    for body in ax.collections[-1:]: body.set_facecolor(col); body.set_edgecolor(col); body.set_alpha(.35)
    ax.scatter([i]*len(vals),vals,s=1.5,alpha=.12,color=col,edgecolors='none')
    ax.scatter(i,np.median(vals),s=22,color=col,edgecolors='white',linewidths=.6,zorder=4)
    ax.text(i,.98,f'n={len(vals):,}',transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=6)
ax.set_xticks([0,1]); ax.set_xticklabels(['Acute\ncholecystitis','Choledocho-\nlithiasis']); ax.set_ylim(.3,1.0); ax.set_ylabel('Risk-standardized completion'); ax.set_title('US hospital-year distributions',fontweight='bold'); ax.yaxis.set_major_formatter(lambda v,pos:f'{v:.0%}'); ax.grid(False)

# Panel b: hospital-year P10-P90 and median.
ax=axs[0,1]; z=hy[hy.outcome=='treatment'].copy(); x=np.arange(len(z)); cols=[blue if c=='acute_cholecystitis' else teal for c in z.cohort]
ax.vlines(x,z.rsr_p10,z.rsr_p90,color=cols,lw=4,alpha=.65); ax.scatter(x,z.rsr_p50,color=cols,s=18,zorder=3)
ax.set_xticks(x); ax.set_xticklabels(['Acute\ncholecystitis','Choledocho-\nlithiasis']); ax.set_ylim(.3,.9); ax.set_ylabel('Risk-standardized completion'); ax.set_title('Hospital-year variation',fontweight='bold'); ax.yaxis.set_major_formatter(lambda v,pos:f'{v:.0%}'); ax.grid(False)

# Panel b: paired hospital-year scatter with a descriptive trend line.
ax=axs[0,1]
rr=rates[(rates.cohort=='acute_cholecystitis') & (rates.eligible_n20.astype(str).str.lower().isin(['true','1']))].copy()
t=rr[rr.outcome=='treatment'][['year','hospital_id','risk_standardized_rate']].rename(columns={'risk_standardized_rate':'completion_rsr'})
r=rr[rr.outcome=='biliary90'][['year','hospital_id','risk_standardized_rate']].rename(columns={'risk_standardized_rate':'readmission_rsr'})
pair=t.merge(r,on=['year','hospital_id'],how='inner'); pair.to_csv(SRC/'figure3_hospital_year_pairs.csv',index=False)
ax.scatter(pair.completion_rsr,pair.readmission_rsr,s=5,alpha=.20,color=blue,edgecolors='none')
q=pair[['completion_rsr','readmission_rsr']].dropna(); slope,inter=np.polyfit(q.completion_rsr,q.readmission_rsr,1); xx=np.linspace(q.completion_rsr.min(),q.completion_rsr.max(),100); ax.plot(xx,slope*xx+inter,color=grey,lw=1.0)
rho=q.corr(method='spearman').iloc[0,1]
ax.text(.98,.97,f'Spearman rho = {rho:.2f}\nPaired hospital-years = {len(q):,}',transform=ax.transAxes,ha='right',va='top',fontsize=6.2,bbox=dict(facecolor='white',edgecolor='#BBBBBB',boxstyle='round,pad=.25',alpha=.9))
ax.set_xlabel('Risk-standardized completion'); ax.set_ylabel('Risk-standardized biliary readmission'); ax.set_title('Ecological hospital-year relation',fontweight='bold'); ax.grid(False)

# Panel c: associated biliary returns under prespecified gap-closure scenarios.
ax=axs[1,0]
ax.fill_between(burden.closure_fraction,burden.readmissions_lcl,burden.readmissions_ucl,color=blue,alpha=.16)
ax.plot(burden.closure_fraction,burden.associated_biliary_readmissions,color=blue,lw=1.0)
ax.scatter(burden.closure_fraction,burden.associated_biliary_readmissions,color=blue,s=16)
ax.set_xlabel('Peer-P75 completion-gap closure'); ax.set_ylabel('Associated biliary returns'); ax.set_title('Associated return burden',fontweight='bold'); ax.set_xticks(burden.closure_fraction); ax.set_xticklabels([f'{v:.0%}' for v in burden.closure_fraction]); ax.grid(False)

# Panel d: the same scenarios expressed as additional recorded completions.
ax=axs[1,1]
ax.fill_between(burden.closure_fraction,burden.newly_completed_lcl,burden.newly_completed_ucl,color=teal,alpha=.16)
ax.plot(burden.closure_fraction,burden.associated_newly_completed,color=teal,lw=1.0)
ax.scatter(burden.closure_fraction,burden.associated_newly_completed,color=teal,s=16)
ax.set_xlabel('Peer-P75 completion-gap closure'); ax.set_ylabel('Associated additional completions'); ax.set_title('Associated completion burden',fontweight='bold'); ax.set_xticks(burden.closure_fraction); ax.set_xticklabels([f'{v:.0%}' for v in burden.closure_fraction]); ax.grid(False)

for ax,tag in zip(axs.flat,['a','b','c','d']): ax.text(-.16,1.06,tag,transform=ax.transAxes,fontweight='bold',fontsize=9,va='top')
for ext,dpi in [('pdf',None),('png',300),('svg',None),('tiff',600)]: fig.savefig(FIG/f'Figure3_CrossSystem_HospitalVariation_and_AssociatedBurden.{ext}',dpi=dpi,bbox_inches='tight')
plt.close(fig)
