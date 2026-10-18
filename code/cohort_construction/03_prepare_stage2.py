# Release adaptation: paths only unless noted in docs/RELEASE_CHANGES.md.
import os as _release_os
def _release_path(env, relative=''):
    root = _release_os.environ.get(env)
    if not root: raise RuntimeError('Set environment variable '+env+'; see docs/REPRODUCIBILITY.md')
    return _release_os.path.join(root, relative)
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
S=Path(_release_path('ODKD_WORK_ROOT', 'Stage2'));P=S.parent/'Stage1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
# Internal historical release-manifest guard omitted; public inputs must be versioned locally.
d=pd.read_parquet(P/'data/cohort.parquet');assert len(d)==5905 and d.SEQN.is_unique
d['group20']=np.select([d.ualb_mgl.ge(20)&d.acr_adjusted.lt(30),d.ualb_mgl.lt(20)&d.acr_adjusted.ge(30),d.ualb_mgl.ge(20)&d.acr_adjusted.ge(30)],['G1','G2','G3'],default='G0')
d['ccNoPIR']=d[['age10','sex','race','cycle','bmi5','hba1c','egfr10','smoking','htn','cvd']].notna().all(axis=1)
d['unknown60']=d.death.eq(0)&d.months.lt(60);d['y60']=np.where(d.unknown60,np.nan,(d.death.eq(1)&d.months.le(60)).astype(int))
d.to_parquet(S/'data/cohort_stage2.parquet',index=False)
rows=[]
for t,v in [(30,'group'),(20,'group20')]:
 for g,z in d.groupby(v):
  q=z[z.year.le(2013)];rows.append(dict(threshold_mgl=t,group=g,n=len(z),deaths=int(z.death.sum()),ccB=int(z.ccB.sum()),death_ccB=int(z.loc[z.ccB,'death'].sum()),early_n=len(q),early_ccB=int(q.ccB.sum()),death60=int(q.y60.eq(1).sum()),unknown60=int(q.unknown60.sum())))
pd.DataFrame(rows).to_csv(S/'audit/group_counts.csv',index=False)
pd.crosstab(d.group,d.group20).to_csv(S/'audit/threshold_switches.csv')
rows=[]
for t,v in [(30,'group'),(20,'group20')]:
 for g,z in d.groupby(v):
  for f in ['ucr_adjusted','ualb_mgl','age','bmi','hba1c','egfr']:
   rows.append(dict(threshold=t,group=g,variable=f,n=int(z[f].notna().sum()),median=z[f].median(),q1=z[f].quantile(.25),q3=z[f].quantile(.75)))
  rows.append(dict(threshold=t,group=g,variable='female_fraction',n=len(z),median=z.sex.eq(2).mean(),q1=np.nan,q3=np.nan))
pd.DataFrame(rows).to_csv(S/'audit/descriptive.csv',index=False)
(S/'audit/locked_stage2_plan.md').write_text('''# Stage2 implementation lock, before Stage2 effect estimation
30 mg/L remains Stage1 primary; only 20 mg/L added post hoc by literature and explicit user request. ACR30, calibrated urine, adult DM eGFRcr>=60 and baseline/followup/design unchanged. Viana 2012 uses ROC-derived ~14 mg/L, not the basis for 20. Historical20 source: PMID10932516 and Diabetes Care2005 28:164. No additional cutpoints.
Cox A/B same Stage1 formulas, Breslow ties; A on B sample separately for matching. Supplemental 30mg/L B + log2(UCr); independent component model B without group + log2(UAlb) + log2(UCr); all strictly positive in frozen cohort. No causal decomposition, no interaction search.
NoPIR: same-B sample model and expanded non-PIR-complete sample, to separate covariate omission from sample change. Limited selection weighting: estimate complete-B probability with survey-weighted logistic using fully observed group30, natural-spline age3, sex,race,cycle; multiply MEC by inverse probability among complete cases; no outcome predictors in selection model. If min probability<.05 stop IPW instead of data-driven trimming. JKn re-estimates selection and Cox to reflect both stages. Assumes selection at random conditional on those covariates, not a correction of unobserved missingness.
Five-year: original1999-2014, N4315, B target N3789 unchanged. Primary logit fits known60-month outcome complete-B cases only, predicts all3789 baseline profiles including 4 with unknown outcome. Uses group + Stage1 B covariates, fixed age spline knots from whole target. Unknown survivors before60 are never observed zeros. Additional hypothetical all-unknown0 and all-unknown1 scenarios, labelled assumptions, not observed outcomes. Exhaustively enumerate16 assignments of the4 unknown in B to find range of refitted adjusted RD; the fifth person is outside B and cannot affect adjusted B. This is outcome-assignment sensitivity, not a nonparametric causal bound or CI.
Design inference: original full positive-MEC survey design before domain subset, same weights, strata,PSU. Primary svyglm(quasibinomial) point estimates checked against weighted IRLS. JKn301 PSU repeats for standardized 4 risks and paired G2-G0/G2-G1 differences, re-fit and standardize each repeat; logit-delta probability CIs, normal RD CIs. Both20 and30 run identical logic. Same-patient paired changes between thresholds retain joint replicates. No ML, outcome expansion, imputation search or further tuning.
New20 Cox standardized60-month risk also uses same Stage1 approach for direct method comparison, with JKn; original30 result read from frozen file. Every original source and Stage1/v3 immutable. Full finite/convergence/pairing and raw classification QC before release.
''',encoding='utf8')
print('Prepared Stage2, frozen Stage1 verified; N',len(d),'noPIR',d.ccNoPIR.sum())
