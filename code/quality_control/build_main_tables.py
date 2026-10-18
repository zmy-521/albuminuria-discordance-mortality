# Release adaptation: paths only unless noted in docs/RELEASE_CHANGES.md.
import os as _release_os
def _release_path(env, relative=''):
    root = _release_os.environ.get(env)
    if not root: raise RuntimeError('Set environment variable '+env+'; see docs/REPRODUCIBILITY.md')
    return _release_os.path.join(root, relative)
from pathlib import Path
import pandas as pd,numpy as np,json
S=Path(_release_path('ODKD_WORK_ROOT', 'Stage3_finalization'));P=Path(_release_path('ODKD_WORK_ROOT', ''));F=Path(_release_path('ODKD_WORK_ROOT', 'raw_reconstruction/v3_submission'));A=S/'audit'
ps=[]
def sheet(name,d,note):
 d=d.copy()
 for c in d:
  if c.lower() in ['p','p_value','pvalue']:
   d[c]=d[c].map(lambda x:'<1e-300 (numerical underflow)' if pd.notna(x) and x==0 else x)
 ps.append(dict(name=name,title=name,note=note,headers=list(d.columns),rows=json.loads(d.to_json(orient='values',double_precision=15))))
def csv(name,path,note):sheet(name,pd.read_csv(path),note)
summary=pd.read_csv(A/'supplemental_Cox.csv');summary=summary[summary.contrast=='G2_G0']
sheet('Final supplements',summary,'All Stage3 additions are post-Stage2. HR/95%CI; main5905 cohort and B5155 remain frozen.')
csv('Classification',A/'weighted_classification.csv','Main mortality cohort5905; numerator=G2. weighted/low/high are survey proportions; denominators explicitly defined.')
sheet('Hospital classification',pd.read_csv(A/'hospital_classification.csv').drop(columns='weighted'),'Unweighted hospital description only; no mortality validation or diagnostic gold standard.')
csv('Group counts',P/'Stage2/audit/group_counts.csv','Same5905 persons at30 and20mg/L; deaths=all follow-up. ccB=complete modelB; early=1999–2014.')
csv('Hospital groups',A/'hospital_group_counts.csv','Each center denominator remains distinct; no combined independent count.')
# Legacy ML-overlap sheet excluded from clinical release.
cox=pd.read_csv(P/'Stage1/audit/cox_results.csv');cox['version']='Stage1 frozen';c2=pd.read_csv(P/'Stage2/audit/cox_results.csv');c2['version']='Stage2 frozen';c2['scenario']='threshold_'+c2.threshold.astype(str);c2['design']='survey_Taylor'
sheet('Cox evidence',pd.concat([cox,c2],ignore_index=True),'All frozen main/sensitivity Cox results, including unweighted comparisons. No selection by significance.')
csv('Stage3 Cox',A/'supplemental_Cox.csv','DomainB+log2albumin;24month landmark resets time; sex-specific B. All contrasts retained.')
csv('Sex counts',A/'sex_counts.csv','sex1 male; sex2 female. Counts for full cohort and B complete cases.')
csv('Sex joint test',A/'sex_interaction.csv','Only one3df joint group×sex survey Wald test; no subgroup selection.')
csv('Landmark flow',A/'landmark_flow.csv','time>24 months required; Surv(time-24,death). Death/censoring≤24 excluded; not original time zero.')
csv('Period Cox',P/'Stage1/audit/PH_period_HR.csv','Frozen0–60/after60 time-varying group model; different estimands, not replication of5year risk.')
r=pd.read_csv(P/'Stage1/audit/risk60.csv');r['method']='Stage1 Cox or KM30';r2=pd.read_csv(P/'Stage2/audit/fixed_risk.csv');r2['method']='Stage2 fixed60 logistic';r3=pd.read_csv(P/'Stage2/audit/cox20_risk.csv');r3['method']='Stage2 Cox20'
sheet('Absolute risk',pd.concat([r,r2,r3],ignore_index=True),'Probability units. known=primary logistic; all0/all1=hypothetical censoring assignment, not observed survival. JKn301.')
csv('Selection sensitivity',P/'Stage2/audit/selection_cox.csv','Frozen complete-case IPW; conditional missingness assumptions. NoPIR contrasts also in Cox evidence.')
csv('UCr common support',P/'Stage2/audit/ucr_support.csv','UCr mg/dL; G1/G2 disjoint support imposed by definitions; adjusted attenuation is not causal explanation.')
(S/'internal/workbook_payload.json').write_text(json.dumps(ps,ensure_ascii=False,allow_nan=False),encoding='utf8')
