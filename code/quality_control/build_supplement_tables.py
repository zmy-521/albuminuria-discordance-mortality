# Release adaptation: paths only unless noted in docs/RELEASE_CHANGES.md.
import os as _release_os
def _release_path(env, relative=''):
    root = _release_os.environ.get(env)
    if not root: raise RuntimeError('Set environment variable '+env+'; see docs/REPRODUCIBILITY.md')
    return _release_os.path.join(root, relative)
from pathlib import Path
import pandas as pd,numpy as np,json
S=Path(_release_path('ODKD_WORK_ROOT', 'Review_supplement_urine_creatinine'));A=S/'audit';ps=[]
def tab(name,file,note):
 d=pd.read_csv(A/file)
 for c in d:
  if c in ['p','F_p']:d[c]=d[c].map(lambda x:'<1e-300' if pd.notna(x) and x==0 else x)
 ps.append(dict(name=name,title=name,note=note,headers=list(d.columns),rows=json.loads(d.to_json(orient='values',double_precision=15))))
tab('Effects','effect_estimates.csv','Original linear eGFR is primary; one fixed RCS3 sensitivity. p=design-covariance Wald chi-square; F_p=finite-df reference. No causal claims.')
tab('Joint covariance','contrast_covariance_summary.csv','Var(sum)=VarAlb+VarUCr+2Cov. Use the full design covariance, not independent coefficient errors.')
tab('A2 counts','A2_counts.csv','A2=30≤UACR<300. G2 versus G3-A2; originalB complete cases only in fitted model. All-domain counts separately shown.')
tab('A2 characteristics','A2_weighted_descriptives.csv','Survey-weighted means with95%CI and descriptive weighted quartiles. female0/1; UCrmg/dL; Albmg/L; ACRmg/g.')
tab('A2 support','A2_common_support.csv','Min/max shared support is descriptive, not proof of sufficient overlap. No post-hoc trimming/refit.')
tab('A2 support bins','A2_support_bins.csv','Fixed descriptive intervals only; not additional model thresholds or subgroup searches.')
tab('Flow groups','flow_counts.csv','2009–2018 valid composite collection subset; deaths are counts only, no new mortality model fitted.')
tab('Flow characteristics','flow_weighted_descriptives.csv','flow mL/min; AERµg/min; CERmg/min; SCrmg/dL. Means95%CI; weighted quartiles descriptive. Not24h excretion.')
tab('G2 AER20','G2_AER20.csv','G2 short-interval AER≥20µg/min; survey logit CI. Not persistent albuminuria or diagnostic gold standard.')
tab('Flow missingness','flow_missingness.csv','Required time missing => unassessable; no imputing missing flow or coding normal.')
tab('Flow selection','flow_selection_description.csv','Valid and missing flow domains compared descriptively; no model of mortality or nonresponse correction.')
tab('Detection limits','LOD_audit.csv','Missing LC fields mean unknown flag availability, not zero below-LOD. Official released values retained unchanged.')
tab('Retained extremes','retained_extremes.csv','True observed positive values retained; no winsorization or outcome-dependent exclusion.')
tab('Original model QC','frozen_model_audit.csv','B5155/1115;148 strata301PSU; original components omit group and reproduce saved coefficients/covariance.')
tab('Independent QC','independent_QC.csv','Independent scalar/matrix, raw hashes, paired rows and rate identities; reported rate discrepancy only rounding.')
tab('eGFR knots','fixed_eGFR_knots.csv','One locked3df restricted/natural cubic basis with4 knots from fullB eGFR5/35/65/95 percentiles; not selected by effects.')
(S/'internal/workbook_payload.json').write_text(json.dumps(ps,ensure_ascii=False,allow_nan=False),encoding='utf8')
