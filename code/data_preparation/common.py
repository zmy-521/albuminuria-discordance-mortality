# Release adaptation: paths only unless noted in docs/RELEASE_CHANGES.md.
import os as _release_os
def _release_path(env, relative=''):
    root = _release_os.environ.get(env)
    if not root: raise RuntimeError('Set environment variable '+env+'; see docs/REPRODUCIBILITY.md')
    return _release_os.path.join(root, relative)
from pathlib import Path
import pandas as pd,numpy as np,hashlib,json,re,os
R=Path(_release_path('ODKD_WORK_ROOT', 'raw_reconstruction'))
S=R/'internal/pseudonym_salt.bin'
if not S.exists():S.write_bytes(os.urandom(32))
salt=S.read_bytes()
def pseud(s):return hashlib.sha256(salt+str(s).encode()).hexdigest()[:24]
def norm(s):return s.fillna('').astype(str).str.strip().str.replace(r'\.0$','',regex=True)
def parse_dt(s):return pd.to_datetime(s,errors='coerce',format='mixed')
def egfr(d):
 female=d.sex.eq(2);k=np.where(female,.7,.9);cr=d.scr_mgdl.where(d.scr_mgdl.gt(0));cys=d.cys_mgl.where(d.cys_mgl.gt(0));q=cr/k
 d['egfr_cr']=142*np.minimum(q,1)**np.where(female,-.241,-.302)*np.maximum(q,1)**-1.2*.9938**d.age*np.where(female,1.012,1)
 d['egfr_cr_cys']=135*np.minimum(q,1)**np.where(female,-.219,-.144)*np.maximum(q,1)**-.544*np.minimum(cys/.8,1)**-.323*np.maximum(cys/.8,1)**-.778*.9961**d.age*np.where(female,.963,1)
 d.loc[~d.sex.isin([1,2])|~d.age.ge(18),['egfr_cr','egfr_cr_cys']]=np.nan
 return d
def labels(d):
 d=egfr(d);d['base_eligible']=d.identity_usable&d.diabetes&d.age.ge(18)&d.egfr_cr.ge(60)
 b=d.base_eligible;a=d.ualb_mgl.notna();u=d.uacr_mgg.notna();c=d.egfr_cr_cys.notna();m=d.uacr_mgg.lt(300);g=d.egfr_cr_cys.ge(60)
 elig=[b&a,b&u,b&u&m,b&u&c&g,b&u&m&c&g];ys=[d.ualb_mgl.ge(30)]+[d.uacr_mgg.ge(30)]*4
 for j,(e,y) in enumerate(zip(elig,ys),1):d[f'D{j}_eligible']=e;d[f'D{j}']=y.astype('Int8').where(e)
 d['common_measurement_complete']=b&a&u&c
 d['all5_common_eligible']=pd.concat(elig,axis=1).all(axis=1)
 d['urine_axis_eligible']=b&a&u
 d['urine_category']=np.select([a&u&d.ualb_mgl.ge(30)&d.uacr_mgg.ge(30),a&u&d.ualb_mgl.ge(30)&d.uacr_mgg.lt(30),a&u&d.ualb_mgl.lt(30)&d.uacr_mgg.ge(30),a&u&d.ualb_mgl.lt(30)&d.uacr_mgg.lt(30)],['both_positive','concentration_only','ACR_only','both_negative'],default='not_evaluable')
 return d
def save_summary(d,center,extra=None):
 d.to_parquet(R/f'data/{center}_canonical.parquet',index=False)
 d[d.base_eligible].to_parquet(R/f'data/{center}_base_analysis.parquet',index=False)
 d[d.all5_common_eligible].to_parquet(R/f'data/{center}_common5.parquet',index=False)
 out={'center':center,'mother_rows':len(d),'unique_keys':d.person_key.nunique(),'identity_usable':int(d.identity_usable.sum()),'adult_diabetes':int((d.age.ge(18)&d.diabetes&d.identity_usable).sum()),'scr_evaluable':int((d.age.ge(18)&d.diabetes&d.identity_usable&d.egfr_cr.notna()).sum()),'base_eligible':int(d.base_eligible.sum()),'common_measurement_complete':int(d.common_measurement_complete.sum()),'all5_common_eligible':int(d.all5_common_eligible.sum()),'definitions':[],'urine_axis':[],'missing':[]}
 for j in range(1,6):
  for pop,mask in [('max',d[f'D{j}_eligible']),('measurement_complete',d.common_measurement_complete&d[f'D{j}_eligible']),('strict_common5',d.all5_common_eligible)]:
   out['definitions'].append({'definition':f'D{j}','population':pop,'n':int(mask.sum()),'positive':int(d.loc[mask,f'D{j}'].sum())})
 for cat,g in d[d.urine_axis_eligible].groupby('urine_category'):
  out['urine_axis'].append({'category':cat,'n':len(g),'age_mean':float(g.age.mean()),'female':int(g.sex.eq(2).sum()),'scr_mean':float(g.scr_mgdl.mean()),'cys_mean':float(g.cys_mgl.mean())})
 g=d[d.identity_usable&d.diabetes&d.age.ge(18)&d.egfr_cr.notna()&d.egfr_cr_cys.notna()]
 out['renal_axis']={'both_measured':len(g),'cr_ge60_cys_ge60':int((g.egfr_cr.ge(60)&g.egfr_cr_cys.ge(60)).sum()),'cr_ge60_cys_lt60':int((g.egfr_cr.ge(60)&g.egfr_cr_cys.lt(60)).sum()),'cr_lt60_cys_ge60':int((g.egfr_cr.lt(60)&g.egfr_cr_cys.ge(60)).sum()),'cr_lt60_cys_lt60':int((g.egfr_cr.lt(60)&g.egfr_cr_cys.lt(60)).sum()),'cr_quantiles':g.egfr_cr.quantile([.05,.25,.5,.75,.95]).to_dict(),'combined_quantiles':g.egfr_cr_cys.quantile([.05,.25,.5,.75,.95]).to_dict(),'near60_cr45_75':int(g.egfr_cr.between(45,75).sum())}
 for col in ['scr_mgdl','cys_mgl','ualb_mgl','ucr_mgdl','uacr_mgg','HGB','UA','HbA1c','ALB','TP','WBC','NEUT','LYMPH','PLT','BMI']:
  if col in d:out['missing'].append({'variable':col,'base_n':int(d.base_eligible.sum()),'missing':int(d.loc[d.base_eligible,col].isna().sum())})
 if extra:out.update(extra)
 (R/f'audit/{center}_summary.json').write_text(json.dumps(out,ensure_ascii=False,indent=2,default=str),encoding='utf-8');print(center,{k:out[k] for k in ['mother_rows','base_eligible','all5_common_eligible']},flush=True)
