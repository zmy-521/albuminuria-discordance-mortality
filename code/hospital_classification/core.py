# Release adaptation: paths only unless noted in docs/RELEASE_CHANGES.md.
import os as _release_os
def _release_path(env, relative=''):
    root = _release_os.environ.get(env)
    if not root: raise RuntimeError('Set environment variable '+env+'; see docs/REPRODUCIBILITY.md')
    return _release_os.path.join(root, relative)
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
R=Path(_release_path('ODKD_WORK_ROOT', 'raw_reconstruction'));V=R/'v2'
ANALYTES=['scr_mgdl','cys_mgl','ualb_mgl','ucr_mgdl','uacr_direct','HGB','UA','HbA1c','ALB','TP','GLB','AG','BUN','GLU','TG','HDL','LDL','CHOL','TBIL','PAB','ALP','ALT','AST','GGT','WBC','NEUT','LYMPH','MONO','BASO','EO','PLT']
URINE=['ualb_mgl','ucr_mgdl','uacr_direct']
def dump(o,p):Path(p).write_text(json.dumps(o,ensure_ascii=False,indent=2,default=str),encoding='utf-8')
def calculate(d):
 d=d.copy();f=d.sex.eq(2);q=d.scr_mgdl/np.where(f,.7,.9);c=d.cys_mgl/.8
 d['egfr_cr']=142*np.minimum(q,1)**np.where(f,-.241,-.302)*np.maximum(q,1)**-1.2*.9938**d.age*np.where(f,1.012,1)
 d['egfr_cr_cys']=135*np.minimum(q,1)**np.where(f,-.219,-.144)*np.maximum(q,1)**-.544*np.minimum(c,1)**-.323*np.maximum(c,1)**-.778*.9961**d.age*np.where(f,.963,1)
 d.loc[~d.sex.isin([1,2])|~d.age.ge(18),['egfr_cr','egfr_cr_cys']]=np.nan
 d['uacr_calculated']=d.ualb_mgl/d.ucr_mgdl.where(d.ucr_mgdl.gt(0))*100
 d['uacr_mgg']=d.uacr_calculated.combine_first(d.uacr_direct)
 d['uacr_source']=np.select([d.uacr_calculated.notna(),d.uacr_direct.notna()],['same_specimen_calculated','direct'],default='missing')
 d['base_eligible']=d.identity_usable&d.diabetes&d.age.ge(18)&d.egfr_cr.ge(60)
 b=d.base_eligible;a=d.ualb_mgl.notna();u=d.uacr_mgg.notna();c=d.egfr_cr_cys.notna();g=d.egfr_cr_cys.ge(60);m=d.uacr_mgg.lt(300)
 for j,(e,y) in enumerate(zip([b&a,b&u,b&u&m,b&u&c&g,b&u&c&g&m],[d.ualb_mgl.ge(30)]+[d.uacr_mgg.ge(30)]*4),1):
  d[f'D{j}_eligible']=e;d[f'D{j}']=y.astype('Int8').where(e)
 d['common_measurement_complete']=b&a&u&c;d['all5_common_eligible']=d[[f'D{j}_eligible' for j in range(1,6)]].all(axis=1)
 d['urine_axis_eligible']=b&a&u
 d['urine_category']=np.select([a&u&d.ualb_mgl.ge(30)&d.uacr_mgg.ge(30),a&u&d.ualb_mgl.ge(30)&d.uacr_mgg.lt(30),a&u&d.ualb_mgl.lt(30)&d.uacr_mgg.ge(30),a&u&d.ualb_mgl.lt(30)&d.uacr_mgg.lt(30)],['both_positive','concentration_only','ACR_only','both_negative'],default='not_evaluable')
 d['qc_extreme_scr']=(d.scr_mgdl.lt(.2)|d.scr_mgdl.gt(20))&d.scr_mgdl.notna()
 d['qc_extreme_cys']=(d.cys_mgl.lt(.2)|d.cys_mgl.gt(10))&d.cys_mgl.notna()
 d['ratio_threshold_disagreement']=d.uacr_calculated.notna()&d.uacr_direct.notna()&d.uacr_calculated.ge(30).ne(d.uacr_direct.ge(30))
 return d
