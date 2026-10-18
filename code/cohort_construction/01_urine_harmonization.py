# Release adaptation: paths only unless noted in docs/RELEASE_CHANGES.md.
import os as _release_os
def _release_path(env, relative=''):
    root = _release_os.environ.get(env)
    if not root: raise RuntimeError('Set environment variable '+env+'; see docs/REPRODUCIBILITY.md')
    return _release_os.path.join(root, relative)
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,pyreadstat,math
S=Path(_release_path('ODKD_WORK_ROOT', 'Stage0'));R=Path(_release_path('ODKD_WORK_ROOT', 'raw_reconstruction'));H=Path(_release_path('NHANES_ROOT', ''))
def out(name,rows):pd.DataFrame(rows).to_csv(S/f'audit/{name}.csv',index=False,encoding='utf-8-sig')
n=pd.read_parquet(S/'data/NHANES_diabetes_mortality.parquet');x=n.ucr_mgdl
adj=np.select([x.lt(75),x.lt(250)],[(1.02*np.sqrt(x)-.36)**2,(1.05*np.sqrt(x)-.74)**2],default=(1.01*np.sqrt(x)-.10)**2)
n['ucr_method_adjusted']=np.where(n.cycle.str[:4].astype(int).lt(2007),adj,x)
n['uacr_method_adjusted']=100*n.ualb_mgl/n.ucr_method_adjusted
n['adjusted_group']=np.select([n.ualb_mgl.ge(30)&n.uacr_method_adjusted.ge(30),n.ualb_mgl.ge(30)&n.uacr_method_adjusted.lt(30),n.ualb_mgl.lt(30)&n.uacr_method_adjusted.ge(30)],['both_positive','concentration_only','ACR_only'],default='both_negative')
rows=[]
for pop,m in [('all_GFR',n.paired&n.linked),('eGFRcr_ge60',n.paired&n.linked&n.egfr_cr.ge(60))]:
 for group,z in n[m].groupby('adjusted_group'):
  rows.append(dict(population=pop,group=group,n=len(z),death=int(z.MORTSTAT.eq(1).sum()),heart_death=int((z.MORTSTAT.eq(1)&z.UCOD_LEADING.eq(1)).sum()),death_5y=int((z.MORTSTAT.eq(1)&z.PERMTH_EXM.le(60)).sum()),censored_before5y=int((z.MORTSTAT.eq(0)&z.PERMTH_EXM.lt(60)).sum())))
out('NHANES_adjusted_groups',rows)
out('NHANES_method_sensitivity',[dict(population=p,n=int(m.sum()),classification_changed=int((m&n.urine_group.ne(n.adjusted_group)).sum()),ACR_up30=int((m&n.uacr_mgg.lt(30)&n.uacr_method_adjusted.ge(30)).sum()),ACR_down30=int((m&n.uacr_mgg.ge(30)&n.uacr_method_adjusted.lt(30)).sum())) for p,m in [('all_GFR',n.paired&n.linked),('eGFRcr_ge60',n.paired&n.linked&n.egfr_cr.ge(60))]])
n.to_parquet(S/'data/NHANES_diabetes_mortality_QC.parquet',index=False)
