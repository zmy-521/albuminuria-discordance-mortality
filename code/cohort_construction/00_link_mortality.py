# Release adaptation: paths only unless noted in docs/RELEASE_CHANGES.md.
import os as _release_os
def _release_path(env, relative=''):
    root = _release_os.environ.get(env)
    if not root: raise RuntimeError('Set environment variable '+env+'; see docs/REPRODUCIBILITY.md')
    return _release_os.path.join(root, relative)
from pathlib import Path
import pandas as pd,numpy as np,json,pyreadstat,re,hashlib
S=Path(_release_path('ODKD_WORK_ROOT', 'Stage0'));R=Path(_release_path('ODKD_WORK_ROOT', 'raw_reconstruction'));H=Path(_release_path('NHANES_ROOT', ''))
tables={}

def save(name,rows):
 d=pd.DataFrame(rows);d.to_csv(S/f'audit/{name}.csv',index=False,encoding='utf-8-sig');tables[name]=d;print(name,len(d),flush=True)


def eg(cr,age,sex,cys=None):
 f=sex.eq(2);q=cr.where(cr.gt(0))/np.where(f,.7,.9)
 if cys is None:return 142*np.minimum(q,1)**np.where(f,-.241,-.302)*np.maximum(q,1)**-1.2*.9938**age*np.where(f,1.012,1)
 t=cys.where(cys.gt(0))/.8
 return 135*np.minimum(q,1)**np.where(f,-.219,-.144)*np.maximum(q,1)**-.544*np.minimum(t,1)**-.323*np.maximum(t,1)**-.778*.9961**age*np.where(f,.963,1)


def xpt(p):return pyreadstat.read_xport(p)[0]


def nhanes():
 n=pd.read_parquet(R/'data/NHANES_canonical.parquet');n=n[n.cycle.str[:4].astype(int).le(2015)].copy()
 folder=H/'archive/xpt/2017-2018';j=xpt(folder/'DEMO_J.xpt')
 for name in ['DIQ_J','GHB_J','BIOPRO_J','ALB_CR_J']:
  z=xpt(folder/(name+'.xpt'));j=j.merge(z,on='SEQN',how='left',validate='one_to_one')
 j['cycle']='2017-2018';j['age']=j.RIDAGEYR;j['sex']=j.RIAGENDR;j['diabetes']=j.DIQ010.eq(1)|j.DIQ050.eq(1)|j.DIQ070.eq(1)|j.LBXGH.ge(6.5)
 j['scr_mgdl']=j.LBXSCR;j['ualb_mgl']=j.URXUMA;j['ucr_mgdl']=j.URXUCR;j['uacr_mgg']=100*j.URXUMA/j.URXUCR.where(j.URXUCR.gt(0));j['egfr_cr']=eg(j.scr_mgdl,j.age,j.sex)
 n=pd.concat([n,j],ignore_index=True);assert not n.SEQN.duplicated().any()
 lm=[]
 for file in sorted((H/'archive/mortality').glob('NHANES_*_MORT_2019_PUBLIC.dat')):
  if not re.fullmatch(r'NHANES_\d{4}_\d{4}_MORT_2019_PUBLIC.dat',file.name):continue
  z=pd.read_fwf(file,colspecs=[(0,6),(14,15),(15,16),(16,19),(19,20),(20,21),(42,45),(45,48)],names=['SEQN','ELIGSTAT','MORTSTAT','UCOD_LEADING','DIABETES_DEATH','HTN_DEATH','PERMTH_INT','PERMTH_EXM'],na_values=['.']);lm.append(z)
 l=pd.concat(lm);assert l.SEQN.is_unique
 n=n.merge(l,on='SEQN',how='left',validate='one_to_one');n=n[n.age.ge(18)&n.diabetes.fillna(False)&~n.RIDEXPRG.eq(1)].copy()
 n['paired']=n.ualb_mgl.notna()&n.uacr_mgg.notna();n['linked']=n.ELIGSTAT.eq(1)&n.MORTSTAT.isin([0,1])&n.PERMTH_EXM.ge(0)
 n['urine_group']=np.select([n.ualb_mgl.ge(30)&n.uacr_mgg.ge(30),n.ualb_mgl.ge(30)&n.uacr_mgg.lt(30),n.ualb_mgl.lt(30)&n.uacr_mgg.ge(30)],['both_positive','concentration_only','ACR_only'],default='both_negative')
 n.loc[~n.paired,'urine_group']='not_evaluable'
 n.to_parquet(S/'data/NHANES_diabetes_mortality.parquet',index=False)
 rows=[]
 for cycle,z in n.groupby('cycle'):
  rows.append(dict(cycle=cycle,adult_diabetes=len(z),linked=int(z.linked.sum()),paired=int(z.paired.sum()),paired_linked=int((z.paired&z.linked).sum()),paired_linked_ge60=int((z.paired&z.linked&z.egfr_cr.ge(60)).sum()),deaths=int((z.paired&z.linked&z.MORTSTAT.eq(1)).sum())))
 save('NHANES_cycles',rows);rows=[]
 for pop,mask in [('all_GFR',n.paired&n.linked),('eGFRcr_ge60',n.paired&n.linked&n.egfr_cr.ge(60))]:
  for group,z in n[mask].groupby('urine_group'):
   rows.append(dict(population=pop,group=group,n=len(z),death=int(z.MORTSTAT.eq(1).sum()),CVD_death=int((z.MORTSTAT.eq(1)&z.UCOD_LEADING.isin([1,5])).sum()),heart_death=int((z.MORTSTAT.eq(1)&z.UCOD_LEADING.eq(1)).sum()),death_unknown_cause=int((z.MORTSTAT.eq(1)&z.UCOD_LEADING.isna()).sum()),median_months=z.PERMTH_EXM.median(),death_5y=int((z.MORTSTAT.eq(1)&z.PERMTH_EXM.le(60)).sum()),known_alive_at5y=int((z.PERMTH_EXM.ge(60)&~(z.MORTSTAT.eq(1)&z.PERMTH_EXM.eq(60))).sum()),censored_before5y=int((z.MORTSTAT.eq(0)&z.PERMTH_EXM.lt(60)).sum())))
 save('NHANES_deaths',rows)


if __name__=="__main__": nhanes()
