# Release adaptation: paths only unless noted in docs/RELEASE_CHANGES.md.
import os as _release_os
def _release_path(env, relative=''):
    root = _release_os.environ.get(env)
    if not root: raise RuntimeError('Set environment variable '+env+'; see docs/REPRODUCIBILITY.md')
    return _release_os.path.join(root, relative)
from pathlib import Path
import pandas as pd,numpy as np,pyreadstat,json,hashlib
S=Path(_release_path('ODKD_WORK_ROOT', 'Stage1')); H=Path(_release_path('NHANES_ROOT', 'archive/xpt'))
n=pd.read_parquet(S.parent/'Stage0/data/NHANES_diabetes_mortality_QC.parquet')
sources=[]; demos=[];cov=[]
for i,y in enumerate(range(1999,2019,2)):
 cycle=f'{y}-{y+1}';suffix='' if i==0 else '_'+chr(65+i)
 parts=[]
 for stem,cols in [('DEMO',['SEQN','SDMVSTRA','SDMVPSU','WTMEC2YR','WTMEC4YR','RIDRETH1','INDFMPIR']),('SMQ',['SEQN','SMQ020','SMQ040']),('BPQ',['SEQN','BPQ020']),('MCQ',['SEQN','MCQ160B','MCQ160C','MCQ160D','MCQ160E','MCQ160F']),('BMX',['SEQN','BMXBMI'])]:
  p=H/cycle/(stem+suffix+'.xpt');d,meta=pyreadstat.read_xport(str(p));assert d.SEQN.is_unique
  sources.append(dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),n=len(d),fields={c:meta.column_names_to_labels.get(c) for c in cols if c in d}))
  d=d.reindex(columns=cols)
  if stem=='DEMO':
   d['cycle']=cycle;d['year']=y;d['weight']=d.WTMEC4YR*.2 if y<2003 else d.WTMEC2YR*.1;demos.append(d.copy())
  parts.append(d)
 c=parts[0]
 for d in parts[1:]:c=c.merge(d,on='SEQN',how='left',validate='one_to_one')
 cov.append(c)
c=pd.concat(cov,ignore_index=True);all_demo=pd.concat(demos,ignore_index=True);assert c.SEQN.is_unique
z=n[n.paired & n.linked & n.egfr_cr.ge(60)].copy();assert z.SEQN.is_unique
z['HbA1c']=z['LBXGH'] # raw cross-cycle field, including 2017-2018 added in Stage0
fields=['SEQN','cycle','age','sex','HbA1c','egfr_cr','ualb_mgl','ucr_mgdl','uacr_mgg','scr_mgdl','PERMTH_EXM','MORTSTAT']
z=z[fields].merge(c.drop(columns='cycle'),on='SEQN',how='left',validate='one_to_one')
x=z.ucr_mgdl
z['ucr_adjusted']=np.where(z.year<2007,np.select([x<75,x<250],[(1.02*np.sqrt(x)-.36)**2,(1.05*np.sqrt(x)-.74)**2],default=(1.01*np.sqrt(x)-.10)**2),x)
z['acr_adjusted']=100*z.ualb_mgl/z.ucr_adjusted
def group(a,u):return np.select([(a>=30)&(u<30),(a<30)&(u>=30),(a>=30)&(u>=30)],['G1','G2','G3'],default='G0')
z['group']=group(z.ualb_mgl,z.acr_adjusted);z['group_raw']=group(z.ualb_mgl,z.uacr_mgg)
z['smoking']=np.select([z.SMQ020.eq(2),z.SMQ020.eq(1)&z.SMQ040.eq(3),z.SMQ020.eq(1)&z.SMQ040.isin([1,2])],['never','former','current'],default=None)
z['htn']=z.BPQ020.map({1:1,2:0});m=z[['MCQ160'+v for v in 'BCDEF']]
z['cvd']=np.select([m.eq(1).any(axis=1),m.eq(2).all(axis=1)],[1.,0.],default=np.nan)
z=z.rename(columns={'BMXBMI':'bmi','INDFMPIR':'pir','RIDRETH1':'race','HbA1c':'hba1c','egfr_cr':'egfr','PERMTH_EXM':'months','MORTSTAT':'death'})
z['time']=np.maximum(z.months,.5);z['age10']=z.age/10;z['egfr10']=z.egfr/10;z['bmi5']=z.bmi/5
A=['age10','sex','race','cycle'];B=A+['bmi5','hba1c','egfr10','smoking','htn','cvd','pir']
z['ccA']=z[A].notna().all(axis=1);z['ccB']=z[B].notna().all(axis=1)
z['domain']=True
for f in ['cohort','source_manifest']:
 pass
z.to_parquet(S/'data/cohort.parquet',index=False);z.to_csv(S/'data/cohort.csv',index=False)
base=all_demo[['SEQN','SDMVSTRA','SDMVPSU','weight','year','cycle']]
full=base.merge(z.drop(columns=['SDMVSTRA','SDMVPSU','weight','year','cycle']),on='SEQN',how='left',validate='one_to_one')
for k in ['domain','ccA','ccB']:full[k]=full[k].fillna(False).astype(int)
full=full[full.weight.gt(0)];full.to_csv(S/'data/survey_population.csv',index=False)
pd.DataFrame([dict(variable=v,missing=int(z[v].isna().sum()),fraction=float(z[v].isna().mean())) for v in B]).to_csv(S/'audit/missing.csv',index=False)
z.groupby(['cycle','group'],observed=True).agg(n=('SEQN','size'),deaths=('death','sum'),median_months=('months','median'),completeB=('ccB','sum')).reset_index().to_csv(S/'audit/counts.csv',index=False)
pd.crosstab(z.group_raw,z.group).to_csv(S/'audit/switches.csv')
rows=[]
for g,d in z.groupby('group'):
 row=dict(group=g,n=len(d),death=int(d.death.sum()),ccB=int(d.ccB.sum()),death_ccB=int(d.loc[d.ccB,'death'].sum()),median_months=d.months.median(),censored_before60=int((d.death.eq(0)&d.months.lt(60)).sum()))
 for v in ['age','bmi','hba1c','egfr','ucr_mgdl','ucr_adjusted','ualb_mgl','acr_adjusted','pir']:
  row[v+'_median']=d[v].median();row[v+'_q25']=d[v].quantile(.25);row[v+'_q75']=d[v].quantile(.75)
 row['female_fraction']=d.sex.eq(2).mean();rows.append(row)
pd.DataFrame(rows).to_csv(S/'audit/group_description.csv',index=False)
(S/'audit/source_manifest.json').write_text(json.dumps(sources,ensure_ascii=False,indent=2),encoding='utf8')
print('N',len(z),'deaths',z.death.sum(),'zero_months',z.months.eq(0).sum(),'ccB',z.ccB.sum())
print(pd.read_csv(S/'audit/missing.csv').to_string(index=False));print(z.groupby('group').agg(n=('SEQN','size'),death=('death','sum')).to_string())
print('design',len(full),full.SDMVSTRA.nunique(),full[['SDMVSTRA','SDMVPSU']].drop_duplicates().shape[0])
