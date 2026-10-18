# Release adaptation: paths only unless noted in docs/RELEASE_CHANGES.md.
import os as _release_os
def _release_path(env, relative=''):
    root = _release_os.environ.get(env)
    if not root: raise RuntimeError('Set environment variable '+env+'; see docs/REPRODUCIBILITY.md')
    return _release_os.path.join(root, relative)
from pathlib import Path
import pandas as pd,numpy as np,pyreadstat,json,hashlib,math
S=Path(_release_path('ODKD_WORK_ROOT', 'Stage1'));H=Path(_release_path('NHANES_ROOT', 'archive'));checks=[];flow=[];raws=[];manifest=[]
def read(c,stem):
 p=H/'xpt'/c/(stem+'.xpt');d=pyreadstat.read_xport(p,encoding='latin1')[0];assert d.SEQN.is_unique
 manifest.append(dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()));return d
for i,y in enumerate(range(1999,2019,2)):
 c=f'{y}-{y+1}';su='' if i==0 else '_'+chr(65+i)
 demo=read(c,'DEMO'+su);q=read(c,'DIQ'+su);gh=read(c,['LAB10','L10_B','L10_C'][i] if i<3 else 'GHB'+su);u=read(c,['LAB16','L16_B','L16_C'][i] if i<3 else 'ALB_CR'+su);cr=read(c,['LAB18','L40_B','L40_C'][i] if i<3 else 'BIOPRO'+su)
 d=demo[['SEQN','RIDAGEYR','RIAGENDR','RIDEXPRG']].merge(q.reindex(columns=['SEQN','DIQ010','DIQ050','DIQ070']),on='SEQN',how='left',validate='one_to_one').merge(gh[['SEQN','LBXGH']],on='SEQN',how='left',validate='one_to_one').merge(u[['SEQN','URXUMA','URXUCR']],on='SEQN',how='left',validate='one_to_one')
 crcol='LBDSCR' if y==2001 else 'LBXSCR';d=d.merge(cr[['SEQN',crcol]].rename(columns={crcol:'cr'}),on='SEQN',how='left',validate='one_to_one')
 if y==1999:d['cr']=d.cr*1.013+.147
 if y==2005:d['cr']=d.cr*.978-.016
 d['egfr']=[142*min(v/(.7 if sx==2 else .9),1)**(-.241 if sx==2 else -.302)*max(v/(.7 if sx==2 else .9),1)**(-1.2)*.9938**a*(1.012 if sx==2 else 1) if pd.notna(v) and v>0 and pd.notna(a) else np.nan for v,sx,a in zip(d.cr,d.RIAGENDR,d.RIDAGEYR)]
 p=H/'mortality'/f'NHANES_{y}_{y+1}_MORT_2019_PUBLIC.dat';mf=pd.read_fwf(p,colspecs=[(0,6),(14,15),(15,16),(45,48)],names=['SEQN','eligible','death','months']);assert mf.SEQN.is_unique
 mf=mf.apply(pd.to_numeric,errors='coerce')
 manifest.append(dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()));d=d.merge(mf,on='SEQN',how='left',validate='one_to_one')
 adult=d.RIDAGEYR.ge(18)&d.RIDEXPRG.ne(1);dm=adult&(d.DIQ010.eq(1)|d.DIQ050.eq(1)|d.DIQ070.eq(1)|d.LBXGH.ge(6.5));paired=d.URXUMA.ge(0)&d.URXUCR.gt(0);linked=d.eligible.eq(1)&d.death.isin([0,1])&d.months.ge(0);g=dm&paired&linked&d.egfr.ge(60)
 flow.append(dict(cycle=c,all_demo=len(d),adult_DM=int(dm.sum()),paired=int((dm&paired).sum()),paired_linked=int((dm&paired&linked).sum()),final=int(g.sum()),deaths=int(d.loc[g,'death'].sum())))
 d=d[g].copy();d['group']=[('G3' if alb>=30 else 'G2') if 100*alb/(((1.02*math.sqrt(uc)-.36)**2 if uc<75 else (1.05*math.sqrt(uc)-.74)**2 if uc<250 else (1.01*math.sqrt(uc)-.10)**2) if y<2007 else uc)>=30 else ('G1' if alb>=30 else 'G0') for alb,uc in zip(d.URXUMA,d.URXUCR)];raws.append(d)
a=pd.concat(raws);b=pd.read_parquet(S/'data/cohort.parquet');assert set(a.SEQN)==set(b.SEQN);q=a.merge(b.drop(columns='group_raw'),on='SEQN',suffixes=('_raw','_new'),validate='one_to_one')
for v in ['group','death','months']:
 assert (q[v+'_raw'].to_numpy()==q[v+'_new'].to_numpy()).all(),v;checks.append(dict(check='raw independent '+v,result='PASS',n=len(q)))
assert np.allclose(q.egfr_raw,q.egfr_new);assert np.allclose(q.LBXGH,q.hba1c,equal_nan=True);checks.append(dict(check='raw independent CKD-EPI and HbA1c',result='PASS',n=len(q)))
pd.DataFrame(flow).to_csv(S/'audit/independent_flow.csv',index=False)
# Unrelated legacy-v3 file-hash check removed from release.
c=pd.read_csv(S/'audit/point_estimate_checks.csv');assert c.max_coef_difference.max()<1e-6;checks.append(dict(check='svycoxph coefficients versus independent weighted Cox',result='PASS',n=len(c)))
for key,m in [('completeB',b.ccB),('excluded_missingB',~b.ccB)]:
 for group,d in b[m].groupby('group'):
  checks.append(dict(check=key+' '+group,result='descriptive',n=len(d),deaths=int(d.death.sum()),age_mean=d.age.mean(),female=d.sex.eq(2).mean(),mean_followup=d.months.mean()))
pd.DataFrame(checks).to_csv(S/'audit/independent_QC.csv',index=False)
(S/'audit/independent_raw_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf8');print('raw independent audit PASS',len(b),'records; source hashes',len(manifest))
