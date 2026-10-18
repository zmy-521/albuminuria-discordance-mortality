# Release adaptation: paths only unless noted in docs/RELEASE_CHANGES.md.
import os as _release_os
def _release_path(env, relative=''):
    root = _release_os.environ.get(env)
    if not root: raise RuntimeError('Set environment variable '+env+'; see docs/REPRODUCIBILITY.md')
    return _release_os.path.join(root, relative)
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,math,pyreadstat
from scipy.stats import norm,chi2
S=Path(_release_path('ODKD_WORK_ROOT', 'Review_supplement_urine_creatinine'));P=Path(_release_path('ODKD_WORK_ROOT', ''));qc=[]
# Internal original-project hash manifest omitted.
src=json.loads((S/'audit/raw_sources.json').read_text(encoding='utf8'))
for r in src:assert hashlib.sha256(Path(r['path']).read_bytes()).hexdigest()==r['sha256']
rounding=[]
for rec in src:
 if Path(rec['path']).name.startswith('ALB_CR'):
  raw,_=pyreadstat.read_xport(rec['path'],encoding='latin1')
  if 'URDACT' in raw:
   valid=raw[['URDACT','URXUMA','URXUCR']].notna().all(axis=1)&raw.URXUCR.gt(0)
   er=(raw.loc[valid,'URDACT']-100*raw.loc[valid,'URXUMA']/raw.loc[valid,'URXUCR']).abs()
   assert er.max()<.00501
   rounding.append(dict(source=Path(rec['path']).name,n=int(valid.sum()),max_rounding=er.max()))
pd.DataFrame(rounding).to_csv(S/'audit/raw_UACR_unit_check.csv',index=False)
q=pd.read_parquet(S/'data/review_cohort.parquet');old=pd.read_parquet(P/'Stage2/data/cohort_stage2.parquet');assert old.equals(q[old.columns]);assert q.SEQN.is_unique
qc.append(dict(check='Original5905 records and every frozen column preserved',n=len(q),max_difference=0,result='PASS'))
b=pd.read_csv(S/'audit/component_coefficients.csv').set_index('term').coefficient;V=pd.read_csv(S/'audit/component_design_covariance.csv',index_col=0);assert not any('group' in k for k in b.index)
e=b.log2alb+b.log2ucr;v=V.loc['log2alb','log2alb']+V.loc['log2ucr','log2ucr']+2*V.loc['log2alb','log2ucr'];se=math.sqrt(v)
z=pd.read_csv(S/'audit/effect_estimates.csv');r=z[(z.scenario=='original_linear_eGFR')&(z.contrast=='UCr_per_doubling_at_fixed_ACR')].iloc[0]
for a,c in [(e,r.logHR),(v,r.variance),(math.exp(e),r.HR),(math.exp(e-1.96*se),r.low),(chi2.sf(e*e/v,1),r.p)]:assert math.isclose(a,c,rel_tol=1e-10,abs_tol=1e-12)
qc.append(dict(check='Independent full covariance sum, Wald and CI',n=5155,max_difference=abs(e-r.logHR),result='PASS'))
# Independent vectorized composite computation, not the original rowwise function.
v=q[['URXVOL1','URXVOL2','URXVOL3']].to_numpy(float);t=q[['URDTIME1','URDTIME2','URDTIME3']].to_numpy(float);last=np.max(np.where(np.isfinite(v)&(v>0),np.arange(1,4),0),axis=1)
inside=np.arange(1,4)[None,:]<=last[:,None];known=((~inside)|(np.isfinite(v)&(v>=0)&np.isfinite(t)&(t>0))).all(axis=1)
amb=((~inside)&~np.isfinite(v)&np.isfinite(t)).any(axis=1);ok=(last>0)&known&~amb&q.year.ge(2009).to_numpy()
vv=np.where(inside,v,0).sum(axis=1);tt=np.where(inside,t,0).sum(axis=1);pred=np.divide(vv,tt,out=np.full(len(q),np.nan),where=ok)
assert np.array_equal(ok,q.flow_valid.eq(1));assert np.allclose(pred[ok],q.loc[ok,'flow'])
assert np.allclose(q.loc[ok,'AER'],q.loc[ok,'ualb_mgl']*pred[ok]);assert np.allclose(q.loc[ok,'CER'],q.loc[ok,'ucr_mgdl']*pred[ok]/100)
qc.append(dict(check='Independent pooled volume/time, AER and CER',n=int(ok.sum()),max_difference=float(np.max(abs(pred[ok]-q.loc[ok,'flow']))),result='PASS'))
err=[];desc=pd.read_csv(S/'audit/flow_weighted_descriptives.csv')
for r in desc.itertuples():
 a=q[q.flow_valid.eq(1)&q.group.eq(r.population)];a=a[a[r.variable].notna()];mu=np.average(a[r.variable],weights=a.weight);err.append(abs(mu-r.mean));assert abs(mu-r.mean)<1e-10
qc.append(dict(check='Weighted means independent of R survey function',n=len(desc),max_difference=max(err),result='PASS'))
g=q[q.flow_valid.eq(1)&q.group.eq('G2')];res=pd.read_csv(S/'audit/G2_AER20.csv').iloc[0];assert len(g)==res.n and g.AER.ge(20).sum()==res['count'];assert abs(np.average(g.AER.ge(20),weights=g.weight)-res.weighted)<1e-12
early=q[q.year<=2013];cc=early[early.ccB];unknown=cc.death.eq(0)&cc.months.lt(60);assert len(early)==4315 and len(cc)==3789 and unknown.sum()==4
qc.append(dict(check='Frozen5year1999-2014 target/known outcomes audit',n=len(cc),max_difference=0,result='PASS'))
single=q[q.flow_valid.eq(1)&q.last_collection.eq(1)&q.URDFLOW1.notna()];dif=(single.flow-single.URDFLOW1).abs()
qc.append(dict(check='First-collection published rate rounding comparison',n=len(single),max_difference=dif.max(),result='reported rounding only; composite rates recomputed'))
pd.DataFrame(qc).to_csv(S/'audit/independent_QC.csv',index=False)
stats=[]
for gr,df in q[q.flow_valid.eq(1)].groupby('group'):
 for col in ['flow','AER','CER']:
  stats.append(dict(group=gr,variable=col,n=len(df),min=df[col].min(),p99=df[col].quantile(.99),max=df[col].max()))
pd.DataFrame(stats).to_csv(S/'audit/retained_extremes.csv',index=False)
pd.DataFrame([dict(variance_albumin=V.loc['log2alb','log2alb'],variance_creatinine=V.loc['log2ucr','log2ucr'],covariance=V.loc['log2alb','log2ucr'],variance_sum=se**2)]).to_csv(S/'audit/contrast_covariance_summary.csv',index=False)
print('Independent QC PASS;',int(ok.sum()),'flow records;',len(hashes),'frozen inputs')
