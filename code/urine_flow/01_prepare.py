# Release adaptation: paths only unless noted in docs/RELEASE_CHANGES.md.
import os as _release_os
def _release_path(env, relative=''):
    root = _release_os.environ.get(env)
    if not root: raise RuntimeError('Set environment variable '+env+'; see docs/REPRODUCIBILITY.md')
    return _release_os.path.join(root, relative)
from pathlib import Path
import pandas as pd,numpy as np,pyreadstat,json,hashlib,re,html,shutil
S=Path(_release_path('ODKD_WORK_ROOT', 'Review_supplement_urine_creatinine'));P=Path(_release_path('ODKD_WORK_ROOT', ''));H=Path(_release_path('NHANES_ROOT', 'archive'))
# Private protocol text and historical document hash are not distributed.
d=pd.read_parquet(P/'Stage2/data/cohort_stage2.parquet');assert len(d)==5905 and d.SEQN.is_unique and d.ccB.sum()==5155
sources=[];documents=[]
def load(p):
 x,m=pyreadstat.read_xport(str(p),encoding="latin1");sources.append(dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),n=len(x),fields=m.column_names_to_labels));return x
drug=load(H/'xpt/1988-2020/RXQ_DRUG.xpt');assert drug.RXDDRGID.is_unique
cn=[c for c in drug if re.match(r'RXD[DI]CN',c)];isdiur=drug[cn].fillna('').apply(lambda x:x.str.contains(r'\bDIURETICS\b',regex=True,case=False)).any(axis=1)
drug['is_diuretic']=isdiur.astype(int);drug[['RXDDRGID','RXDDRUG','is_diuretic']].to_csv(S/'audit/drug_classification_dictionary.csv',index=False)
dm=drug.set_index('RXDDRGID').is_diuretic
flows=[];rxs=[];hfs=[];lod=[]
for i,y in enumerate(range(1999,2019,2)):
 cyc=f'{y}-{y+1}';s='' if i==0 else '_'+chr(65+i);stem=['LAB16','L16_B','L16_C','ALB_CR_D','ALB_CR_E','ALB_CR_F','ALB_CR_G','ALB_CR_H','ALB_CR_I','ALB_CR_J'][i]
 a=load(H/'xpt'/cyc/(stem+'.xpt'));m=d[d.year==y][['SEQN','ualb_mgl','ucr_mgdl','group','ccB']].merge(a,on='SEQN',how='left',validate='one_to_one')
 assert np.allclose(m.ualb_mgl,m.URXUMA) and np.allclose(m.ucr_mgdl,m.URXUCR)
 lc=[c for c in a if c.endswith('LC')];al=next((c for c in lc if 'UMA' in c),None);cr=next((c for c in lc if 'UCR' in c),None)
 cp=H/'codebooks'/cyc/(stem+'.htm');text=re.sub(r'\s+',' ',html.unescape(re.sub('<[^>]*>',' ',cp.read_text(encoding='utf8')))) if cp.exists() else ''
 pos=text.lower().find('detection limits',text.lower().find('analytic notes'))
 documents.append(dict(cycle=cyc,file=stem,kind='urine_LLOD',path=str(cp),excerpt=text[pos:pos+1700] if pos>=0 else 'not found'))
 for g,q in m.groupby('group'):
  lod.append(dict(cycle=cyc,group=g,n=len(q),alb_min=q.ualb_mgl.min(),ucr_min=q.ucr_mgdl.min(),albumin_LC_field=al,creatinine_LC_field=cr,albumin_below=int(q[al].eq(1).sum()) if al else None,creatinine_below=int(q[cr].eq(1).sum()) if cr else None))
 if y<2009:continue
 f=load(H/'xpt'/cyc/('UCFLOW'+s+'.xpt'));assert f.SEQN.is_unique;f['flow_file_present']=1;flows.append(f)
 cp=H/'codebooks'/cyc/('UCFLOW'+s+'.htm');tx=re.sub(r'\s+',' ',html.unescape(re.sub('<[^>]*>',' ',cp.read_text(encoding='utf8'))));assert 'composit' in tx.lower() and 'total volume' in tx.lower()
 ix=tx.find('Collected samples');documents.append(dict(cycle=cyc,file='UCFLOW'+s,kind='pooled_specimen',path=str(cp),excerpt=tx[ix:ix+800]))
 rx=load(H/'xpt'/cyc/('RXQ_RX'+s+'.xpt'));rx['mapped_diuretic']=rx.RXDDRGID.map(dm)
 for seq,q in rx.groupby('SEQN'):
  if q.mapped_diuretic.eq(1).any():val=1.;status='identified_diuretic'
  elif q.RXDUSE.eq(2).all():val=0.;status='explicit_no_prescription'
  elif q.RXDUSE.eq(1).all() and q.mapped_diuretic.notna().all():val=0.;status='all_reported_drugs_mapped_non_diuretic'
  else:val=np.nan;status='unresolved_or_unknown_use'
  rxs.append(dict(SEQN=seq,diuretic=val,rx_status=status))
 hf=load(H/'xpt'/cyc/('MCQ'+s+'.xpt'))[['SEQN','MCQ160B']];hf['heart_failure']=hf.MCQ160B.map({1:1,2:0});hfs.append(hf[['SEQN','heart_failure']])
u=pd.concat(flows,ignore_index=True);assert u.SEQN.is_unique
q=d.merge(u,on='SEQN',how='left',validate='one_to_one').merge(pd.DataFrame(rxs),on='SEQN',how='left',validate='one_to_one').merge(pd.concat(hfs),on='SEQN',how='left',validate='one_to_one')
def flow(r):
 if r.year<2009:return (np.nan,np.nan,np.nan,0,'out_of_flow_cycles')
 if pd.isna(r.flow_file_present):return(np.nan,np.nan,np.nan,0,'no_flow_record')
 v=np.array([getattr(r,'URXVOL'+str(j)) for j in [1,2,3]],float);t=np.array([getattr(r,'URDTIME'+str(j)) for j in [1,2,3]],float)
 if np.any(np.isfinite(v)&(v<0)):return(np.nan,np.nan,np.nan,0,'negative_volume')
 pos=np.where(np.isfinite(v)&(v>0))[0]
 if not len(pos):return(np.nan,np.nan,np.nan,0,'no_positive_volume')
 last=pos[-1]+1
 # A nonmissing duration after a missing volume cannot be assumed absent from the pooled specimen.
 if any(np.isnan(v[j]) and np.isfinite(t[j]) for j in range(last,3)):return(np.nan,np.nan,np.nan,last,'ambiguous_later_collection')
 if not np.isfinite(v[:last]).all():return(np.nan,np.nan,np.nan,last,'missing_volume_inside_pool')
 if not np.isfinite(t[:last]).all():return(np.nan,np.nan,np.nan,last,'missing_time_inside_pool')
 if (t[:last]<=0).any():return(np.nan,np.nan,np.nan,last,'nonpositive_time')
 vol=v[:last].sum();dur=t[:last].sum();return(vol/dur,vol,dur,last,'valid')
out=[flow(r) for r in q.itertuples()];q[['flow','total_volume','total_time','last_collection','flow_status']]=pd.DataFrame(out,index=q.index)
q['AER']=q.ualb_mgl*q.flow;q['CER']=q.ucr_adjusted*q.flow/100;q['female']=(q.sex==2).astype(int);q['flow_valid']=q.flow_status.eq('valid').astype(int);q['AER20']=np.where(q.flow_valid.eq(1),q.AER.ge(20).astype(float),np.nan)
q.to_parquet(S/'data/review_cohort.parquet',index=False);q[['SEQN','flow','AER','CER','female','heart_failure','diuretic','flow_valid','AER20','total_volume','total_time']].to_csv(S/'data/flow_to_R.csv',index=False)
q.groupby(['cycle','flow_status'],dropna=False).size().reset_index(name='n').to_csv(S/'audit/flow_missingness.csv',index=False)
q[q.year>=2009].groupby(['group','flow_status']).agg(n=('SEQN','size'),deaths=('death','sum')).reset_index().to_csv(S/'audit/flow_group_availability.csv',index=False)
pd.DataFrame(lod).to_csv(S/'audit/LOD_audit.csv',index=False)
(S/'audit/raw_sources.json').write_text(json.dumps(sources,ensure_ascii=False,indent=2),encoding='utf8');(S/'audit/documentation_evidence.json').write_text(json.dumps(documents,ensure_ascii=False,indent=2),encoding='utf8')
print('Flow eligible',q.year.ge(2009).sum(),'valid',q.flow_valid.sum(),'groups',q[q.flow_valid.eq(1)].groupby('group').size().to_dict(),flush=True)

