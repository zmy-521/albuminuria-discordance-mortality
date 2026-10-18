# Release adaptation: paths only unless noted in docs/RELEASE_CHANGES.md.
import os as _release_os
def _release_path(env, relative=''):
    root = _release_os.environ.get(env)
    if not root: raise RuntimeError('Set environment variable '+env+'; see docs/REPRODUCIBILITY.md')
    return _release_os.path.join(root, relative)
from core import *
import sys,gc
center=sys.argv[1]
if center=='NHANES':
 d=pd.read_parquet(R/'data/NHANES_canonical.parquet')
 # These are newly reconstructed v1 raw-XPT variables, not historical paper data.
 d=calculate(d);d['assessment_id']=d.person_key;d['index_rule']='single_MEC_exam';d['window_days']=np.nan
 valid=d.identity_usable&d.diabetes&d.age.ge(18)&d.egfr_cr.notna()
 u=d.uacr_mgg.notna();a=d.ualb_mgl.notna();c=d.cys_mgl.notna()
 masks={'urine':valid&a&u,'renal':valid&c,'D1':valid&a,'D2':valid&u,'D4':valid&u&c,'complete':valid&a&u&c}
 for axis,mask in masks.items():d[mask].to_parquet(V/f'data/{center}_{axis}_indexed.parquet',index=False)
 dump({'source_rows':len(d),'adult_diabetes':int((d.identity_usable&d.diabetes&d.age.ge(18)).sum()),'eligible_assessments':{k:int(v.sum()) for k,v in masks.items()},'independent_axes':True,'window':'one_MEC_exam; no exact dates inferred'},V/f'audit/{center}_construction.json')
 print(center,{k:int(v.sum()) for k,v in masks.items()});sys.exit()

d=pd.read_parquet(R/f'data/{center}_record_audit.parquet')
base=d[d.identity_usable&d.day.notna()].copy()
# Neither pid nor the LIS serial is a verified encounter identifier. Same code + day only.
keys=['person_key','pid','day']
status=base.groupby(keys,dropna=False).agg(diabetes=('diabetes','max'),age=('age','first'),sex=('sex','first')).reset_index()
status=status[status.diabetes&status.age.ge(18)&status.sex.isin([1,2])]
raw=base[base.analyte.ne('')&base['sample'].ne('')].merge(status[keys],on=keys,validate='many_to_one')
raw=raw.sort_values(keys+['time','sample','raw_index'])
# Patient-code / calendar-day / specimen / analyte is the deterministic measurement unit.
sk=keys+['sample']
ag=raw.groupby(sk+['analyte'],dropna=False).agg(n_values=('value','nunique'),value=('value','first'),time=('time','min'),raw_indices=('raw_index',lambda s:'|'.join(sorted(set(s)))),raw_rows=('raw_index','size')).reset_index()
ag.loc[ag.n_values.ne(1),'value']=np.nan
ag.to_parquet(V/f'data/{center}_specimen_measurement_trace.parquet',index=False)
spec=ag.pivot(index=sk,columns='analyte',values='value').reset_index()
for col in ANALYTES:
 if col not in spec:spec[col]=np.nan
times=ag.groupby(sk).time.min().rename('index_time').reset_index();spec=spec.merge(times,on=sk,validate='one_to_one').merge(status,on=keys,validate='many_to_one')
spec['assessment_id']=[hashlib.sha256((center+'|'+str(p)+'|'+str(i)+'|'+str(day)+'|'+str(s)).encode()).hexdigest()[:24] for p,i,day,s in spec[sk].itertuples(index=False,name=None)]
spec['uacr_calc']=100*spec.ualb_mgl/spec.ucr_mgdl.where(spec.ucr_mgdl.gt(0));spec['uacr_available']=spec.uacr_calc.notna()|spec.uacr_direct.notna()

def attach(anchor,analyte):
 # Prefer the exact specimen, then closest within the same code and calendar day.
 z=ag[ag.analyte.eq(analyte)&ag.value.notna()][keys+['sample','value','time','raw_indices']].rename(columns={'sample':'match_sample','time':'match_time','value':'match_value'})
 x=anchor[['assessment_id']+keys+['sample','index_time']].merge(z,on=keys,how='inner')
 x['same_specimen']=x['sample'].eq(x.match_sample);x['delta_h']=(x.match_time-x.index_time).dt.total_seconds()/3600
 x=x.sort_values(['assessment_id','same_specimen','delta_abs','match_time','match_sample']) if 'delta_abs' in x else x.assign(delta_abs=x.delta_h.abs()).sort_values(['assessment_id','same_specimen','delta_abs','match_time','match_sample'],ascending=[True,False,True,True,True])
 return x.drop_duplicates('assessment_id').set_index('assessment_id')

ur=spec[spec.ualb_mgl.notna()|spec.uacr_available].copy()
ren=spec[spec.scr_mgdl.notna()].copy()
for cand,required in [(ur,['scr_mgdl','cys_mgl']),(ren,['cys_mgl'])]:
 for a in required:
  m=attach(cand,a);cand[a]=cand.assessment_id.map(m.match_value);cand[a+'_delta_h']=cand.assessment_id.map(m.delta_h);cand[a+'_matched_sample']=cand.assessment_id.map(m.match_sample)
ur=ur[ur.scr_mgdl.notna()];ren=ren[ren.cys_mgl.notna()]
axes={'urine':ur[ur.ualb_mgl.notna()&ur.uacr_available], 'D1':ur[ur.ualb_mgl.notna()], 'D2':ur[ur.uacr_available], 'D4':ur[ur.uacr_available&ur.cys_mgl.notna()], 'complete':ur[ur.ualb_mgl.notna()&ur.uacr_available&ur.cys_mgl.notna()], 'renal':ren}
counts={};chosen={}
for axis,cand in axes.items():
 # Required measurements, not eGFR threshold or any positive result, select the index.
 ix=cand.sort_values(['person_key','index_time','sample','assessment_id']).drop_duplicates('person_key').copy()
 counts[axis]={'candidate_specimens':len(cand),'people_before_egfr_restriction':len(ix)}
 chosen[axis]=ix
 allc=cand[['assessment_id','person_key','index_time','scr_mgdl','cys_mgl','ualb_mgl','uacr_calc']].copy();allc['selected']=allc.assessment_id.isin(ix.assessment_id)
 allc.to_parquet(V/f'data/{center}_{axis}_candidate_selection.parquet',index=False)
del d,base,raw;gc.collect()

trace=[]
for axis,out in chosen.items():
 out=out.copy()
 # Urine values remain anchored to the one urine specimen. Renal axis has no urine borrowing.
 for a in ANALYTES:
  if a in URINE:continue
  m=attach(out,a)
  out[a]=out.assessment_id.map(m.match_value)
  out[a+'_delta_h']=out.assessment_id.map(m.delta_h)
  z=m.reset_index()[['assessment_id','match_sample','match_time','raw_indices','delta_h','match_value']];z['axis']=axis;z['analyte']=a;trace.append(z)
 if axis=='renal':
  for a in URINE:out[a]=np.nan
 out['identity_usable']=True;out['center']=center;out['index_rule']='first_required_measurements_same_pid_day';out['window_days']=0
 out=calculate(out)
 counts[axis].update({'base_ge60':int(out.base_eligible.sum()),'extreme_scr':int(out.qc_extreme_scr.sum()),'extreme_cys':int(out.qc_extreme_cys.sum()),'ratio_disagreements':int(out.ratio_threshold_disagreement.sum())})
 out.drop(columns=['pid']).to_parquet(V/f'data/{center}_{axis}_indexed.parquet',index=False)
 print(center,axis,counts[axis],flush=True)
pd.concat(trace,ignore_index=True).to_parquet(V/f'data/{center}_axis_selected_trace.parquet',index=False)
dump({'eligible_adult_dm_people':status.person_key.nunique(),'same_day_code_assessments':len(status),'conflicting_specimen_analytes':int(ag.n_values.gt(1).sum()),'extra_raw_duplicate_measurements':int(ag.raw_rows.sum()-len(ag)),'axis_counts':counts,'window_rule':'No verified encounter key; retain same patient code and same calendar day for nominal 0/3/7-day analyses; no cross-day expansion'},V/f'audit/{center}_construction.json')
