# Release adaptation: paths only unless noted in docs/RELEASE_CHANGES.md.
import os as _release_os
def _release_path(env, relative=''):
    root = _release_os.environ.get(env)
    if not root: raise RuntimeError('Set environment variable '+env+'; see docs/REPRODUCIBILITY.md')
    return _release_os.path.join(root, relative)
"""Export frozen plot fields only; no model estimation or new survival curves."""
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,shutil,datetime
R=Path(_release_path('ODKD_WORK_ROOT', '')); V=Path(_release_path('ODKD_WORK_ROOT', 'raw_reconstruction/v2')); O=R/'Figure_v3_PlotData_export'; P=O/'package';I=O/'internal';P.mkdir(exist_ok=True);I.mkdir(exist_ok=True)
sources=[R/'Stage1/data/cohort.parquet',R/'Stage1/data/cox_models.rds',R/'Stage1/data/survey_design.rds',V/'data/Wanbei_urine_indexed.parquet',R/'Review_supplement_urine_creatinine/data/review_cohort.parquet',R/'Review_supplement_urine_creatinine/data/review_survey_design.rds',R/'Review_supplement_urine_creatinine/models/supplemental_models.rds',R/'Review_supplement_urine_creatinine/audit/flow_weighted_descriptives.csv',R/'Review_supplement_urine_creatinine/02_supplement_statistical_results.xlsx',R/'Stage3_finalization/19_final_clinical_results.xlsx',R/'Stage1/code/01_prepare.py',R/'Stage1/code/02_survival.R',R/'Stage1/code/03_absolute_risk.R',V/'code/01_reconstruct.py',V/'code/core.py',V.parent/'code/04_hospitals.py',R/'Review_supplement_urine_creatinine/code/01_prepare.py',R/'Review_supplement_urine_creatinine/code/02_analysis.R']
for n in ['component_coefficients.csv','component_design_covariance.csv']:sources.append(R/'Review_supplement_urine_creatinine/audit'/n)
sources+=list(pd.read_csv(I/'RDS_inventory.csv')['file'].map(Path))
sources=[p for p in dict.fromkeys(sources) if p.suffix.lower() not in ['.py','.r']];hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources};(I/'input_hashes.json').write_text(json.dumps(hashes,ensure_ascii=False,indent=2),encoding='utf-8')
nh=pd.read_parquet(sources[0]); wb=pd.read_parquet(V/'data/Wanbei_urine_indexed.parquet');fl=pd.read_parquet(R/'Review_supplement_urine_creatinine/data/review_cohort.parquet')
mb=pd.read_csv(I/'actual_ModelB_members.csv');a2=pd.read_csv(I/'actual_A2_members.csv')
assert len(nh)==5905 and nh.SEQN.is_unique and mb.SEQN.is_unique and set(mb.SEQN).issubset(set(nh.SEQN))
nh=nh.sample(frac=1,random_state=20261014).reset_index(drop=True);nh['anonymous_row_id']=[f'N{i:05d}' for i in range(1,len(nh)+1)]
nh['model_B_actual']=nh.SEQN.isin(mb.SEQN).astype(int);assert nh.model_B_actual.sum()==5155 and nh.loc[nh.model_B_actual.eq(1),'death'].sum()==1115
assert np.array_equal(nh.model_B_actual,nh.ccB.astype(int))
nm=nh.merge(mb,on='SEQN',suffixes=('','_model'),validate='one_to_one');assert (nm.group==nm.group_model).all() and (nm.cycle==nm.cycle_model).all() and np.allclose(nm.weight,nm.weight_model,rtol=1e-13)
nh['A2_ModelB_actual']=nh.SEQN.isin(a2.SEQN).astype(int);assert (nh.A2_ModelB_actual.eq(1)==(nh.model_B_actual.eq(1)&nh.acr_adjusted.ge(30)&nh.acr_adjusted.lt(300))).all()
assert nh.groupby('group').size().to_dict()==dict(G0=3953,G1=314,G2=251,G3=1387)
assert nh[nh.A2_ModelB_actual.eq(1)].group.value_counts().to_dict()==dict(G3=898,G2=222)
assert np.allclose(nh.weight,np.where(nh.year<2003,nh.WTMEC4YR*.2,nh.WTMEC2YR*.1))
x=nh.ucr_mgdl;adj=np.where(nh.year<2007,np.select([x<75,x<250],[(1.02*np.sqrt(x)-.36)**2,(1.05*np.sqrt(x)-.74)**2],default=(1.01*np.sqrt(x)-.10)**2),x);assert np.allclose(adj,nh.ucr_adjusted,rtol=1e-13)
cols={'ualb_mgl':'urine_albumin_mg_L','ucr_adjusted':'urine_creatinine_mg_dL_CDC','acr_adjusted':'UACR_mg_g_CDC','group':'original_group','weight':'original_survey_weight','cycle':'survey_cycle'}
def save(d,p):d.to_csv(P/p,index=False,encoding='utf-8-sig',float_format='%.17g',na_rep='')
save(nh[['anonymous_row_id']+list(cols)+['model_B_actual','A2_ModelB_actual']].rename(columns=cols),'01_NHANES_5905_plot.csv')
nh[['anonymous_row_id','SEQN']].to_csv(I/'NHANES_local_key_mapping.csv',index=False)
wb=wb[wb.base_eligible&wb.D1_eligible&wb.D2_eligible].copy().sample(frac=1,random_state=20261014).reset_index(drop=True)
assert len(wb)==2261 and wb.person_key.is_unique
wb['anonymous_row_id']=[f'W{i:05d}' for i in range(1,len(wb)+1)]
wb['original_group']=wb.urine_category.map(dict(both_negative='G0',concentration_only='G1',ACR_only='G2',both_positive='G3'))
assert wb.original_group.value_counts().to_dict()==dict(G0=1533,G3=475,G2=152,G1=101)
assert wb.uacr_source.eq('same_specimen_calculated').all()
save(wb[['anonymous_row_id','ualb_mgl','ucr_mgdl','uacr_mgg','original_group']].rename(columns={'ualb_mgl':'urine_albumin_mg_L','ucr_mgdl':'urine_creatinine_mg_dL','uacr_mgg':'UACR_mg_g'}),'02_Wanbei_2261_plot.csv')
wb[['anonymous_row_id','person_key']].to_csv(I/'Wanbei_local_key_mapping.csv',index=False)
fl=fl[fl.year>=2009].copy();fl['anonymous_row_id']=fl.SEQN.map(nh.set_index('SEQN').anonymous_row_id);assert fl.anonymous_row_id.notna().all() and fl.SEQN.is_unique
assert len(fl)==3539 and fl.flow_valid.sum()==3261
assert fl[fl.flow_valid.eq(1)].group.value_counts().to_dict()==dict(G0=2276,G3=686,G1=176,G2=123)
check=fl.merge(nh,on='SEQN',suffixes=('','_nh'),validate='one_to_one')
for c in ['ualb_mgl','ucr_adjusted','acr_adjusted','weight']:assert np.allclose(check[c],check[c+'_nh'],rtol=1e-13)
assert (check.group==check.group_nh).all()
extra={'flow':'onsite_urine_flow_mL_min','CER':'short_CER_mg_min','AER':'short_AER_ug_min','flow_valid':'flow_evaluable','flow_status':'existing_flow_status','total_volume':'pooled_volume_mL','total_time':'pooled_time_min'}
save(fl[['anonymous_row_id']+list(cols)+list(extra)].rename(columns=cols|extra).sort_values('anonymous_row_id'),'03_Flow_plot.csv')
checks=[]
def ck(name,value,limit=1e-10):checks.append(dict(check=name,value=float(value),tolerance=limit,pass_check=bool(value<=limit)));assert value<=limit,name
ck('NHANES UACR formula absolute max error',abs(nh.acr_adjusted-100*nh.ualb_mgl/nh.ucr_adjusted).max())
ck('Wanbei UACR formula absolute max error',abs(wb.uacr_mgg-100*wb.ualb_mgl/wb.ucr_mgdl).max())
ok=fl.flow_valid.eq(1)
ck('Flow = frozen pooled volume/time',abs(fl.loc[ok,'flow']-fl.loc[ok,'total_volume']/fl.loc[ok,'total_time']).max())
ck('AER formula',abs(fl.loc[ok,'AER']-fl.loc[ok,'ualb_mgl']*fl.loc[ok,'flow']).max())
ck('CER formula',abs(fl.loc[ok,'CER']-fl.loc[ok,'ucr_adjusted']*fl.loc[ok,'flow']/100).max())
frozen=pd.read_csv(R/'Review_supplement_urine_creatinine/audit/flow_weighted_descriptives.csv');quant=[]
for g,d in fl[ok].groupby('group'):
 for var in ['flow','AER','CER']:
  z=frozen[(frozen.population==g)&frozen.variable.eq(var)].iloc[0];v=d[var].to_numpy();w=d.weight.to_numpy();idx=np.argsort(v);q=v[idx][np.searchsorted(np.cumsum(w[idx])/sum(w),[.25,.5,.75])]
  mean=np.average(v,weights=w);ck(f'{g}/{var} weighted mean',abs(mean-z['mean']))
  for key,val in zip(['q1','median','q3'],q):ck(f'{g}/{var} weighted {key}',abs(val-z[key]))
  quant.append([g,var,len(d),mean,*q])
pd.DataFrame(quant,columns=['group','variable','n','mean','q1','median','q3']).to_csv(I/'flow_QC_summary.csv',index=False)
pd.DataFrame(checks).to_csv(I/'formula_and_summary_checks.csv',index=False)
profiles=[]
for label,d,m in [('NHANES',nh,{'albumin':'ualb_mgl','creatinine':'ucr_adjusted','UACR':'acr_adjusted'}),('Wanbei',wb,{'albumin':'ualb_mgl','creatinine':'ucr_mgdl','UACR':'uacr_mgg'}),('Flow all3539',fl,{'flow':'flow','AER':'AER','CER':'CER'})]:
 for name,c in m.items():profiles.append([label,name,len(d),int(d[c].isna().sum()),int(d[c].le(0).sum()),float(d[c].min()),float(d[c].max())])
pd.DataFrame(profiles,columns=['cohort','variable','n','missing','nonpositive','min','max']).to_csv(I/'range_QC.csv',index=False)
# Copy existing aggregate coefficients/covariance verbatim; no recomputation or surface generation.
for orig,new in [('component_coefficients.csv','07_Frozen_component_coefficients.csv'),('component_design_covariance.csv','08_Frozen_component_covariance.csv')]:shutil.copy2(R/'Review_supplement_urine_creatinine/audit'/orig,P/new)
for f in (O/'code').glob('*'):
 if f.is_file(): (P/'scripts').mkdir(exist_ok=True);shutil.copy2(f,P/'scripts'/f.name)
print('Export complete: NHANES5905; Wanbei2261; flow3539/3261; actual ModelB5155/1115; A2 1120. No fits/curves created.')
