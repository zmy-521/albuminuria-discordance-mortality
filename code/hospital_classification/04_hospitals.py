# Release adaptation: paths only unless noted in docs/RELEASE_CHANGES.md.
import os as _release_os
def _release_path(env, relative=''):
    root = _release_os.environ.get(env)
    if not root: raise RuntimeError('Set environment variable '+env+'; see docs/REPRODUCIBILITY.md')
    return _release_os.path.join(root, relative)
from common import *
import sys
center=sys.argv[1];wb=center=='Wanbei'
raw=pd.read_parquet(R/f'data/{"wanbei" if wb else "anyi"}_raw_records.parquet')
if wb:
 rename={'患者编码':'pid','患者身份证号':'stable_id','患者性别':'sex_raw','患者年龄':'age_raw','项目名称':'project','项目编码':'code','结果':'result','单位':'unit','标本类型':'specimen','条码号':'sample','检验日期':'date_raw','采集时间':'collection','审核时间':'review','诊断':'diagnosis'}
else:
 rename={'登记号':'pid','病案号':'stable_id','性别':'sex_raw','年龄':'age_raw','检测项目':'project','结果':'result','单位':'unit','标本类型':'specimen','检验号':'sample','审核日期':'date_raw','接收日期':'collection','诊断':'diagnosis'}
d=raw.rename(columns=rename);d['pid']=norm(d.pid);d['stable_id']=norm(d.stable_id).str.upper();d['sample']=norm(d['sample']);d['source_row']=d.source_row.astype(str)
d['raw_index']=d.source_file+':'+d.source_row
if not wb:d['review']=d.date_raw;d['code']='not_provided'
# Identity linking uses source stable identifiers only, never names or labels.
valid=d.stable_id.str.fullmatch(r'\d{17}[\dX]') if wb else d.stable_id.str.fullmatch(r'\d{4,}')
links=d.loc[valid,['pid','stable_id']].drop_duplicates();conf=links.groupby('pid').stable_id.nunique();bad=set(conf[conf.gt(1)].index)
link=links[~links.pid.isin(bad)].drop_duplicates('pid').set_index('pid').stable_id
d['stable_link']=d.pid.map(link);d['identity_usable']=d.stable_link.notna()&~d.pid.isin(bad)
d['person_key']=[pseud(center+':ID:'+a) if isinstance(a,str) else pseud(center+':unresolved:'+b) for a,b in zip(d.stable_link,d.pid)]
d['sex']=d.sex_raw.map({'男':1,'女':2});sx=d.groupby('person_key').sex.nunique();badsex=set(sx[sx.gt(1)].index);d.loc[d.person_key.isin(badsex),'identity_usable']=False
d['age']=pd.to_numeric(d.age_raw.str.extract(r'^(\d+(?:\.\d+)?)岁$')[0],errors='coerce');d['diabetes']=d.diagnosis.str.contains('糖尿病',na=False)&~d.diagnosis.str.contains('妊娠|孕',na=False)
col=parse_dt(d.collection);rev=parse_dt(d.review);date=parse_dt(d.date_raw)
# Use specimen collection/receipt if present and not after review; retain fallback.
good=col.notna()&col.dt.year.between(1990,2030)&(rev.isna()|col.le(rev))
d['time']=col.where(good,rev.where(rev.notna(),date));d['time_source']=np.where(good,'collection' if wb else 'receipt','review_or_test_date');d['day']=d.time.dt.normalize()
d['num']=pd.to_numeric(d.result.str.strip(),errors='coerce');d['result_qualified']=d.result.str.contains(r'[<>≤≥]',regex=True,na=False);d.loc[d.result_qualified|d.num.lt(0),'num']=np.nan
d['unit_norm']=d.unit.str.replace('μ','u').str.replace('µ','u').str.replace(' ','').str.lower();d['analyte']='';d['value']=np.nan
units={
 'scr_mgdl':{'umol/l':1/88.4,'mg/dl':1},'cys_mgl':{'mg/l':1},'ualb_mgl':{'mg/l':1,'ug/ml':1},'ucr_mgdl':{'mmol/l':11.312,'mg/dl':1,'umol/l':.011312},'uacr_direct':{'mg/g':1,'mg/gcr':1,'mg/mmol.cr':8.84,'mg/mmol':8.84},
 'HGB':{'g/l':1,'g/dl':10},'UA':{'umol/l':1,'mg/dl':59.48},'HbA1c':{'%':1},'ALB':{'g/l':1},'TP':{'g/l':1},'GLB':{'g/l':1},'AG':{'':1,'无':1},'BUN':{'mmol/l':1},'GLU':{'mmol/l':1},'TG':{'mmol/l':1},'HDL':{'mmol/l':1},'LDL':{'mmol/l':1},'CHOL':{'mmol/l':1},'TBIL':{'umol/l':1},'PAB':{'mg/l':1},'ALP':{'u/l':1},'ALT':{'u/l':1},'AST':{'u/l':1},'GGT':{'u/l':1}}
for v in ['WBC','NEUT','LYMPH','MONO','BASO','EO','PLT']:units[v]={'x10^9/l':1,'10^9/l':1,'10*9/l':1,'×10^9/l':1}
serum=d.specimen.isin(['血清','血','血浆']);blood=d.specimen.isin(['全血','抗凝血']);urine=d.specimen.isin(['随机尿','尿','尿液'])&~d.project.str.contains('24|小时',na=False)
if wb:
 amap={'CREA':'scr_mgdl','CREA1':'scr_mgdl','Cys-C':'cys_mgl','UMALB':'ualb_mgl','UCR':'ucr_mgdl','ACR':'uacr_direct','HGB':'HGB','HB':'HGB','UA':'UA','HbA1c':'HbA1c','ALB':'ALB','TP':'TP','GLB':'GLB','A/G':'AG','BUN':'BUN','GLU':'GLU','TG':'TG','HDL':'HDL','LDL-C':'LDL','CHOL':'CHOL','TBIL':'TBIL','PAB':'PAB','ALP':'ALP','ALT':'ALT','AST':'AST','GGT':'GGT','WBC':'WBC','NEUT#':'NEUT','LYMPH#':'LYMPH','MONO#':'MONO','BASO#':'BASO','EO#':'EO','PLT':'PLT'}
 mapped=d.code.map(amap)
else:
 names=d.project.str.replace(r'^[△*]+','',regex=True)
 amap={'肌酐':'scr_mgdl','胱抑素C':'cys_mgl','尿微量白蛋白':'ualb_mgl','尿微量白蛋白/肌酐':'uacr_direct','Alb/肌酐':'uacr_direct','血红蛋白':'HGB','尿酸':'UA','糖化血红蛋白':'HbA1c','白蛋白':'ALB','总蛋白':'TP','球蛋白':'GLB','白球比':'AG','尿素':'BUN','葡萄糖':'GLU','甘油三酯':'TG','高密度脂蛋白胆固醇':'HDL','低密度脂蛋白胆固醇':'LDL','总胆固醇':'CHOL','总胆红素':'TBIL','前白蛋白':'PAB','碱性磷酸酶':'ALP','丙氨酸氨基转移酶':'ALT','天门冬氨酸氨基转移酶':'AST','谷氨酰转移酶':'GGT','白细胞计数':'WBC','中性粒细胞绝对数':'NEUT','淋巴细胞绝对数':'LYMPH','单核细胞绝对数':'MONO','嗜碱性粒细胞绝对数':'BASO','嗜酸性粒细胞绝对数':'EO','血小板计数':'PLT'}
 amap['Alb']='ualb_mgl';amap['γ-谷氨酰转移酶']='GGT'
 mapped=names.map(amap);mapped.loc[urine&names.eq('肌酐')]='ucr_mgdl'
for a,conversion in units.items():
 group=urine if a in ['ualb_mgl','ucr_mgdl','uacr_direct'] else blood if a in ['HGB','HbA1c','WBC','NEUT','LYMPH','MONO','BASO','EO','PLT'] else serum
 hit=mapped.eq(a)&group;d.loc[hit,'analyte']=a;d.loc[hit,'value']=d.loc[hit,'num']*d.loc[hit,'unit_norm'].map(conversion)
d.loc[d.analyte.isin(['scr_mgdl','cys_mgl','ucr_mgdl'])&d.value.le(0),'value']=np.nan
d['unit_supported']=d.analyte.ne('')&d.value.notna()
# Preserve raw record trace, including unmapped units and qualified results.
d.drop(columns=['stable_link','stable_id'],errors='ignore').to_parquet(R/f'data/{center}_record_audit.parquet',index=False)
d[['person_key','pid','stable_link','identity_usable']].drop_duplicates().to_parquet(R/f'internal/{center}_identity_bridge.parquet',index=False)
prof=d[d.analyte.ne('')].groupby(['project','code','unit','specimen','analyte']).agg(records=('raw_index','size'),numeric_usable=('value','count')).reset_index();prof.to_csv(R/f'audit/{center}_mapping.csv',index=False,encoding='utf-8-sig')
# First serum-creatinine assessment, independent of urine and CysC completeness/eGFR result.
s=d[d.analyte.eq('scr_mgdl')&d.value.notna()&d.age.ge(18)&d.diabetes&d.day.notna()].sort_values(['person_key','time','sample','raw_index'])
idx=s.drop_duplicates('person_key')[['person_key','identity_usable','time','day','time_source','age','sex','raw_index','pid']].rename(columns={'time':'index_time','raw_index':'index_source'})
idx['diabetes']=True;idx['center']=center
# Enforce same-day matching. Never borrow future visits to fill missing values.
t=d[d.analyte.ne('')&d['sample'].ne('')].merge(idx[['person_key','day']],on=['person_key','day'],how='inner',validate='many_to_one')
u=t[t.analyte.isin(['ualb_mgl','ucr_mgdl','uacr_direct'])].sort_values(['person_key','time','sample','raw_index'])
first=u.drop_duplicates('person_key')[['person_key','sample']];u=u.merge(first,on=['person_key','sample'],validate='many_to_one')
v=t[~t.analyte.isin(['ualb_mgl','ucr_mgdl','uacr_direct'])].sort_values(['person_key','analyte','time','sample','raw_index'])
first=v.drop_duplicates(['person_key','analyte'])[['person_key','analyte','sample']];v=v.merge(first,on=['person_key','analyte','sample'],validate='many_to_one')
sel=pd.concat([u,v]);ag=sel.groupby(['person_key','analyte']).agg(n_values=('value','nunique'),value=('value','first'),raw_indices=('raw_index',lambda x:'|'.join(sorted(set(x)))),time=('time','min'),sample=('sample','first')).reset_index();ag.loc[ag.n_values.gt(1),'value']=np.nan
ag.to_parquet(R/f'data/{center}_selected_measurement_trace.parquet',index=False)
wide=ag.pivot(index='person_key',columns='analyte',values='value');out=idx.merge(wide,on='person_key',how='left',validate='one_to_one')
for a in units:
 if a not in out:out[a]=np.nan
out['uacr_calculated']=out.ualb_mgl/out.ucr_mgdl*100;out['uacr_mgg']=out.uacr_calculated.combine_first(out.uacr_direct);out['uacr_source']=np.select([out.uacr_calculated.notna(),out.uacr_direct.notna()],['same_specimen_calculated','direct_report'],default='missing')
out['ratio_threshold_disagreement']=out.uacr_calculated.notna()&out.uacr_direct.notna()&out.uacr_calculated.ge(30).ne(out.uacr_direct.ge(30))
out['AG']=out.AG.combine_first(out.ALB/out.GLB.where(out.GLB.gt(0)));out['NLR']=out.NEUT/out.LYMPH.where(out.LYMPH.gt(0));out['PLR']=out.PLT/out.LYMPH.where(out.LYMPH.gt(0));out['SII']=out.NEUT*out.PLT/out.LYMPH.where(out.LYMPH.gt(0))
# Diagnostic text is not a validated negative history. Do not manufacture HTN/CVD/MicroVCs zeros.
out['HTN']=np.nan;out['CVD']=np.nan;out['MicroVCs']=np.nan;out['BMI']=np.nan
out=labels(out);out.drop(columns='pid').to_parquet(R/f'data/{center}_canonical.parquet',index=False)
paired=out[out.base_eligible&out.uacr_calculated.notna()&out.uacr_direct.notna()]
extra={'raw_rows':len(raw),'raw_patient_codes':d.pid.nunique(),'stable_ids':links.stable_id.nunique(),'ambiguous_patient_codes':len(bad),'sex_conflict_keys':len(badsex),'selected_conflicting_analytes':int(ag.n_values.gt(1).sum()),'direct_recalculation':{'n':len(paired),'median_absolute_difference':float((paired.uacr_direct-paired.uacr_calculated).abs().median()),'threshold30_disagreements':int(paired.ratio_threshold_disagreement.sum())},'index_time_min':str(out.index_time.min()),'index_time_max':str(out.index_time.max()),'index_time_sources':out.time_source.value_counts().to_dict(),'max_time_difference_hours':float((ag.time-ag.person_key.map(out.set_index('person_key').index_time)).dt.total_seconds().abs().max()/3600)}
save_summary(out.drop(columns='pid'),center,extra)
