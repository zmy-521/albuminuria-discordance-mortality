# Release adaptation: paths only unless noted in docs/RELEASE_CHANGES.md.
import os as _release_os
def _release_path(env, relative=''):
    root = _release_os.environ.get(env)
    if not root: raise RuntimeError('Set environment variable '+env+'; see docs/REPRODUCIBILITY.md')
    return _release_os.path.join(root, relative)
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,xlrd
R=Path(_release_path('ODKD_WORK_ROOT', 'raw_reconstruction'));W=Path(_release_path('WANBEI_RAW_ROOT', ''))
out=[];parts=[]
for p in W.glob('*.xls'):
 if p.name.startswith('凝血'):continue
 b=xlrd.open_workbook(p,on_demand=True);s=b.sheet_by_index(0);hi=next(i for i in range(5) if '项目名称' in s.row_values(i));headers=[str(x).strip() for x in s.row_values(hi)]
 d=pd.DataFrame([s.row_values(i) for i in range(hi+1,s.nrows)],columns=headers);d['source_file']=p.name;d['source_row']=np.arange(hi+2,s.nrows+1)
 # Preserve identifiers locally for audited links, never print identifiers.
 if '患者编码' in d:
  keep=['患者类型','患者编码','患者性别','患者年龄','检验日期','标本号','项目编码','项目名称','结果','单位','标本类型','诊断','条码号','患者出生日期','患者身份证号','采集时间','审核时间','source_file','source_row']
  d=d[keep].astype(str);parts.append(d)
  out.append({'file':p.name,'rows':len(d),'patient_codes':d['患者编码'].nunique(),'projects':d['项目名称'].nunique(),'date_min':d['检验日期'].min(),'date_max':d['检验日期'].max()})
 else:
  d.drop(columns=['病人姓名','医生','床号'],errors='ignore').astype(str).to_parquet(R/'data/wanbei_diagnosis_raw.parquet',index=False);out.append({'file':p.name,'rows':len(d),'admission_keys':d['住院号'].nunique()})
 b.release_resources();print(p.name,len(d),flush=True)
d=pd.concat(parts,ignore_index=True);d.to_parquet(R/'data/wanbei_raw_records.parquet',index=False)
d.groupby(['项目名称','项目编码','单位','标本类型'],dropna=False).size().reset_index(name='records').to_csv(R/'audit/wanbei_project_units.csv',index=False,encoding='utf-8-sig')
(R/'audit/wanbei_source_profile.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
