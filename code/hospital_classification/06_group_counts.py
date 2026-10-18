# Release adaptation: paths only unless noted in docs/RELEASE_CHANGES.md.
import os as _release_os
def _release_path(env, relative=''):
    root = _release_os.environ.get(env)
    if not root: raise RuntimeError('Set environment variable '+env+'; see docs/REPRODUCIBILITY.md')
    return _release_os.path.join(root, relative)
from pathlib import Path
import pandas as pd
S=Path(_release_path('ODKD_WORK_ROOT', 'Stage3_finalization'));V=Path(_release_path('ODKD_WORK_ROOT', 'raw_reconstruction/v2'))
counts=[]
for center in ['Wanbei','Anyi']:
 q=pd.read_parquet(V/f'data/{center}_urine_indexed.parquet');q=q[q.base_eligible&q.D1_eligible&q.D2_eligible]
 for t in [30,20]:
  for a,b in [(0,0),(1,0),(0,1),(1,1)]:counts.append(dict(source=center,threshold=t,group={(0,0):'G0',(1,0):'G1',(0,1):'G2',(1,1):'G3'}[(a,b)],n=int(((q.ualb_mgl.ge(t)==bool(a))&(q.uacr_mgg.ge(30)==bool(b))).sum())))
pd.DataFrame(counts).to_csv(S/'audit/hospital_group_counts.csv',index=False)
