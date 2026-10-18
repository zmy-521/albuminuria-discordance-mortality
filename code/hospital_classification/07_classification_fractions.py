# Release adaptation: paths only unless noted in docs/RELEASE_CHANGES.md.
import os as _release_os
def _release_path(env, relative=''):
    root = _release_os.environ.get(env)
    if not root: raise RuntimeError('Set environment variable '+env+'; see docs/REPRODUCIBILITY.md')
    return _release_os.path.join(root, relative)
from pathlib import Path
import pandas as pd,numpy as np
S=Path(_release_path('ODKD_WORK_ROOT', 'Stage3_finalization'));V=Path(_release_path('ODKD_WORK_ROOT', 'raw_reconstruction/v2'));rows=[]
classrows=[]
for center in ['Wanbei','Anyi']:
 a=pd.read_parquet(V/f'data/{center}_urine_indexed.parquet');a=a[a.base_eligible&a.D1_eligible&a.D2_eligible].copy();rows.append(dict(dataset=center+'_urine',n=len(a),overlap_mortality=0))
 for t in [30,20]:
  g=np.select([(a.ualb_mgl>=t)&(a.uacr_mgg<30),(a.ualb_mgl<t)&(a.uacr_mgg>=30),(a.ualb_mgl>=t)&(a.uacr_mgg>=30)],['G1','G2','G3'],default='G0')
  for den,mask in [('concentration_negative',a.ualb_mgl<t),('ACR_positive',a.uacr_mgg>=30)]:classrows.append(dict(source=center,threshold=t,denominator=den,n=int(mask.sum()),numerator=int((g=='G2').sum()),fraction=float((g=='G2').sum()/mask.sum()),weighted=False))
pd.DataFrame(classrows).to_csv(S/'audit/hospital_classification.csv',index=False)
