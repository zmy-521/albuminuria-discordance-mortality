"""New packaging check: static syntax/privacy checks only; never imports analysis scripts."""
from pathlib import Path
import ast,re,zipfile,sys
repo=Path(__file__).resolve().parents[1]
blocked={'.rds','.rdata','.rda','.parquet','.feather','.pkl','.pickle','.joblib','.xpt','.dta','.xls','.dat','.bin','.pem','.key'}
patterns={
 'absolute_windows_path':r'(?<![A-Za-z])[A-Za-z]:[\\/]',
 'personal_home':r'/(?:Users|home)/[^/\s]+',
 'private_key':r'-----BEGIN [A-Z ]*PRIVATE KEY-----',
 'github_token':r'(?:gh[pousr]_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{30,})',
 'cloud_access_key':r'AKIA[0-9A-Z]{16}',
 'credential_assignment':r'(?i)(?:api[_-]?key|password|access[_-]?token|cookie)\s*[:=]\s*[\x22\x27][A-Za-z0-9_+/=.-]{12,}',
}
errors=[];n=0
for p in repo.rglob('*'):
 if not p.is_file() or any(x in {'.git','__pycache__'} for x in p.parts):continue
 rel=p.relative_to(repo).as_posix();n+=1
 if rel.startswith(('work/','data/')) and rel!='data/README.md':errors.append((rel,'unexpected local-data file'))
 if p.suffix.lower() in blocked:errors.append((rel,'restricted extension'))
 if p.stat().st_size>10_000_000:errors.append((rel,'unexpected size >10 MB'))
 if p.suffix=='.py':ast.parse(p.read_text(encoding='utf-8-sig'),filename=rel)
 if p.suffix=='.xlsx':
  with zipfile.ZipFile(p) as z:text='\n'.join(z.read(x).decode('utf-8','replace') for x in z.namelist() if x.endswith('.xml'))
 else:
  try:text=p.read_text(encoding='utf-8-sig')
  except UnicodeError:errors.append((rel,'unreviewed binary'));continue
 for name,pattern in patterns.items():
  if re.search(pattern,text):errors.append((rel,name))
 if re.search(r'(?<![.\d])\d{17}[\dXx](?![.\d])',text):errors.append((rel,'possible 18-character ID; inspect'))
print('Scanned',n,'files; findings',len(errors))
for rel,name in errors:print(rel,name)
sys.exit(bool(errors))
