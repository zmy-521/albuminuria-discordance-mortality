"""New release setup utility, not original statistical analysis code.

Creates empty writable directories and copies configuration/aggregate inputs only.
It never fits a model, downloads records or copies hospital participant data.
"""
from pathlib import Path
import os,shutil

repo=Path(__file__).resolve().parents[1]
def root(env):
    value=os.environ.get(env)
    if not value:raise SystemExit('Set '+env+' before running setup.')
    return Path(value).resolve()
def copy_once(src,dst):
    dst.parent.mkdir(parents=True,exist_ok=True)
    if dst.exists():
        if dst.read_bytes()!=src.read_bytes():raise SystemExit('Refusing to overwrite differing file: '+str(dst))
    else:shutil.copy2(src,dst)
work=root('ODKD_WORK_ROOT');nh=root('NHANES_ROOT');fig=root('ODKD_FIGURE_ROOT')
for stage in ['raw_reconstruction','raw_reconstruction/v2','Stage0','Stage1','Stage2','Stage3_finalization','Review_supplement_urine_creatinine','Figure_v3_PlotData_export','Table1_flow_completion']:
    for folder in ['data','audit','models','internal','code','package']:(work/stage/folder).mkdir(parents=True,exist_ok=True)
copy_once(repo/'config/raw_reconstruction_original.json',work/'raw_reconstruction/code/config.json')
copy_once(repo/'config/nhanes_variable_index.csv',nh/'index/variables.csv')
for folder in ['Inputs','Data','internal','Main_candidates','Supplement_candidates','render_project/Code','release_output']:(fig/folder).mkdir(parents=True,exist_ok=True)
for p in (repo/'code/figures/source').glob('*.py'):copy_once(p,fig/'render_project/Code'/p.name)
for p in (repo/'results/figure_inputs').iterdir():copy_once(p,fig/'Inputs'/p.name)
print('Setup complete. No participant files copied; no models run.')
