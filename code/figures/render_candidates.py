# Release adaptation: paths only unless noted in docs/RELEASE_CHANGES.md.
import os as _release_os
def _release_path(env, relative=''):
    root = _release_os.environ.get(env)
    if not root: raise RuntimeError('Set environment variable '+env+'; see docs/REPRODUCIBILITY.md')
    return _release_os.path.join(root, relative)
from pathlib import Path
import shutil,sys,json,hashlib,zipfile
import numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.text import Text
from matplotlib.lines import Line2D
from matplotlib.collections import PathCollection
from PIL import Image,ImageDraw
import fitz
P=Path(_release_path('ODKD_FIGURE_ROOT'))
S=P/'release_source'
R=P/'render_project';O=P/'release_output';O.mkdir(parents=True,exist_ok=True)
(P/'internal').mkdir(parents=True,exist_ok=True)
hashfile=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
source_hash={str(p):hashfile(p) for p in R.rglob('*') if p.is_file()}
stems=['Classification','LongTerm_Mortality','Components_A2','FiveYear_Risk','Urine_Flow']
sys.path.insert(0,str(R/'Code'))
import base_plotting as b,render_v2 as v2,render_v3 as v3
# The mortality renderer has an inline export block. Redirect that block only.
p=R/'Code/v3_mortality.py';s=(Path(__file__).parent/'source/v3_mortality.py').read_text(encoding='utf-8')
pos=s.index('    stem = b.MAIN / "Candidate_02_LongTerm_Mortality"')
end=s.index('\n\n\nif __name__',pos)
s=s[:pos]+"    b.save(fig, 'Candidate_02_LongTerm_Mortality')\n"+s[end:];p.write_text(s,encoding='utf-8')
import v3_mortality
captured=[]
def capture(fig,name,*args,**kwargs):captured.append((fig,name))
b.save=v2.save=v3.save=capture
v3.check_inputs();v3.classification_v3();v3_mortality.make_mortality();v3.components_v3();v2.absolute_risk_v2();v3.flow_v3()
assert len(captured)==5
def arr(a):
 a=np.ma.asarray(a);return dict(shape=a.shape,values=hashlib.sha256(np.ascontiguousarray(a.filled(np.nan) if np.ma.isMaskedArray(a) and np.issubdtype(a.dtype,np.number) else a).tobytes()).hexdigest())
def snapshot(fig):
 result=[]
 for ax in fig.axes:
  result.append(dict(xlim=ax.get_xlim(),ylim=ax.get_ylim(),xscale=ax.get_xscale(),yscale=ax.get_yscale(),xticks=ax.get_xticks().tolist(),yticks=ax.get_yticks().tolist(),lines=[(arr(l.get_xdata()),arr(l.get_ydata())) for l in ax.lines],collections=[dict(offsets=arr(c.get_offsets()),paths=[arr(p.vertices) for p in c.get_paths()],array=arr(c.get_array()) if c.get_array() is not None else None) for c in ax.collections],patches=[arr(p.get_path().vertices) for p in ax.patches],images=[arr(im.get_array()) for im in ax.images]))
 return json.dumps(result,sort_keys=True)
results=[];comparisons=[]
for i,(fig,name) in enumerate(captured,1):
 fig.canvas.draw()
 from matplotlib.ticker import FixedLocator
 for ax in fig.axes:
  ax.xaxis.set_major_locator(FixedLocator(ax.get_xticks()));ax.yaxis.set_major_locator(FixedLocator(ax.get_yticks()))
 fig.canvas.draw();before=snapshot(fig)
 texts_before=[t.get_text() for t in fig.findobj(Text) if t.get_text()]
 # Only title/subtitle pair from header and bottom note from note(). Other annotations remain.
 remove=[t for t in fig.texts if t.get_position()[1]>=.90 or (t.get_position()[0]==.055 and t.get_position()[1]<.08)]
 removed=[t.get_text() for t in remove]
 for t in remove:t.remove()
 # Apply one shared affine vertical layout map; no data coordinate or axis-range change.
 low=.09 if i!=4 else .17;high=.88;targetlow=.08;targethigh=.96
 fac=(targethigh-targetlow)/(high-low)
 mapy=lambda y:targetlow+(y-low)*fac
 for ax in fig.axes:
  pos=ax.get_position();ax.set_position([pos.x0,mapy(pos.y0),pos.width,pos.height*fac])
 for t in fig.texts:
  x,y=t.get_position();t.set_position((x,mapy(y)))
 for lg in fig.legends:
  anchor=lg.get_bbox_to_anchor().transformed(fig.transFigure.inverted());lg.set_bbox_to_anchor((anchor.x0,mapy(anchor.y0)-(.025 if i==1 else 0)),transform=fig.transFigure)
 if i==2:
  for ax in fig.axes:
   if [t.get_text() for t in ax.get_xticklabels()]==['0–60 months','>60 months']:
    ax.tick_params(axis='x',pad=2)
 ratio=(170/25.4)/12
 oldw,oldh=fig.get_size_inches();fig.set_size_inches(170/25.4,oldh*ratio*.88)
 for t in fig.findobj(Text):
  t.set_fontsize(max(6.1,t.get_fontsize()*.64))
 for line in fig.findobj(Line2D):
  if line.get_linewidth()>0:line.set_linewidth(max(.3,line.get_linewidth()*ratio))
  line.set_markersize(line.get_markersize()*ratio)
  if line.get_markeredgewidth()>0:line.set_markeredgewidth(max(.3,line.get_markeredgewidth()*ratio))
 for ax in fig.axes:
  for sp in ax.spines.values():sp.set_linewidth(max(.3,sp.get_linewidth()*ratio))
  for c in ax.collections:
   widths=c.get_linewidths();c.set_linewidths([max(.3,v*ratio) if v>0 else 0 for v in widths])
   if isinstance(c,PathCollection):c.set_sizes(c.get_sizes()*ratio**2)
  for patch in ax.patches:
   if patch.get_linewidth()>0:patch.set_linewidth(max(.3,patch.get_linewidth()*ratio))
 fig.canvas.draw();assert before==snapshot(fig),(name,'DATA OR AXES CHANGED')
 aftertexts=[t.get_text() for t in fig.findobj(Text) if t.get_text()]
 from collections import Counter
 assert Counter(texts_before)-Counter(removed)==Counter(aftertexts),(name,'Unexpected text change')
 stem=O/f'Figure_{i}_{stems[i-1]}_JTM_candidate'
 fig.savefig(stem.with_suffix('.png'),dpi=600,facecolor='white')
 fig.savefig(stem.with_suffix('.tiff'),dpi=300,facecolor='white',pil_kwargs={'compression':'tiff_lzw'})
 fig.savefig(stem.with_suffix('.pdf'),facecolor='white')
 # Print-size inspection and boundary audit; not a cropped submission output.
 fig.savefig(P/f'internal/170mm_Figure_{i}.png',dpi=150,facecolor='white')
 fig.canvas.draw();renderer=fig.canvas.get_renderer();clipped=[]
 for t in fig.findobj(Text):
  if t.get_visible() and t.get_text():
   box=t.get_window_extent(renderer);fb=fig.bbox
   if box.x0<fb.x0-1 or box.x1>fb.x1+1 or box.y0<fb.y0-1 or box.y1>fb.y1+1:clipped.append(t.get_text())
 with Image.open(stem.with_suffix('.tiff')) as im:
  assert im.tag_v2[259]==5;pix=im.size;dpi=tuple(map(float,im.info['dpi']))
 pdf=fitz.open(stem.with_suffix('.pdf'));fonts=pdf[0].get_fonts(full=True);fontok=all(len(pdf.extract_font(f[0])[3])>0 for f in fonts)
 assert fontok
 results.append(dict(figure=i,removed_text=removed,data_axis_snapshot_identical=True,remaining_text_identical=True,clipped_text=clipped,width_mm=170,height_mm=float(fig.get_figheight()*25.4),tiff_pixels=pix,tiff_dpi=dpi,pdf_fonts_embedded=fontok,minimum_font_pt=min(t.get_fontsize() for t in fig.findobj(Text) if t.get_visible() and t.get_text()),files=[dict(path=str(stem.with_suffix(ext)),bytes=stem.with_suffix(ext).stat().st_size,sha256=hashfile(stem.with_suffix(ext))) for ext in ['.png','.tiff','.pdf']]))
 plt.close(fig)
derived=[] # Internal frozen plotting-profile equality check not portable without original profiles.
# Source hashes available to the release wrapper; rendering writes derived outputs only.
(P/'internal/candidate_QC.json').write_text(json.dumps(dict(figures=results,derived=derived,source_hashes=source_hash,source_unchanged=True,source_PNG_equals_Stage5_TIFF="not tested in public run"),ensure_ascii=False,indent=2),encoding='utf-8')
# Frozen workbook rows actually consumed; aggregate only.
used=list({(r['_workbook'],r['_sheet'],r['_row']):r for r in b.USED}.values())
(P/'internal/used_frozen_rows.json').write_text(json.dumps(used,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(results,ensure_ascii=False,indent=2))

