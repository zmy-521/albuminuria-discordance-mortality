# Release adaptation: paths only unless noted in docs/RELEASE_CHANGES.md.
import os as _release_os
def _release_path(env, relative=''):
    root = _release_os.environ.get(env)
    if not root: raise RuntimeError('Set environment variable '+env+'; see docs/REPRODUCIBILITY.md')
    return _release_os.path.join(root, relative)
from pathlib import Path
import csv, hashlib, json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LogNorm
from matplotlib.patches import Rectangle, PathPatch
from matplotlib.path import Path as MPath
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
import base_plotting as b

ROOT=b.ROOT
C=b.C; LIGHT=b.LIGHT; INK=b.INK; MUTED=b.MUTED; GROUPS=b.GROUPS
one=b.one; panel=b.panel; header=b.header; note=b.note

def save(fig,name,supp=False):
    target=b.SUPP if supp else b.MAIN
    fig.savefig(target/(name+'.png'),dpi=300)
    fig.savefig(target/(name+'.pdf'))
    plt.close(fig)

def log_axis(ax,axis,ticks):
    if axis=='x':ax.set_xscale('log',base=2); a=ax.xaxis
    else:ax.set_yscale('log',base=2); a=ax.yaxis
    a.set_major_locator(FixedLocator(ticks));a.set_major_formatter(FixedFormatter([f'{v:g}' for v in ticks]));a.set_minor_locator(NullLocator())

def ribbon(ax,x0,x1,y0,y1,h,color,alpha=.65,z=1):
    path=MPath([(x0,y0),(x0+.34,y0),(x1-.34,y1),(x1,y1),(x1,y1+h),
                (x1-.34,y1+h),(x0+.34,y0+h),(x0,y0+h),(x0,y0)],
               [MPath.MOVETO,MPath.CURVE4,MPath.CURVE4,MPath.CURVE4,MPath.LINETO,
                MPath.CURVE4,MPath.CURVE4,MPath.CURVE4,MPath.CLOSEPOLY])
    ax.add_patch(PathPatch(path,fc=color,ec='white',lw=.45,alpha=alpha,zorder=z))

def alluvial(ax,counts):
    n=sum(counts.values()); scale=.76/n; gap=.12; bottom=.04
    left={'G0':bottom,'G2':bottom+counts['G0']*scale}
    left['G1']=bottom+(counts['G0']+counts['G2'])*scale+gap
    left['G3']=left['G1']+counts['G1']*scale
    right={'G0':bottom,'G1':bottom+counts['G0']*scale}
    right['G2']=bottom+(counts['G0']+counts['G1'])*scale+gap
    right['G3']=right['G2']+counts['G2']*scale
    for g in ['G0','G3','G1','G2']:
        h=counts[g]*scale
        ribbon(ax,.10,.90,left[g],right[g],h,C[g],.72 if g=='G2' else .42,3 if g=='G2' else 1)
        for x,ys in [(.065,left),(.90,right)]:
            ax.add_patch(Rectangle((x,ys[g]),.035,h,fc=C[g],lw=0,zorder=4))
    ax.text(.045,1.04,'Albumin',ha='left',fontsize=10,weight='bold')
    ax.text(.95,1.04,'UACR',ha='right',fontsize=10,weight='bold')
    ax.text(.00,.25,'−',ha='right',va='center',fontsize=13)
    ax.text(.00,.87,'+',ha='right',va='center',fontsize=13)
    ax.text(1,.25,'−',ha='left',va='center',fontsize=13)
    ax.text(1,.87,'+',ha='left',va='center',fontsize=13)
    ax.set_xlim(-.06,1.06);ax.set_ylim(0,1.15);ax.axis('off')

def classification_v2():
    fig=plt.figure(figsize=(12,8.8))
    header(fig,'Classification geometry and observed discordance',
           'Albumin concentration threshold 30 mg/L · UACR threshold 30 mg/g')
    panel(fig,.055,.842,'A','How the two definitions partition the plane')
    ax=fig.add_axes([.08,.30,.355,.445])
    x=np.geomspace(25,400,301);y=np.geomspace(7.5,120,301);xx,yy=np.meshgrid(x,y)
    groups=(yy>=30).astype(int)+2*(100*yy/xx>=30).astype(int)
    ax.pcolormesh(xx,yy,groups,cmap=ListedColormap([LIGHT[g] for g in GROUPS]),vmin=-.5,vmax=3.5,shading='nearest',rasterized=True)
    ax.axhline(30,color=INK,lw=1.6)
    ax.plot(x,.3*x,color=INK,lw=1.6,ls=(0,(5,3)))
    for g,xy,desc in [('G0',(175,11),'Alb− / UACR−'),('G1',(265,39),'Alb+ / UACR−'),('G2',(38,22),'Alb− / UACR+'),('G3',(50,82),'Alb+ / UACR+')]:
        ax.text(*xy,g,fontsize=15,weight='bold',ha='center',color=C[g])
        ax.annotate(desc,xy,xytext=(0,-17),textcoords='offset points',ha='center',fontsize=9,color=C[g])
    log_axis(ax,'x',[25,50,100,200,400]);log_axis(ax,'y',[7.5,15,30,60,120])
    ax.set_xlim(25,400);ax.set_ylim(7.5,120)
    ax.set_xlabel('Urinary creatinine (mg/dL; log scale)',labelpad=10)
    ax.set_ylabel('Urinary albumin (mg/L; log scale)',labelpad=8)
    ax.plot(100,30,'o',ms=5,color=INK)
    ax.annotate('(100, 30)',(100,30),xytext=(9,10),textcoords='offset points',fontsize=9)
    fig.text(.08,.21,'Solid: albumin = 30 mg/L\nDashed: UACR = 30 mg/g\nClassification schematic; no patient points.',fontsize=10,color=MUTED,linespacing=1.7,va='top')
    for ytop,source,total,letter in [(.842,'NHANES',5905,'B'),(.455,'Wanbei',2261,'C')]:
        panel(fig,.515,ytop,letter,f'{source} · n = {total:,}')
        if source=='NHANES':
            fig.text(.545,ytop-.042,'Diabetes; eGFRcr ≥60 · 1999–2018',fontsize=9.5,color=MUTED)
        else:fig.text(.545,ytop-.042,'Clinical paired urine measurements',fontsize=9.5,color=MUTED)
        counts={g:(one('main','Group counts',threshold_mgl=30,group=g) if source=='NHANES' else one('main','Hospital groups',source='Wanbei',threshold=30,group=g))['n'] for g in GROUPS}
        ax=fig.add_axes([.55,ytop-.255,.37,.19]);alluvial(ax,counts)
        for j,g in enumerate(GROUPS):
            fig.text(.545+j*.1,ytop-.282,f'{g}: {counts[g]:,}',fontsize=10,color=C[g],weight='bold' if g=='G2' else 'normal')
        if source=='NHANES':
            r1=one('main','Classification',source='NHANES5905',threshold=30,denominator='concentration_negative')
            r2=one('main','Classification',source='NHANES5905',threshold=30,denominator='ACR_positive')
            st=f"G2 among Alb−: {r1['weighted']*100:.2f}% ({r1['low']*100:.2f}–{r1['high']*100:.2f})\nG2 among UACR+: {r2['weighted']*100:.2f}% ({r2['low']*100:.2f}–{r2['high']*100:.2f})"
        else:
            r1=one('main','Hospital classification',source='Wanbei',threshold=30,denominator='concentration_negative')
            r2=one('main','Hospital classification',source='Wanbei',threshold=30,denominator='ACR_positive')
            st=f"G2 among Alb−: {r1['fraction']*100:.2f}% ({r1['numerator']}/{r1['n']:,})\nG2 among UACR+: {r2['fraction']*100:.2f}% ({r2['numerator']}/{r2['n']:,})"
        fig.text(.545,ytop-.313,st,fontsize=9.5,color=MUTED,linespacing=1.5,va='top')
    note(fig,'B–C: within each source, ribbon widths are proportional to unweighted counts; each source is normalized separately.\nRibbons do not represent transitions over time. NHANES proportions use survey weights and 95% CIs;\nWanbei proportions are unweighted. Wanbei supplies no mortality validation.',y=.021)
    save(fig,'Candidate_01_Classification')

def components_v2():
    ba=one('extra','Effects',scenario='original_linear_eGFR',contrast='Alb_per_doubling')
    bc=one('extra','Effects',scenario='original_linear_eGFR',contrast='UCr_per_doubling_at_fixed_Alb')
    fixed=one('extra','Effects',scenario='original_linear_eGFR',contrast='UCr_per_doubling_at_fixed_ACR')
    cov=one('extra','Joint covariance')
    B=np.array([ba['logHR'],bc['logHR']]);V=np.array([[cov['variance_albumin'],cov['covariance']],[cov['covariance'],cov['variance_creatinine']]])
    assert np.allclose(B.sum(),fixed['logHR'],atol=1e-14)
    assert np.allclose(V.sum(),fixed['variance'],atol=2e-15)
    x=np.geomspace(25,400,241);y=np.geomspace(7.5,120,241);xx,yy=np.meshgrid(x,y)
    da=np.log2(yy/30);dc=np.log2(xx/100);eta=B[0]*da+B[1]*dc
    var=V[0,0]*da**2+V[1,1]*dc**2+2*V[0,1]*da*dc
    hr=np.exp(eta);lo=np.exp(eta-1.96*np.sqrt(np.maximum(var,0)));hi=np.exp(eta+1.96*np.sqrt(np.maximum(var,0)))
    with (ROOT/'Data/component_surface_grid.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['UCr_mg_dL','Alb_mg_L','conditional_HR','pointwise95_low','pointwise95_high','grid_type'])
        for vals in zip(xx.ravel(),yy.ravel(),hr.ravel(),lo.ravel(),hi.ravel()):w.writerow([*vals,'model_display_grid_not_patient_data'])
    fig=plt.figure(figsize=(12,8.4))
    header(fig,'Conditional associations in the urine-component plane',
           'Post hoc component Cox model · n = 5,155; deaths = 1,115 · Other covariates held fixed')
    panel(fig,.055,.842,'A','Model-derived conditional HR surface')
    ax=fig.add_axes([.082,.32,.40,.45])
    mesh=ax.pcolormesh(xx,yy,hr,norm=LogNorm(vmin=.4,vmax=2.5),cmap='RdBu_r',shading='auto',rasterized=True)
    contours=ax.contour(xx,yy,hr,levels=[.5,.75,1,1.5,2],colors=INK,linewidths=.8,alpha=.8)
    ax.clabel(contours,fmt=lambda h:f'HR {h:g}',fontsize=9,inline=True)
    for acr in [10,30,100]:
        Y=x*acr/100;valid=(Y>=7.5)&(Y<=120)
        ax.plot(x[valid],Y[valid],color='white',lw=2.4,ls=(0,(5,3)))
        ax.plot(x[valid],Y[valid],color='#303B43',lw=.75,ls=(0,(5,3)))
    log_axis(ax,'x',[25,50,100,200,400]);log_axis(ax,'y',[7.5,15,30,60,120])
    ax.set_xlim(25,400);ax.set_ylim(7.5,120)
    ax.set_xlabel('Urinary creatinine (mg/dL; log scale)',labelpad=9)
    ax.set_ylabel('Urinary albumin (mg/L; log scale)',labelpad=8)
    ax.plot(100,30,'o',mfc='white',mec=INK,ms=7,mew=1.5,zorder=10)
    ax.annotate('Reference',(100,30),xytext=(8,-20),textcoords='offset points',fontsize=9,weight='bold')
    for acr,px in [(10,230),(30,270),(100,70)]:
        ax.text(px,px*acr/100,f'UACR {acr}',rotation=36,ha='center',va='bottom',fontsize=8.5,color=INK,bbox=dict(fc='white',ec='none',alpha=.8,pad=1))
    cb=fig.colorbar(mesh,cax=fig.add_axes([.505,.32,.015,.45]),ticks=[.5,.75,1,1.5,2])
    cb.ax.set_yticklabels(['0.5','0.75','1','1.5','2']);cb.ax.minorticks_off()
    fig.text(.082,.246,'Reference coordinate: Alb 30 mg/L; UCr 100 mg/dL.\nSolid: equal HR. Dashed: equal UACR (mg/g).\nIllustrative window; observed joint density is not shown.',fontsize=9.5,color=MUTED,linespacing=1.6,va='top')
    panel(fig,.62,.842,'B','Along a fixed-UACR line')
    ax2=fig.add_axes([.67,.49,.28,.235])
    k=np.geomspace(.5,2,201);t=np.log2(k);e=np.exp(t*B.sum());s=np.abs(t)*np.sqrt(V.sum())
    low=np.exp(t*B.sum()-1.96*s);high=np.exp(t*B.sum()+1.96*s)
    ax2.fill_between(k,low,high,color=C['G2'],alpha=.16,lw=0)
    ax2.plot(k,e,color=C['G2'],lw=2.3)
    ax2.axhline(1,color=MUTED,lw=.8,ls=':');ax2.plot([1,2],[1,fixed['HR']],'o',color=C['G2'],ms=5)
    log_axis(ax2,'x',[.5,1,2]);ax2.set_ylim(.77,1.29);ax2.set_xlim(.5,2)
    ax2.set_xlabel('Both concentrations multiplied by k',fontsize=10)
    ax2.set_ylabel('HR vs reference coordinate',fontsize=10);b.tidy(ax2,'y')
    fig.text(.65,.393,f"Per doubling at fixed UACR\nHR {b.fmt(fixed,3)}\nWald P = {fixed['p']:.6f}",fontsize=11,color=C['G2'],linespacing=1.55,va='top',weight='bold')
    fig.text(.65,.267,'Band: pointwise 95% CI from the full\ncoefficient covariance. The surface and\ncurve display the same fitted model.',fontsize=9.5,color=MUTED,linespacing=1.5,va='top')
    a2=one('extra','Effects',scenario='original_linear_eGFR',contrast='G2_vs_G3_A2_adjust_log2ACR')
    fig.text(.055,.126,f"Within A2, G2 vs G3-A2: HR {b.fmt(a2,3)}; P = {a2['p']:.3f}.",fontsize=10.5,weight='bold')
    note(fig,'The reference coordinate is not the G0 reference group. Surface colors do not represent group-average HRs.\nThese are conditional associations under the specified log-linear model, not causal effects or a new clinical score.\nThe A2 comparison is a different contrast; its support distribution is retained in Supplement 02.',y=.028)
    save(fig,'Candidate_03_Components_A2')
    with (ROOT/'Data/fixed_UACR_curve.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['concentration_multiplier','conditional_HR','pointwise95_low','pointwise95_high']);w.writerows(zip(k,e,low,high))
    assert np.allclose([e[-1],low[-1],high[-1]],[fixed['HR'],fixed['low'],fixed['high']],atol=1e-13)

def absolute_risk_v2():
    fig=plt.figure(figsize=(12,7.8))
    header(fig,'Five-year mortality: absolute risks and risk differences',
           'NHANES 1999–2014 · Standardization target n = 3,789 · Estimates with 95% CIs')
    panel(fig,.055,.837,'A','Standardized five-year risk')
    panel(fig,.60,.837,'B','Between-group risk differences')
    ax=fig.add_axes([.075,.24,.43,.46]);ax.set_ylim(0,17.2);ax.set_xlim(-.65,3.65)
    ax.set_xticks(range(4),GROUPS);ax.set_ylabel('Five-year mortality risk (%)');b.tidy(ax,'y')
    ax2=fig.add_axes([.72,.24,.23,.46]);ax2.set_ylim(-.5,2.9);ax2.invert_yaxis();ax2.set_xlim(-6.3,9)
    ax2.set_yticks([.5,2],['G2 − G0','G2 − G1']);ax2.tick_params(axis='y',pad=10)
    ax2.set_xlabel('Risk difference (percentage points)',fontsize=10);b.tidy(ax2)
    ax2.axvline(0,color=MUTED,lw=.8,ls=(0,(3,3)))
    styles=[('logistic',-.17,INK,'o',False),('cox',.17,'#8296A3','s',True)]
    for method,offset,col,marker,opened in styles:
        for j,g in enumerate(GROUPS):
            key=f'30_known_risk_{g}' if method=='logistic' else f'AdjustedB_{g}'
            r=one('main','Absolute risk',estimand=key);e,l,h=[100*r[k] for k in ['estimate','low','high']]
            ax.bar(j+offset,e,width=.29,color=col,alpha=.9,edgecolor='white',lw=.7)
            ax.errorbar(j+offset,e,yerr=[[e-l],[h-e]],fmt='none',ecolor=col,elinewidth=1.25,capsize=3)
            ax.text(j+offset,h+.37,f'{e:.2f}',ha='center',fontsize=10,color=col)
        for y,contrast in [(.5,'G2_G0'),(2,'G2_G1')]:
            key=f'30_known_RD_{contrast}' if method=='logistic' else f'AdjustedB_RD_{contrast}.AdjustedB_G2'
            r=one('main','Absolute risk',estimand=key)
            b.interval(ax2,y+offset,r,col,marker,opened,('estimate','low','high'),100)
            ax2.text(.02,y+offset+.19,b.fmt(r,2,('estimate','low','high'),100),transform=ax2.get_yaxis_transform(),fontsize=9,va='top',color=col)
    fig.legend(handles=[Rectangle((0,0),1,1,fc=INK,label='Fixed-horizon logistic (primary)'),Rectangle((0,0),1,1,fc='#8296A3',label='Cox standardization (sensitivity)')],loc='upper left',bbox_to_anchor=(.068,.788),frameon=False,ncol=2,fontsize=10)
    note(fig,'Logistic fitting: 3,785 participants with known 60-month outcome; standardization target: 3,789.\nCox uses the frozen 1999–2014 complete-case analysis. Risk-difference CIs are the frozen joint estimates.\nG0: Alb−/UACR−; G1: Alb+/UACR−; G2: Alb−/UACR+; G3: Alb+/UACR+. Thresholds: 30 mg/L and 30 mg/g.',y=.045)
    save(fig,'Candidate_04_FiveYear_Risk')

def quartile_box(ax,j,r,g):
    q1,m,q3=r['q1'],r['median'],r['q3'];assert q1<=m<=q3
    ax.add_patch(Rectangle((j-.25,q1),.5,q3-q1,fc=LIGHT[g],ec=C[g],lw=1.5))
    ax.plot([j-.25,j+.25],[m,m],lw=2.7,color=C[g])
    ax.annotate(f'{m:.2f}' if m<10 else f'{m:.1f}',(j,q3),xytext=(0,8),textcoords='offset points',ha='center',fontsize=10.5,color=C[g])

def flow_v2():
    fig=plt.figure(figsize=(12,8.8))
    header(fig,'Urinary creatinine, flow and short-interval excretion',
           'NHANES 2009–2018 · Evaluable subset n = 3,261 / 3,539 · Survey-weighted medians and interquartile ranges')
    settings=[('ucr_adjusted','Urinary creatinine','mg/dL',(0,290),[0,100,200],.055,.84,'A'),('flow','Urine flow','mL/min',(0,2.45),[0,1,2],.55,.84,'B'),('CER','Estimated creatinine excretion','mg/min',(0,1.52),[0,.5,1,1.5],.055,.48,'C'),('AER','Estimated albumin excretion','µg/min · log scale',(2,480),[3,10,30,100,300],.55,.48,'D')]
    for variable,title,unit,limits,ticks,x,top,letter in settings:
        panel(fig,x,top,letter,title)
        ax=fig.add_axes([x+.04,top-.28,.365,.215]);ax.set_xlim(-.65,3.65);ax.set_ylim(*limits)
        if variable=='AER':log_axis(ax,'y',ticks)
        else:ax.set_yticks(ticks)
        for j,g in enumerate(GROUPS):quartile_box(ax,j,one('extra','Flow characteristics',population=g,variable=variable),g)
        ax.set_xticks(range(4),GROUPS);ax.set_ylabel(unit);b.tidy(ax,'y')
        for t,g in zip(ax.get_xticklabels(),GROUPS):t.set_color(C[g]);t.set_weight('bold' if g=='G2' else 'normal')
    fig.text(.095,.143,'G0: n = 2,276       G1: n = 176       G2: n = 123       G3: n = 686',fontsize=10.5)
    note(fig,'Boxes span Q1–Q3; central lines and labels are medians. These are distribution summaries, not confidence intervals.\nWhiskers, outliers and density are not shown; no patient values were simulated. AER/CER are short-interval estimates\nusing the same urine analytes and flow. They are not independent diagnostic reference measures or 24-hour excretion.',y=.036)
    save(fig,'Candidate_05_Urine_Flow')

def a2_supp():
    fig=plt.figure(figsize=(12,7.5))
    header(fig,'Within-A2 comparison and observed UACR support',
           'Supplement candidate · A2: UACR 30–<300 mg/g · Complete-case comparison n = 1,120; deaths = 350')
    r=one('extra','Effects',scenario='original_linear_eGFR',contrast='G2_vs_G3_A2_adjust_log2ACR')
    panel(fig,.055,.832,'A','Adjusted G2 vs G3-A2 comparison')
    fig.text(.085,.755,f"HR {b.fmt(r,3)}\nP = {r['p']:.3f}",fontsize=18,weight='bold',linespacing=1.7,va='top',color=C['G2'])
    fig.text(.085,.57,'Model B covariates + log₂ UACR\nG2: 222 participants / 76 deaths\nG3-A2: 898 participants / 274 deaths',fontsize=10.5,color=MUTED,linespacing=1.8,va='top')
    panel(fig,.53,.832,'B','UACR distribution in the A2 model')
    ax=fig.add_axes([.6,.46,.34,.28]);cols=['#D5E2E7','#91B6C3','#4D879E','#225266']
    bins=['30-<60','60-<100','100-<200','200-<300']
    for j,g,n in [(0,'G2',222),(1,'G3',898)]:
        left=0
        for k,inter in enumerate(bins):
            r=one('extra','A2 support bins',population='A2_B',group=g,interval=inter);w=r['n']/n*100
            ax.barh(j,w,left=left,height=.46,color=cols[k],edgecolor='white',lw=1)
            if w>8:ax.text(left+w/2,j,f'{w:.1f}%',ha='center',va='center',color=INK if k<2 else 'white',fontsize=10)
            left+=w
        assert np.isclose(left,100)
    ax.set_yticks([0,1],['G2','G3-A2']);ax.invert_yaxis();ax.set_xlim(0,100)
    ax.set_xlabel('Within-group proportion (unweighted, %)');ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)
    fig.legend(handles=[Rectangle((0,0),1,1,fc=co,label=bi.replace('-<','–<')) for co,bi in zip(cols,bins)],bbox_to_anchor=(.60,.40),loc='upper left',ncol=2,frameon=False,title='UACR (mg/g)',fontsize=10)
    # Exact bin counts retain small and zero segments without invented density.
    fig.text(.085,.315,'Bin counts (30–<60 / 60–<100 / 100–<200 / 200–<300 mg/g)\nG2: 171 / 39 / 12 / 0        G3-A2: 337 / 253 / 224 / 84',fontsize=10.5,color=MUTED,linespacing=1.8,va='top')
    note(fig,'This contrast is imprecise and adjusts for log₂ UACR despite limited distributional overlap. It does not independently\nvalidate the component model. A non-significant result neither establishes equivalence nor excludes a distinct phenotype.',y=.07)
    save(fig,'Supplement_02_A2_Support',True)

def main():
    # Keep the primary HR plot and sensitivity estimates; change representation only where data support it.
    b.save=save
    classification_v2();b.main_mortality();components_v2();absolute_risk_v2()
    b.flow()
    for ext in ['png','pdf']:
        (b.MAIN/f'Candidate_05_Urine_Flow.{ext}').replace(b.SUPP/f'Supplement_03_Flow_Means_CI.{ext}')
    flow_v2();b.sensitivity();a2_supp()
    unique={(r['_workbook'],r['_sheet'],r['_row']):r for r in b.USED}
    obj={'version':'v2 review candidates; no new model fits','input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [b.MAINFILE,b.EXTRAFILE]},'source_rows':list(unique.values())}
    (ROOT/'Data/figure_source_rows.json').write_text(json.dumps(obj,ensure_ascii=False,indent=2))
    with (ROOT/'Data/figure_source_cells.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.writer(f);w.writerow(['Workbook','Sheet','Excel row','Field','Frozen value'])
        for r in unique.values():
            for k,v in r.items():
                if not k.startswith('_'):w.writerow([r['_workbook'],r['_sheet'],r['_row'],k,v])
    print(json.dumps({'pngs':len(list(ROOT.rglob('*.png'))),'source_rows':len(unique),'output':str(ROOT)}))

if __name__=='__main__':main()
