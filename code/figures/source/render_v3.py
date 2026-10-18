# Release adaptation: paths only unless noted in docs/RELEASE_CHANGES.md.
import os as _release_os
def _release_path(env, relative=''):
    root = _release_os.environ.get(env)
    if not root: raise RuntimeError('Set environment variable '+env+'; see docs/REPRODUCIBILITY.md')
    return _release_os.path.join(root, relative)
from pathlib import Path
import csv, hashlib, json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
from matplotlib.lines import Line2D
import base_plotting as b
import render_v2 as v2

ROOT=b.ROOT; C=b.C; LIGHT=b.LIGHT; INK=b.INK; MUTED=b.MUTED; GROUPS=b.GROUPS
one=b.one; header=b.header; panel=b.panel; note=b.note; save=v2.save
N=pd.read_csv(ROOT/'Inputs/01_NHANES_5905_plot.csv')
W=pd.read_csv(ROOT/'Inputs/02_Wanbei_2261_plot.csv')
F=pd.read_csv(ROOT/'Inputs/03_Flow_plot.csv')
ALB='urine_albumin_mg_L'; CR='urine_creatinine_mg_dL_CDC'; ACR='UACR_mg_g_CDC'
WT='original_survey_weight'; G='original_group'
AUDIT={}

def logaxis(ax,axis,ticks):
    if axis=='x':ax.set_xscale('log');a=ax.xaxis
    else:ax.set_yscale('log');a=ax.yaxis
    a.set_major_locator(FixedLocator(ticks));a.set_major_formatter(FixedFormatter([f'{t:g}' for t in ticks]));a.set_minor_locator(NullLocator())

def empirical_cdf(values,weights):
    u,inv=np.unique(values,return_inverse=True)
    w=np.bincount(inv,weights=np.asarray(weights,dtype=float));cdf=np.cumsum(w)/w.sum()
    return u,cdf

def check_inputs():
    for name,d,expected in [('NHANES',N,5905),('Wanbei',W,2261),('Flow',F,3539)]:
        assert len(d)==expected and d.anonymous_row_id.is_unique
        AUDIT[name]={'n':len(d),'groups':d[G].value_counts().sort_index().to_dict()}
    assert N[G].value_counts().to_dict()=={'G0':3953,'G1':314,'G2':251,'G3':1387}
    assert W[G].value_counts().to_dict()=={'G0':1533,'G1':101,'G2':152,'G3':475}
    assert N.model_B_actual.sum()==5155 and N.A2_ModelB_actual.sum()==1120
    assert np.array_equal(N.A2_ModelB_actual.astype(bool),N.model_B_actual.astype(bool)&N[ACR].ge(30)&N[ACR].lt(300))
    for d,cr,acr in [(N,CR,ACR),(W,'urine_creatinine_mg_dL','UACR_mg_g')]:
        assert np.allclose(d[acr],100*d[ALB]/d[cr],rtol=1e-12)
        groups=np.select([d[ALB].lt(30)&d[acr].lt(30),d[ALB].ge(30)&d[acr].lt(30),d[ALB].lt(30)&d[acr].ge(30)],['G0','G1','G2'],'G3')
        assert np.array_equal(groups,d[G])
    f=F[F.flow_evaluable.eq(1)];assert len(f)==3261
    assert np.allclose(f.short_AER_ug_min,f[ALB]*f.onsite_urine_flow_mL_min)
    assert np.allclose(f.short_CER_mg_min,f[CR]*f.onsite_urine_flow_mL_min/100)
    AUDIT['limits']={'curves_available':False,'individual_followup_available':False,'no_models_fitted':True}

def scatter_plane(ax,data,creat,groups=True):
    assert data[creat].between(3,1000).all() and data[ALB].between(.15,20000).all()
    if groups:
        for g in ['G0','G3','G1','G2']:
            d=data[data[G].eq(g)]
            ax.scatter(d[creat],d[ALB],s=9 if g in ['G0','G3'] else 12,c=C[g],alpha=.22 if g in ['G0','G3'] else .53,lw=0,rasterized=True)
    else:ax.scatter(data[creat],data[ALB],s=8,c='#53636D',alpha=.20,lw=0,rasterized=True)
    x=np.geomspace(3,1000,301)
    ax.axhline(30,color=INK,lw=1.1)
    ax.plot(x,.3*x,color=INK,lw=1.1,ls=(0,(5,3)))
    ax.set_xlim(3,1000);ax.set_ylim(.15,20000)
    logaxis(ax,'x',[3,10,30,100,300,1000]);logaxis(ax,'y',[.3,3,30,300,3000])
    ax.set_xlabel('Urinary creatinine (mg/dL; log scale)',fontsize=10)
    ax.set_ylabel('Urinary albumin (mg/L; log scale)',fontsize=10)
    ax.tick_params(labelsize=9)

def classification_v3():
    fig=plt.figure(figsize=(12,10.4))
    header(fig,'Observed urine measurements and classification discordance',
           'Real participant coordinates · Albumin threshold 30 mg/L; UACR threshold 30 mg/g')
    for x,d,cr,name,letter in [(.055,N,CR,'NHANES 1999–2018 · n = 5,905','A'),(.535,W,'urine_creatinine_mg_dL','Wanbei clinical sample · n = 2,261','B')]:
        panel(fig,x,.842,letter,name)
        ax=fig.add_axes([x+.037,.478,.385,.30]);scatter_plane(ax,d,cr)
        ax.text(4,39,'Albumin 30',fontsize=8.5,color=INK,bbox=dict(fc='white',ec='none',alpha=.85,pad=1.5))
        ax.text(160,120,'UACR 30',fontsize=8.5,color=INK,rotation=20,bbox=dict(fc='white',ec='none',alpha=.85,pad=1.5))
    fig.legend(handles=[Line2D([0],[0],marker='o',color=C[g],lw=0,label=g,ms=6) for g in GROUPS],loc='center',bbox_to_anchor=(.52,.419),frameon=False,ncol=4,fontsize=10)
    for x,d,source,letter in [(.055,N,'NHANES','C'),(.535,W,'Wanbei','D')]:
        panel(fig,x,.373,letter,f'{source}: paired classifications')
        counts=d[G].value_counts().to_dict()
        for g in GROUPS:
            r=one('main','Group counts',threshold_mgl=30,group=g) if source=='NHANES' else one('main','Hospital groups',source='Wanbei',threshold=30,group=g)
            assert r['n']==counts[g]
        ax=fig.add_axes([x+.047,.174,.37,.151]);v2.alluvial(ax,counts)
        for j,g in enumerate(GROUPS):fig.text(x+.037+j*.1,.154,f'{g}: {counts[g]:,}',color=C[g],fontsize=10,weight='bold' if g=='G2' else 'normal')
        if source=='NHANES':
            r1=one('main','Classification',source='NHANES5905',threshold=30,denominator='concentration_negative');r2=one('main','Classification',source='NHANES5905',threshold=30,denominator='ACR_positive')
            st=f"G2 among Alb−: {100*r1['weighted']:.2f}% ({100*r1['low']:.2f}–{100*r1['high']:.2f})\nG2 among UACR+: {100*r2['weighted']:.2f}% ({100*r2['low']:.2f}–{100*r2['high']:.2f})"
        else:
            r1=one('main','Hospital classification',source='Wanbei',threshold=30,denominator='concentration_negative');r2=one('main','Hospital classification',source='Wanbei',threshold=30,denominator='ACR_positive')
            st=f"G2 among Alb−: {100*r1['fraction']:.2f}% ({r1['numerator']}/{r1['n']:,})\nG2 among UACR+: {100*r2['fraction']:.2f}% ({r2['numerator']}/{r2['n']:,})"
        fig.text(x+.037,.123,st,fontsize=9.5,color=MUTED,va='top',linespacing=1.5)
    note(fig,'All observed coordinates are retained; overlapping points are possible. Points and ribbons are unweighted.\nRibbon totals are normalized separately by source and do not represent time transitions. NHANES proportions\nuse survey weights and 95% CIs; Wanbei proportions are unweighted. Wanbei provides no mortality validation.',y=.014)
    save(fig,'Candidate_01_Classification')

def get_component_params():
    coeff=pd.read_csv(ROOT/'Inputs/07_Frozen_component_coefficients.csv').set_index('term')['coefficient']
    cov=pd.read_csv(ROOT/'Inputs/08_Frozen_component_covariance.csv',index_col=0)
    terms=['log2alb','log2ucr'];B=coeff.loc[terms].to_numpy();V=cov.loc[terms,terms].to_numpy()
    ba=one('extra','Effects',scenario='original_linear_eGFR',contrast='Alb_per_doubling');bc=one('extra','Effects',scenario='original_linear_eGFR',contrast='UCr_per_doubling_at_fixed_Alb')
    f=one('extra','Effects',scenario='original_linear_eGFR',contrast='UCr_per_doubling_at_fixed_ACR')
    q=one('extra','Joint covariance')
    assert np.allclose(B,[ba['logHR'],bc['logHR']],rtol=0,atol=1e-14)
    assert np.allclose(V,[[q['variance_albumin'],q['covariance']],[q['covariance'],q['variance_creatinine']]],rtol=0,atol=1e-14)
    ci=np.exp(B.sum()+np.array([-1,1])*1.96*np.sqrt(V.sum()))
    assert np.allclose([np.exp(B.sum()),*ci],[f['HR'],f['low'],f['high']],atol=1e-13)
    return B,V,f

def components_v3():
    B,V,fixed=get_component_params();model=N[N.model_B_actual.eq(1)]
    in_window=model[CR].between(25,400)&model[ALB].between(7.5,120)
    AUDIT['component_model_display_window']={'n_full':len(model),'n_in_window':int(in_window.sum()),'ucr_mg_dL':[25,400],'alb_mg_L':[7.5,120],'window_set_in_v2_before_patient_data':True}
    fig=plt.figure(figsize=(12,10.7))
    header(fig,'Component associations within observed urine measurements',
           'Post hoc explanatory analyses · Frozen log-linear Cox model · Other covariates held fixed')
    panel(fig,.055,.848,'A','Observed Model B measurements · n = 5,155')
    ax=fig.add_axes([.088,.548,.348,.246]);scatter_plane(ax,model,CR,False)
    ax.add_patch(Rectangle((25,7.5),375,112.5,fill=False,ec='#AA7C28',lw=2))
    ax.text(.025,.95,f'Outlined window: {int(in_window.sum()):,} participants',transform=ax.transAxes,va='top',fontsize=9,color='#815C1F',bbox=dict(fc='white',ec='none',alpha=.88,pad=2))
    panel(fig,.535,.848,'B','Conditional HR in the outlined window')
    ax2=fig.add_axes([.575,.548,.322,.246])
    x=np.geomspace(25,400,241);y=np.geomspace(7.5,120,241);xx,yy=np.meshgrid(x,y)
    da=np.log2(yy/30);dc=np.log2(xx/100);eta=B[0]*da+B[1]*dc;hr=np.exp(eta)
    variance=V[0,0]*da**2+V[1,1]*dc**2+2*V[0,1]*da*dc
    mesh=ax2.pcolormesh(xx,yy,hr,norm=LogNorm(vmin=.4,vmax=2.5),cmap='RdBu_r',shading='auto',rasterized=True)
    cs=ax2.contour(xx,yy,hr,levels=[.5,.75,1,1.5,2],colors=INK,linewidths=.7,alpha=.8)
    ax2.clabel(cs,fmt=lambda a:f'{a:g}',fontsize=8.5,inline=True)
    for acr in [10,30,100]:
        z=x*acr/100;valid=(z>=7.5)&(z<=120)
        ax2.plot(x[valid],z[valid],color='white',lw=2.6,ls=(0,(5,3)))
        ax2.plot(x[valid],z[valid],color='#303B43',lw=.7,ls=(0,(5,3)))
    logaxis(ax2,'x',[25,50,100,200,400]);logaxis(ax2,'y',[7.5,15,30,60,120])
    ax2.set_xlim(25,400);ax2.set_ylim(7.5,120)
    ax2.set_xlabel('Urinary creatinine (mg/dL; log scale)',fontsize=10)
    ax2.set_ylabel('Urinary albumin (mg/L; log scale)',fontsize=10)
    ax2.plot(100,30,'o',mfc='white',mec=INK,ms=6,mew=1.3)
    cb=fig.colorbar(mesh,cax=fig.add_axes([.916,.548,.013,.246]),ticks=[.5,1,2]);cb.ax.set_yticklabels(['0.5','1','2']);cb.ax.minorticks_off();cb.ax.set_title('HR',fontsize=9,pad=8)
    fig.text(.575,.494,'Solid: equal HR. Dashed: UACR 10, 30, 100 mg/g.\nReference: Alb 30 mg/L, UCr 100 mg/dL (white point).',fontsize=9,color=MUTED,va='top',linespacing=1.4)
    panel(fig,.055,.423,'C','Fixed-UACR contrast with pointwise 95% CI')
    ax3=fig.add_axes([.088,.172,.348,.176])
    k=np.geomspace(.5,2,201);t=np.log2(k);e=np.exp(t*B.sum());s=np.abs(t)*np.sqrt(V.sum());lo=np.exp(t*B.sum()-1.96*s);hi=np.exp(t*B.sum()+1.96*s)
    ax3.fill_between(k,lo,hi,color=C['G2'],alpha=.17,lw=0);ax3.plot(k,e,color=C['G2'],lw=2)
    ax3.axhline(1,color=MUTED,lw=.7,ls=':');ax3.plot([1,2],[1,fixed['HR']],'o',color=C['G2'],ms=4)
    logaxis(ax3,'x',[.5,1,2]);ax3.set_xlim(.5,2);ax3.set_ylim(.77,1.29)
    ax3.set_xlabel('Both concentrations multiplied by k',fontsize=10);ax3.set_ylabel('HR vs reference coordinate',fontsize=10);b.tidy(ax3,'y')
    fig.text(.088,.378,f"Per doubling: {b.fmt(fixed,3)}; P = {fixed['p']:.6f}",color=C['G2'],fontsize=10,weight='bold')
    panel(fig,.535,.423,'D','Within A2: estimate and UACR distribution')
    a2=one('extra','Effects',scenario='original_linear_eGFR',contrast='G2_vs_G3_A2_adjust_log2ACR')
    fig.text(.575,.378,f"G2 vs G3-A2: {b.fmt(a2,3)}; P = {a2['p']:.3f}",fontsize=10,weight='bold')
    ax4=fig.add_axes([.59,.172,.34,.176]);a2data=N[N.A2_ModelB_actual.eq(1)]
    curves=[]
    for g,label in [('G2','G2 (n = 222)'),('G3','G3-A2 (n = 898)')]:
        d=a2data[a2data[G].eq(g)];u,cdf=empirical_cdf(d[ACR],d[WT])
        ax4.step(np.r_[30,u,300],np.r_[0,100*cdf,100],where='post',color=C[g],lw=1.8,label=label)
        curves += [{'group':g,'UACR_mg_g':float(a),'weighted_cumulative_percent':float(p*100)} for a,p in zip(u,cdf)]
    ax4.set_xlim(30,300);ax4.set_ylim(0,103);ax4.set_yticks([0,25,50,75,100])
    logaxis(ax4,'x',[30,60,100,200,300]);ax4.set_xlabel('UACR (mg/g; log scale)',fontsize=10);ax4.set_ylabel('Weighted cumulative %',fontsize=10)
    ax4.legend(frameon=False,fontsize=8.5,loc='lower right');b.tidy(ax4,'y')
    note(fig,'A: all observed coordinates, unweighted display. B–C: the same frozen model, not separate evidence; the reference is\na coordinate, not G0. The window does not guarantee joint support at every location. D: weighted empirical distribution,\nnot survival; A2 adjusts for log₂ UACR and is not independent validation. All associations remain model-dependent.',y=.027)
    save(fig,'Candidate_03_Components_A2')
    pd.DataFrame(curves).to_csv(ROOT/'Data/A2_weighted_empirical_CDF.csv',index=False)
    pd.DataFrame({'multiplier':k,'HR':e,'pointwise95_low':lo,'pointwise95_high':hi}).to_csv(ROOT/'Data/fixed_UACR_curve.csv',index=False)
    pd.DataFrame({'UCr_mg_dL':xx.ravel(),'Alb_mg_L':yy.ravel(),'conditional_HR':hr.ravel(),'pointwise95_low':np.exp(eta-1.96*np.sqrt(np.maximum(variance,0))).ravel(),'pointwise95_high':np.exp(eta+1.96*np.sqrt(np.maximum(variance,0))).ravel()}).to_csv(ROOT/'Data/component_surface_display_grid.csv',index=False)

def flow_v3():
    f=F[F.flow_evaluable.eq(1)].copy();fig=plt.figure(figsize=(12,10))
    header(fig,'Observed distributions of urine concentration and excretion',
           'NHANES 2009–2018 · Flow-evaluable n = 3,261 / 3,539 · Survey-weighted distribution profiles')
    settings=[(CR,'ucr_adjusted','Urinary creatinine','mg/dL',(3,800),[3,10,30,100,300],.055,.845,'A'),
              ('onsite_urine_flow_mL_min','flow','Urine flow','mL/min',(.004,80),[.01,.1,1,10,50],.55,.845,'B'),
              ('short_CER_mg_min','CER','Short-interval creatinine excretion','mg/min',(.007,20),[.01,.1,1,10],.055,.464,'C'),
              ('short_AER_ug_min','AER','Short-interval albumin excretion','µg/min',(.15,11000),[1,10,100,1000,10000],.55,.464,'D')]
    profiles=[]
    for col,key,title,unit,limits,ticks,x,top,letter in settings:
        panel(fig,x,top,letter,title)
        ax=fig.add_axes([x+.04,top-.281,.365,.220]);ax.set_xlim(-.58,3.58);ax.set_ylim(*limits)
        logaxis(ax,'y',ticks);ax.set_ylabel(unit+' · log scale',fontsize=10)
        edges=np.linspace(np.log10(limits[0]),np.log10(limits[1]),37);width=edges[1]-edges[0]
        densities={}
        for g in GROUPS:
            d=f[f[G].eq(g)];vals=d[col].to_numpy();weights=d[WT].to_numpy()
            assert np.all((vals>=limits[0])&(vals<=limits[1]))
            counts,_=np.histogram(np.log10(vals),bins=edges);mass,_=np.histogram(np.log10(vals),bins=edges,weights=weights)
            assert counts.sum()==len(d) and np.isclose(mass.sum(),weights.sum())
            den=mass/mass.sum()/width;densities[g]=den
            for k in range(len(den)):profiles.append({'variable':key,'group':g,'bin_low':10**edges[k],'bin_high':10**edges[k+1],'unweighted_n':int(counts[k]),'weighted_fraction':float(mass[k]/mass.sum()),'log10_density':float(den[k])})
        peak=max(v.max() for v in densities.values())
        for j,g in enumerate(GROUPS):
            d=f[f[G].eq(g)];r=one('extra','Flow characteristics',population=g,variable=key)
            assert r['n']==len(d)
            # Symmetric histogram outline: no KDE, no invented patient coordinates.
            yy=10**np.repeat(edges,2)[1:-1];half=np.repeat(densities[g]/peak*.34,2)
            ax.fill_betweenx(yy,j-half,j+half,facecolor=C[g],alpha=.30,edgecolor='none',lw=0)
            ax.plot([j,j],[d[col].min(),d[col].max()],color=C[g],lw=.8,alpha=.75)
            ax.add_patch(Rectangle((j-.066,r['q1']),.132,r['q3']-r['q1'],fc='white',ec=C[g],lw=1.15,zorder=4))
            ax.plot([j-.075,j+.075],[r['median'],r['median']],color=C[g],lw=2.2,zorder=5)
            md=r['median'];label=f'{md:.2f}' if md<10 else f'{md:.1f}'
            ax.text(j,.985,label,transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=10.5,color=C[g],weight='bold' if g=='G2' else 'normal')
        ax.set_xticks(range(4),[f'{g}\nn = {int((f[G]==g).sum()):,}' for g in GROUPS]);ax.tick_params(axis='x',labelsize=9)
        for t,g in zip(ax.get_xticklabels(),GROUPS):t.set_color(C[g])
        b.tidy(ax,'y')
    note(fig,'Filled profiles: weighted histograms in 36 equal log₁₀ bins, normalized within each group; widths share one scale\nwithin each panel. Boxes: frozen weighted IQR; central lines and top labels: medians. Thin lines: observed min–max.\nAll evaluable records and extremes are retained. These are distributions, not confidence intervals. Short-interval\nAER/CER use the same urine analytes and flow; they are not independent diagnostic standards or 24-hour excretion.',y=.030)
    save(fig,'Candidate_05_Urine_Flow')
    pd.DataFrame(profiles).to_csv(ROOT/'Data/Flow_weighted_histogram_profiles.csv',index=False)

def main():
    check_inputs();b.save=save
    classification_v3()
    from v3_mortality import make_mortality
    make_mortality()
    components_v3();v2.absolute_risk_v2();b.sensitivity();v2.a2_supp()
    b.flow()
    for ext in ['png','pdf']:(b.MAIN/f'Candidate_05_Urine_Flow.{ext}').replace(b.SUPP/f'Supplement_03_Flow_Means_CI.{ext}')
    flow_v3()
    unique={(r['_workbook'],r['_sheet'],r['_row']):r for r in b.USED}
    obj={'version':'v3 review candidates; no new model fits','input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'Inputs').iterdir() if p.is_file()},'source_rows':list(unique.values()),'audit':AUDIT}
    (ROOT/'Data/figure_source_rows.json').write_text(json.dumps(obj,ensure_ascii=False,indent=2))
    with (ROOT/'Data/figure_source_cells.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.writer(f);w.writerow(['Workbook','Sheet','Excel row','Field','Frozen value'])
        for r in unique.values():
            for k,v in r.items():
                if not k.startswith('_'):w.writerow([r['_workbook'],r['_sheet'],r['_row'],k,v])
    print(json.dumps({'output':str(ROOT),'png_count':len(list(ROOT.rglob('*.png'))),'audit':AUDIT},ensure_ascii=False))

if __name__=='__main__':main()
