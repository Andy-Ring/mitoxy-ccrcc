"""Manuscript figures rendered from numerical outputs, without PDF overlays."""
import numpy as np
import pandas as pd
from scipy.stats import t
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .data import ROOT
from .plot_style import style, forest, C, box, ptext, letter

def finish(fig,number,sources):
    out=ROOT/'figures';out.mkdir(exist_ok=True)
    for ax,label in zip([a for a in fig.axes if a.get_label()!='<colorbar>'],'ABCDEFGHIJKLMNOPQRSTUVWXYZ'):
        letter(ax,label)
    for ext in ['pdf','png']:fig.savefig(out/f'Figure_{number}.{ext}',dpi=300,bbox_inches='tight')
    plt.close(fig)
    import json
    (out/f'Figure_{number}.sources.json').write_text(json.dumps(sources,indent=2)+'\n')


def source():
    style();d=ROOT/'data/processed/figure_7_clean';fig,axs=plt.subplots(3,2,figsize=(12,12),layout='constrained')
    p=pd.read_csv(d/'panel_b_pca.csv');groups=['Media','LPS','CS','LPS+CS','LPS+CS+Torin']
    for i,g in enumerate(groups):
        x=p[p.condition.eq(g)];axs[0,0].scatter(x.pc1,x.pc2,label=g,color=[C[0],C[1],C[2],'#CC6677',C[3]][i],s=22)
    axs[0,0].set(xlabel=f'PC1 ({p.pc1_variance_pct.iloc[0]:.1f}%)',ylabel=f'PC2 ({p.pc2_variance_pct.iloc[0]:.1f}%)');axs[0,0].legend(fontsize=8)
    means=pd.read_csv(d/'panel_c_score_trajectories.csv')
    for j,name in enumerate(['Equal4','Core3','SLC7A11']):
        q=means[means.score.eq(name)].set_index('condition').loc[groups]
        axs[0,1].errorbar(np.arange(5)+.08*(j-1),q['mean'],yerr=q['sem']*t.ppf(.975,q.n-1),marker='o',color=[C[0],C[2],C[1]][j],label=name,capsize=3,markersize=4)
    axs[0,1].set(xticks=range(5),xticklabels=['Media','LPS','CS','LPS + CS','LPS + CS\n+ Torin'],ylabel='Standardized score');axs[0,1].legend(fontsize=8,ncol=3)
    effects=pd.read_csv(d/'panel_d_contrasts.csv').rename(columns={'effect':'beta'})
    effects=effects.assign(_order=~effects.contrast.str.contains('interaction',case=False)).sort_values('_order',kind='stable').drop(columns='_order')
    forest(axs[1,0],effects,[r.contrast.replace('LPS-by-CS interaction','LPS × CS')+': '+r.score for r in effects.itertuples()],colors=[C[0] if 'interaction' in r.contrast.lower() else C[1] for r in effects.itertuples()],xlabel='Effect on standardized score (95% CI)')
    null=pd.read_csv(d/'panel_e_specificity.csv');axs[1,1].barh(null.contrast,null.percentile,color=C[:2])
    axs[1,1].set(xlim=(0,112),xlabel='Matched-module percentile');axs[1,1].invert_yaxis()
    for i,r in null.iterrows():box(axs[1,1],r.percentile+2,i,f'{r.percentile:.2f}%',C[i])
    axs[1,1].text(.02,-.20,f'Joint exceedance: {null.joint_empirical_p.iloc[0]:.2%}',transform=axs[1,1].transAxes,va='top')
    g=pd.read_csv(ROOT/'results/source_de/GSEA_all.csv');g=g[g.exclude_core.eq(False)].copy()
    for col in ['NES','FDR q-val']:g[col]=pd.to_numeric(g[col])
    names=['HALLMARK_INTERFERON_GAMMA_RESPONSE','HALLMARK_INFLAMMATORY_RESPONSE',
           'HALLMARK_TNFA_SIGNALING_VIA_NFKB','HALLMARK_APOPTOSIS','REACTOME_PYROPTOSIS',
           'HALLMARK_MTORC1_SIGNALING','HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY',
           'HALLMARK_OXIDATIVE_PHOSPHORYLATION','HALLMARK_HYPOXIA']
    labels=['IFN-γ response','Inflammatory response','TNFα/NF-κB','Apoptosis','Pyroptosis','mTORC1','ROS','OXPHOS','Hypoxia']
    contrasts=['LPS_vs_Media','Interaction','Torin_vs_Combined']
    selected=g[g.Term.isin(names)&g.contrast.isin(contrasts)]
    matrix=selected.pivot(index='Term',columns='contrast',values='NES').reindex(index=names,columns=contrasts)
    from matplotlib.colors import LinearSegmentedColormap
    heatmap_cmap=LinearSegmentedColormap.from_list('reference_diverging',['#5C83B6','#F2F4F7','#B35C61'])
    ax=axs[2,0];values=matrix.to_numpy();lim=np.nanmax(np.abs(values));im=ax.imshow(values,cmap=heatmap_cmap,vmin=-lim,vmax=lim,aspect='auto')
    ax.grid(False);ax.set(yticks=range(len(names)),yticklabels=labels,xticks=range(3),xticklabels=['LPS','LPS × CS','Torin − combined'])
    ax.tick_params(length=0);fig.colorbar(im,ax=ax,shrink=.7,label='Normalized enrichment score')
    for i,n in enumerate(names):
        for j,c in enumerate(contrasts):
            q=selected.loc[selected.Term.eq(n)&selected.contrast.eq(c),'FDR q-val']
            if len(q):ax.text(j,i,f'{values[i,j]:.2f}'+('*' if q.iloc[0]<.05 else ''),ha='center',va='center',fontsize=7,color='white' if abs(values[i,j])>lim*.55 else C[4])
    q=g[g.contrast.eq('Torin_vs_Combined')&g.NES.lt(0)&g['FDR q-val'].lt(.05)].sort_values('FDR q-val').head(8).copy();q['NES']=-q.NES;q=q.sort_values('NES',ascending=False)
    short={'GOBP_FERROPTOSIS':'Ferroptosis','HALLMARK_MTORC1_SIGNALING':'mTORC1','HALLMARK_OXIDATIVE_PHOSPHORYLATION':'OXPHOS','HALLMARK_PI3K_AKT_MTOR_SIGNALING':'PI3K–AKT–mTOR','HALLMARK_UNFOLDED_PROTEIN_RESPONSE':'Endoplasmic-reticulum UPR','HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY':'ROS'}
    import textwrap
    labels=[textwrap.fill(short.get(n,n.replace('HALLMARK_','').replace('_',' ').capitalize()),25) for n in q.Term]
    ax=axs[2,1];ax.barh(labels,q.NES,color=C[3]);ax.invert_yaxis();ax.set_xlabel('Attenuation NES (combined − Torin)');ax.margins(x=.30)
    for i,(_,r) in enumerate(q.iterrows()):ax.text(r.NES+.03,i,'q < resolution' if r['FDR q-val']==0 else f"q={ptext(r['FDR q-val'])}",va='center',fontsize=7)
    finish(fig,'1',{'scores':str(d.relative_to(ROOT)),'pathways':'results/source_de/GSEA_all.csv','direction':'Panel E reports Torin minus combined; panel F reverses the NES sign to display attenuation'})

def architecture():
    style();fig,axs=plt.subplots(2,2,figsize=(12,9),layout='constrained')
    load=pd.read_csv(ROOT/'results/TRACERx_PC1_loadings.csv');bars=axs[0,0].bar(load.gene,load.loading,color=[C[1] if g=='SLC7A11' else C[0] for g in load.gene]);axs[0,0].set_ylabel('PC1 loading');axs[0,0].bar_label(bars,fmt='%.2f',padding=4,fontsize=9);axs[0,0].margins(y=.15)
    m=pd.read_csv(ROOT/'results/TRACERx_models.csv').set_index('model_id')
    names=['Equal4','PC1_4','BAK1_centered','Core3','BAK1_SLC7A11','BAK1_only']
    for color,(outcome,label) in zip([C[0],C[1],C[3]],[('grade','Grade'),('stage','Stage'),('punctuated.full','Punctuated evolution')]):
        d=m.loc[[f'TRACERx.{outcome}.{n}' for n in names]]
        axs[0,1].plot(range(6),d.beta,'o-',color=color,label=label)
    axs[0,1].set(xticks=range(6),xticklabels=['Equal4','PC1','BAK1-centered','Core3','BAK1+SLC7A11','BAK1']);axs[0,1].tick_params(axis='x',rotation=30);axs[0,1].set_ylabel('Adjusted coefficient');axs[0,1].legend(fontsize=8)
    genes=['BAK1','BAX','BID','SLC7A11','Equal4'];d=m.loc[[f'TRACERx.grade.{g}' for g in genes]]
    forest(axs[1,0],d,genes,colors=[C[0]]*4+[C[2]],xlabel='Grade coefficient per SD (95% CI)')
    keys=[f'TRACERx.joint.{y}.{g}' for y in ['grade','punctuated'] for g in ['BAK1_joint','Rest3']]
    forest(axs[1,1],m.loc[keys],['Grade: BAK1','Grade: Rest3','Punctuated: BAK1','Punctuated: Rest3'],colors=[C[0],C[3],C[0],C[3]],xlabel='Joint coefficient per predictor unit (95% CI)')
    finish(fig,'4',{'models':'results/TRACERx_models.csv','loadings':'results/TRACERx_PC1_loadings.csv','Rest3':'mean of complete-case gene z scores; not subsequently standardized'})

def specificity():
    style();fig,axs=plt.subplots(3,2,figsize=(12,12),layout='constrained');d=ROOT/'data/processed/figure_2_clean'
    p=pd.read_csv(d/'panels_a_b_purity_models.csv')
    for ax,(_,group) in zip(axs[0],p.groupby('analysis',sort=False)):
        forest(ax,group,[f'{r.endpoint}   n={r.n}' for r in group.itertuples()],colors=C[:2],xlabel='Adjusted effect per SD Mitoxy (β)')
    p=pd.read_csv(d/'panel_c_locked_matched_null.csv');axs[1,0].bar(p.endpoint,p.percentile,color=C[:len(p)]);axs[1,0].set(ylabel='Observed-score percentile',ylim=(90,102));axs[1,0].axhline(97.5,color='#ABB5BE',ls='--',lw=.8)
    for i,r in p.iterrows():box(axs[1,0],i,r.percentile+.7,f'{r.percentile:.2f}%\nP = {ptext(r.empirical_p)}',C[i]).set_ha('center')
    p=pd.read_csv(d/'panel_d_global_pc_robustness.csv');axs[1,1].plot(p.adjustment,p.percentile,'o-',color=C[3]);axs[1,1].set(ylabel='Observed-score percentile',ylim=(90,103))
    for i,r in p.iterrows():box(axs[1,1],i,r.percentile+1,f'{r.percentile:.2f}%\nP = {ptext(r.empirical_p)}',C[3]).set_ha('center')
    for ax,outcome in zip(axs[2],['grade','stage']):
        labels=[];rows=[]
        for cohort in ['TCGA','CPTAC']:
            m=pd.read_csv(ROOT/f'results/bulk/{cohort}_models.csv')
            for variant,comparison,label in [('main','base','Base'),('main','plus_8','Eight programs'),('main','plus_ComparatorPC1','Comparator PC1'),('no_overlap','plus_8','Eight, no overlap'),('no_overlap','plus_ComparatorPC1','PC1, no overlap')]:
                rows.append(m[m.variant.eq(variant)&m.comparison.eq(comparison)&m.outcome.eq(outcome)].iloc[0]);labels.append(f'{cohort}: {label}')
        forest(ax,pd.DataFrame(rows),labels,colors=([C[0],C[1],C[2],C[3],C[3]])*2,xlabel=f'{outcome.title()} coefficient per cohort SD (95% CI)')
    finish(fig,'3',{'purity_null':str(d.relative_to(ROOT)),'comparators':'results/bulk','scaling':'Clinical input and score reference populations are specified in docs/METHODS.md'})

def evolution():
    from matplotlib.patches import FancyBboxPatch, Circle
    import json
    style();fig=plt.figure(figsize=(12,14),layout='constrained')
    gs=fig.add_gridspec(5,2,height_ratios=[1.25,2.2,1.4,2,1.5])
    axa=fig.add_subplot(gs[0,:]);axb=fig.add_subplot(gs[1,:]);axc=fig.add_subplot(gs[2,0]);axd=fig.add_subplot(gs[2,1]);axe=fig.add_subplot(gs[3,:]);axf=fig.add_subplot(gs[4,:])
    counts=json.loads((ROOT/'results/TRACERx_sample_summary.json').read_text())
    axa.axis('off');axa.set(xlim=(0,1),ylim=(0,1))
    axa.text(.16,.91,'TRACERx multiregion expression',ha='center',fontsize=9)
    for x,y,label in [(.05,.57,'R1'),(.13,.72,'R2'),(.21,.52,'R3'),(.29,.66,'Rn')]:
        axa.scatter([x],[y],s=420,facecolors='#E3F3EF',edgecolors=C[2],linewidths=.8);axa.text(x,y,label,ha='center',va='center',fontsize=8)
    for x,title,detail,c in [(.40,'Region-level Mitoxy','Mean of component gene z scores',C[0]),(.73,'Patient-level mean','Primary unit for clinical models',C[3])]:
        axa.add_patch(FancyBboxPatch((x,.38),.24,.45,boxstyle='round,pad=.01',fc='#F3F6F8',ec=c,lw=.8))
        axa.text(x+.12,.68,title,ha='center',weight='bold',fontsize=9);axa.text(x+.12,.50,detail,ha='center',fontsize=7)
    for x1,x2 in [(.33,.39),(.65,.72)]:axa.annotate('',xy=(x2,.6),xytext=(x1,.6),arrowprops=dict(arrowstyle='->',color='#9BA7B4'))
    for x,value,label,c in [(.16,counts['complete_case_regions'],'Complete-case regions',C[2]),(.52,counts['complete_case_patients'],'Complete-case patients',C[0]),(.85,counts['multiregion_complete_case_patients'],'Patients with ≥2 regions',C[3])]:
        axa.text(x,.19,str(value),ha='center',weight='bold',color=c,fontsize=13);axa.text(x,.04,label,ha='center',fontsize=8)
    r=pd.read_csv(ROOT/'results/evolution/regions_complete_case.csv');order=r.groupby('Patient').Equal4.mean().sort_values().index
    for i,patient in enumerate(order):
        a=r[r.Patient.eq(patient)];axb.plot([i,i],[a.Equal4.min(),a.Equal4.max()],color='#A8B5C0',lw=.7)
        axb.scatter(np.repeat(i,len(a)),a.Equal4,s=11,color=C[2],alpha=.65,edgecolor='white',lw=.2,zorder=3);axb.scatter(i,a.Equal4.mean(),color=C[3],s=14,zorder=4)
    axb.axhspan(r.Equal4.min(),r.Equal4.max(),color='#F3F6F8',zorder=0)
    axb.set(xlabel='Complete-case patients ordered by mean Mitoxy',ylabel='Regional Mitoxy score',xticks=[])
    axb.text(.99,.97,f"{len(r)} regions / {len(order)} complete-case patients",ha='right',va='top',transform=axb.transAxes,fontsize=8,color=C[3])
    v=pd.read_csv(ROOT/'data/processed/figure_4_clean/panel_c_variance_partition.csv');left=0
    for i,row in v.iterrows():
        axc.barh(0,row.variance_fraction,left=left,color=[C[2],C[1]][i],height=.5)
        axc.text(left+row.variance_fraction/2,0,f'{row.variance_fraction:.0%}\n{row.component}',ha='center',va='center',fontsize=8,color='white');left+=row.variance_fraction
    axc.set(yticks=[],xlabel='Variance share',xlim=(0,1),ylim=(-.6,.8));box(axc,.54,.49,f'ICC = {v.variance_fraction.iloc[0]:.2f}',C[2]);axc.spines['left'].set_visible(False)
    m=pd.read_csv(ROOT/'results/TRACERx_models.csv').set_index('model_id');labels=['Clinical','+ WGI','+ ITH','+ WGI and ITH'];keys=['clinical','clinical_wgii','clinical_ITH','full']
    forest(axd,m.loc[[f'TRACERx.punctuated.{k}.Equal4' for k in keys]],labels,color=C[3],xlabel='Punctuated-versus-other coefficient (95% CI)')
    d=pd.read_csv(ROOT/'results/evolution/ITH_dispersion.csv');d=d[d.genomic_feature.eq('ITH')]
    metrics=list(dict.fromkeys(d.metric));display={'range':'Range','sd':'SD','mad':'MAD','mean_pairwise':'Mean pairwise difference','median_pairwise':'Median pairwise difference'}
    for j,(model,c) in enumerate([('raw',C[0]),('adjusted',C[1])]):
        q=d[d.model.eq(model)].set_index('metric').loc[metrics]
        for i,(_,row) in enumerate(q.iterrows()):axe.errorbar(row.rho,i+(j-.5)*.20,xerr=[[row.rho-row.ci_low],[row.ci_high-row.rho]],fmt='o',color=c,ms=4,capsize=2,label=model.title() if i==0 else None)
    for i in range(len(metrics)):
        if i%2==0:axe.axhspan(i-.45,i+.45,color='#F0F5F7',zorder=0)
    axe.set(yticks=range(len(metrics)),yticklabels=[display.get(k,k.replace('_',' ')) for k in metrics],xlabel='Association with genomic ITH (Spearman ρ)',ylim=(len(metrics)-.5,-.5));axe.axvline(0,color='#AAB5C1',ls='--',lw=.8);axe.legend(ncol=2,loc='upper left');axe.grid(axis='y',visible=False)
    b=pd.read_csv(ROOT/'results/bulk/TRACERx_models.csv');b=b[b.outcome.eq('Mitoxy')&b.comparison.eq('comparator_pc1')]
    forest(axf,b,b.variant.str.replace('_',' ').tolist(),colors=[C[2],C[3]][:len(b)],xlabel='Punctuated effect with comparator PC1 (95% CI)')
    finish(fig,'5',{'regional':'results/evolution','models':'results/TRACERx_models.csv','comparators':'results/bulk/TRACERx_models.csv'})

def clinical():
    """Four cohort cards and four association panels; no score-definition block."""
    import json
    from matplotlib.patches import Rectangle
    style();fig=plt.figure(figsize=(13,11),layout='constrained')
    gs=fig.add_gridspec(3,2,height_ratios=[1.15,1.4,1.3])
    a=fig.add_subplot(gs[0,:]);b=fig.add_subplot(gs[1,0]);c=fig.add_subplot(gs[1,1]);d=fig.add_subplot(gs[2,0]);e=fig.add_subplot(gs[2,1])
    root=ROOT/'data/processed/figure_1_clean'
    grade=pd.read_csv(root/'panel_b_grade_summary.csv');stage=pd.read_csv(root/'panel_c_stage_summary.csv')
    patients=pd.concat([pd.read_csv(root/'panel_b_grade_data.csv'),pd.read_csv(root/'panel_c_stage_data.csv')])
    counts=patients.groupby('cohort').sample_id.nunique()
    tx=json.loads((ROOT/'results/TRACERx_sample_summary.json').read_text())
    a.axis('off');a.set(xlim=(0,1),ylim=(0,1))
    cards=[('TCGA-KIRC','TCGA-KIRC',f"{counts['TCGA-KIRC']} primary tumors",'Bulk RNA-seq',C[0]),
           ('CPTAC','CPTAC ccRCC',f"{counts['CPTAC']} clinically matched tumors",'RNA-seq + proteomics',C[1]),
           ('E-MTAB-1980','E-MTAB-1980',f"{counts['E-MTAB-1980']} tumors",'Expression microarray',C[2]),
           ('TRACERx','TRACERx Renal',f"{tx['primary_regions']} primary regions · {tx['expression_patients']} patients",'Multiregion RNA-seq',C[3])]
    for i,(_,title,count,tag,color) in enumerate(cards):
        x=.02+.51*(i%2);y=.54 if i<2 else .04
        a.add_patch(Rectangle((x,y),.45,.40,fc='white',ec='#D3DCE3',lw=.8))
        a.add_patch(Rectangle((x+.015,y+.04),.009,.32,fc=color,ec='none'))
        a.text(x+.04,y+.29,title,color=color,weight='bold',fontsize=12,va='center')
        a.text(x+.04,y+.18,count,weight='bold',fontsize=10,va='center')
        a.text(x+.04,y+.075,tag,color='#778290',fontsize=8,va='center')
    order=['TCGA-KIRC','CPTAC','E-MTAB-1980','TRACERx']
    for ax,q,label in [(b,grade,'histologic grade'),(c,stage,'stage / T category')]:
        q=q.set_index('cohort').loc[order].reset_index()
        forest(ax,q,[f'{r.cohort}   n={r.n}' for r in q.itertuples()],value='rho',symbol='ρ',colors=C[:4],xlabel=f'Spearman correlation with {label} (ρ)')
    q=pd.read_csv(root/'panel_d_cptac_coherence.csv')
    forest(d,q,q.comparison,value='rho',symbol='ρ',color=C[1],xlabel='Spearman correlation (ρ)')
    q=pd.read_csv(root/'panel_e_tracerx_summary.csv').rename(columns={'q025':'ci_low','q975':'ci_high'})
    forest(e,q,['Histologic grade','Pathologic stage'],value='median',symbol='median ρ',color=C[3],annotate=False,xlabel='Spearman correlation across one-region resamples (ρ)')
    e.set_xlim(min(-.04,q.ci_low.min()-.04),max(.78,q.ci_high.max()+.3))
    for i,r in q.iterrows():box(e,.73,i,f'median ρ = {r["median"]:.2f}\n{100*r.direction_positive_frac:.0f}% positive',C[3],e.get_yaxis_transform())
    finish(fig,'2',{'clinical':str(root.relative_to(ROOT)),'cohort_counts':'union of clinical input patients; TRACERx_sample_summary.json','panel_A':'Cohort cards only'})
