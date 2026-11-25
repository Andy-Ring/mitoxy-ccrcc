"""Supplementary figures computed from the current numerical output tables."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from .data import ROOT
from .figures import forest,finish,C
from .plot_style import style, box, ptext, letter
from .cellular import exact_signed_rank

def comparators():
    style();fig,axs=plt.subplots(2,2,figsize=(13,10),layout='constrained');root=ROOT/'results/bulk'
    d=pd.read_csv(root/'TCGA_main_patient_scores.csv')
    names=['Apoptosis','Ferroptosis','Necroptosis','Pyroptosis','UPRmt','ROS','OxidativeStress','Hypoxia'];rho=[spearmanr(d[['Mitoxy',n]].dropna()).statistic for n in names]
    bars=axs[0,0].barh(names,rho,color=C[0]);axs[0,0].invert_yaxis();axs[0,0].bar_label(bars,fmt='%.2f',padding=3,fontsize=8);axs[0,0].margins(x=.15);axs[0,0].set_xlabel('Spearman correlation with Equal4')
    for ax,cohort in [(axs[0,1],'TCGA'),(axs[1,0],'CPTAC')]:
        m=pd.read_csv(root/f'{cohort}_models.csv');m=m[m.variant.eq('main')&m.outcome.eq('grade')&m.comparison.isin(['plus_'+n for n in names])]
        forest(ax,m,m.comparison.str.removeprefix('plus_').tolist());ax.set_title(cohort+' grade')
    v=pd.read_csv(root/'TCGA_VIF.csv');v=v[v.variant.eq('main')&v.outcome.eq('grade')&v.term.ne('const')]
    axs[1,1].barh(v.term,v.VIF,color='#78899A');axs[1,1].invert_yaxis();axs[1,1].axvline(5,ls='--',color='gray');axs[1,1].set_xlabel('Variance inflation factor')
    finish(fig,'S2',{'tables':'results/bulk','scaling':'cohort-wide gene and score population SD'})

def evolution():
    style();fig,axs=plt.subplots(2,2,figsize=(13,11),layout='constrained')
    root=ROOT/'data/processed/supplementary_figure_s4';d=pd.read_csv(root/'panel_a_gene_regional_variation.csv')
    x=np.arange(len(d));axs[0,0].bar(x-.18,d.icc,.36,label='ICC',color=C[0]);axs[0,0].bar(x+.18,d.median_within_range,.36,label='Median within-patient range',color=C[1]);axs[0,0].set(xticks=x,xticklabels=d.gene);axs[0,0].legend(fontsize=8)
    d=pd.read_csv(root/'panel_b_evolution_omnibus.csv').rename(columns={'effect':'beta'})
    forest(axs[0,1],d,[f'{r.classification}: {r.level}' for r in d.itertuples()],colors=[C[3] if 'condition' in r.classification.lower() else C[2] for r in d.itertuples()]);axs[0,1].tick_params(axis='y',labelsize=7)
    m=pd.read_csv(ROOT/'results/bulk/TRACERx_models.csv')
    d=m[m.variant.eq('main')&m.comparison.eq('clinical_full')&m.outcome.ne('Mitoxy')]
    if d.empty:raise ValueError('Comparator-specific evolutionary models missing')
    forest(axs[1,0],d,d.outcome.tolist(),xlabel='Punctuated − other, adjusted coefficient')
    d=pd.read_csv(ROOT/'results/evolution/reverse_PC1.csv')
    if d.empty:raise ValueError('Reverse-logistic models missing')
    d=d.assign(effect=np.exp(d.beta),ci_low=np.exp(d.ci_low),ci_high=np.exp(d.ci_high));forest(axs[1,1],d,d.variant.str.replace('_',' '),value='effect',null=1,symbol='OR',colors=[C[2],C[3]],xlabel='Odds ratio per SD Mitoxy (95% CI)')
    finish(fig,'S4',{'regional':str(root.relative_to(ROOT)),'comparators':'results/bulk/TRACERx_models.csv'})

def paired(ax,d,group,value):
    w=d.pivot(index='patient',columns=group,values=value)[['Macro_MitoxyLow','Macro_MitoxyHigh']].dropna()
    for _,r in w.iterrows():ax.plot([0,1],r,color='#AAB5C1',lw=1)
    ax.scatter(np.zeros(len(w)),w.iloc[:,0],color=C[0]);ax.scatter(np.ones(len(w)),w.iloc[:,1],color=C[1])
    ax.set(xticks=[0,1],xticklabels=['Low','High'],ylabel={'median_UMI':'Median UMI per macrophage','n_detected':'Detected interactions','sum_lr_means':'Summed interaction magnitude'}[value],xlim=(-.2,1.2))
    ax.set_title(f'Exact paired P={exact_signed_rank(w.iloc[:,1],w.iloc[:,0]):.3g}',fontsize=10)

def communication():
    style();root=ROOT/'results/LIANA';fig=plt.figure(figsize=(14,14),layout='constrained')
    gs=fig.add_gridspec(3,6,height_ratios=[1,1.8,1.1]);ax=[fig.add_subplot(gs[0,0:2]),fig.add_subplot(gs[0,2:4]),fig.add_subplot(gs[0,4:6]),fig.add_subplot(gs[1,0:3]),fig.add_subplot(gs[1,3:6]),fig.add_subplot(gs[2,0:3]),fig.add_subplot(gs[2,3:6])]
    paired(ax[0],pd.read_csv(root/'main/macrophage_subset_qc.csv'),'ccc_group','median_UMI')
    summary=pd.read_csv(root/'main/macro_to_T_per_patient_summary.csv')
    paired(ax[1],summary,'source','n_detected');paired(ax[2],summary,'source','sum_lr_means')
    for a in ax[:3]:a.set_xticklabels(['Mitoxy-low','Mitoxy-high'])
    d=pd.read_csv(root/'main/macro_to_T_interaction_consistency.csv');adj=pd.read_csv(root/'umi_adjusted/macro_to_T_interaction_consistency.csv')
    d=d[d.all_patients_same_direction.astype(bool)].copy()
    d=pd.concat([d[d.median_diff.lt(0)].nsmallest(10,'median_diff'),d[d.median_diff.gt(0)].nlargest(10,'median_diff')]).sort_values('median_diff',ascending=False)
    retained=set(adj.loc[adj.all_patients_same_direction.astype(bool),'pair'])
    bars=ax[3].barh(d.pair.str.replace('^',' → ',regex=False),d.median_diff,color=[C[1] if v>0 else C[0] for v in d.median_diff]);ax[3].tick_params(axis='y',labelsize=7);ax[3].invert_yaxis()
    for bar,name in zip(bars,d.pair):
        if name not in retained:bar.set_hatch('///')
    ax[3].set_xlabel('Median within-patient difference in LR magnitude\n(Mitoxy-high − Mitoxy-low)')
    names=['CXCL16','CXCL10','CXCL9','TGFB1','CD86','LGALS9','CD274','PDCD1LG2','HLA-DRA','SPP1','IL10','APOE','C1QA']
    main=pd.read_csv(root/'main/macrophage_ligands_High_vs_Low.csv').set_index('gene').loc[names]
    adjusted=pd.read_csv(root/'umi_adjusted/macrophage_ligands_High_vs_Low.csv').set_index('gene').loc[names]
    y=np.arange(len(names));ax[4].scatter(main.median_diff,y-.10,color=C[4],s=16,label='Main split');ax[4].scatter(adjusted.median_diff,y+.10,facecolors='white',edgecolors='#8B98A6',s=16,label='Library-size adjusted')
    lo=min(main.median_diff.min(),adjusted.median_diff.min());hi=max(main.median_diff.max(),adjusted.median_diff.max());span=hi-lo
    for i,r in enumerate(main.itertuples()):
        if i%2==0:ax[4].axhspan(i-.45,i+.45,color='#F0F5F7',zorder=0)
        ax[4].text(hi+.07*span,i,f'{r.n_high_greater}/{r.n}; P={ptext(r.p_exact)}',va='center',fontsize=7)
    ax[4].set(yticks=y,yticklabels=names,ylim=(len(names)-.5,-1.5),xlim=(lo-.1*span,hi+.8*span),xlabel='Median difference, log-normalized\n(Mitoxy-high − Mitoxy-low)');ax[4].legend(loc='upper left',bbox_to_anchor=(0,1.12),fontsize=7);ax[4].axvline(0,ls='--',color='#AAB5C1',lw=.8);ax[4].grid(axis='y',visible=False)
    d=pd.read_csv(root/'main/patient_T_programs_vs_macrophage_Equal4.csv')
    for a,col in zip([ax[5],ax[6]],['T_cytotoxic','T_exhaustion']):
        a.scatter(d.macrophage_Equal4,d[col],s=30+120*d.n_T_cells/d.n_T_cells.max(),color=C[3],alpha=.65,edgecolor='white',lw=.5)
        for r in d.itertuples():a.annotate(str(r[1]).removeprefix('SS_'),(r.macrophage_Equal4,getattr(r,col)),xytext=(4,4),textcoords='offset points',fontsize=7)
        r,p=spearmanr(d.macrophage_Equal4,d[col]);a.set(xlabel='Mean macrophage Equal4 (per patient)',ylabel=col.replace('_',' ')+' program',title=f'Spearman ρ={r:.2f}; P={ptext(p)} (n={len(d)})');a.margins(.18)
    finish(fig,'S8',{'tables':'results/LIANA/main','sensitivity':'results/LIANA/umi_adjusted','note':'Up to ten consistently signed interactions per direction, retaining the current display selection; hatching denotes loss of sign consistency after UMI adjustment'})


def treatment():
    style();fig,axs=plt.subplots(2,2,figsize=(13,9),layout='constrained');root=ROOT/'data/processed/figure_6_clean'
    specs=[('panel_a_2d_subgroups.csv','subgroup'),('panel_b_formal_interactions.csv','score'),('panel_c_context_interactions.csv','context_stratum')]
    for ax,(file,label),color in zip(axs.flat,specs,[C[0],C[3],C[2]]):
        d=pd.read_csv(root/file);ratio_forest(ax,d,'hr',d[label].str.replace(' / ',' /\n',regex=False).tolist(),color=color)
    d=pd.read_csv(ROOT/'results/CheckMate025/interaction_models.csv');d=d[d.score.eq('Equal4')]
    ratio_forest(axs[1,1],d,'HR',[f'{r.endpoint} ({r.ratio_type}, n={r.n})' for r in d.itertuples()],color=C[3],symbol='ratio')
    for ax,label in zip(axs.flat,['Avelumab + axitinib vs sunitinib PFS hazard ratio','JAVELIN treatment-interaction hazard ratio','Mitoxy × treatment interaction hazard ratio','Equal4 × treatment interaction ratio (95% CI)']):ax.set_xlabel(label)
    finish(fig,'7',{'JAVELIN':str(root.relative_to(ROOT)),'CheckMate':'results/CheckMate025/interaction_models.csv','source':'Motzer 2020 publisher workbook; recalculated models'})

def ratio_forest(ax,d,value,labels,color=C[0],symbol='HR',colors=None):
    forest(ax,d,labels,value=value,null=1,symbol=symbol,color=color,colors=colors,
           xlabel='Hazard / odds ratio (95% CI)')


def immune_context():
    from scipy.stats import kruskal
    from statsmodels.stats.multitest import multipletests
    style();fig,axs=plt.subplots(3,2,figsize=(13,13),layout='constrained');ax=axs.flat
    d=pd.read_csv(ROOT/'results/GSE159115/patient_celltype_scores.csv')
    w=d.pivot(index='patient',columns='cell_type',values='Equal4')[['Tumor','Macro']].dropna()
    for _,r in w.iterrows():ax[0].plot([0,1],r,color='#B7C1CB',lw=.8,zorder=1)
    ax[0].scatter(np.zeros(len(w)),w.Tumor,color=C[0],s=22,zorder=2);ax[0].scatter(np.ones(len(w)),w.Macro,color=C[1],s=22,zorder=2)
    tests=pd.read_csv(ROOT/'results/GSE159115/macrophage_tumor_paired.csv').set_index('score')
    ax[0].set(xticks=[0,1],xticklabels=['Tumor cells','Macrophages'],ylabel='Equal4 score (pseudobulk)',xlim=(-.2,1.2),title=f"Exact paired P={tests.loc['Equal4','p']:.4f}; BH q={tests.loc['Equal4','q']:.3f}")
    d=pd.read_csv(ROOT/'results/GSE159115/component_logCPM.csv');genes=['BAX','BAK1','BID','SLC7A11'];diffs=[];ps=[]
    for g in genes:
        w=d.pivot(index='patient',columns='cell_type',values=g+'_log2CPM1')[['Macro','Tumor']].dropna();diffs.append(w.Macro-w.Tumor);ps.append(exact_signed_rank(w.Macro,w.Tumor))
    qs=multipletests([tests.loc['Equal4','p']]+ps,method='fdr_bh')[1][1:]
    for i,(g,v,p,q) in enumerate(zip(genes,diffs,ps,qs)):
        ax[1].scatter(i+np.linspace(-.10,.10,len(v)),v,color=C[0],s=22,alpha=.9)
        ax[1].plot([i-.18,i+.18],[v.mean()]*2,color=C[4],lw=1.5)
        ax[1].text(i,v.max()+.25,f'P={p:.3f}\nq={q:.3f}',ha='center',fontsize=7)
    ax[1].set(xticks=range(4),xticklabels=genes,ylabel='Macrophage − tumor log2(CPM+1)');ax[1].margins(y=.25);ax[1].axhline(0,color='#AAB5C1',ls='--',lw=.8)
    d=pd.read_csv(ROOT/'results/immune/TCGA_immune_models.csv')
    labels={'base':'Base model','relative_macrophage':'+ Relative macrophage','tissue_macrophage':'+ Tissue-scaled macrophage','eight_tissue_proxies':'+ Eight tissue-scaled\nimmune proxies','purity':'+ Purity','relative_macrophage_plus_purity':'+ Relative macrophage\n+ purity','tissue_macrophage_plus_purity':'+ Tissue-scaled macrophage\n+ purity'}
    for a,sel,title,color in [(ax[2],~d.model.str.startswith('myeloid_'),'Grade coefficient (95% CI)',C[0]),(ax[3],d.model.str.startswith('myeloid_'),'MYD88/TNF/TLR2/TLR4\ncoefficient (95% CI)',C[1])]:
        q=d[sel];forest(a,q,[f"{labels.get(r.model.removeprefix('myeloid_'),r.model)} (n={r.n})" for r in q.itertuples()],color=color,xlabel=title)
    j=pd.read_csv(ROOT/'results/JAVELIN/analysis_patients.csv');groups=['m1','m2','m3','m4','Other'];data=[j.loc[j.tcga_cluster.eq(g),'Equal4_mean_gene_z'].dropna().values for g in groups]
    v=ax[4].violinplot(data,positions=range(5),showextrema=False,showmedians=True)
    for body,color in zip(v['bodies'],['#A5C7DC','#5288B4','#EBA546','#72B568','#AAAAAA']):body.set_facecolor(color);body.set_edgecolor('white');body.set_alpha(.9)
    v['cmedians'].set_color(C[4]);v['cmedians'].set_linewidth(1.2)
    ax[4].set(xticks=range(5),xticklabels=[f'{g}\nn={len(v)}' for g,v in zip(groups,data)],ylabel='Equal4 score (mean gene-z units)',title=f'Kruskal–Wallis P={ptext(kruskal(*data).pvalue)}')
    d=pd.read_csv(ROOT/'results/JAVELIN/context_correlations.csv');d=d[d.comparison.ne('Myeloid trigger')]
    labels={'TEFF':'T-effector','MCDERMOTT_MYELOID':'Myeloid inflammation','ANGIO':'Angiogenesis'}
    forest(ax[5],d,d.comparison.map(labels),value='rho',symbol='ρ',color=C[3],xlabel='Spearman correlation with Equal4 (95% CI)')
    finish(fig,'6',{'cellular':'results/GSE159115','immune':'results/immune','JAVELIN':'results/JAVELIN'})


def clinical_survival():
    style();fig=plt.figure(figsize=(13,10),layout='constrained');gs=fig.add_gridspec(2,2);a=fig.add_subplot(gs[0,0]);sub=gs[0,1].subgridspec(1,2);b=fig.add_subplot(sub[0,0]);b2=fig.add_subplot(sub[0,1]);c=fig.add_subplot(gs[1,0]);d=fig.add_subplot(gs[1,1])
    root=ROOT/'data/processed/supplementary_figure_s1'
    q=pd.read_csv(root/'panel_a_tracerx_clinical.csv');forest(a,q,q.endpoint.tolist(),value='rho',symbol='ρ',colors=C[:len(q)],xlabel='Spearman correlation')
    q=pd.read_csv(root/'panel_b_emtab_adjusted.csv');v=q[q.effect_type.eq('beta')];forest(b,v,v.endpoint.tolist(),value='effect',colors=C[:len(v)],xlabel='Coefficient per SD')
    v=q[~q.effect_type.eq('beta')];ratio_forest(b2,v,'effect',v.endpoint.tolist(),symbol='OR',colors=C[2:4]);b2.set_xlabel('Odds ratio per SD')
    for ax,name,title in [(c,'emtab','E-MTAB-1980 OS'),(d,'tcga','TCGA-KIRC OS')]:
        q=pd.read_csv(root/f'panel_{"c" if name=="emtab" else "d"}_{name}_survival.csv');ratio_forest(ax,q,'effect',[f'{r.model} (n={r.n})' for r in q.itertuples()],colors=C[:len(q)]);ax.set_title(title);ax.set_xlabel('Hazard ratio per SD')
    # Panel B has two independent axes because beta and OR have different nulls.
    for ax,letter in zip([a,b,b2,c,d],['A','B','','C','D']):ax.text(-.10,1.04,letter,transform=ax.transAxes,weight='bold',fontsize=15)
    out=ROOT/'figures'
    for ext in ['pdf','png']:fig.savefig(out/f'Figure_S1.{ext}',dpi=300,bbox_inches='tight')
    plt.close(fig)

def architecture():
    style();fig=plt.figure(figsize=(13,12),layout='constrained');gs=fig.add_gridspec(3,2,height_ratios=[1.5,1.8,2]);ax=[fig.add_subplot(gs[0,0]),fig.add_subplot(gs[0,1]),fig.add_subplot(gs[1,0]),fig.add_subplot(gs[1,1]),fig.add_subplot(gs[2,:])];root=ROOT/'results/architecture'
    for i,file,label in [(0,'pc1_concordance.csv','comparison'),(1,'emtab_components.csv','comparison'),(2,'inflammatory_coupling.csv','comparison')]:
        q=pd.read_csv(root/file);labels=q[label].tolist() if i!=2 else (q.cohort+' / '+q.score).tolist()
        colors=[C[0]]*len(q) if i==0 else [C[0] if r.score in ['Equal4','Mitoxy'] else C[2] if r.score=='Core3' else C[1] for r in q.itertuples()]
        forest(ax[i],q,labels,value='rho',symbol='ρ',colors=colors,xlabel='Spearman correlation (95% CI)')
    q=pd.read_csv(ROOT/'results/TRACERx_models.csv').set_index('model_id');names=['Equal4','PC1_4','BAK1_centered','Core3','BAK1_SLC7A11','BAK1_only']
    for outcome,color in [('grade',C[0]),('stage',C[1])]:
        d=q.loc[[f'TRACERx.{outcome}.{n}' for n in names]];ax[3].plot(range(6),d.beta,'o-',color=color,ms=4,label=outcome.title())
    ax[3].set(xticks=range(6),xticklabels=['Equal4','PC1','BAK1-centered','Core3','BAK1 + SLC7A11','BAK1 only'],ylabel='Adjusted coefficient');ax[3].tick_params(axis='x',rotation=35);ax[3].legend(ncol=2)
    q=pd.read_csv(ROOT/'results/TCGA_sensitivities/models.csv');q=q[q.analysis.eq('joint_four_gene')]
    forest(ax[4],q,(q.outcome+' — '+q.signature).tolist(),colors=[C[0] if r.outcome.lower()=='grade' else C[3] for r in q.itertuples()],xlabel='Conditional effect per gene SD (95% CI)')
    finish(fig,'S3',{'architecture':'results/architecture','TRACERx':'results/TRACERx_models.csv','TCGA':'results/TCGA_sensitivities/models.csv'})


def spatial():
    import tarfile,json,io
    from PIL import Image
    from matplotlib.colors import TwoSlopeNorm
    style();fig,axs=plt.subplots(3,2,figsize=(12,13),gridspec_kw={'height_ratios':[2,2,1]},layout='constrained');root=ROOT/'results/GSE210041';ds={s:pd.read_csv(root/f'{s}_spatial_scores.csv') for s in ['T1','T7']}
    tissue={}
    for name,prefix in [('T1','GSM6415705_GLMF1'),('T7','GSM6415706_GLMF2')]:
        path=ROOT/'data/raw/GSE210041'/f'{prefix}_spatial.tar.gz'
        if not path.exists():raise FileNotFoundError(f'Missing public histology archive: {path}. Download the GSE210041 inputs first.')
        with tarfile.open(path) as z:
            tissue[name]=(np.asarray(Image.open(io.BytesIO(z.extractfile('spatial/tissue_hires_image.png').read()))),json.load(z.extractfile('spatial/scalefactors_json.json')))
    for j,col in enumerate(['Equal4','Myeloid4']):
        lim=np.quantile(np.abs(np.concatenate([d[col].values for d in ds.values()])),.98)
        for i,(name,d) in enumerate(ds.items()):
            a=axs[i,j];background,scale=tissue[name];factor=scale['tissue_hires_scalef'];a.imshow(background,alpha=.60)
            im=a.scatter(d.pxl_col*factor,d.pxl_row*factor,c=d[col],s=5,cmap='coolwarm',norm=TwoSlopeNorm(vmin=-lim,vcenter=0,vmax=lim),alpha=.70,linewidths=0,rasterized=True)
            a.set_aspect('equal');a.axis('off');fig.colorbar(im,ax=a,shrink=.85,label=name+' / '+col,extend='both')
    q=pd.read_csv(root/'definition_ranges.csv');q=q[q.family.eq('grid')].rename(columns={'median_rho':'beta','q025':'ci_low','q975':'ci_high'});forest(axs[2,0],q,q['sample'].tolist(),symbol='ρ',colors=C[:2],xlabel='Spearman ρ (95% definition range)')
    old=pd.read_csv(root/'archived_summary.csv');new=pd.read_csv(root/'hypoxia_coupling.csv');new=new[new.k.eq(6)];rows=[]
    for name in ['T1','T7']:
        for label,d in [('Archived epithelial',old),('Hypoxia targets',new)]:
            r=d[d['sample'].eq(name)].iloc[0];rows.append([name+' / '+label,f'{int(r.n_tumor)} / {int(r.n_myeloid)}',f'{r.rho:.3f}',ptext(r.p_two_sided)])
    a=axs[2,1];a.axis('off');table=a.table(cellText=rows,colLabels=['Tumor / definition','Tumor / myeloid','ρ','Permutation P'],colWidths=[.41,.24,.11,.24],loc='center',cellLoc='left');table.auto_set_font_size(False);table.set_fontsize(7);table.scale(1,1.7)
    for (i,j),cell in table.get_celld().items():cell.set_edgecolor('white');cell.set_facecolor('#F0F5F7' if i%2 else 'white');cell.set_text_props(color=C[4],weight='bold' if i==0 else 'normal')
    finish(fig,'S6',{'tables':'results/GSE210041','tissue':'data/raw/GSE210041/*_spatial.tar.gz','panel_F':'Archived epithelial versus documented hypoxia-target definition; six nearest myeloid spots'})
