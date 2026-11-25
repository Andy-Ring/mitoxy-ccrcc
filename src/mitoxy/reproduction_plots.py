"""New panels render exclusively from computed patient/spot-level results."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from .data import ROOT
from .cellular import GENES,exact_signed_rank
from statsmodels.stats.multitest import multipletests

from .plot_style import style, BLUE, ORANGE, TEAL, PURPLE, INK, letter

def save(fig,name):
    out=ROOT/'figures';out.mkdir(parents=True,exist_ok=True)
    for ext in ['pdf','png']:
        fig.savefig(out/(name+'.'+ext),dpi=300,bbox_inches='tight')
    plt.close(fig)


def cellular():
    style();out=ROOT/'results/GSE159115'
    d=pd.read_csv(out/'patient_celltype_scores.csv');tests=pd.read_csv(out/'macrophage_tumor_paired.csv')
    components=pd.read_csv(out/'component_logCPM.csv')
    fig,axs=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    w=d.pivot(index='patient',columns='cell_type',values='Equal4')[['Tumor','Macro']].dropna()
    for _,r in w.iterrows():axs[0,0].plot([0,1],r,color='#AAB5C1',lw=1,zorder=1)
    axs[0,0].scatter(np.zeros(len(w)),w.Tumor,color=BLUE,zorder=2)
    axs[0,0].scatter(np.ones(len(w)),w.Macro,color=ORANGE,zorder=2)
    t=tests.loc[tests.score.eq('Equal4')].iloc[0]
    axs[0,0].set(xticks=[0,1],xticklabels=['Tumor cells','Macrophages'],ylabel='Equal4 pseudobulk score',xlim=(-.25,1.25))
    ordered=w.Macro.sort_values();positions=ordered.to_numpy().copy();gap=(w.to_numpy().max()-w.to_numpy().min())*.045
    for i in range(1,len(positions)):positions[i]=max(positions[i],positions[i-1]+gap)
    for (patient,value),y in zip(ordered.items(),positions):axs[0,0].annotate(patient,xy=(1,value),xytext=(1.10,y),fontsize=7,va='center',arrowprops=dict(arrowstyle='-',color='#AAB5C1',lw=.5))
    axs[0,0].set_xlim(-.2,1.45)
    diff=[]
    for g in GENES:
        a=components.pivot(index='patient',columns='cell_type',values=g+'_log2CPM1')[['Macro','Tumor']].dropna()
        diff.append(dict(gene=g,mean=(a.Macro-a.Tumor).mean(),p=exact_signed_rank(a.Macro,a.Tumor),n=len(a)))
    diff=pd.DataFrame(diff);diff['q']=multipletests([t.p]+diff.p.tolist(),method='fdr_bh')[1][1:]
    diff=diff.set_index('gene').loc[['BAX','BAK1','BID','SLC7A11']].reset_index()
    colors=[BLUE,BLUE,BLUE,ORANGE]
    axs[0,1].bar(diff.gene,diff['mean'],color=colors,width=.65)
    for i,r in diff.iterrows():axs[0,1].text(i,r['mean']+.08,f'P={r.p:.3f}\nq={r.q:.3f}',ha='center',fontsize=8)
    axs[0,1].set(ylabel='Macrophage − tumor, log2(CPM+1)');axs[0,1].margins(y=.3)
    labels=['Equal4']+diff.gene.tolist();qs=np.r_[t.q,diff.q.values]
    bars=axs[1,0].barh(labels,-np.log10(qs),color=PURPLE,height=.7);axs[1,0].invert_yaxis();axs[1,0].margins(x=.3)
    for i,(v,q) in enumerate(zip(-np.log10(qs),qs)):axs[1,0].text(v+.03,i,f'q={q:.3f}',va='center',fontsize=8)
    axs[1,0].axvline(-np.log10(.05),ls='--',lw=.8,color='#929CA6');axs[1,0].set(xlabel='−log10(BH q)')
    means=d.groupby('cell_type').Equal4.mean().sort_values(ascending=False)
    axs[1,1].barh(means.index,means.values,color=[ORANGE if n=='Macro' else BLUE if n=='Tumor' else '#78899A' for n in means.index],height=.65);axs[1,1].invert_yaxis()
    axs[1,1].set(xlabel='Mean patient-level Equal4')
    for ax,label in zip(axs.flat,'ABCD'):letter(ax,label)
    save(fig,'Figure_S5')


def umap():
    style();d=pd.read_csv(ROOT/'results/GSE159115/UMAP_coordinates_scores.csv.gz')
    fig,axs=plt.subplots(1,3,figsize=(14,4.5),layout='constrained')
    for ax,col in zip(axs[:2],['anno','patient']):
        values=sorted(d[col].unique());cmap=plt.get_cmap('tab20')
        for i,label in enumerate(values):
            a=d[d[col].eq(label)];ax.scatter(a.UMAP1,a.UMAP2,s=1.3,color=cmap(i),label=label,linewidths=0,rasterized=True)
        ax.legend(frameon=False,fontsize=6,markerscale=3,loc='upper left',ncol=2 if col=='anno' else 1)
    lim=np.quantile(np.abs(d.Equal4),.98)
    im=axs[2].scatter(d.UMAP1,d.UMAP2,c=d.Equal4,s=1.3,cmap='coolwarm',linewidths=0,norm=TwoSlopeNorm(vmin=-lim,vcenter=0,vmax=lim),rasterized=True)
    fig.colorbar(im,ax=axs[2],shrink=.65,label='Equal4',extend='both')
    for ax,label in zip(axs,'ABC'):
        ax.set(xlabel='UMAP1',ylabel='UMAP2',xticks=[],yticks=[]);ax.grid(False);letter(ax,label)
    save(fig,'Figure_S7')
