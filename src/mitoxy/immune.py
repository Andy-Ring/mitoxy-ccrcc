"""TCGA immune adjustment from public Xena expression and GDC immune fractions.

Reconstructed from the supplied HC3/t-inference analysis. Missing cell fractions
remain missing; they are never filled with zero. All source inputs are public.
"""
import numpy as np
import pandas as pd
from scipy import stats
from .data import ROOT
from .analysis import zscore

def ols_hc3(df, outcome, predictors, focal='Mitoxy'):
    z = df[[outcome] + predictors].dropna().copy()
    y = z[outcome].to_numpy(float)
    X = np.column_stack([np.ones(len(z)), z[predictors].to_numpy(float)])
    inv = np.linalg.inv(X.T @ X)
    beta = inv @ X.T @ y
    resid = y - X @ beta
    h = np.sum((X @ inv) * X, axis=1)
    meat = X.T @ (((resid / (1 - h)) ** 2)[:, None] * X)
    cov = inv @ meat @ inv
    se = np.sqrt(np.diag(cov))
    j = predictors.index(focal) + 1
    df_resid = len(z) - X.shape[1]
    t = beta[j] / se[j]
    p = 2 * stats.t.sf(abs(t), df_resid)
    crit = stats.t.ppf(0.975, df_resid)
    ci = beta[j] + np.array([-1, 1]) * crit * se[j]
    return {'n':len(z),'outcome':outcome,'predictors':'+'.join(predictors),
            'term':focal,'beta':beta[j],'se_HC3':se[j],'p':p,
            'ci_low':ci[0],'ci_high':ci[1],'df_resid':df_resid}

def vif_values(df, columns):
    z=df[columns].dropna().to_numpy(float)
    out=[]
    for j,name in enumerate(columns):
        y=z[:,j]
        X=np.column_stack([np.ones(len(z)),np.delete(z,j,axis=1)])
        coef=np.linalg.lstsq(X,y,rcond=None)[0]
        ssr=np.sum((y-X@coef)**2); sst=np.sum((y-y.mean())**2)
        r2=1-ssr/sst
        out.append({'term':name,'VIF':1/(1-r2),'n':len(z)})
    return pd.DataFrame(out)

def analyze():
    out=ROOT/'results/immune';out.mkdir(parents=True,exist_ok=True)
    d=pd.read_csv(ROOT/'results/bulk/TCGA_main_patient_scores.csv').rename(columns={'grade':'grade_ord'})
    e=pd.read_csv(ROOT/'data/raw/tcga_kirc_public/TCGA.KIRC.sampleMap_HiSeqV2.gz',sep='\t',index_col=0)
    e=e.loc[:,[c for c in e if c[13:15]=='01']];e.columns=e.columns.str[:12]
    if e.columns.duplicated().any():raise ValueError('Duplicate TCGA patients')
    genes=['MYD88','TNF','TLR2','TLR4'];gz=zscore(e.loc[genes].T)
    my=zscore(gz.mean(axis=1)).rename('MyeloidTrigger4')
    d=d.merge(my.rename_axis('patient').reset_index(),on='patient',how='left',validate='one_to_one')
    # Keep primary tumor aliquots; collapse technical/aliquot duplicates to patient mean.
    c=pd.read_csv(ROOT/'data/raw/TCGA-immune/cibersort_relative.tsv',sep='\t')
    c=c[c['CancerType'].eq('KIRC')].copy()
    c['sample']=c['SampleID'].astype(str).str.replace('.','-',regex=False)
    c=c[c['sample'].str.split('-').str[3].str[:2].eq('01')]
    c['patient']=c['sample'].str[:12]
    cellcols=[c for c in c.columns if c.startswith(('B.cells','Plasma.cells','T.cells','NK.cells','Monocytes','Macrophages','Dendritic.cells','Neutrophils'))]
    c=c.groupby('patient')[cellcols].mean().reset_index()
    leuk=pd.read_csv(ROOT/'data/raw/TCGA-immune/estimate_leukocyte.tsv',sep='\t',header=None,names=['CancerType','barcode','leukocyte_fraction'])
    leuk=leuk[leuk['CancerType'].eq('KIRC')].copy()
    leuk=leuk[leuk['barcode'].astype(str).str.split('-').str[3].str[:2].eq('01')]
    leuk['patient']=leuk['barcode'].astype(str).str[:12]
    leuk=leuk.groupby('patient')['leukocyte_fraction'].mean().reset_index()
    d=d.merge(c,on='patient',how='left').merge(leuk,on='patient',how='left')
    d['macrophage_relative']=d[['Macrophages.M0','Macrophages.M1','Macrophages.M2']].sum(axis=1, min_count=3)
    groups={
      'macrophage_tissue':['Macrophages.M0','Macrophages.M1','Macrophages.M2'],
      'CD8_tissue':['T.cells.CD8'],
      'CD4_tissue':['T.cells.CD4.naive','T.cells.CD4.memory.resting','T.cells.CD4.memory.activated','T.cells.follicular.helper','T.cells.regulatory..Tregs.','T.cells.gamma.delta'],
      'B_tissue':['B.cells.naive','B.cells.memory','Plasma.cells'],
      'NK_tissue':['NK.cells.resting','NK.cells.activated'],
      'monocyte_tissue':['Monocytes'],
      'dendritic_tissue':['Dendritic.cells.resting','Dendritic.cells.activated'],
      'neutrophil_tissue':['Neutrophils']}
    for name,cols in groups.items():
        d[name]=d[cols].sum(axis=1,min_count=len(cols))*d['leukocyte_fraction']
    d.to_csv(out/'TCGA_immune_patient_level.csv',index=False)
    covars=['age','male','purity','PC1']
    models=[]
    for label,extra in [('base',[]),('relative_macrophage',['macrophage_relative']),('tissue_macrophage',['macrophage_tissue']),('eight_tissue_proxies',list(groups))]:
        models.append({'model':label,**ols_hc3(d,'grade_ord',['Mitoxy']+covars+extra)})
    # Myeloid models are restricted to the same purity-complete 498-patient set.
    common=d.dropna(subset=['MyeloidTrigger4','purity'])
    for label,extra in [('base',[]),('relative_macrophage',['macrophage_relative']),('purity',['purity']),('relative_macrophage_plus_purity',['macrophage_relative','purity']),('tissue_macrophage',['macrophage_tissue']),('tissue_macrophage_plus_purity',['macrophage_tissue','purity'])]:
        models.append({'model':'myeloid_'+label,**ols_hc3(common,'MyeloidTrigger4',['Mitoxy']+extra)})
    pd.DataFrame(models).to_csv(out/'TCGA_immune_models.csv',index=False)
    vif=[]
    for label,cols in [
      ('myeloid_tissue_macrophage_plus_purity',['Mitoxy','macrophage_tissue','purity']),
      ('grade_eight_tissue_proxies',['Mitoxy']+covars+list(groups))]:
        v=vif_values(d,cols); v.insert(0,'model',label); vif.append(v)
    pd.concat(vif,ignore_index=True).to_csv(out/'TCGA_immune_VIF.csv',index=False)
