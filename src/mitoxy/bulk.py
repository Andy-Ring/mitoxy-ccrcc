"""Canonical death/stress comparisons on patient-level public expression.

Fixed definitions and cohort-wide gene/score standardization follow the supplied
revision code. Base and adjusted models use identical complete-case patients.
"""
import json
import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.decomposition import PCA
from statsmodels.stats.outliers_influence import variance_inflation_factor
from .data import ROOT
from .analysis import GENES,zscore,model

def read_gmt(p):
    return {x[0]:x[2:] for x in (line.rstrip().split('\t') for line in p.read_text().splitlines()) if len(x)>2}

def definitions():
    d=ROOT/'data/raw/genesets';hall=read_gmt(d/'h.all.v2024.1.Hs.symbols.gmt')
    go=read_gmt(d/'c5.go.bp.v2024.1.Hs.symbols.gmt');reactome=read_gmt(d/'c2.cp.reactome.v2024.1.Hs.symbols.gmt')
    custom=json.loads((ROOT/'config/custom_comparators.json').read_text())
    return dict(Apoptosis=hall['HALLMARK_APOPTOSIS'],Ferroptosis=custom['Ferroptosis'],
      Necroptosis=go['GOBP_NECROPTOTIC_PROCESS'],Pyroptosis=reactome['REACTOME_PYROPTOSIS'],UPRmt=custom['UPRmt'],
      ROS=hall['HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY'],OxidativeStress=go['GOBP_RESPONSE_TO_OXIDATIVE_STRESS'],Hypoxia=hall['HALLMARK_HYPOXIA'])

def prepare(cohort):
    if cohort=='TCGA':
        e=pd.read_csv(ROOT/'data/raw/tcga_kirc_public/TCGA.KIRC.sampleMap_HiSeqV2.gz',sep='\t',index_col=0)
        e=e.loc[:,[c for c in e if c[13:15]=='01']];e.columns=e.columns.str[:12]
        d=pd.read_csv(ROOT/'data/processed/figure_2_clean/tcga_kirc_figure_2_analysis_table.csv').set_index('sample_id')
        # Original public notebook standardizes the top 5,000 variable genes before PCA.
        top=e.var(axis=1,ddof=0).nlargest(5000).index
        pcs=PCA(n_components=1,svd_solver='full').fit_transform(zscore(e.loc[top].T))[:,0]
        d['PC1']=pd.Series(pcs,index=e.columns)
        cov=['age','male','purity','PC1'];e=e.T
    elif cohort=='CPTAC':
        e=pd.read_csv(ROOT/'data/raw/CPTAC/HS_CPTAC_CCRCC_RNAseq_fpkm_log2_Tumor.cct',sep='\t',index_col=0).groupby(level=0).mean().T
        d=pd.read_csv(ROOT/'data/derived/CPTAC_analysis_patients.csv',index_col='Case_ID').rename(columns={'grade_num':'grade','stage_num':'stage'})
        cov=['Age_lower_bound','Male']
    elif cohort=='TRACERx':
        import pyreadr
        e=pyreadr.read_r(str(ROOT/'data/raw/TRACERx/tumour_tpm.rds'))[None]
        a=pd.read_csv(ROOT/'data/raw/TRACERx/tx_annotation.tsv',sep='\t');a=a[a.type_collapsed.eq('PRIMARY')].set_index('sample')
        e=e.loc[:,a.index.intersection(e.columns)].T
        d=pd.read_csv(ROOT/'data/derived/TRACERx_all_patients.csv',index_col=0)
        cov=['Age','purity']
        return e,d,cov,a.Patient
    else:raise ValueError(cohort)
    if not e.index.is_unique:raise ValueError('Duplicate patients')
    return e,d,cov,None

def analyze(cohort):
    sets=definitions();e,meta,cov,patient=prepare(cohort)
    if not e.columns.is_unique:raise ValueError('Duplicate expression gene symbols')
    gz=(e-e.mean())/e.std(ddof=0).replace(0,np.nan)
    out=ROOT/'results/bulk';out.mkdir(parents=True,exist_ok=True)
    models=[];coverage=[];vifs=[]
    for variant,remove in [('main',False),('no_overlap',True)]:
        scores=pd.DataFrame(index=e.index)
        scores['Mitoxy']=gz[GENES].mean(axis=1)
        for name,defined in sets.items():
            use=sorted((set(defined)-set(GENES) if remove else set(defined))&set(e.columns))
            if not use:raise ValueError(f'Unmeasured signature: {name}')
            scores[name]=gz[use].mean(axis=1)
            coverage.append(dict(cohort=cohort,variant=variant,signature=name,n_defined=len(set(defined)),n_measured=len(use),genes=';'.join(use)))
        if patient is not None:scores=scores.groupby(patient.reindex(scores.index)).mean()
        scores=zscore(scores)
        # Select metadata columns only so prior derived score columns cannot overwrite new calculations.
        d=meta[[c for c in ['grade','stage','punctuated','wgii','ITH']+cov if c in meta]].join(scores)
        pcdata=d[list(sets)].dropna();pc=PCA(n_components=1,svd_solver='full').fit(pcdata)
        d['Comparator_PC1']=zscore(pd.Series(pc.transform(pcdata)[:,0],index=pcdata.index))
        d.to_csv(out/f'{cohort}_{variant}_patient_scores.csv',index_label='patient')
        for outcome in ['grade','stage']:
            common=d.dropna(subset=[outcome,'Mitoxy']+cov+list(sets)+['Comparator_PC1'])
            for label,extra in [('base',[]),*[(f'plus_{name}',[name]) for name in sets],('plus_8',list(sets)),('plus_ComparatorPC1',['Comparator_PC1'])]:
                models.append(dict(cohort=cohort,variant=variant,comparison=label,**model(common,outcome,['Mitoxy']+cov+extra,'Mitoxy',f'{cohort}.{variant}.{outcome}.{label}')[0]))
            for name in ['Mitoxy']+list(sets):
                models.append(dict(cohort=cohort,variant=variant,comparison='standalone_'+name,**model(common,outcome,[name]+cov,name,f'{cohort}.{variant}.{outcome}.standalone.{name}')[0]))
            x=sm.add_constant(common[['Mitoxy']+list(sets)+cov])
            for j,name in enumerate(x):vifs.append(dict(cohort=cohort,variant=variant,outcome=outcome,term=name,VIF=variance_inflation_factor(x.to_numpy(),j)))
        if cohort=='TRACERx':
            covariates=['Age','purity','grade','stage','wgii','ITH']
            for name in ['Mitoxy']+list(sets):
                for extra,label in [([], 'clinical_full'),(['Comparator_PC1'],'comparator_pc1')]:
                    models.append(dict(cohort=cohort,variant=variant,comparison=label,**model(d,name,['punctuated']+covariates+extra,'punctuated',f'{cohort}.{variant}.evolution.{name}.{label}')[0]))
    pd.DataFrame(models).to_csv(out/f'{cohort}_models.csv',index=False)
    pd.DataFrame(coverage).to_csv(out/f'{cohort}_gene_coverage.csv',index=False)
    pd.DataFrame(vifs).to_csv(out/f'{cohort}_VIF.csv',index=False)
    (out/'definitions.json').write_text(json.dumps(sets,indent=2)+'\n')
