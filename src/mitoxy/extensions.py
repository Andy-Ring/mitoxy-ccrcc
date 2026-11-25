"""Clinical, signature and inflammatory-context analyses from reconstructed inputs."""
import json
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from .data import ROOT
from .analysis import GENES,zscore,model
from .bulk import prepare,read_gmt
from .survival import cox,bh

def correlation(x,y,label):
    d=pd.concat([x.rename('x'),y.rename('y')],axis=1).dropna();n=len(d)
    if n<4:raise ValueError('Insufficient correlation observations: '+label)
    r,p=stats.spearmanr(d.x,d.y);ci=np.tanh(np.arctanh(np.clip(r,-.999999,.999999))+np.array([-1,1])*stats.norm.ppf(.975)/np.sqrt(n-3))
    return dict(comparison=label,n=n,rho=r,p=p,ci_low=ci[0],ci_high=ci[1])

def javelin():
    out=ROOT/'results/JAVELIN';out.mkdir(parents=True,exist_ok=True)
    d=pd.read_csv(ROOT/'data/processed/figure_6_clean/JAVELIN_all_patients.csv',keep_default_na=False)
    for name in ['Mitoxy','Core3','SLC7A11','Myeloid trigger','TEFF','MCDERMOTT_MYELOID','ANGIO','age','pfs','pfs_event','treatment']:
        d[name]=pd.to_numeric(d[name],errors='coerce')
    rows=[correlation(d.Mitoxy,d[n],n) for n in ['Myeloid trigger','TEFF','MCDERMOTT_MYELOID','ANGIO']]
    pd.DataFrame(rows).to_csv(out/'context_correlations.csv',index=False)
    known=d[~d.tcga_cluster.isin(['','NA','nan','Unknown'])].copy()
    g=known.groupby('tcga_cluster').Equal4_mean_gene_z.agg(n='size',median='median',q25=lambda x:x.quantile(.25),q75=lambda x:x.quantile(.75))
    kw=stats.kruskal(*[v.Mitoxy for _,v in known.groupby('tcga_cluster')]);g['omnibus_p']=kw.pvalue;g.to_csv(out/'cluster_summary.csv')
    x=pd.get_dummies(known[['sex','pdl1','tcga_cluster']].astype(str),drop_first=True,dtype=float)
    x['treatment']=known.treatment;x['age']=known.age.rank();x['ANGIO']=known.ANGIO.rank();x=sm.add_constant(x)
    a=sm.OLS(known.Mitoxy.rank(),x).fit().resid;b=sm.OLS(known['Myeloid trigger'].rank(),x).fit().resid
    rho=np.corrcoef(a,b)[0,1];df=len(known)-x.shape[1]-1;p=2*stats.t.sf(abs(rho)*np.sqrt(df/(1-rho**2)),df)
    pd.DataFrame([dict(n=len(known),rho=rho,p=p,df=df,covariates='age, sex, treatment, PD-L1, TCGA cluster, ANGIO')]).to_csv(out/'partial_coupling.csv',index=False)
    # Within-arm Cox models, pooled population-SD scaling and a pooled median.
    d['Equal4_SD']=zscore(d.Mitoxy);d['Equal4_high']=(d.Mitoxy>d.Mitoxy.median()).astype(int);rows=[]
    for label,frame in [('clinical',d),('clinical_cluster',known.assign(Equal4_SD=d.Equal4_SD,Equal4_high=d.Equal4_high))]:
        for arm,sub in frame.groupby('treatment'):
            for focal in ['Equal4_SD','Equal4_high']:
                cats=['sex','pdl1']+(['tcga_cluster'] if label.endswith('cluster') else [])
                a=sub.dropna(subset=['pfs','pfs_event','age',focal]).copy();x=pd.get_dummies(a[cats].astype(str),drop_first=True,dtype=float);x.insert(0,focal,a[focal]);x['Age_per_10y']=a.age/10
                x=x.loc[:,x.nunique()>1];fit=cox(a.pfs,a.pfs_event,x)
                rows.append(dict(model=label,arm='Avelumab + axitinib' if arm==1 else 'Sunitinib',contrast=focal,n=len(a),events=int(a.pfs_event.sum()),HR=fit['HR'][0],ci_low=fit['ci_low'][0],ci_high=fit['ci_high'][0],p=fit['p'][0]))
    m=pd.DataFrame(rows)
    for label,ix in m.groupby('model').groups.items():m.loc[ix,'q']=bh(m.loc[ix,'p'])
    m.to_csv(out/'within_arm_PFS.csv',index=False)
    d.to_csv(out/'analysis_patients.csv',index=False)

def signatures():
    out=ROOT/'results/signatures';out.mkdir(parents=True,exist_ok=True)
    genes=json.loads((ROOT/'config/signature45.json').read_text())['genes'];rows=[];cor=[];coverage=[]
    for cohort in ['TCGA','CPTAC','TRACERx']:
        e,meta,cov,patient=prepare(cohort);gz=(e-e.mean())/e.std(ddof=0).replace(0,np.nan);scores=pd.read_csv(ROOT/f'results/bulk/{cohort}_main_patient_scores.csv',index_col=0)
        for label,gs in [('Equal45',genes),('Equal45_no_overlap',[g for g in genes if g not in GENES])]:
            use=[g for g in gs if g in gz];missing=sorted(set(gs)-set(use))
            if missing:raise ValueError(f'{cohort}: missing signature genes {missing}')
            s=gz[use].mean(axis=1)
            if patient is not None:s=s.groupby(patient.reindex(s.index)).mean()
            s=zscore(s);d=scores.join(s.rename(label));cor.append(dict(cohort=cohort,**correlation(d.Mitoxy,d[label],label)))
            coverage.append(dict(cohort=cohort,signature=label,genes=';'.join(use),n_genes=len(use)))
            for y in ['grade','stage']:
                common=d.dropna(subset=[y,label,'Mitoxy']+cov).copy();common[label]=zscore(common[label])
                for kind,preds,focal in [('signature_only',[label]+cov,label),('joint_Mitoxy',['Mitoxy',label]+cov,'Mitoxy')]:
                    rows.append(dict(cohort=cohort,signature=label,analysis=kind,**model(common,y,preds,focal,f'{cohort}.{y}.{label}.{kind}')[0]))
        d.to_csv(out/f'{cohort}_scores.csv')
    pd.DataFrame(rows).to_csv(out/'equal45_models.csv',index=False);pd.DataFrame(cor).to_csv(out/'equal45_correlations.csv',index=False);pd.DataFrame(coverage).to_csv(out/'equal45_coverage.csv',index=False)

def tcga_sensitivities():
    out=ROOT/'results/TCGA_sensitivities';out.mkdir(parents=True,exist_ok=True)
    e,meta,cov,_=prepare('TCGA');d=pd.read_csv(ROOT/'results/bulk/TCGA_main_patient_scores.csv',index_col=0);gz=(e-e.mean())/e.std(ddof=0).replace(0,np.nan);rows=[]
    gmt=read_gmt(ROOT/'data/raw/genesets/h.all.v2024.1.Hs.symbols.gmt')
    for name,gs in [('E2F',gmt['HALLMARK_E2F_TARGETS']),('G2M',gmt['HALLMARK_G2M_CHECKPOINT']),('E2F_G2M_union',sorted(set(gmt['HALLMARK_E2F_TARGETS'])|set(gmt['HALLMARK_G2M_CHECKPOINT'])))]:
        use=sorted((set(gs)-set(GENES))&set(e));d[name]=zscore(gz[use].mean(axis=1))
    for g in GENES:d[g]=gz[g];d['without_'+g]=zscore(gz[[n for n in GENES if n!=g]].mean(axis=1))
    for y in ['grade','stage']:
        common=d.dropna(subset=[y,'Mitoxy']+cov+['E2F','G2M','E2F_G2M_union'])
        for sig in ['E2F','G2M','E2F_G2M_union']:
            rows.append(dict(analysis='proliferation',signature=sig,**model(common,y,['Mitoxy']+cov+[sig],'Mitoxy',f'TCGA.{y}.plus_{sig}')[0]))
        for g in GENES:
            rows.append(dict(analysis='leave_one_out',signature='without_'+g,**model(d,y,['without_'+g]+cov,'without_'+g,f'TCGA.{y}.without_{g}')[0]))
            rows.append(dict(analysis='joint_four_gene',signature=g,**model(d,y,GENES+['age','male'],g,f'TCGA.{y}.joint.{g}')[0]))
    pd.DataFrame(rows).to_csv(out/'models.csv',index=False)
    from .clearcode import predict,xn
    missing=set(xn.index)-set(e)
    if missing:raise ValueError('Unmapped ClearCode34 features: '+str(missing))
    calls=predict(e.loc[:,xn.index].T);calls.to_csv(out/'ClearCode34_calls.csv',index_label='patient')
    d=d.join(calls);d['ccB']=(d.ClearCode34=='ccB').astype(int)
    rows=[]
    for y in ['grade','stage']:
        for label,adjust in [('base',cov),('subtype',cov+['ccB']),('subtype_replaces_PC1',[v for v in cov if v!='PC1']+['ccB'])]:
            rows.append(dict(analysis=label,**model(d,y,['Mitoxy']+adjust,'Mitoxy',f'TCGA.{y}.{label}')[0]))
    pd.DataFrame(rows).to_csv(out/'ClearCode34_adjusted_models.csv',index=False)
    d['Equal4_mean_gene_z']=gz[GENES].mean(axis=1)
    summary=d.groupby('ClearCode34').Mitoxy.agg(n='size',median='median')
    a=d.loc[d.ClearCode34.eq('ccB'),'Mitoxy'];b=d.loc[d.ClearCode34.eq('ccA'),'Mitoxy'];test=stats.mannwhitneyu(a,b,alternative='two-sided');summary['p_MannWhitney']=test.pvalue;summary['rank_biserial_ccB_minus_ccA']=2*test.statistic/(len(a)*len(b))-1
    summary.to_csv(out/'ClearCode34_score_summary.csv')
    from .clearcode import benchmark
    benchmark(calls).to_csv(out/'ClearCode34_benchmark.csv',index=False)
    d.to_csv(out/'analysis_patients.csv')

def architecture():
    import re
    out=ROOT/'results/architecture';out.mkdir(parents=True,exist_ok=True)
    d=pd.read_csv(ROOT/'data/derived/TRACERx_analysis_patients.csv',index_col=0)
    pd.DataFrame([correlation(d.Equal4,d.PC1_4,'Equal4 vs PC1')]).to_csv(out/'pc1_concordance.csv',index=False)
    base=ROOT/'data/raw/e_mtab_1980_public';e=pd.read_csv(base/'ccRCC_exp_log_quantile_normalized.txt',sep='\t',index_col=0)
    probes={'BAK1':['A_23_P145357'],'BAX':['A_33_P3341144'],'BID':['A_24_P187948','A_33_P3227041'],'SLC7A11':['A_33_P3242623']}
    samples=[c for c in e if re.fullmatch(r'ccRCC-\d+',str(c),re.I)];g=pd.DataFrame({gene:e.loc[ids,samples].mean(axis=0) for gene,ids in probes.items()});gz=zscore(g)
    d=pd.DataFrame(dict(Equal4=gz.mean(axis=1),Core3=gz[GENES[:3]].mean(axis=1),SLC7A11=gz.SLC7A11))
    meta=pd.read_excel(base/'Sato_2013_Supplementary_Table_1_patient_characteristics.xlsx',header=1,dtype=str);meta.columns=meta.columns.str.strip();meta=meta.set_index('sample ID')
    d['grade']=pd.to_numeric(meta['Fuhrman grade'],errors='coerce');d['T category']=pd.to_numeric(meta['Stage at diagnosis'].str.extract(r'T([1-4])',expand=False),errors='coerce')
    rows=[]
    for sig in ['Equal4','Core3','SLC7A11']:
        for endpoint in ['grade','T category']:rows.append(dict(score=sig,endpoint=endpoint,**correlation(d[sig],d[endpoint],sig+' / '+endpoint)))
    pd.DataFrame(rows).to_csv(out/'emtab_components.csv',index=False)
    j=pd.read_csv(ROOT/'results/JAVELIN/analysis_patients.csv');rows=[]
    for sig in ['Mitoxy','Core3','SLC7A11']:rows.append(dict(cohort='JAVELIN',score=sig,**correlation(j[sig],j['Myeloid trigger'],sig)))
    e=pd.read_csv(ROOT/'data/raw/CPTAC/HS_CPTAC_CCRCC_RNAseq_fpkm_log2_Tumor.cct',sep='\t',index_col=0).groupby(level=0).mean().T
    gz=zscore(e[GENES+['MYD88','TNF','TLR2','TLR4']]);trigger=gz[['MYD88','TNF','TLR2','TLR4']].mean(axis=1)
    for sig,genes in [('Mitoxy',GENES),('Core3',GENES[:3]),('SLC7A11',['SLC7A11'])]:rows.append(dict(cohort='CPTAC',score=sig,**correlation(gz[genes].mean(axis=1),trigger,sig)))
    pd.DataFrame(rows).to_csv(out/'inflammatory_coupling.csv',index=False)

def membership():
    import gseapy
    out=ROOT/'results/signatures';out.mkdir(parents=True,exist_ok=True)
    genes=json.loads((ROOT/'config/signature45.json').read_text())['genes'];rows=[]
    for cohort in ['TCGA','CPTAC','TRACERx']:
        e,meta,cov,patient=prepare(cohort);cov=[c for c in cov if c!='PC1']
        gs={'Published45':genes,'Published45_no_overlap':[g for g in genes if g not in GENES]}
        res=gseapy.ssgsea(data=e.T,gene_sets=gs,outdir=None,min_size=5,max_size=500,permutation_num=0,sample_norm_method='rank',correl_norm_type='rank',threads=1,seed=20261110,no_plot=True)
        score=res.res2d.pivot(index='Name',columns='Term',values='ES').astype(float)
        if patient is not None:score=score.groupby(patient.reindex(score.index)).mean()
        score=zscore(score);d=pd.read_csv(ROOT/f'results/bulk/{cohort}_main_patient_scores.csv',index_col=0).join(score)
        d.to_csv(out/f'{cohort}_ssGSEA_patient_scores.csv')
        for y in ['grade','stage']:
            common=d.dropna(subset=[y,'Mitoxy']+cov+list(gs)).copy();common['Mitoxy']=zscore(common.Mitoxy)
            for sig in gs:
                common[sig]=zscore(common[sig])
                for label,extra in [('base',[]),('joint',[sig])]:rows.append(dict(cohort=cohort,signature=sig,analysis=label,**model(common,y,['Mitoxy']+cov+extra,'Mitoxy',f'{cohort}.{y}.{sig}.{label}')[0]))
    pd.DataFrame(rows).to_csv(out/'ssGSEA45_models.csv',index=False)

def immune_subtypes():
    out=ROOT/'results/subtypes';out.mkdir(parents=True,exist_ok=True);source=ROOT/'data/raw/CPTAC/Clark2019_TableS7.xlsx';rows=[]
    tcga=pd.read_excel(source,sheet_name='TCGA classification',header=4)
    tcga['patient']=tcga.submitter_id.astype(str).str.replace('.','-',regex=False).str[:12]
    labels={'TCGA':tcga.set_index('patient').CPTAC_Immune_subtype}
    mapping=pd.read_csv(next((ROOT/'data/raw/CPTAC-package').glob('*aliquot_to_patient*')),sep='\t').set_index('aliquot_ID').patient_ID
    cptac=pd.read_excel(source,sheet_name='xCell Signatures',header=3,index_col=0).loc['Immune Group'];cptac.index=cptac.index.map(mapping);cptac=cptac[cptac.index.notna()&~cptac.index.str.endswith(('.N','.C'),na=False)]
    if cptac.groupby(level=0).nunique().max()>1:raise ValueError('Conflicting immune subtype assignments')
    labels['CPTAC']=cptac.groupby(level=0).first()
    for cohort,subtype in labels.items():
        e,meta,cov,_=prepare(cohort);cov=[c for c in cov if c!='PC1'];d=pd.read_csv(ROOT/f'results/bulk/{cohort}_main_patient_scores.csv',index_col=0).join(subtype.rename('immune_subtype'),how='inner').dropna(subset=['immune_subtype'])
        x=pd.get_dummies(d.immune_subtype,drop_first=True,dtype=float).add_prefix('subtype_');d=d.join(x);d.to_csv(out/f'{cohort}_patient_subtypes.csv')
        for y in ['grade','stage']:
            cc=d.dropna(subset=[y,'Mitoxy']+cov+list(x)).copy();cc.Mitoxy=zscore(cc.Mitoxy)
            for label,extra in [('base',[]),('subtype_adjusted',list(x))]:rows.append(dict(cohort=cohort,analysis=label,**model(cc,y,['Mitoxy']+cov+extra,'Mitoxy',cohort+'.'+y+'.'+label)[0]))
    pd.DataFrame(rows).to_csv(out/'immune_subtype_models.csv',index=False)

def component_models():
    out=ROOT/'results/architecture';out.mkdir(parents=True,exist_ok=True);rows=[]
    for cohort in ['TCGA','CPTAC','TRACERx']:
        e,meta,cov,patient=prepare(cohort);cov=[c for c in cov if c!='PC1'];gz=zscore(e[GENES]);sc=gz.copy();sc['Core3']=gz[GENES[:3]].mean(axis=1);sc['Equal4']=gz.mean(axis=1)
        if patient is not None:sc=sc.groupby(patient.reindex(sc.index)).mean()
        d=meta[[c for c in ['grade','stage']+cov if c in meta]].join(zscore(sc))
        for y in ['grade','stage']:
            cc=d.dropna(subset=[y]+cov+list(sc))
            for sig in sc:rows.append(dict(cohort=cohort,score=sig,**model(cc,y,[sig]+cov,sig,f'{cohort}.{y}.{sig}')[0]))
    base=ROOT/'data/raw/e_mtab_1980_public';raw=pd.read_csv(base/'ccRCC_exp_log_quantile_normalized.txt',sep='\t',index_col=0)
    probes={'BAK1':['A_23_P145357'],'BAX':['A_33_P3341144'],'BID':['A_24_P187948','A_33_P3227041'],'SLC7A11':['A_33_P3242623']}
    m=pd.read_excel(base/'Sato_2013_Supplementary_Table_1_patient_characteristics.xlsx',header=1,dtype=str).set_index('sample ID');samples=raw.columns.intersection(m.index)
    e=pd.DataFrame({g:raw.loc[ids,samples].mean(axis=0) for g,ids in probes.items()});gz=zscore(e);sc=gz.copy();sc['Core3']=gz[GENES[:3]].mean(axis=1);sc['Equal4']=gz.mean(axis=1);d=zscore(sc)
    d['age']=pd.to_numeric(m.Age,errors='coerce');d['male']=m.Sex.map({'M':1.,'F':0.});d['grade']=pd.to_numeric(m['Fuhrman grade'],errors='coerce');d['T_category']=pd.to_numeric(m['Stage at diagnosis'].str.extract(r'T([1-4])',expand=False),errors='coerce')
    for y in ['grade','T_category']:
        cc=d.dropna(subset=[y,'age','male']+list(sc))
        for sig in sc:rows.append(dict(cohort='E-MTAB-1980',score=sig,**model(cc,y,[sig,'age','male'],sig,f'EMTAB.{y}.{sig}')[0]))
    pd.DataFrame(rows).to_csv(out/'common_patient_components.csv',index=False)
