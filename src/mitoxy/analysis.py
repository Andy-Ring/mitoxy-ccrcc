"""Mitoxy cohort analyses. Model statistics are computed, never encoded in figures."""
import hashlib
import json
import numpy as np
import pandas as pd
import pyreadr
import statsmodels.api as sm
from sklearn.decomposition import PCA
from .data import ROOT, sha256

GENES = ['BAK1', 'BAX', 'BID', 'SLC7A11']
CLINICAL = ['purity', 'Age', 'grade', 'stage']
FULL = CLINICAL + ['wgii', 'ITH']

def zscore(x):
    sd = x.std(ddof=0)
    if np.any(~np.isfinite(sd)) or np.any(sd <= 0):
        raise ValueError('Cannot standardize a constant or invalid variable')
    return (x-x.mean())/sd

def model(data, outcome, predictors, term, model_id, use_t=False):
    d = data[list(dict.fromkeys([outcome]+predictors))].dropna()
    x = sm.add_constant(d[predictors], has_constant='add')
    if np.linalg.matrix_rank(x) != x.shape[1]:
        raise ValueError(f'Rank-deficient model: {model_id}')
    fit = sm.OLS(d[outcome], x).fit(cov_type='HC3', use_t=use_t)
    ci = fit.conf_int().loc[term]
    row = dict(model_id=model_id, outcome=outcome, term=term, covariates=';'.join(predictors),
               n=len(d), beta=float(fit.params[term]), se_HC3=float(fit.bse[term]),
               ci_low=float(ci.iloc[0]), ci_high=float(ci.iloc[1]), p=float(fit.pvalues[term]),
               inference='HC3 t' if use_t else 'HC3 normal', df_resid=float(fit.df_resid),
               patient_sha256=hashlib.sha256('\n'.join(sorted(map(str, d.index))).encode()).hexdigest())
    return row, fit, d

def tracerx():
    raw = ROOT/'data/raw/TRACERx'
    expression = pyreadr.read_r(str(raw/'tumour_tpm.rds'))[None]
    annotations = pd.read_csv(raw/'tx_annotation.tsv', sep='\t')
    clinical = pd.read_csv(raw/'TRACERx_s1_1_clinical.txt', sep='\t', skiprows=1, index_col='Subject')
    ann = annotations.loc[annotations.type_collapsed.eq('PRIMARY')].set_index('sample')
    ann = ann.loc[ann.index.intersection(expression.columns)]
    if not ann.index.is_unique or not clinical.index.is_unique:
        raise ValueError('Duplicate sample or clinical patient identifier')
    if ann.groupby('Patient').Condition.nunique().max() != 1:
        raise ValueError('Inconsistent evolutionary state across patient regions')
    # Canonical harmonized score: gene z across primary regions, then region score,
    # patient mean, final patient z across all expression-eligible patients.
    regional_z = zscore(expression.loc[GENES, ann.index].T)
    regions = regional_z.copy()
    regions['Equal4'] = regional_z.mean(axis=1)
    patient_unscaled = regions.groupby(ann.Patient).mean()
    patient = zscore(patient_unscaled)
    patient = patient.join(ann.groupby('Patient')[['purity', 'wgii', 'ITH']].mean())
    patient = patient.join(clinical[['Age', 'Fuhrman Grade', 'Overall Stage']])
    patient['Age'] = pd.to_numeric(patient.Age, errors='coerce')
    patient['grade'] = pd.to_numeric(patient['Fuhrman Grade'], errors='coerce')
    patient['stage'] = patient['Overall Stage'].map({'I':1, 'II':2, 'III':3, 'IV':4})
    patient['punctuated'] = ann.groupby('Patient').Condition.first().eq('Punctuated').astype(int)
    d = patient.dropna(subset=FULL).copy()
    # Retain historical architecture alternatives as explicitly named sensitivities.
    # Their patient-component weights are distinct from canonical regional Equal4.
    pc = PCA(n_components=1).fit(d[GENES])
    sign = 1 if pc.components_[0].sum() > 0 else -1
    d['PC1_4'] = zscore(pd.Series(pc.transform(d[GENES])[:,0]*sign, index=d.index))
    d['Equal4_patient_components'] = zscore(d[GENES].mean(axis=1))
    d['BAK1_centered'] = zscore((d.BAK1+d[GENES[1:]].mean(axis=1))/2)
    d['Core3'] = zscore(d[GENES[:3]].mean(axis=1))
    d['BAK1_SLC7A11'] = zscore(d[['BAK1','SLC7A11']].mean(axis=1))
    d['BAK1_only'] = zscore(d.BAK1)
    # Historical Figure 4D Rest3 is the un-restandardized mean of the
    # three patient-standardized components (axis is coefficient, not per SD).
    d['Rest3'] = zscore(d[GENES[1:]]).mean(axis=1)
    d['BAK1_joint'] = zscore(d.BAK1)  # explicit complete-case component scaling for Figure 4D
    rows, influence, loo = [], [], []
    for g in GENES+['Equal4']:
        row, fit, used = model(d, 'grade', [g,'purity','Age'], g, f'TRACERx.grade.{g}')
        rows.append(row)
        inf = fit.get_influence()
        for idx, h, cook in zip(used.index, inf.hat_matrix_diag, inf.cooks_distance[0]):
            influence.append(dict(model_id=row['model_id'], patient=idx, leverage=float(h), cooks_distance=float(cook)))
        if g == 'SLC7A11':
            for idx in used.index:
                rr, _, _ = model(d.drop(idx), 'grade', [g,'purity','Age'], g, 'TRACERx.grade.SLC7A11.leave_one_out')
                loo.append({**rr, 'omitted_patient':idx})
    variants = ['Equal4','PC1_4','BAK1_centered','Core3','BAK1_SLC7A11','BAK1_only','Equal4_patient_components']
    for v in variants:
        for outcome in ['grade','stage']:
            mid=f'TRACERx.{outcome}.{v}'
            if not any(r['model_id']==mid for r in rows):
                rows.append(model(d, outcome, [v,'purity','Age'], v, mid)[0])
        rows.append(model(d, v, ['punctuated']+FULL, 'punctuated', f'TRACERx.punctuated.full.{v}')[0])
    for label, covariates in [('clinical',CLINICAL),('clinical_wgii',CLINICAL+['wgii']),('clinical_ITH',CLINICAL+['ITH'])]:
        rows.append(model(d, 'Equal4', ['punctuated']+covariates, 'punctuated', f'TRACERx.punctuated.{label}.Equal4')[0])
    for outcome,predictors in [('grade',['BAK1_joint','Rest3','purity','Age']),
                               ('punctuated',['BAK1_joint','Rest3']+FULL)]:
        used=d[[outcome]+predictors].dropna()
        if outcome=='grade':
            for term in ['BAK1_joint','Rest3']:
                rows.append(model(d,outcome,predictors,term,f'TRACERx.joint.{outcome}.{term}')[0])
        else:
            x=sm.add_constant(used[predictors]);fit=sm.Logit(used[outcome],x).fit(disp=False)
            for term in ['BAK1_joint','Rest3']:
                ci=fit.conf_int().loc[term]
                rows.append(dict(model_id=f'TRACERx.joint.{outcome}.{term}',outcome=outcome,term=term,
                    covariates=';'.join(predictors),n=len(used),beta=float(fit.params[term]),
                    se_HC3=np.nan,ci_low=float(ci.iloc[0]),ci_high=float(ci.iloc[1]),p=float(fit.pvalues[term]),
                    inference='logistic model covariance',df_resid=float(fit.df_resid),
                    patient_sha256=hashlib.sha256('\n'.join(sorted(map(str,used.index))).encode()).hexdigest()))
    pd.DataFrame(rows).to_csv(ROOT/'results/TRACERx_models.csv',index=False)
    pd.DataFrame(influence).to_csv(ROOT/'results/TRACERx_influence.csv',index=False)
    pd.DataFrame(loo).to_csv(ROOT/'results/TRACERx_SLC7A11_leave_one_out.csv',index=False)
    pd.DataFrame({'gene':GENES,'loading':pc.components_[0]*sign}).to_csv(ROOT/'results/TRACERx_PC1_loadings.csv',index=False)
    patient.to_csv(ROOT/'data/derived/TRACERx_all_patients.csv')
    d.to_csv(ROOT/'data/derived/TRACERx_analysis_patients.csv')
    regions.join(ann[['Patient','purity']]).to_csv(ROOT/'data/derived/TRACERx_regions.csv')
    pd.DataFrame({'variable':GENES+['Equal4'], 'SD_all_patients':patient[GENES+['Equal4']].std(ddof=0),
                  'SD_model_patients':d[GENES+['Equal4']].std(ddof=0)}).to_csv(ROOT/'results/TRACERx_scaling.csv',index=False)
    counts={'primary_regions':len(ann),'expression_patients':len(patient),'complete_case_patients':len(d),
            'complete_case_regions':int(ann.Patient.isin(d.index).sum()),
            'multiregion_complete_case_patients':int((ann[ann.Patient.isin(d.index)].groupby('Patient').size()>1).sum()),
            'score_definition':'raw TPM gene-z on primary regions; equal mean; patient mean; final z on all expression-eligible patients',
            'architecture_alternatives':'patient-component scaling and final z on complete-case patients; explicitly secondary'}
    (ROOT/'results/TRACERx_sample_summary.json').write_text(json.dumps(counts,indent=2)+'\n')
    return pd.DataFrame(rows)

def cptac():
    raw=ROOT/'data/raw/CPTAC'
    expr=pd.read_csv(raw/'HS_CPTAC_CCRCC_RNAseq_fpkm_log2_Tumor.cct',sep='\t',index_col=0)
    expr=expr.groupby(level=0).mean()
    gz=(expr.T-expr.T.mean())/expr.T.std(ddof=0).replace(0,np.nan)
    from .clinical import cptac_clinical
    clinical=cptac_clinical(expr.columns)
    clinical_provenance='public study Table S1, clear-cell histology only'
    # The public portal supplies a lower bound for top-coded ages.
    portal=pd.read_excel(raw/'S044_CPTAC_CCRCC_Discovery_Cohort_Clinical_Data_r1_Sept2018.xlsx',header=2,index_col='Case_ID')
    clinical['Age_lower_bound']=pd.to_numeric(clinical.Age,errors='coerce')
    missing=clinical.index[clinical.Age_lower_bound.isna()]
    for idx in missing:
        value=str(portal.loc[idx,'Age'])
        if value.strip() in ['>90','> 90','90+','≥90','>=90']:
            clinical.loc[idx,'Age_lower_bound']=90
        else:
            raise ValueError(f'Unresolved age for {idx}: {value}')
    clinical['Male']=clinical.Gender.map({'Male':1,'Female':0})
    clinical['Mitoxy']=zscore(gz[GENES].mean(axis=1))
    gmt={}
    for line in (ROOT/'data/raw/genesets/h.all.v2024.1.Hs.symbols.gmt').read_text().splitlines():
        f=line.split('\t');gmt[f[0]]=set(f[2:])
    defs={'E2F_targets':gmt['HALLMARK_E2F_TARGETS'],'G2M_checkpoint':gmt['HALLMARK_G2M_CHECKPOINT'],
          'E2F_G2M_union':gmt['HALLMARK_E2F_TARGETS']|gmt['HALLMARK_G2M_CHECKPOINT']}
    coverage=[]
    for label,gs in defs.items():
        mapped=sorted((gs-set(GENES))&set(gz.columns))
        clinical[label]=zscore(gz[mapped].mean(axis=1))
        coverage.append(dict(signature=label,n_defined=len(gs-set(GENES)),n_measured=len(mapped),genes=';'.join(mapped)))
    rows=[]
    for outcome in ['grade_num','stage_num']:
        # One common complete-case frame for the base and every proliferation model.
        d=clinical.dropna(subset=[outcome,'Mitoxy','Age_lower_bound','Male']+list(defs))
        for label,extra in [('base',[])]+[(k,[k]) for k in defs]:
            for term in ['Mitoxy']+extra:
                row,_,_=model(d,outcome,['Mitoxy','Age_lower_bound','Male']+extra,term,f'CPTAC.{outcome}.{label}.{term}')
                rows.append({**row,'cohort':'CPTAC','signature':label,'clinical_provenance':clinical_provenance,'expression_source_sha256':sha256(raw/'HS_CPTAC_CCRCC_RNAseq_fpkm_log2_Tumor.cct')})
    pd.DataFrame(coverage).to_csv(ROOT/'results/CPTAC_proliferation_gene_coverage.csv',index=False)
    pd.DataFrame(rows).to_csv(ROOT/'results/CPTAC_proliferation_models.csv',index=False)
    clinical.to_csv(ROOT/'data/derived/CPTAC_analysis_patients.csv')
    return pd.DataFrame(rows)
