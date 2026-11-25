"""Reconstruct CheckMate 025 from Braun et al. public supplementary workbooks."""
import json
import numpy as np
import pandas as pd
import statsmodels.api as sm
from openpyxl import load_workbook
from .data import ROOT,sha256
from .survival import cox,bh,km_median

GENES=['BAK1','BAX','BID','SLC7A11'];MY4=['MYD88','TNF','TLR2','TLR4']

def event_diagnostics(d):
    rows=[]
    for arm,g in d.groupby('Arm'):
        for endpoint in ['PFS','OS']:
            a=g.dropna(subset=[endpoint,endpoint+'_CNSR'])
            for mapping in ['1=event','0=event']:
                events=a[endpoint+'_CNSR'].to_numpy(int)
                if mapping=='0=event':events=1-events
                rows.append(dict(arm=arm,endpoint=endpoint,mapping=mapping,n=len(a),
                    events=int(events.sum()),km_median=km_median(a[endpoint].to_numpy(float),events)))
    out=ROOT/'results/CheckMate025';out.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(rows).to_csv(out/'event_coding_diagnostics.csv',index=False)
    return pd.DataFrame(rows)


def extract_inputs():
    source=ROOT/'data/raw/CheckMate025/Braun2020_Supplementary_Tables.xlsx'
    out=ROOT/'data/derived/CheckMate025';out.mkdir(parents=True,exist_ok=True)
    wb=load_workbook(source,read_only=True,data_only=True)
    clinical=None;expression=None;identified=[]
    targets=set(GENES+MY4+['RICTOR'])
    for sheet in wb:
        rows=iter(sheet.values)
        for _ in range(20):
            row=next(rows,None)
            if row is None:break
            if 'SUBJID' in row and 'RNA_ID' in row:
                if clinical is not None:raise ValueError('Ambiguous clinical worksheet')
                clinical=pd.DataFrame(rows,columns=row).replace('NA',None).dropna(subset=['SUBJID'])
                identified.append(dict(role='clinical',sheet=sheet.title));break
            if 'gene_name' in row:
                if expression is not None:raise ValueError('Ambiguous expression worksheet')
                col=list(row).index('gene_name');selected=[r for r in rows if r[col] in targets]
                expression=pd.DataFrame(selected,columns=row)
                if set(expression.gene_name)!=targets or expression.gene_name.duplicated().any():
                    raise ValueError('Missing or duplicated target gene in public expression')
                identified.append(dict(role='expression',sheet=sheet.title));break
    wb.close()
    if clinical is None or expression is None:
        raise ValueError('Required clinical/expression sheets not identified')
    clinical.to_csv(out/'clinical_all.csv',index=False)
    expression.to_csv(out/'target_gene_expression.csv',index=False)
    (out/'source.json').write_text(json.dumps(dict(public_file=str(source.relative_to(ROOT)),sha256=sha256(source),sheets=identified),indent=2)+'\n')
    return clinical,expression


def analyze(expression_transform='as_released'):
    clinical,rna=extract_inputs();rna=rna.set_index('gene_name').T.apply(pd.to_numeric,errors='raise')
    cc=clinical.loc[clinical.Cohort.eq('CM-025') & clinical.RNA_ID.notna()].copy()
    if cc.SUBJID.duplicated().any() or cc.RNA_ID.duplicated().any() or not rna.index.is_unique:
        raise ValueError('Nonunique CheckMate join keys')
    d=cc.merge(rna,left_on='RNA_ID',right_index=True,how='inner',validate='one_to_one')
    for col in ['Age','PFS','OS','PFS_CNSR','OS_CNSR']:d[col]=pd.to_numeric(d[col],errors='raise')
    v=d[GENES+MY4+['RICTOR']].astype(float)
    if expression_transform=='log2p1':v=np.log2(v+1)
    elif expression_transform!='as_released':raise ValueError('Unknown expression transform')
    if not np.isfinite(v.to_numpy()).all():raise ValueError('Nonfinite expression')
    z=(v-v.mean())/v.std(ddof=0)
    scores={'Equal4':z[GENES].mean(axis=1),'Core3':z[GENES[:3]].mean(axis=1),
            'Myeloid':z[MY4].mean(axis=1),**{g:z[g] for g in GENES+['RICTOR']}}
    for name,values in scores.items():d[name+'_SD']=(values-values.mean())/values.std(ddof=0)
    d['age_decade']=d.Age/10;d['male']=d.Sex.str.upper().map({'M':1,'MALE':1,'F':0,'FEMALE':0})
    d['high_arm_median']=(d.Equal4_SD>=d.groupby('Arm').Equal4_SD.transform('median')).astype(int)
    d['high_pooled_median']=(d.Equal4_SD>=d.Equal4_SD.median()).astype(int)
    d['nivolumab']=d.Arm.eq('NIVOLUMAB').astype(int)
    results=[]
    for family,contrasts,covs in [('primary',['continuous_per_pooled_SD','within_arm_median'],[]),
        ('age_sex_sensitivity',['continuous_per_pooled_SD','within_arm_median'],['age_decade','male']),
        ('pooled_median_sensitivity',['pooled_median'],[])]:
        for arm in ['NIVOLUMAB','EVEROLIMUS']:
            for endpoint in ['PFS','OS']:
                for contrast in contrasts:
                    focal={'continuous_per_pooled_SD':'Equal4_SD','within_arm_median':'high_arm_median','pooled_median':'high_pooled_median'}[contrast]
                    a=d.loc[d.Arm.eq(arm)].dropna(subset=[endpoint,endpoint+'_CNSR',focal]+covs)
                    fit=cox(a[endpoint],a[endpoint+'_CNSR'],a[[focal]+covs])
                    results.append(dict(family=family,arm=arm,endpoint=endpoint,contrast=contrast,
                       n=len(a),events=int(a[endpoint+'_CNSR'].sum()),covariates=';'.join(covs),
                       **{k:float(fit[k][0]) for k in ['beta','se','HR','ci_low','ci_high','p']}))
    results=pd.DataFrame(results)
    results['q']=results.groupby('family').p.transform(lambda v:bh(v.to_numpy()))
    interactions=[]
    for name in scores:
        d['score']=d[name+'_SD'];d['interaction']=d.score*d.nivolumab
        for endpoint in ['ORR','PFS','OS']:
            xs=['score','nivolumab','interaction']
            if endpoint=='ORR':
                response=d.ORR.astype(str).str.upper().map({'CR':1,'PR':1,'CRPR':1,'SD':0,'PD':0})
                use=response.notna();f=sm.Logit(response[use],sm.add_constant(d.loc[use,xs])).fit(disp=False)
                ci=f.conf_int().loc['interaction'];row=dict(beta=f.params['interaction'],se=f.bse['interaction'],
                    HR=np.exp(f.params['interaction']),ci_low=np.exp(ci.iloc[0]),ci_high=np.exp(ci.iloc[1]),p=f.pvalues['interaction'])
                n=int(use.sum());events=int(response[use].sum())
            else:
                a=d.dropna(subset=[endpoint,endpoint+'_CNSR']+xs)
                fit=cox(a[endpoint],a[endpoint+'_CNSR'],a[xs]);row={k:float(fit[k][2]) for k in ['beta','se','HR','ci_low','ci_high','p']}
                n=len(a);events=int(a[endpoint+'_CNSR'].sum())
            interactions.append(dict(score=name,endpoint=endpoint,n=n,events_or_responses=events,
                    treatment_reference='everolimus',ratio_type='OR' if endpoint=='ORR' else 'HR',**row))
    interactions=pd.DataFrame(interactions);interactions['q_BH_24']=bh(interactions.p)
    out=ROOT/'results/CheckMate025';out.mkdir(parents=True,exist_ok=True)
    event_diagnostics(d)
    results.to_csv(out/'within_arm_models.csv',index=False);interactions.to_csv(out/'interaction_models.csv',index=False)
    d.drop(columns=['score','interaction']).to_csv(ROOT/'data/derived/CheckMate025/patient_scores.csv',index=False)
    unmatched=cc.loc[~cc.RNA_ID.isin(rna.index),['SUBJID','RNA_ID','Arm']]
    unmatched.to_csv(out/'unmatched_clinical_RNA_IDs.csv',index=False)
    (out/'methods.json').write_text(json.dumps(dict(expression_transform=expression_transform,
        gene_z_ddof=0,score_z_ddof=0,event_mapping='1=event; 0=censored; retained study specification',
        clinical_n=len(cc),matched_n=len(d),matched_by_arm=d.Arm.value_counts().to_dict(),
        ties='Breslow',primary_family='8 within-arm contrasts',interaction_family='8 scores × 3 endpoints',
        historical_236_case_mask='not available; exact public RNA_ID join is used'),indent=2)+'\n')
    return results,interactions
