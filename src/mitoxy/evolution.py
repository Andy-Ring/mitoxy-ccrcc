"""TRACERx regional dispersion with the manuscript's complete-case patient mask."""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr,rankdata,norm
from .data import ROOT

def analyze():
    r=pd.read_csv(ROOT/'data/derived/TRACERx_regions.csv',index_col=0)
    patients=pd.read_csv(ROOT/'data/derived/TRACERx_analysis_patients.csv',index_col=0)
    r=r[r.Patient.isin(patients.index)].copy()
    sizes=r.groupby('Patient').size();eligible=sizes[sizes>=2].index;r=r[r.Patient.isin(eligible)]
    summary=r.groupby('Patient').Equal4.agg(mean='mean',sd='std',range=lambda x:x.max()-x.min(),mad=lambda x:np.median(np.abs(x-np.median(x))))
    summary['n_regions']=sizes;summary=summary.join(patients[['ITH','wgii','purity']])
    rows=[]
    for feature in ['ITH','wgii']:
        for metric in ['range','sd','mad']:
            d=summary[[metric,feature,'purity','n_regions']].dropna();n=len(d)
            rho,p=spearmanr(d[metric],d[feature])
            cov=np.column_stack([np.ones(n),rankdata(d.n_regions),rankdata(d.purity)])
            x=rankdata(d[metric]);y=rankdata(d[feature])
            xr=x-cov@np.linalg.lstsq(cov,x,rcond=None)[0];yr=y-cov@np.linalg.lstsq(cov,y,rcond=None)[0]
            # Partial correlation of ranked variables; two nuisance covariates.
            from scipy.stats import t
            adjusted=np.corrcoef(xr,yr)[0,1];ap=2*t.sf(abs(adjusted)*np.sqrt((n-4)/(1-adjusted**2)),n-4)
            for kind,value,pv,k in [('raw',rho,p,0),('adjusted',adjusted,ap,2)]:
                ci=np.tanh(np.arctanh(np.clip(value,-.999999,.999999))+np.array([-1,1])*norm.ppf(.975)/np.sqrt(n-k-3))
                rows.append(dict(genomic_feature=feature,metric=metric,model=kind,n=n,rho=value,ci_low=ci[0],ci_high=ci[1],p=pv))
    out=ROOT/'results/evolution';out.mkdir(parents=True,exist_ok=True)
    r.to_csv(out/'regions_complete_case.csv');summary.to_csv(out/'patient_dispersion.csv');pd.DataFrame(rows).to_csv(out/'ITH_dispersion.csv',index=False)

    # Reverse-logistic sensitivity uses the same complete-case covariates.
    import statsmodels.api as sm
    rows=[]
    for variant in ['main','no_overlap']:
        d=pd.read_csv(ROOT/f'results/bulk/TRACERx_{variant}_patient_scores.csv',index_col=0)
        cols=['Mitoxy','Age','purity','grade','stage','wgii','ITH','Comparator_PC1']
        d=d.dropna(subset=['punctuated']+cols)
        fit=sm.Logit(d.punctuated,sm.add_constant(d[cols])).fit(disp=False)
        if not fit.mle_retvals['converged']:raise RuntimeError('Reverse model did not converge')
        lo,hi=fit.conf_int().loc['Mitoxy'];rows.append(dict(variant=variant,n=len(d),beta=fit.params.Mitoxy,ci_low=lo,ci_high=hi,p=fit.pvalues.Mitoxy))
    pd.DataFrame(rows).to_csv(out/'reverse_PC1.csv',index=False)
