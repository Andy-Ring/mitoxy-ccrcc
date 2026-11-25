"""Patient-level Hallmark inference, reconstructed from supplied permutation code."""
import json
import numpy as np
import pandas as pd
from .data import ROOT
from .bulk import prepare,read_gmt
from .analysis import GENES
from .survival import bh

def inputs():
    e,meta,cov,_=prepare('TCGA')
    scores=pd.read_csv(ROOT/'results/bulk/TCGA_main_patient_scores.csv',index_col=0).rename(columns={'grade':'grade_ord'})
    scores['high']=(scores.Mitoxy>scores.Mitoxy.median()).astype(int)
    d=scores.dropna(subset=['high','grade_ord','age','male','purity','PC1'])
    return e,d

def competitive():
    wide,D=inputs();wide=wide.loc[D.index];output_dir=ROOT/"results/pathways";output_dir.mkdir(parents=True,exist_ok=True)
    Y=wide.to_numpy(float);genes=np.array(wide.columns.astype(str));valid=~np.isin(genes,['BAK1','BAX','BID','SLC7A11']);genes=genes[valid];Y=Y[:,valid];v=Y.var(0)>0;genes=genes[v];Y=Y[:,v]
    # Use one common model for each gene: expression ~ group + age + sex + grade + purity + PC1.
    base_cov=D[['age','male','grade_ord','purity','PC1']].to_numpy(float)

    def gene_t(y):
     X=np.column_stack([np.ones(len(y)),y.astype(float),base_cov])
     bread=np.linalg.inv(X.T@X);B=bread@X.T@Y;res=Y-X@B;h=np.einsum('ij,jk,ik->i',X,bread,X);c=bread[1,:]@X.T
     var=np.sum((c[:,None]*res/(1-h)[:,None])**2,axis=0)
     return np.nan_to_num(B[1,:]/np.sqrt(var),nan=0.0,posinf=0.0,neginf=0.0)
    obs_t=gene_t(D.high.to_numpy());
    # Hallmark sets and set-member positions.
    sets=[]
    for line in (ROOT/'data/raw/genesets/h.all.v2024.1.Hs.symbols.gmt').read_text().splitlines():
     f=line.split('\t'); members=set(f[2:])&set(genes)
     if 10<=len(members)<=500: sets.append((f[0].replace('HALLMARK_',''),members))
    lookup={g:i for i,g in enumerate(genes)}
    # Weighted running-sum ES (|t| weighting) using the ranked gene list.
    def es_for_set(rank_order,stat,set_idx):
     inv=np.empty_like(rank_order);inv[rank_order]=np.arange(len(rank_order));loc=np.sort(inv[np.asarray(set_idx,dtype=int)])
     m=len(loc); N=len(stat)
     if not m or m==N:return 0.
     w=np.abs(stat[rank_order[loc]]); total=w.sum()
     if total<=0:return 0.
     hit=np.cumsum(w)/total
     # After each hit, count misses since previous hit.
     miss=(loc-np.arange(m))/(N-m)
     vals=hit-miss
     before=np.r_[0.0,hit[:-1]]-miss
     mx=vals.max(); mn=before.min()
     return float(mx if mx>=abs(mn) else mn)
    indices=[np.array([lookup[g] for g in ms],dtype=int) for _,ms in sets]
    def all_es(stat):
     order=np.argsort(-stat,kind='stable')
     return np.array([es_for_set(order,stat,idx) for idx in indices])
    obs=all_es(obs_t)
    # Grade-stratified shuffles of high/low labels preserve the observed counts in each grade stratum.
    rng=np.random.default_rng(20261110); y0=D.high.to_numpy(); strata=[np.flatnonzero(D.grade_ord.to_numpy()==g) for g in np.sort(D.grade_ord.unique())]
    null=np.empty((1000,len(sets)),dtype=np.float32)
    for b in range(1000):
     y=y0.copy()
     for ix in strata:y[ix]=rng.permutation(y[ix])
     null[b]=all_es(gene_t(y))
     if (b+1)%250==0:print('permutations',b+1,flush=True)
    # Direction-specific ES normalization; empirical nominal P and BH q values.
    rows=[]
    for j,(name,members) in enumerate(sets):
     x=obs[j]; ns=null[:,j]
     if x>=0:
      same=ns[ns>=0];den=float(same.mean()) if len(same) else 1.;p=(1+int(np.sum(ns>=x)))/(1+len(same))
     else:
      same=np.abs(ns[ns<0]);den=float(same.mean()) if len(same) else 1.;p=(1+int(np.sum(ns<=x)))/(1+len(same))
     rows.append([name,len(members),x,x/den,p,'higher_in_Equal4_high' if x>0 else 'higher_in_Equal4_low'])
    out=pd.DataFrame(rows,columns=['Hallmark','set_size','ES','NES','permutation_p','direction'])
    pv=out.permutation_p.to_numpy();ordr=np.argsort(pv);adj=np.minimum.accumulate((pv[ordr]*len(pv)/np.arange(1,len(pv)+1))[::-1])[::-1];q=np.empty_like(pv);q[ordr]=np.minimum(adj,1);out['BH_q']=q
    out=out.sort_values('permutation_p');out.to_csv(output_dir/'TCGA_Hallmark_competitive_sample_permutation.csv',index=False)
    np.savez_compressed(output_dir/'TCGA_Hallmark_competitive_null.npz',null_es=null,observed_es=obs,sets=np.array([x[0] for x in sets]))
    summary={'n_patients':len(D),'high_n':int(D.high.sum()),'low_n':int((1-D.high).sum()),'n_genes':len(genes),'n_sets':len(sets),'permutations':1000,'grade_stratified':True,'covariates':['age','sex','grade','ABSOLUTE purity','global PC1'],'seed':20261110,'n_q05':int((out.BH_q<.05).sum()),'n_q25':int((out.BH_q<.25).sum()),'q25_directions':out.loc[out.BH_q<.25,'direction'].value_counts().to_dict()}
    (output_dir/'TCGA_Hallmark_competitive_methods.json').write_text(json.dumps(summary,indent=2));print(summary);print(out.head(15).to_string(index=False))

def patient_scores():
    patient,scores=inputs();output_dir=ROOT/'results/pathways';output_dir.mkdir(parents=True,exist_ok=True)
    gene_z=(patient-patient.mean())/patient.std(ddof=0).replace(0,np.nan)
    gmt={k:set(v) for k,v in read_gmt(ROOT/'data/raw/genesets/h.all.v2024.1.Hs.symbols.gmt').items()}
    all_gene_sets={k:v-{'BAK1','BAX','BID','SLC7A11'} for k,v in gmt.items()}
    pathways={}
    for name,genes in all_gene_sets.items():
        mapped=sorted(genes & set(gene_z.columns))
        if 10<=len(mapped)<=500: pathways[name]=gene_z[mapped].mean(axis=1)
    P=pd.DataFrame(pathways).join(scores[['Mitoxy','grade_ord','age','male','purity','PC1']])
    P['Equal4_high']=(P.Mitoxy>scores.Mitoxy.median()).astype(float)
    required=['Equal4_high','grade_ord','age','male','purity','PC1']
    P=P.dropna(subset=required)
    path_names=list(pathways); Y=P[path_names].to_numpy(float)
    covars=P[['age','male','grade_ord','purity','PC1']].to_numpy(float)
    X0=np.column_stack([np.ones(len(P)),covars])
    X=np.column_stack([np.ones(len(P)),P.Equal4_high.to_numpy(float),covars])
    if np.linalg.matrix_rank(X)!=X.shape[1]: raise ValueError('Rank-deficient Hallmark design')
    br0=np.linalg.inv(X0.T@X0); fitted0=X0@(br0@X0.T@Y); resid0=Y-fitted0
    bread=np.linalg.inv(X.T@X); c=(bread@X.T)[1,:]; h=np.sum((X@bread)*X,axis=1)
    beta=bread@X.T@Y; resid=Y-X@beta
    se=np.sqrt(np.sum((c[:,None]*resid/(1-h)[:,None])**2,axis=0))
    t_obs=beta[1,:]/se
    rng=np.random.default_rng(20261110); exceed=np.zeros(len(path_names),dtype=np.int64)
    for _ in range(10000):
        Ystar=fitted0+resid0[rng.permutation(len(P)),:]
        bstar=bread@X.T@Ystar; rstar=Ystar-X@bstar
        ses=np.sqrt(np.sum((c[:,None]*rstar/(1-h)[:,None])**2,axis=0))
        tstar=bstar[1,:]/ses
        exceed+=(np.abs(tstar)>=np.abs(t_obs))
    pvals=(exceed+1)/(10000+1)
    order=np.argsort(pvals); q=np.empty(len(pvals)); q_sorted=np.minimum.accumulate((pvals[order]*len(pvals)/np.arange(1,len(pvals)+1))[::-1])[::-1]; q[order]=np.minimum(q_sorted,1)
    out=pd.DataFrame({'pathway':path_names,'n_genes':[(len(all_gene_sets[n]&set(gene_z.columns))) for n in path_names],
                      'n_patients':len(P),'beta_high_vs_low':beta[1,:],'se_HC3':se,
                      'ci95_low':beta[1,:]-1.96*se,'ci95_high':beta[1,:]+1.96*se,
                      'p_Freedman_Lane':pvals,'BH_q_50_sets':q,
                      'direction':np.where(beta[1,:]>0,'higher in Equal4-high','higher in Equal4-low'),
                      'permutations':10000,'permutation_unit':'patient residual rows'})
    out=out.sort_values('p_Freedman_Lane')
    out.to_csv(output_dir/'TCGA_Hallmark_sample_permutation.csv',index=False)
