"""Public GSE210041 counts → spatial scores, definitions and permutations.

Archived, grid and cross-family specifications remain separate. The original
submitted marker sets are not recoverable from the available source scripts.
"""
import itertools
import json
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
from scipy.stats import spearmanr
from sklearn.neighbors import NearestNeighbors
from .data import ROOT, sha256
from .tenx import read_counts, positions

SLIDES={'T1':'GSM6415705_GLMF1','T7':'GSM6415706_GLMF2'}
TUMOR=['KRT8','KRT18','KRT19','KRT7','EPCAM','CA9','VIM']
MYELOID=['C1QB','C1QC','TYROBP','C1QA','FCER1G','CTSS','AIF1','CSF1R','PTPRC','LST1','CD68']
EQUAL=['BAK1','BAX','BID','SLC7A11']; MY4=['MYD88','TNF','TLR2','TLR4']
TP=['CA9','NDUFA4L2','ANGPTL4','EGLN3','VEGFA','PAX8','PAX2','KRT18','KRT8','NNMT','SLC2A1','HILPDA']
MP=['CD68','CD163','CD14','C1QA','C1QB','C1QC','LYZ','AIF1','CSF1R','TYROBP','FCER1G','MRC1']


def public_inputs():
    out=ROOT/'data/derived/GSE210041';out.mkdir(parents=True,exist_ok=True)
    records=[];slides={}
    for name,prefix in SLIDES.items():
        raw=ROOT/'data/raw/GSE210041';h5=raw/(prefix+'_filtered_feature_bc_matrix.h5')
        coords=raw/(prefix+'_spatial.tar.gz')
        x,bars,genes,_=read_counts(h5);pos=positions(coords)
        keep=bars.isin(pos.index[pos.in_tissue.eq(1)]);x=x[keep];bars=bars[keep];pos=pos.loc[bars]
        available=[g for g in dict.fromkeys(TUMOR+MYELOID+EQUAL+MY4+TP+MP) if g in genes]
        if any(np.sum(genes==g)!=1 for g in available):
            raise ValueError('Duplicated required spatial marker requires explicit resolution')
        required=set(EQUAL+MY4+TP+MP)
        if not required <= set(available):
            raise ValueError(f'Missing required spatial markers: {required-set(available)}')
        ix=[int(np.flatnonzero(genes==g)[0]) for g in available]
        d=pd.DataFrame(x[:,ix].toarray(),index=bars,columns=available)
        d['total_UMI']=np.asarray(x.sum(axis=1)).ravel()
        d=d.join(pos[['pxl_col','pxl_row']]);slides[name]=(x,bars,genes,pos,d)
        d.to_csv(out/f'{name}_raw_marker_counts_and_library_totals.csv.gz')
        records.append(dict(slide=name,n_spots=len(bars),n_genes=len(genes),
                            counts_sha256=sha256(h5),coordinates_sha256=sha256(coords)))
    pd.DataFrame(records).to_csv(out/'public_inputs.csv',index=False)
    return slides


def zmarkers(d,genes,round_to=None):
    y=np.log1p(d[genes].to_numpy(float)*1e4/d.total_UMI.to_numpy()[:,None])
    sd=y.std(axis=0,ddof=0)
    if np.any(sd<=0):
        raise ValueError('Constant spatial marker')
    z=pd.DataFrame((y-y.mean(axis=0))/sd,index=d.index,columns=genes)
    return z.round(round_to) if round_to is not None else z


def archived_score(x, genes, selected):
    """Retain sparse normalization and per-gene reduction order of source code.

    Spatial ranks contain ties; vectorizing these reductions can break ties at
    machine precision and change a Spearman statistic or permutation count.
    """
    lib=np.asarray(x.sum(axis=1)).ravel()
    y=x.multiply(1e4/np.maximum(lib,1)[:,None]).tocsr();y.data=np.log1p(y.data)
    idx={g:i for i,g in enumerate(genes)};values=[]
    for g in selected:
        if g in idx:
            v=y[:,idx[g]].toarray().ravel();sd=v.std(ddof=0)
            values.append((v-v.mean())/sd if sd>0 else np.zeros_like(v))
    return np.sum(values,axis=0)/len(selected)


def evaluate(eq,my4,tumor,myeloid,xy,k=6,nperm=0,seed=20261110,neighbor='ckdtree',quantiles=(.75,.75,.5),round_to=None):
    tq,mq,ex=quantiles
    mal=(tumor>np.quantile(tumor,tq))&(myeloid<np.quantile(myeloid,ex))
    mye=(myeloid>np.quantile(myeloid,mq))&(tumor<np.quantile(tumor,ex))
    r=dict(n_tumor=int(mal.sum()),n_myeloid=int(mye.sum()),valid=False,rho=np.nan)
    if mal.sum()<10 or mye.sum()<k:return r,None,None
    if neighbor=='stable_barcode':
        idx=np.argsort(cdist(xy[mal],xy[mye],metric='sqeuclidean'),axis=1,kind='stable')[:,:k]
    elif neighbor=='sklearn':
        _,idx=NearestNeighbors(n_neighbors=k).fit(xy[mye]).kneighbors(xy[mal])
    else:
        _,idx=cKDTree(xy[mye]).query(xy[mal],k=k)
    vals=my4[mye];nb=vals[idx].mean(axis=1)
    if round_to is not None:nb=np.round(nb,round_to)
    rho=float(spearmanr(eq[mal],nb).statistic);r.update(valid=bool(np.isfinite(rho)),rho=rho)
    if nperm:
        rng=np.random.default_rng(seed)
        null=np.asarray([spearmanr(eq[mal],rng.permutation(vals)[idx].mean(axis=1)).statistic for _ in range(nperm)])
        r.update(nperm=nperm,seed=seed,p_two_sided=float((1+np.sum(np.abs(null)>=abs(rho)))/(nperm+1)))
    else:null=None
    return r,dict(malignant=mal,myeloid=mye,neighbor_index=idx,neighbor_mean=nb),null


def analyze(nperm=5000,n_random=2000):
    slides=public_inputs();out=ROOT/'results/GSE210041';out.mkdir(parents=True,exist_ok=True)
    summary=[];grid=[];random=[];pre=[];rng=np.random.default_rng(20261110)
    shared=set.intersection(*(set(v[4].columns) for v in slides.values()))
    tp=[g for g in TUMOR if g in shared];mp=[g for g in MYELOID if g in shared]
    panels=lambda gs:[('full',gs)]+[('omit_'+g,[x for x in gs if x!=g]) for g in gs]
    combos=lambda gs:[c for k in range(1,6) for c in itertools.combinations(gs,k)]
    tc,mc=combos(TP),combos(MP)
    maps=[]
    for name,(x,bars,genes,pos,d) in slides.items():
        # Archived alternative: float64 normalization, original H5 ordering,
        # fixed nominal marker denominators and sklearn neighbor tie handling.
        cols=list(dict.fromkeys(tp+mp+EQUAL+MY4));z=zmarkers(d,cols)
        eq=archived_score(x,genes,EQUAL);my4=archived_score(x,genes,MY4)
        tumor=archived_score(x,genes,TUMOR);myel=archived_score(x,genes,MYELOID)
        xy=d[['pxl_col','pxl_row']].to_numpy(float)
        r,masks,null=evaluate(eq,my4,tumor,myel,xy,nperm=nperm,seed=20260913,neighbor='sklearn')
        summary.append(dict(sample=name,analysis='archived_alternative',**r))
        np.save(out/f'{name}_archived_permutation_null.npy',null)
        mapping=pd.DataFrame(dict(barcode=bars,Equal4=eq,Myeloid4=my4,tumor_score=tumor,
                 myeloid_abundance=myel,malignant_enriched=masks['malignant'],myeloid_enriched=masks['myeloid']),index=bars)
        mapping=mapping.join(pos[['pxl_col','pxl_row']]);mapping['sample']=name
        mapping.to_csv(out/f'{name}_spatial_scores.csv',index=False);maps.append(mapping)
        # New grid: sort barcodes and use explicit 12-decimal numerical ties.
        gd=d.sort_index();gz=zmarkers(gd,cols,round_to=12)
        geq=gz[EQUAL].mean(axis=1).round(12).to_numpy();gm4=gz[MY4].mean(axis=1).round(12).to_numpy()
        gxy=gd[['pxl_col','pxl_row']].to_numpy(float)
        for (tn,tgs),(mn,mgs) in itertools.product(panels(tp),panels(mp)):
            t=gz[tgs].mean(axis=1).round(12).to_numpy();m=gz[mgs].mean(axis=1).round(12).to_numpy()
            for tq,mq,ex in itertools.product([.70,.75,.80],[.70,.75,.80],[.40,.50]):
                r,_,_=evaluate(geq,gm4,t,m,gxy,neighbor='stable_barcode',quantiles=(tq,mq,ex),round_to=12)
                grid.append(dict(sample=name,definition_id=f'{tn}__{mn}__t{tq:.2f}_m{mq:.2f}_x{ex:.2f}',
                     tumor_panel=tn,myeloid_panel=mn,tumor_genes=';'.join(tgs),myeloid_genes=';'.join(mgs),
                     tumor_quantile=tq,myeloid_quantile=mq,reciprocal_exclusion_quantile=ex,**r))
        # Cross-family scripts used scanpy float32 normalization and cKDTree.
        # Retain that numerical convention separately from the archived/grid runs.
        import scanpy as sc
        import anndata as ad
        a=ad.AnnData(x.astype(np.float32),obs=pos.copy(),var=pd.DataFrame(index=genes))
        a.var_names_make_unique()
        sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a)
        cgenes=list(dict.fromkeys(TP+MP+EQUAL+MY4));v=pd.DataFrame(a[:,cgenes].X.toarray(),index=bars,columns=cgenes)
        cz=(v-v.mean())/v.std(ddof=0);ce=cz[EQUAL].mean(axis=1).to_numpy();cm=cz[MY4].mean(axis=1).to_numpy()
        cxy=pos[['pxl_row','pxl_col']].to_numpy(float)
        hyp=['CA9','NDUFA4L2','ANGPTL4','EGLN3'];my=['CD68','CD163','CD14','C1QA','C1QB','C1QC','LYZ','AIF1']
        for k in [6,4,8,10]:
            r,_,null=evaluate(ce,cm,cz[hyp].mean(axis=1).to_numpy(),cz[my].mean(axis=1).to_numpy(),cxy,k=k,nperm=nperm)
            pre.append(dict(sample=name,k=k,**r));np.save(out/f'{name}_hypoxia_k{k}_permutation_null.npy',null)
        for b in range(n_random):
            ts=list(tc[rng.integers(len(tc))]);ms=list(mc[rng.integers(len(mc))])
            r,_,_=evaluate(ce,cm,cz[ts].mean(axis=1).to_numpy(),cz[ms].mean(axis=1).to_numpy(),cxy)
            random.append(dict(sample=name,definition=b,tumor_genes=';'.join(ts),myeloid_genes=';'.join(ms),**r))
        print(f'{name}: archived, grid and cross-family analyses computed',flush=True)
    for filename,rows in [('archived_summary.csv',summary),('spatial_grid_all_definitions.csv.gz',grid),
                          ('cross_family_all_definitions.csv.gz',random),('hypoxia_coupling.csv',pre)]:
        pd.DataFrame(rows).to_csv(out/filename,index=False)
    ranges=[]
    for family,rows in [('grid',grid),('random',random)]:
        for sample,dd in pd.DataFrame(rows).groupby('sample'):
            r=dd.loc[dd.valid,'rho'];ranges.append(dict(family=family,sample=sample,n_attempted=len(dd),n_valid=len(r),
               median_rho=r.median(),q025=r.quantile(.025),q975=r.quantile(.975),min_rho=r.min(),max_rho=r.max(),fraction_positive=(r>0).mean()))
    pd.DataFrame(ranges).to_csv(out/'definition_ranges.csv',index=False)
    (out/'methods.json').write_text(json.dumps(dict(nperm=nperm,n_random=n_random,
       original_submitted='Not regenerated: original marker definitions are not recorded in available source code',
       archived_seed=20260913,cross_family_seed=20261110,grid_rounding_decimals=12,
       missing_archived_tumor_markers=sorted(set(TUMOR)-set(tp))),indent=2)+'\n')
    return pd.DataFrame(summary),pd.DataFrame(ranges)
