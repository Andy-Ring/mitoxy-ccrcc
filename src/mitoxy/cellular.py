"""GSE159115 public counts → author annotations → patient-level pseudobulks."""
import itertools
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests
from .data import ROOT, sha256
from .tenx import read_counts

GENES = ['BAK1','BAX','BID','SLC7A11']


def exact_signed_rank(a, b):
    d = np.asarray(a, float) - np.asarray(b, float)
    if not np.isfinite(d).all():
        raise ValueError('Nonfinite paired values')
    d = d[d != 0]
    if len(d) == 0:
        return 1.0
    ranks = stats.rankdata(np.abs(d), method='average')
    observed = abs(np.dot(np.sign(d), ranks))
    null = np.asarray(list(itertools.product([-1,1], repeat=len(d)))) @ ranks
    return float(np.mean(np.abs(null) >= observed - 1e-12))


def build_counts():
    raw = ROOT/'data/raw/GSE159115'
    annotations = pd.read_csv(raw/'GSE159115_ccRCC_anno.csv.gz').set_index('cell')
    if not annotations.index.is_unique:
        raise ValueError('Duplicate author cell identifiers')
    metadata = json.loads((ROOT/'config/GSE159115_sample_metadata.json').read_text())
    by_sample = {s['title']:s for s in metadata}
    matrices = []; observations = []; gene_ids = None; symbols = None; evidence = []
    for sample in sorted(annotations['sample'].unique()):
        s = by_sample[sample]
        files = [u.rsplit('/',1)[-1] for u in s['urls'] if u.endswith('.h5')]
        if len(files) != 1:
            raise ValueError(f'Ambiguous matrix for {sample}')
        path = raw/files[0]; x,bars,syms,ids = read_counts(path)
        if gene_ids is None:
            gene_ids, symbols = ids,syms
        elif not np.array_equal(ids, gene_ids) or not np.array_equal(syms, symbols):
            raise ValueError('Gene order differs across samples; explicit alignment required')
        cells = pd.Index(sample+'_'+bars)
        wanted = annotations.loc[annotations['sample'].eq(sample)].copy()
        idx = cells.get_indexer(wanted.index)
        if np.any(idx < 0):
            raise ValueError(f'Author annotated cells absent from public matrix: {sample}')
        part = x[idx]; actual = np.asarray(part.sum(axis=1)).ravel()
        evidence.append(dict(sample=sample, n_cells=len(wanted), raw_sha256=sha256(path),
                             total_UMI_mismatches=int(np.sum(actual != wanted.total_UMI.to_numpy()))))
        wanted['public_total_UMI'] = actual
        matrices.append(part); observations.append(wanted)
    x = sparse.vstack(matrices, format='csr'); obs = pd.concat(observations)
    # Preserve author cell order, rather than relying on filesystem ordering.
    order = obs.index.get_indexer(annotations.index); x=x[order]; obs=obs.iloc[order]
    out=ROOT/'data/derived/GSE159115'; out.mkdir(parents=True,exist_ok=True)
    sparse.save_npz(out/'GSE159115_ccRCC_cell_gene_counts.npz',x)
    obs.to_csv(out/'GSE159115_ccRCC_cell_metadata.csv.gz')
    pd.DataFrame({'gene_id':gene_ids,'gene_symbol':symbols}).to_csv(out/'GSE159115_gene_order.csv',index=False)
    pd.DataFrame(evidence).to_csv(out/'public_matrix_annotation_checks.csv',index=False)
    if (obs.public_total_UMI <= 0).any():
        raise ValueError('Zero-count annotated cells')
    return x,obs,symbols


def pseudobulk(min_cells=20, ddof=1):
    x,obs,symbols = build_counts()
    records=[]; arrays=[]
    for (patient,celltype), ids in obs.groupby(['patient','anno'],sort=True).groups.items():
        loc = obs.index.get_indexer(ids)
        if len(loc) < min_cells:
            continue
        arrays.append(np.asarray(x[loc].sum(axis=0)).ravel())
        records.append(dict(patient=patient,cell_type=celltype,n_cells=len(loc)))
    counts = np.asarray(arrays); m=pd.DataFrame(records)
    # Source matrices may contain multiple stable IDs for a symbol. Sum counts
    # before CPM normalization; do not average normalized duplicate features.
    expr=pd.DataFrame(counts,columns=symbols).T.groupby(level=0,sort=False).sum().T
    libs=expr.sum(axis=1); log=np.log2(expr.div(libs,axis=0)*1e6+1)
    if not set(GENES) <= set(log.columns):
        raise ValueError('Required score genes missing')
    # The retained pseudobulk means use sample SD (ddof=1). Export the
    # population-SD alternative as a transparent scaling sensitivity below.
    z=(log[GENES]-log[GENES].mean())/log[GENES].std(ddof=ddof)
    if not np.isfinite(z.to_numpy()).all():
        raise ValueError('Cannot standardize pseudobulk score components')
    scored=pd.concat([m,z],axis=1); scored['Equal4']=z.mean(axis=1)
    out=ROOT/'results/GSE159115';out.mkdir(parents=True,exist_ok=True)
    scored.to_csv(out/'patient_celltype_scores.csv',index=False)
    pd.concat([m,log[GENES].add_suffix('_log2CPM1')],axis=1).to_csv(out/'component_logCPM.csv',index=False)
    population_z=(log[GENES]-log[GENES].mean())/log[GENES].std(ddof=0)
    population=pd.concat([m,population_z],axis=1);population['Equal4']=population_z.mean(axis=1)
    population.to_csv(out/'patient_celltype_scores_population_SD.csv',index=False)
    (out/'methods.json').write_text(json.dumps(dict(min_cells=min_cells,gene_z_ddof=ddof,
        normalization='log2(CPM+1), summed raw counts per patient/cell type',
        standardization_population='all eligible patient-by-author-cell-type pseudobulks',
        scoring_note='ddof=1 reconstructs retained means; ddof=0 sensitivity exported. General manuscript population-SD statement needs this exception.',
        n_pseudobulks=len(m),n_author_cells=len(obs),n_patients=int(obs.patient.nunique())),indent=2)+'\n')
    summary=scored.groupby('cell_type').agg(n_patients=('patient','nunique'),n_cells=('n_cells','sum'),
               mean_Equal4=('Equal4','mean'),median_Equal4=('Equal4','median')).reset_index()
    summary.to_csv(out/'celltype_summary.csv',index=False)
    rows=[]
    for g in ['Equal4']+GENES:
        pair=scored.pivot(index='patient',columns='cell_type',values=g)[['Macro','Tumor']].dropna()
        delta=pair.Macro-pair.Tumor
        rows.append(dict(score=g,n_patients=len(pair),mean_macrophage=float(pair.Macro.mean()),
            mean_tumor=float(pair.Tumor.mean()),mean_difference=float(delta.mean()),
            median_difference=float(delta.median()),p=exact_signed_rank(pair.Macro,pair.Tumor)))
    result=pd.DataFrame(rows);result['q']=multipletests(result.p,method='fdr_bh')[1]
    result.to_csv(out/'macrophage_tumor_paired.csv',index=False)
    return summary,result


def embedding(seed=20261110):
    """Recalculate the manuscript-described UMAP from public UMI counts."""
    from sklearn.decomposition import PCA
    from umap import UMAP
    raw=ROOT/'data/derived/GSE159115'
    x=sparse.load_npz(raw/'GSE159115_ccRCC_cell_gene_counts.npz').astype(float)
    obs=pd.read_csv(raw/'GSE159115_ccRCC_cell_metadata.csv.gz').set_index('cell')
    genes=pd.read_csv(raw/'GSE159115_gene_order.csv').gene_symbol
    lib=np.asarray(x.sum(axis=1)).ravel()
    log=x.multiply((1e6/lib)[:,None]).tocsr();log.data=np.log2(log.data+1)
    mean=np.asarray(log.mean(axis=0)).ravel()
    variance=np.asarray(log.power(2).mean(axis=0)).ravel()-mean**2
    selected=np.argsort(-variance,kind='stable')[:2000]
    a=log[:,selected].toarray();a=(a-a.mean(axis=0))/a.std(axis=0,ddof=0)
    pc=PCA(n_components=30,svd_solver='randomized',random_state=seed).fit_transform(a)
    coords=UMAP(n_neighbors=15,min_dist=.1,n_components=2,n_epochs=300,
                random_state=seed,n_jobs=1,compatibility_layout=True).fit_transform(pc)
    out=ROOT/'results/GSE159115'
    z=[]
    for g in GENES:
        ix=np.flatnonzero(genes.eq(g))
        if len(ix)!=1:raise ValueError(f'Ambiguous score gene: {g}')
        v=log[:,ix[0]].toarray().ravel();obs[g+'_log2CPM1']=v
        z.append((v-v.mean())/v.std(ddof=0))
    obs['Equal4']=np.mean(z,axis=0);obs['UMAP1']=coords[:,0];obs['UMAP2']=coords[:,1]
    obs.to_csv(out/'UMAP_coordinates_scores.csv.gz')
    pd.DataFrame({'gene_symbol':genes.iloc[selected].values,'variance_log2CPM1':variance[selected]}).to_csv(out/'UMAP_variable_genes.csv',index=False)
    np.save(out/'UMAP_PC_scores.npy',pc)
    (out/'UMAP_methods.json').write_text(json.dumps(dict(seed=seed,n_cells=len(obs),n_variable_genes=len(selected),
       normalization='log2(CPM+1)',variance_selection='all cells, descending variance, stable ties',
       gene_scaling='population SD',pca_solver='randomized',n_pcs=30,n_neighbors=15,min_dist=.1,n_epochs=300,
       compatibility_layout=True,status='new computation of manuscript-described embedding; original generating implementation absent'),indent=2)+'\n')
    return obs

def compartments():
    """Tumor-only broad compartments; no inferred normal-tissue labels."""
    root=ROOT/'data/derived/GSE159115';x=sparse.load_npz(root/'GSE159115_ccRCC_cell_gene_counts.npz');obs=pd.read_csv(root/'GSE159115_ccRCC_cell_metadata.csv.gz',index_col=0);symbols=pd.read_csv(root/'GSE159115_gene_order.csv').gene_symbol
    mapping={'Macro':'Macrophages','Macro_MKI67':'Macrophages','Tcell':'T cells','Tcell_CD8':'T cells','Tumor':'Tumor epithelial','Bcell':'B cells','Endo_ACKR1':'Endothelial','Endo_PLVAP':'Endothelial'}
    obs['compartment']=obs.anno.replace(mapping);keep=~obs.cellfilt.astype(bool)&~obs.doublet.astype(bool);records=[];vectors=[]
    for (patient,label),group in obs[keep].groupby(['patient','compartment']):
        if len(group)<30:continue
        v=np.asarray(x[obs.index.get_indexer(group.index)].sum(axis=0)).ravel();vectors.append(v);records.append(dict(patient=patient,compartment=label,n_cells=len(group)))
    e=pd.DataFrame(vectors,columns=symbols).T.groupby(level=0).sum().T;log=np.log2(e.div(e.sum(axis=1),axis=0)*1e6+1);z=(log[GENES]-log[GENES].mean())/log[GENES].std(ddof=0)
    d=pd.concat([pd.DataFrame(records),z],axis=1);d['Equal4']=z.mean(axis=1);out=ROOT/'results/GSE159115';d.to_csv(out/'broad_compartment_scores.csv',index=False);rows=[];rng=np.random.default_rng(20261110)
    for score in ['Equal4']+GENES:
        w=d.pivot(index='patient',columns='compartment',values=score)
        for a,b in itertools.combinations(w,2):
            pair=w[[a,b]].dropna()
            if len(pair)<5:continue
            delta=pair[a]-pair[b];boot=rng.choice(delta.to_numpy(),size=(5000,len(delta))).mean(axis=1)
            rows.append(dict(score=score,left=a,right=b,n=len(pair),mean_difference=delta.mean(),ci_low=np.quantile(boot,.025),ci_high=np.quantile(boot,.975),p=exact_signed_rank(pair[a],pair[b])))
    d=pd.DataFrame(rows);d['q']=d.groupby('score').p.transform(lambda v:multipletests(v,method='fdr_bh')[1]);d.to_csv(out/'broad_compartment_tests.csv',index=False)
