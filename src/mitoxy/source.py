"""GSE235046 count-based differential expression and pathway enrichment.

Ported from the revision's 02_GSE235046_Contrasts notebook. This is the
rounded-RSEM DE workflow; source-score inference is kept in its original notebook.
"""
import json, os
import numpy as np
import pandas as pd
from .data import ROOT, check_source

CONTRASTS = {'LPS_vs_Media':[-1,1,0,0,0], 'CS_vs_Media':[-1,0,1,0,0],
 'Combined_vs_LPS':[0,-1,0,1,0], 'Combined_vs_CS':[0,0,-1,1,0],
 'Interaction':[1,-1,-1,1,0], 'Torin_vs_Combined':[0,0,0,-1,1]}

def differential_expression():
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    from pydeseq2.default_inference import DefaultInference
    p=ROOT/'data/raw/GSE235046/GSE235046_Count_table.txt.gz'
    spec=next(s for s in json.loads((ROOT/'config/public_sources.json').read_text()) if s['path']=='GSE235046/GSE235046_Count_table.txt.gz')
    check_source(p,spec)
    meta=pd.read_csv(ROOT/'config/GSE235046_samples.csv',index_col='sample_id')
    raw=pd.read_csv(p,sep='\t')
    counts=raw.dropna(subset=['geneSymbol']).groupby('geneSymbol')[list(meta.index)].sum().T
    if not np.isfinite(counts.to_numpy()).all() or (counts<0).any().any():
        raise ValueError('Expected finite nonnegative RSEM counts')
    fractional=float((np.abs(counts.to_numpy()-counts.to_numpy().round())>1e-6).mean())
    counts=counts.round().astype(int)
    counts=counts.loc[:,counts.ge(5).sum().ge(3)]
    groups=['Media','LPS','CS','LPS_CS','LPS_CS_Torin']
    if set(meta.condition)!=set(groups) or not meta.condition.value_counts().eq(3).all():
        raise ValueError('Expected five conditions with three independent libraries each')
    design=pd.get_dummies(pd.Categorical(meta.condition,categories=groups),dtype=float)
    design.index=meta.index
    threads=int(os.environ.get('MITOXY_THREADS','1'))
    inference=DefaultInference(n_cpus=threads)
    dds=DeseqDataSet(counts=counts,metadata=meta,design=design,refit_cooks=True,inference=inference,quiet=True)
    dds.deseq2()
    out=ROOT/'results/source_de';out.mkdir(parents=True,exist_ok=True)
    for label,vector in CONTRASTS.items():
        ds=DeseqStats(dds,contrast=np.array(vector,dtype=float),inference=inference,quiet=True)
        ds.summary();d=ds.results_df.copy();d.index.name='gene'
        d['CI_low']=d.log2FoldChange-1.96*d.lfcSE;d['CI_high']=d.log2FoldChange+1.96*d.lfcSE
        d.to_csv(out/f'DE_{label}.csv')
    pd.DataFrame(dds.layers['normed_counts'],index=counts.index,columns=counts.columns).to_csv(out/'normalized_counts.csv.gz')
    (out/'methods.json').write_text(json.dumps(dict(method='PyDESeq2',count_source='GEO RSEM expected counts',
       rounding='sum by original symbol then nearest integer',fraction_fractional=fractional,
       filter='at least 5 counts in at least 3 libraries',n_genes=counts.shape[1],n_samples=counts.shape[0],
       design='five independent groups, no pairing or inferred batches',contrasts=CONTRASTS),indent=2)+'\n')

def pathways():
    import gseapy
    definitions={}
    for line in (ROOT/'config/genesets/enrichment_mouse.gmt').read_text().splitlines():
        f=line.split('\t');definitions[f[0]]=f[2:]
    out=ROOT/'results/source_de';rows=[]
    for contrast in CONTRASTS:
        d=pd.read_csv(out/f'DE_{contrast}.csv',index_col='gene')
        for exclude in [False,True]:
            rank=d.stat.replace([np.inf,-np.inf],np.nan).dropna()
            if exclude:rank=rank.drop(['Bak1','Bax','Bid','Slc7a11'],errors='ignore')
            rank=rank.sort_values(ascending=False,kind='stable')
            fit=gseapy.prerank(rnk=rank,gene_sets=definitions,min_size=5,max_size=500,
                permutation_num=1000,seed=20261109,threads=int(os.environ.get('MITOXY_THREADS','1')),
                outdir=None,verbose=False)
            table=fit.res2d.assign(contrast=contrast,exclude_core=exclude)
            table.to_csv(out/f'GSEA_{contrast}_exclude_core_{exclude}.csv',index=False);rows.append(table)
    pd.concat(rows,ignore_index=True).to_csv(out/'GSEA_all.csv',index=False)
