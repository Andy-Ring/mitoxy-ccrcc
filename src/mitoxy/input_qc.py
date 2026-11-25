"""Validate GSE235046 public counts against independently downloaded GEO metadata."""
import json
import re
import numpy as np
import pandas as pd
from .data import ROOT,sha256

def validate_gse235046():
    raw=ROOT/'data/raw/GSE235046'
    soft=(raw/'GSE235046_family.soft').read_text()
    rows=[]
    for block in soft.split('^SAMPLE = ')[1:]:
        def get(key):
            match=re.search(r'^!Sample_'+key+r' = (.+)$',block,re.M)
            if match is None:raise ValueError(f'Missing GEO field {key}')
            return match.group(1).strip()
        treatment=re.search(r'^!Sample_characteristics_ch1 = treatment: (.+)$',block,re.M).group(1).strip()
        rows.append({'gsm':block.splitlines()[0].strip(),'sample_id':get('description'),
                     'title':get('title'),'treatment':treatment,'organism':get('organism_ch1')})
    metadata=pd.DataFrame(rows).set_index('sample_id')
    if not metadata.index.is_unique:raise ValueError('Duplicate source samples')
    count_file=raw/'GSE235046_Count_table.txt.gz'
    counts=pd.read_csv(count_file,sep='\t')
    matrix=counts.loc[:,metadata.index]
    if not np.isfinite(matrix.to_numpy()).all() or (matrix.to_numpy()<0).any():
        raise ValueError('Invalid public gene counts')
    missing=set(metadata.index)-set(counts.columns)
    if missing:raise ValueError(f'GEO sample missing from count table: {missing}')
    metadata['library_total_counts']=matrix.sum(axis=0)
    metadata['genes_detected']=(matrix>0).sum(axis=0)
    metadata.to_csv(ROOT/'data/derived/GSE235046_public_sample_QC.csv')
    summary={'gene_rows':len(counts),'samples':len(metadata),
             'condition_counts':metadata.treatment.value_counts().to_dict(),
             'count_sha256':sha256(count_file),
             'scope':'Public count matrix and sample-identity validation. The frozen Figure 1 DESeq2/score/matched-module workflow is not included.'}
    (ROOT/'results/GSE235046_input_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return metadata
