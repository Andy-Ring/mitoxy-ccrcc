"""EPIC TRef deconvolution from documented linear TPM and matched clinical models."""
import gzip,re,json,os,subprocess
import numpy as np
import pandas as pd
from .data import ROOT
from .analysis import GENES,zscore,model
from .bulk import prepare,read_gmt
from .survival import bh

def matrices():
    out=ROOT/'data/derived/EPIC';out.mkdir(parents=True,exist_ok=True)
    root=ROOT/'data/raw/TCGA-TPM';info=json.loads((root/'TCGA-KIRC.star_tpm.tsv.json').read_text())
    if info['unit']!='log2(tpm+1)':raise ValueError('Unexpected TCGA TPM units')
    e=pd.read_csv(root/'TCGA-KIRC.star_tpm.tsv.gz',sep='\t',index_col=0)
    e=e.loc[:,[c for c in e if c[13:15]=='01']];e=np.exp2(e)-1
    mapping=pd.read_csv(root/'gencode.v36.annotation.gtf.gene.probemap',sep='\t').set_index('id').gene
    e.index=e.index.map(mapping);e=e.loc[e.index.notna()].groupby(level=0).sum();e.columns=e.columns.str[:12];e=e.T.groupby(level=0).mean().T
    e.to_csv(out/'TCGA_TPM.tsv',sep='\t')
    root=ROOT/'data/raw/CPTAC-package';desc=pd.read_csv(next(root.glob('*sample_descriptions*')),sep='\t');desc=desc[desc.cohort.eq('CCRCC')&desc.tissue_type.eq('Tumor')].set_index('sample_id')
    t=pd.read_csv(next(root.glob('*rsem_transcripts_tpm*')),sep='\t');cols=t.columns.intersection(desc.index);mapping={}
    with gzip.open(next(root.glob('*gencode*')),'rt') as f:
        for line in f:
            if line.startswith('#'):continue
            gene=re.search('gene_id "([^"]+)"',line);symbol=re.search('gene_name "([^"]+)"',line)
            if gene and symbol:mapping[gene.group(1).split('.')[0]]=symbol.group(1)
    t['symbol']=t.gene_id.str.split('.').str[0].map(mapping);e=t.dropna(subset=['symbol']).groupby('symbol')[list(cols)].sum();e.columns=desc.loc[cols].GDC_id.values;e=e.T.groupby(level=0).mean().T;e.to_csv(out/'CPTAC_TPM.tsv',sep='\t')
    return out

def analyze():
    source=matrices();out=ROOT/'results/EPIC';out.mkdir(parents=True,exist_ok=True);rows=[];coverage=[]
    hall=read_gmt(ROOT/'data/raw/genesets/h.all.v2024.1.Hs.symbols.gmt')
    for cohort in ['TCGA','CPTAC']:
        subprocess.run([os.environ.get('MITOXY_RSCRIPT','Rscript'),str(ROOT/'scripts/run_epic.R'),str(source/f'{cohort}_TPM.tsv'),str(out/f'{cohort}_fractions.csv'),str(ROOT/'data/raw/EPIC/EPIC.tar.gz')],check=True,cwd=ROOT)
        frac=pd.read_csv(out/f'{cohort}_fractions.csv',index_col=0);qc=pd.read_csv(out/f'{cohort}_fractions.csv.fit.csv',index_col=0);good=qc.index[pd.to_numeric(qc.convergeCode,errors='coerce').eq(0)];frac=frac.loc[frac.index.intersection(good)]
        e,meta,cov,_=prepare(cohort);cov=[c for c in cov if c!='PC1'];d=pd.read_csv(ROOT/f'results/bulk/{cohort}_main_patient_scores.csv',index_col=0).join(frac,how='inner');gz=(e-e.mean())/e.std(ddof=0).replace(0,np.nan)
        for label in ['HALLMARK_INFLAMMATORY_RESPONSE','HALLMARK_TNFA_SIGNALING_VIA_NFKB']:d[label]=gz[sorted((set(hall[label])-set(GENES))&set(gz))].mean(axis=1)
        cells=[c for c in frac if c!='otherCells'];coverage.append(dict(cohort=cohort,estimated=len(qc),converged=len(good),matched=len(d)))
        d.to_csv(out/f'{cohort}_patient_scores.csv')
        for y in ['grade','stage','HALLMARK_INFLAMMATORY_RESPONSE','HALLMARK_TNFA_SIGNALING_VIA_NFKB']:
            cc=d.dropna(subset=[y,'Mitoxy']+cov+cells).copy();cc['Mitoxy']=zscore(cc.Mitoxy)
            if y.startswith('HALLMARK'):cc[y]=zscore(cc[y])
            for label,extra in [('base',[]),('macrophage',['Macrophages']),('major_fractions',cells)]:rows.append(dict(cohort=cohort,analysis=label,**model(cc,y,['Mitoxy']+cov+extra,'Mitoxy',cohort+'.'+y+'.'+label)[0]))
    d=pd.DataFrame(rows);d['q']=d.groupby(['cohort','outcome']).p.transform(lambda x:bh(x));d.to_csv(out/'adjusted_models.csv',index=False);pd.DataFrame(coverage).to_csv(out/'coverage.csv',index=False)
