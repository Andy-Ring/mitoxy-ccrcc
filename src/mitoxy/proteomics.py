"""CPTAC pathway context and peptide-verified phosphosites from public deposits."""
import re,json
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from .data import ROOT
from .analysis import GENES,zscore,model
from .bulk import prepare
from .extensions import correlation
from .survival import bh

TARGETS=['RHEB','RPS6','MTOR','RICTOR','MAPKAP1','AKT1','AKT2','SGK1','NDRG1','RPS6KB1','EIF4EBP1']
SITES=json.loads((ROOT/'config/phosphosites.json').read_text())['sites']
ACCESSIONS={'NDRG1':'Q92597','RICTOR':'Q6R327','PRKCA':'P17252','MTOR':'P42345'}

def phospho():
    root=ROOT/'data/raw/CPTAC-package';raw=pd.read_csv(next(root.glob('*phosphoproteomics*')),sep='\t')
    fields=raw.Index.str.split('|',expand=True);site=fields[7].str.split('_',expand=True)
    meta=pd.DataFrame(dict(gene=fields[6],export_site=site[5],peptide=raw.Peptide,protein_id=fields[0],detected=pd.to_numeric(site[3]),localized=pd.to_numeric(site[4])))
    duplicate=meta.duplicated(['gene','export_site','peptide','protein_id'],keep=False)
    keep=meta.export_site.notna()&~(duplicate&meta.detected.ne(meta.localized));meta=meta[keep]
    mapping=pd.read_csv(next(root.glob('*aliquot_to_patient*')),sep='\t').set_index('aliquot_ID').patient_ID
    cols=[c for c in raw if c in mapping.index];mat=raw.loc[keep,cols].subtract(raw.loc[keep,'ReferenceIntensity'],axis=0).T
    mat.index=mapping.loc[cols].values
    mat=mat.loc[~mat.index.str.endswith(('.N','.C'))].groupby(level=0).mean()
    return mat,meta

def analyze():
    out=ROOT/'results/proteomics';out.mkdir(parents=True,exist_ok=True)
    e,meta,cov,_=prepare('CPTAC');ids=meta.index;e=e.loc[ids]
    protein=pd.read_csv(ROOT/'data/raw/cptac_ccrcc_public/HS_CPTAC_CCRCC_proteome_Tumor.txt',sep='\t',index_col=0).groupby(level=0).mean().T.reindex(ids)
    # Fixed clear-cell reference population for the correlations, matching manuscript layer analyses.
    equal=zscore(e[GENES]).mean(axis=1);rows=[]
    for layer,values in [('RNA',e),('protein',protein)]:
        for gene in TARGETS:rows.append(dict(layer=layer,gene=gene,**correlation(equal,values[gene],gene)))
    d=pd.DataFrame(rows);d['q']=d.groupby('layer').p.transform(lambda x:bh(x));d.to_csv(out/'mTOR_RNA_protein.csv',index=False)
    pd.DataFrame([correlation(zscore(e[GENES[:3]]).mean(axis=1),zscore(protein[GENES[:3]]).mean(axis=1),'Core3 RNA versus protein')]).to_csv(out/'Core3_RNA_protein.csv',index=False)
    mat,pm=phospho();seq={g:''.join((ROOT/f'data/raw/UniProt/{acc}.fasta').read_text().splitlines()[1:]) for g,acc in ACCESSIONS.items()}
    selected={};identities=[]
    for i,r in pm[pm.gene.isin(SITES)].iterrows():
        if not re.fullmatch('[STY][0-9]+',r.export_site):continue
        positions=[]
        for pep in r.peptide.split(';'):
            start=seq[r.gene].find(pep.upper())
            if start<0 or seq[r.gene].find(pep.upper(),start+1)>=0:continue
            positions.extend(pep[j].upper()+str(start+j+1) for j in range(len(pep)) if pep[j].islower())
        positions=sorted(set(positions));identities.append(dict(gene=r.gene,export_site=r.export_site,peptide=r.peptide,canonical_sites=';'.join(positions),accession=ACCESSIONS[r.gene],unique_site=len(positions)==1))
        if len(positions)==1 and positions[0] in SITES[r.gene]:selected.setdefault((r.gene,positions[0]),[]).append(mat[i])
    pd.DataFrame(identities).to_csv(out/'phosphosite_identity_mapping.csv',index=False);rows=[];measure={}
    for gene,sites in SITES.items():
        for site in sites:
            if (gene,site) not in selected:raise ValueError('No sequence-verified measurement for '+gene+' '+site)
            value=pd.concat(selected[gene,site],axis=1).mean(axis=1);measure[gene+'_'+site]=value
            rows.append(dict(layer='phospho',site=gene+'_'+site,**correlation(equal,value,gene+'_'+site)))
            cc=pd.DataFrame(dict(score=equal,phospho=value,protein=protein[gene])).dropna().rank();x=sm.add_constant(cc.protein)
            a=sm.OLS(cc.score,x).fit().resid;b=sm.OLS(cc.phospho,x).fit().resid;r=np.corrcoef(a,b)[0,1];df=len(cc)-3;p=2*stats.t.sf(abs(r)*np.sqrt(df/(1-r*r)),df)
            ci=np.tanh(np.arctanh(r)+np.array([-1,1])*1.95996398454/np.sqrt(len(cc)-4))
            rows.append(dict(layer='protein_adjusted_partial_rank',site=gene+'_'+site,n=len(cc),rho=r,p=p,ci_low=ci[0],ci_high=ci[1]))
    d=pd.DataFrame(rows);d['q']=d.groupby('layer').p.transform(lambda x:bh(x));d.to_csv(out/'seven_phosphosites.csv',index=False)
    pd.DataFrame(measure).reindex(ids).to_csv(out/'phosphosite_patient_values.csv')
    rows=[]
    for cohort in ['TCGA','CPTAC']:
        expression,meta,cov,_=prepare(cohort);cov=[c for c in cov if c!='PC1'];d=pd.read_csv(ROOT/f'results/bulk/{cohort}_main_patient_scores.csv',index_col=0);gz=(expression-expression.mean())/expression.std(ddof=0).replace(0,np.nan)
        for gene in ['MTOR','RICTOR','MAPKAP1','RPTOR','AKT1','AKT2','AKT3','SGK1','NDRG1']:
            d['target']=gz[gene];rows.append(dict(cohort=cohort,gene=gene,**model(d,'target',['Mitoxy']+cov,'Mitoxy',cohort+'.'+gene)[0]))
    d=pd.DataFrame(rows);d['q']=d.groupby('cohort').p.transform(lambda x:bh(x));d.to_csv(out/'adjusted_mTOR_RNA.csv',index=False)

def kinase_activity():
    out=ROOT/'results/proteomics';out.mkdir(parents=True,exist_ok=True);root=ROOT/'data/raw/kinases'
    net=pd.read_csv(root/'substrates.tsv',sep='\t');kinases=['MTOR','AKT1','AKT2','AKT3','SGK1','RPS6KB1'];net=net[net.enzyme_genesymbol.isin(kinases)&net.modification.eq('phosphorylation')&net.references.notna()].copy();net['site']=net.residue_type+net.residue_offset.astype(str)
    seq={}
    for f in sorted(root.glob('sequences_*.fasta')):
        for block in f.read_text().split('>')[1:]:
            lines=block.splitlines();seq[lines[0].split('|')[1]]=''.join(lines[1:])
    required=set(net.substrate)
    if not required.issubset(seq):raise ValueError('Missing reference sequences: '+str(sorted(required-set(seq))))
    mapping=net.groupby('substrate_genesymbol').substrate.first().to_dict();mat,pm=phospho();measure={};identity=[]
    for i,r in pm[pm.gene.isin(mapping)].iterrows():
        if not re.fullmatch('[STY][0-9]+',r.export_site):continue
        acc=mapping[r.gene];positions=[]
        for peptide in r.peptide.split(';'):
            start=seq.get(acc,'').find(peptide.upper())
            if start<0 or seq[acc].find(peptide.upper(),start+1)>=0:continue
            positions.extend(peptide[j].upper()+str(start+j+1) for j in range(len(peptide)) if peptide[j].islower())
        positions=sorted(set(positions));identity.append(dict(gene=r.gene,export_site=r.export_site,peptide=r.peptide,accession=acc,canonical_sites=';'.join(positions)))
        if len(positions)==1:measure.setdefault((r.gene,positions[0]),[]).append(mat[i])
    protein=pd.read_csv(ROOT/'data/raw/cptac_ccrcc_public/HS_CPTAC_CCRCC_proteome_Tumor.txt',sep='\t',index_col=0).groupby(level=0).mean().T
    e,meta,cov,_=prepare('CPTAC');d=pd.read_csv(ROOT/'results/bulk/CPTAC_main_patient_scores.csv',index_col=0);residuals={}
    for (gene,site),vs in measure.items():
        if gene not in protein:continue
        cc=pd.DataFrame(dict(phospho=pd.concat(vs,axis=1).mean(axis=1),protein=protein[gene])).dropna()
        if len(cc)<30:continue
        residuals[gene,site]=zscore(sm.OLS(cc.phospho,sm.add_constant(cc.protein)).fit().resid)
    rows=[];coverage=[];used=[]
    for kinase in kinases:
        edges=net[net.enzyme_genesymbol.eq(kinase)];keys=sorted(set(zip(edges.substrate_genesymbol,edges.site))&set(residuals));keys=[k for k in keys if len(residuals[k].index.intersection(e.index))>=.8*len(e)];genes=sorted(set(k[0] for k in keys));eligible=len(keys)>=5 and len(genes)>=3
        coverage.append(dict(kinase=kinase,sites=len(keys),substrate_proteins=len(genes),eligible=eligible))
        if not eligible:continue
        v=pd.DataFrame({g:pd.concat([residuals[k] for k in keys if k[0]==g],axis=1).mean(axis=1,skipna=False) for g in genes}).mean(axis=1,skipna=False)
        cc=d.join(v.rename('activity')).dropna(subset=['Mitoxy','activity']+cov).copy();cc.Mitoxy=zscore(cc.Mitoxy);cc.activity=zscore(cc.activity)
        rows.append(dict(kinase=kinase,**model(cc,'activity',['Mitoxy']+cov,'Mitoxy','CPTAC.'+kinase)[0]))
        for gene,site in keys:used.append(dict(kinase=kinase,gene=gene,site=site))
    pd.DataFrame(identity).to_csv(out/'kinase_peptide_mapping.csv',index=False);pd.DataFrame(coverage).to_csv(out/'kinase_coverage.csv',index=False);pd.DataFrame(used).to_csv(out/'kinase_used_sites.csv',index=False)
    if not rows:raise ValueError('No kinase met prespecified coverage thresholds')
    d=pd.DataFrame(rows);d['q']=bh(d.p);d.to_csv(out/'kinase_activity_models.csv',index=False)
