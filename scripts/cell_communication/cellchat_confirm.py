"""Method-sensitivity run of the CellChat method (LIANA implementation) on the prespecified Mitoxy-high/low macrophage split."""
import warnings, json, argparse, sys, numpy as np, pandas as pd, scipy.stats as st
from paired_exact import exact_signed_rank_p
from pathlib import Path
warnings.filterwarnings('ignore')
import liana as li
parser = argparse.ArgumentParser(description='CellChat-method sensitivity run through LIANA.')
parser.add_argument('--input-dir', required=True)
parser.add_argument('--output-dir', default='outputs')
args = parser.parse_args()
_here = Path(__file__).resolve().parent
_old_argv = sys.argv[:]
sys.argv = ['liana_gse159115.py', '--input-dir', args.input_dir, '--output-dir', args.output_dir]
exec((_here / 'liana_gse159115.py').read_text().split('# LIANA per patient')[0])  # identical data preparation and split
sys.argv = _old_argv
res=[]
for p in sorted(A.obs.patient.unique()):
    S=A[A.obs.patient==p].copy()
    li.mt.cellchat(S, groupby='ccc_group', resource=resource, expr_prop=0.1, min_cells=10, n_perms=1000, seed=SEED, use_raw=False, verbose=False)
    r=S.uns['liana_res'].copy(); r['patient']=p; res.append(r)
R=pd.concat(res,ignore_index=True); R.to_csv(OUT/'cellchat_all_patients.csv.gz',index=False)
M=R[R.source.isin(['Macro_MitoxyHigh','Macro_MitoxyLow'])&(R.target=='T cell')].copy(); M['pair']=M.ligand_complex+' -> '+M.receptor_complex
M['sig']=M.cellchat_pvals<0.05
s=M.groupby(['patient','source']).agg(n_sig=('sig','sum'),sum_prob=('lr_probs','sum')).reset_index(); s.to_csv(OUT/'cellchat_macro_to_T_summary.csv',index=False)
w=s.pivot(index='patient',columns='source')
def pr(a,b):
    d=(a-b).values; p_exact=exact_signed_rank_p(a,b)
    return dict(n_high_greater=int((d>0).sum()),n_low_greater=int((d<0).sum()),median_high=float(np.median(a)),median_low=float(np.median(b)),p_exact=p_exact)
out={k:pr(w[k]['Macro_MitoxyHigh'],w[k]['Macro_MitoxyLow']) for k in ['n_sig','sum_prob']}
P=M.pivot_table(index=['pair','patient'],columns='source',values='lr_probs').dropna(); P['d']=P.Macro_MitoxyHigh-P.Macro_MitoxyLow
c=P.groupby('pair').d.agg(n='size',hi=lambda x:int((x>0).sum()),lo=lambda x:int((x<0).sum()),med='median').reset_index()
c=c[c.n>=5]; c.to_csv(OUT/'cellchat_interaction_consistency.csv',index=False)
key=['C1QB -> C1QBP','C3 -> CD81','C3 -> IFITM1','CD14 -> ITGB2','CXCL12 -> CXCR4','MRC1 -> PTPRC','CXCL16 -> CXCR6','ICOSLG -> ICOS','ICOSLG -> CD28','ICAM2 -> ITGAL_ITGB2','CD55 -> ADGRE5']
out['key_axes']=c[c.pair.isin(key)].set_index('pair')[['n','hi','lo','med']].to_dict('index')
json.dump(out,open(OUT/'cellchat_run_log.json','w'),indent=1,default=float); print(json.dumps(out,indent=1,default=float))
