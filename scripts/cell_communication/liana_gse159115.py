#!/usr/bin/env python3
"""Exploratory ligand-receptor analysis of Mitoxy-high vs Mitoxy-low macrophages in GSE159115 (Reviewer 2, comment 4).

Analysis specification (exploratory; not preregistered):
  * Cells: 20,748 author-annotated ccRCC cells, 7 patients. Author labels merged into major groups:
    Tumor; Macrophage (Macro); T cell (Tcell + Tcell_CD8); Endothelial (Endo_PLVAP + Endo_ACKR1);
    Mural (vSMC + Peri); B/Plasma (Bcell + Plasma); Mast. Excluded: 'ua' (unassigned) and Macro_MKI67 (37 cells).
  * Expression: counts normalized to 10,000 per cell and log1p-transformed.
  * Cell-level Equal4: mean of the four genes' z scores (BAK1, BAX, BID, SLC7A11), z computed within patient
    across all retained cells. Within each patient, macrophages are split at the median cell-level Equal4
    into Macro_MitoxyHigh and Macro_MitoxyLow.
  * LIANA rank_aggregate (consensus resource; 1,000 permutations; expr_prop 0.1; min_cells 10) run separately in
    each patient. Interactions with consensus specificity rank < 0.05 are called "detected".
  * Primary read-outs, source macrophage subset -> target T cells, per patient:
      (i) number of detected interactions; (ii) summed LR magnitude (lr_means) over the union of interactions.
    Compared High vs Low with exact two-sided paired Wilcoxon tests across patients (n=7).
  * Interaction-level: High minus Low lr_means per patient; sign consistency across patients reported
    descriptively (with n=7 the minimum exact paired P is 0.0156; directional consistency is descriptive,
    and multiple-testing correction requires the full tested family).
  * Secondary: macrophage immunoregulatory ligands and T-cell activation/exhaustion programs (pseudobulk per
    patient); library size (UMI) of the two macrophage subsets as a technical control.
Interpretation: exploratory and descriptive; patient is the unit of inference.
"""
import json, warnings, argparse
from paired_exact import exact_signed_rank_p
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp, scipy.stats as st
import anndata as ad, scanpy as sc, liana as li

warnings.filterwarnings('ignore')
parser = argparse.ArgumentParser(description='Patient-stratified LIANA analysis for GSE159115.')
parser.add_argument('--input-dir', required=True, help='Directory containing GSE159115_ccRCC_cell_gene_counts.npz, metadata CSV.gz, and gene-order CSV.')
parser.add_argument('--output-dir', default=None, help='Output directory. Defaults to outputs/ or outputs_umi_adjusted/.')
parser.add_argument('--umi-adjusted', action='store_true', help='Residualize cell-level Equal4 against log10 total UMI within patient before the macrophage split.')
args = parser.parse_args()
U = Path(args.input_dir).expanduser().resolve()
SENS = bool(args.umi_adjusted)
OUT = Path(args.output_dir) if args.output_dir else Path('outputs_umi_adjusted' if SENS else 'outputs')
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20261110

z = np.load(U / 'GSE159115_ccRCC_cell_gene_counts.npz', allow_pickle=False)
X = sp.csr_matrix((z['data'], z['indices'], z['indptr']), shape=tuple(z['shape'])).astype(np.float32)
meta = pd.read_csv(U / 'GSE159115_ccRCC_cell_metadata.csv.gz').set_index('cell')
genes = pd.read_csv(U / 'GSE159115_gene_order.csv')
A = ad.AnnData(X=X, obs=meta, var=pd.DataFrame(index=genes.gene_symbol.astype(str).values))
A.var_names_make_unique()
if A.n_obs != len(meta) or A.n_vars != len(genes): raise ValueError('Count/annotation dimensions differ')

merge = {'Tumor': 'Tumor', 'Macro': 'Macrophage', 'Tcell': 'T cell', 'Tcell_CD8': 'T cell', 'Endo_PLVAP': 'Endothelial',
         'Endo_ACKR1': 'Endothelial', 'vSMC': 'Mural', 'Peri': 'Mural', 'Bcell': 'B/Plasma', 'Plasma': 'B/Plasma', 'Mast': 'Mast'}
A = A[A.obs.anno.isin(merge)].copy(); A.obs['group'] = A.obs.anno.map(merge)
A.layers['counts'] = A.X.copy()
sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)

# cell-level Equal4 (z within patient across all retained cells)
E4 = ['BAK1', 'BAX', 'BID', 'SLC7A11']
G = pd.DataFrame(A[:, E4].X.toarray(), index=A.obs_names, columns=E4)
Gz = G.groupby(A.obs.patient).transform(lambda s: (s - s.mean()) / s.std(ddof=0))
A.obs['Equal4_cell'] = Gz.mean(1).values
if SENS:  # sensitivity: remove library-size dependence of cell-level Equal4 within each patient
    lu = np.log10(A.obs.total_UMI.astype(float))
    for p_, idx in A.obs.groupby('patient').groups.items():
        y_ = A.obs.loc[idx, 'Equal4_cell']; x_ = lu.loc[idx]; b_ = np.polyfit(x_, y_, 1)
        A.obs.loc[idx, 'Equal4_cell'] = y_ - np.polyval(b_, x_)

lab = A.obs.group.astype(str).copy()
for p, idx in A.obs[A.obs.group == 'Macrophage'].groupby('patient').groups.items():
    v = A.obs.loc[idx, 'Equal4_cell']; r = v.rank(method='first')
    lab.loc[idx] = np.where(r > len(v) / 2, 'Macro_MitoxyHigh', 'Macro_MitoxyLow')
A.obs['ccc_group'] = pd.Categorical(lab)

# technical control: library size of the two macrophage subsets
mac = A.obs[A.obs.ccc_group.astype(str).str.startswith('Macro_')]
umi = mac.groupby(['patient', 'ccc_group'], observed=True).agg(n=('total_UMI', 'size'), median_UMI=('total_UMI', 'median'),
                                                              median_genes=('no_genes', 'median'), mean_Equal4=('Equal4_cell', 'mean')).reset_index()
umi.to_csv(OUT / 'macrophage_subset_qc.csv', index=False)

# Pin and verify the packaged interaction resource before inference.
resource_path = Path(__file__).resolve().parents[2] / 'config/liana_consensus_1.10.0.csv'
resource = pd.read_csv(resource_path)
if not resource.equals(li.resource.select_resource('consensus')):
    raise RuntimeError('Installed LIANA resource differs from the recorded snapshot')

# LIANA per patient
res = []
for p in sorted(A.obs.patient.unique()):
    S = A[A.obs.patient == p].copy()
    li.mt.rank_aggregate(S, groupby='ccc_group', resource=resource, expr_prop=0.1, min_cells=10,
                         n_perms=1000, seed=SEED, use_raw=False, verbose=False)
    r = S.uns['liana_res'].copy(); r['patient'] = p; res.append(r)
R = pd.concat(res, ignore_index=True); R.to_csv(OUT / 'liana_rank_aggregate_all_patients.csv.gz', index=False)

src = ['Macro_MitoxyHigh', 'Macro_MitoxyLow']
M = R[R.source.isin(src) & (R.target == 'T cell')].copy()
M['detected'] = M.specificity_rank < 0.05
M['pair'] = M.ligand_complex + ' -> ' + M.receptor_complex

# (i)/(ii) per-patient summaries
summ = M.groupby(['patient', 'source']).agg(n_detected=('detected', 'sum'), sum_lr_means=('lr_means', 'sum'),
                                            n_tested=('pair', 'size')).reset_index()
summ.to_csv(OUT / 'macro_to_T_per_patient_summary.csv', index=False)
wide = summ.pivot(index='patient', columns='source')


def paired(a, b):
    d = (a - b).values
    p_exact = exact_signed_rank_p(a,b)
    return dict(n=len(d), n_high_greater=int((d > 0).sum()), n_low_greater=int((d < 0).sum()),
                median_high=float(np.median(a)), median_low=float(np.median(b)), median_diff=float(np.median(d)),
                p_exact=p_exact)


primary = {k: paired(wide[k]['Macro_MitoxyHigh'], wide[k]['Macro_MitoxyLow']) for k in ['n_detected', 'sum_lr_means']}

# interaction-level sign consistency (pairs tested in both subsets in >=5 patients)
P = M.pivot_table(index=['pair', 'patient'], columns='source', values='lr_means').dropna()
P['diff'] = P['Macro_MitoxyHigh'] - P['Macro_MitoxyLow']
cons = P.groupby('pair')['diff'].agg(n_patients='size', n_higher_in_High=lambda s: int((s > 0).sum()),
                                     n_higher_in_Low=lambda s: int((s < 0).sum()), median_diff='median').reset_index()
cons = cons[cons.n_patients >= 5]
det = M.groupby(['pair', 'source']).detected.sum().unstack(fill_value=0)
cons = cons.merge(det.rename(columns=lambda c: 'n_patients_detected_' + c.replace('Macro_Mitoxy', '')).reset_index(), on='pair', how='left')
cons['all_patients_same_direction'] = (cons.n_higher_in_High == cons.n_patients) | (cons.n_higher_in_Low == cons.n_patients)
cons = cons.sort_values(['all_patients_same_direction', 'median_diff'], ascending=[False, False])
cons.to_csv(OUT / 'macro_to_T_interaction_consistency.csv', index=False)

# secondary: pseudobulk marker programs per patient
def pb(mask, gl):
    S = A[mask]; gl = [g for g in gl if g in A.var_names]
    return pd.DataFrame(S[:, gl].X.toarray(), index=S.obs_names, columns=gl).groupby(S.obs.patient.values).mean()

LIG = ['CD274', 'PDCD1LG2', 'LGALS9', 'CD86', 'CD80', 'HLA-DRA', 'CXCL9', 'CXCL10', 'CXCL16', 'IL10', 'TGFB1', 'SPP1', 'C1QA', 'APOE', 'TNF', 'IL1B']
TAC = ['GZMB', 'PRF1', 'IFNG', 'GZMK', 'NKG7']; TEX = ['PDCD1', 'HAVCR2', 'LAG3', 'TIGIT', 'CTLA4', 'TOX']
mh = A.obs.ccc_group == 'Macro_MitoxyHigh'; ml = A.obs.ccc_group == 'Macro_MitoxyLow'
H, L = pb(mh, LIG), pb(ml, LIG)
lig = pd.DataFrame([dict(gene=g, **paired(H[g], L[g])) for g in H.columns])
lig.to_csv(OUT / 'macrophage_ligands_High_vs_Low.csv', index=False)

Tm = A.obs.group == 'T cell'
Tpb = pd.concat([pb(Tm, TAC).mean(1).rename('T_cytotoxic'), pb(Tm, TEX).mean(1).rename('T_exhaustion')], axis=1)
macE = A.obs[A.obs.group == 'Macrophage'].groupby('patient').Equal4_cell.mean().rename('macrophage_Equal4')
fracH = None
pt = Tpb.join(macE)
pt['n_T_cells'] = A.obs[Tm].groupby('patient').size()
pt.to_csv(OUT / 'patient_T_programs_vs_macrophage_Equal4.csv')
cor = {k: dict(zip(['rho', 'p'], st.spearmanr(pt.macrophage_Equal4, pt[k]))) for k in ['T_cytotoxic', 'T_exhaustion']}

ucheck = umi.pivot(index='patient', columns='ccc_group', values='median_UMI')
log = dict(n_cells=int(A.n_obs), n_patients=int(A.obs.patient.nunique()), primary=primary,
           umi_control=paired(ucheck['Macro_MitoxyHigh'], ucheck['Macro_MitoxyLow']),
           n_pairs_consistency=int(len(cons)), n_pairs_all_same_direction=int(cons.all_patients_same_direction.sum()),
           n_pairs_all_higher_in_High=int((cons.n_higher_in_High == cons.n_patients).sum()),
           n_pairs_all_higher_in_Low=int((cons.n_higher_in_Low == cons.n_patients).sum()),
           patient_level_T_program_correlations=cor,
           versions=dict(liana=li.__version__, scanpy=sc.__version__, anndata=ad.__version__))
json.dump(log, open(OUT / 'run_log.json', 'w'), indent=2, default=float)
print(json.dumps(log, indent=2, default=float))
