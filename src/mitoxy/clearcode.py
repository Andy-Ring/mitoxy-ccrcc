#!/usr/bin/env python3
"""ClearCode34 ccA/ccB classification of TCGA-KIRC, reimplementing the authors' classifier.

Source: KidneyCancerSubtypes (A. Culhane; Brooks et al., Eur Urol 2014), R/ClearCode34.R:
  training  : pamr.train(list(x = scale(xn), y = classes))      # 34 genes x 40 reference samples
  new data  : log2(x + 1) - Median533 (per gene); then scale() per sample
  prediction: pamr.predict(model, scale(newdata), threshold = 0) # nearest centroid, no shrinkage
pamr at threshold 0 (nsc/diag.disc): x* = (x - overall)/(sd + s0), delta_k = (centroid_k - overall)/(sd + s0),
score_k = x*.delta_k - |delta_k|^2/2 + log prior_k, posterior = softmax(score);
sd = pooled within-class SD, s0 = median(sd), prior = training class proportions.

The public training matrix and published median-scaling factors are downloaded by the workflow.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, pyreadr

from .data import ROOT
D = ROOT/'data/raw/ClearCode34'
xn = pyreadr.read_r(str(D / 'xn.rda'))['xn']
y = pyreadr.read_r(str(D / 'classes.rda'))['classes']['classes'].astype(str).values
fd = pyreadr.read_r(str(D / 'fData.model.Rda'))['fData.model']
med = pd.Series(np.array((D / 'Median533_scalingFactor.txt').read_text().split(), float), index=xn.index)


def rscale(m):  # R scale() on a genes x samples matrix: centre/scale each column (sample), SD with n-1
    return (m - m.mean(0)) / m.std(0, ddof=1)


# ---- train (pamr.train at threshold 0)
X = rscale(xn.to_numpy(float)); K = ['ccA', 'ccB']
cent = np.column_stack([X[:, y == k].mean(1) for k in K]); n = X.shape[1]
resid = X - cent[:, [K.index(c) for c in y]]
sd = np.sqrt((resid ** 2).sum(1) / (n - len(K))); s0 = np.median(sd); sdp = sd + s0
overall = X.mean(1); prior = np.array([(y == k).mean() for k in K])
delta = (cent - overall[:, None]) / sdp[:, None]


def predict(newlog2):  # newlog2: genes x samples, log2(expr+1), rows in xn order
    Z = rscale(newlog2.sub(med, axis=0).to_numpy(float))
    xs = (Z - overall[:, None]) / sdp[:, None]
    sc = xs.T @ delta - 0.5 * (delta ** 2).sum(0) + np.log(prior)
    sc -= sc.max(1, keepdims=True); p = np.exp(sc); p /= p.sum(1, keepdims=True)
    return pd.DataFrame({'ccA_prob': p[:, 0], 'ccB_prob': p[:, 1],
                         'ClearCode34': np.where(p[:, 0] >= .5, 'ccA', 'ccB')}, index=newlog2.columns)




def benchmark(calls):
    """Extract the fixed publisher DOC table's ASCII table-cell records.

    This is deliberately restricted to the checksum-locked source document;
    the record count, unique identifiers and probability columns are validated.
    """
    import re
    raw=(D/'Brooks2014_TableS3.doc').read_bytes().decode('latin1')
    pattern=r'(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})\x07([0-9.Ee+-]+)\x07([0-9.Ee+-]+)\x07(cc[AB])\x07'
    ref=pd.DataFrame(re.findall(pattern,raw),columns=['patient','pA','pB','published_class'])
    if len(ref)!=380 or ref.patient.duplicated().any():raise ValueError('Published ClearCode34 table extraction is incomplete or ambiguous')
    ref[['pA','pB']]=ref[['pA','pB']].astype(float)
    if not np.allclose(ref.pA+ref.pB,1,atol=.002):raise ValueError('Invalid published posterior probabilities')
    joined=calls.join(ref.set_index('patient'),how='inner',validate='one_to_one')
    joined.to_csv(ROOT/'results/TCGA_sensitivities/ClearCode34_reference_comparison.csv')
    return pd.DataFrame([dict(n_overlap=len(joined),n_concordant=int(joined.ClearCode34.eq(joined.published_class).sum()),fraction_concordant=float(joined.ClearCode34.eq(joined.published_class).mean()),ccA_probability_r=float(np.corrcoef(joined.ccA_prob,joined.pA)[0,1]))])
