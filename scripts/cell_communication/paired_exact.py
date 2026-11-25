import itertools
import numpy as np
from scipy.stats import rankdata

def exact_signed_rank_p(a,b):
    """Two-sided conditional signed-rank randomization; zeros dropped, ties midranked."""
    d=np.asarray(a,dtype=float)-np.asarray(b,dtype=float)
    d=d[np.isfinite(d) & (d!=0)]
    if not len(d):return 1.0
    if len(d)>20:raise ValueError('Exact enumeration restricted to <=20 nonzero pairs')
    ranks=rankdata(np.abs(d),method='average')
    observed=abs(np.dot(np.sign(d),ranks))
    null=np.array([abs(np.dot(signs,ranks)) for signs in itertools.product([-1,1],repeat=len(d))])
    return float(np.mean(null>=observed-1e-12))
