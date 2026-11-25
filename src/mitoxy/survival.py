"""Breslow Cox likelihood recovered from the study analysis code."""
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp
from scipy.stats import norm

def cox(time,event,X):
    time=np.asarray(time,float);event=np.asarray(event,float);X=np.asarray(X,float)
    if X.ndim==1:X=X[:,None]
    if len(time)!=len(event) or len(time)!=len(X):raise ValueError('Incompatible survival dimensions')
    if not (np.isfinite(X).all() and np.isfinite(time).all() and (time>0).all()):
        raise ValueError('Survival times must be positive and predictors finite')
    if not set(event).issubset({0,1}) or event.sum()==0:
        raise ValueError('Events must be binary, with at least one observed event')
    event=event.astype(int)
    groups=[(np.flatnonzero((time==t)&(event==1)),np.flatnonzero(time>=t)) for t in np.unique(time[event==1])]
    def calc(b, want_h=False):
        eta=X@b;ll=0.;g=np.zeros(X.shape[1]);h=np.zeros((X.shape[1],X.shape[1]))
        for es,rs in groups:
            k=len(es);v=eta[rs];ls=logsumexp(v);w=np.exp(v-ls);rx=X[rs];mu=w@rx
            ll+=eta[es].sum()-k*ls;g+=X[es].sum(axis=0)-k*mu
            h+=k*((rx.T*w)@rx-np.outer(mu,mu))
        return (-ll,-g,h) if want_h else (-ll,-g)
    fit=minimize(calc,np.zeros(X.shape[1]),jac=True,method='BFGS',options={'gtol':1e-8,'maxiter':1000})
    _,g,h=calc(fit.x,True)
    if not np.isfinite(g).all() or np.max(np.abs(g))>1e-5:raise RuntimeError('Cox score equation failed '+str(g))
    se=np.sqrt(np.diag(np.linalg.inv(h)));b=fit.x
    return {'beta':b,'se':se,'HR':np.exp(b),'ci_low':np.exp(b-1.95996398454*se),'ci_high':np.exp(b+1.95996398454*se),'p':2*norm.sf(np.abs(b/se)),'max_abs_score':float(np.max(np.abs(g)))}

def bh(p):
    p=np.asarray(p);o=np.argsort(p);q=np.empty(len(p));q[o]=np.minimum(1,np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return q

def km_median(t,e):
    s=1.
    for u in np.unique(t[e==1]):
        s*=1-np.sum((t==u)&(e==1))/np.sum(t>=u)
        if s<=.5:return float(u)
    return None
