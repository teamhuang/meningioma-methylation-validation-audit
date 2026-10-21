import warnings
import numpy as np
from lifelines import CoxPHFitter
from lifelines.exceptions import ConvergenceWarning
from lifelines.utils.concordance import _concordance_summary_statistics

def fit(d,cols,pen=.05):
    active=[x for x in cols if d[x].nunique()>1]
    with warnings.catch_warnings():
        warnings.simplefilter('error',ConvergenceWarning)
        m=CoxPHFitter(penalizer=pen)
        m.fit(d[['T','E']+active],'T','E')
    if not np.all(np.isfinite(m.params_)):raise ValueError('nonfinite coefficients')
    return m

def concord(d,score):
    a,b,n=_concordance_summary_statistics(d['T'].to_numpy(),-np.asarray(score),d.E.to_numpy())
    if n==0:raise ValueError('no comparable pairs')
    return float((a+.5*b)/n),float(a+.5*b),int(n)
