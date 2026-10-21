"""Small independent hand-worked checks of censoring weights and concordance."""
import json
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.utils import concordance_index
from run_strengthening import weights,concord,horizon_metrics,OUT

checks=[]
def check(name,ok):
    checks.append({'check':name,'passed':bool(ok)})
    assert ok,name
d=pd.DataFrame({'T':[1.,2.,3.,4.],'E':[1,0,1,0]})
w,y,ct,g=weights(d,d,3.)
check('hand calculated censor survival',np.isclose(g,2/3))
check('hand calculated IPCW weights',np.allclose(w,[1,0,1.5,1.5]))
check('IPCW risk equals KM risk without event censor ties',np.isclose(np.sum(w*y)/4,.625) and np.isclose(1-KaplanMeierFitter().fit(d['T'],d.E).predict(3),.625))
lp=np.array([3.,2.,1.,0.])
check('comparable pair accumulation equals lifelines C',np.isclose(concord(d,lp)[0],concordance_index(d['T'],-lp,d.E)))
unc=pd.DataFrame({'T':[1.,2.,4.,5.],'E':[1,1,1,1]})
r=np.array([.9,.7,.2,.1]); mm=horizon_metrics(unc,unc,r,r,3.,'test')
check('Brier reduces to ordinary squared error with no censoring',np.isclose(mm['brier_molecular'],np.mean((np.array([1,1,0,0])-r)**2)))
check('AUC equals 1 for hand worked perfect ranking',np.isclose(mm['auc_molecular'],1.))
tie=pd.DataFrame({'T':[1.,1.,2.,3.],'E':[1,0,1,0]})
w,y,ct,g=weights(tie,tie,2.)
check('events use G(T minus) at tied censor times',np.isclose(w[0],1.))
check('control uses horizon censor survival',np.isclose(w[-1],1/g))
(OUT/'statistical_checks.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
print(json.dumps(checks,indent=2))
