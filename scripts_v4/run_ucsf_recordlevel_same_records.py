"""Record-level counterpart of the UCSF donor-grouped analysis, on the identical
198-record analysis set.

Purpose: the deployment-audit counterfactual table compares estimates with and without
donor grouping. For Baylor and Tuebingen both estimates already exist on the same
records in the archive; for UCSF the record-level archive uses 201 records whereas the
donor-grouped analysis uses 198. Running the record-level variant on the same 198
records isolates the effect of donor grouping from the effect of the identity-exclusion
ledger.
"""
from pathlib import Path
import sys, json, time
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold

BASE = Path(__file__).resolve().parent.parent
IDENT = BASE / 'identity_preprocessing_20261019'
sys.path.insert(0, str(IDENT))
import donor_model_utilities as legacy

OUT = BASE / 'results_v4'
OUT.mkdir(exist_ok=True)
SEED = 2026101902
CLIN = ['grade', 'gtr', 'male']
MOL = ['hm', 'ie']

d = pd.read_csv(OUT / 'input_ucsf_published_overlap_guard.csv')
repeats, predictions, failures = [], [], []
for rep in range(20):
    cv = StratifiedKFold(5, shuffle=True, random_state=SEED + rep)
    tot = np.zeros(3)
    for fold, (tr, te) in enumerate(cv.split(d, d.E)):
        a, b = d.iloc[tr], d.iloc[te]
        try:
            m0 = legacy.fit(a, CLIN)
            m1 = legacy.fit(a, CLIN + MOL)
            x = m0.predict_log_partial_hazard(b).to_numpy()
            y = m1.predict_log_partial_hazard(b).to_numpy()
            _, c, n = legacy.concord(b, x)
            _, cc, nn = legacy.concord(b, y)
            assert n == nn
            tot += np.array([c, cc, n])
            for (_, r), s0, s1 in zip(b.iterrows(), x, y):
                predictions.append({'repeat': rep, 'fold': fold, 'id': r.id,
                                    'LP_clinical': s0, 'LP_molecular': s1})
        except Exception as exc:
            failures.append({'stage': 'CV', 'repeat': rep, 'fold': fold, 'error': str(exc)})
    repeats.append({'repeat': rep, 'C_clinical': tot[0] / tot[2],
                    'C_molecular': tot[1] / tot[2],
                    'deltaC': (tot[1] - tot[0]) / tot[2]})

rp = pd.DataFrame(repeats)
rp.to_csv(OUT / 'cv_ucsf_published_recordlevel_same198.csv', index=False)
pd.DataFrame(predictions).to_csv(OUT / 'oof_ucsf_published_recordlevel_same198.csv', index=False)
pd.DataFrame(failures).to_csv(OUT / 'failures_ucsf_published_recordlevel_same198.csv', index=False)

r = {'comparison': 'ucsf_published_recordlevel_same198', 'records': len(d),
     'events': int(d.E.sum()), 'donor_clusters': d.donor.nunique(),
     'C_cv_clinical': rp.C_clinical.mean(), 'C_cv_molecular': rp.C_molecular.mean(),
     'delta_cv': rp.deltaC.mean(),
     'delta_cv_split_p025': rp.deltaC.quantile(.025),
     'delta_cv_split_p975': rp.deltaC.quantile(.975),
     'positive_repeats': int((rp.deltaC > 0).sum()),
     'n_failures': len(failures)}
pd.DataFrame([r]).to_csv(OUT / 'summary_ucsf_published_recordlevel_same198.csv', index=False)
print('DONE', json.dumps(r, indent=2), flush=True)