"""Validate the local comparable-pair concordance implementation against the archived
concordance function used for every reported C-index in this study.

The local implementation must reproduce lifelines' (and therefore the archive's)
C exactly when the IPCW weights are all one, and must equal 0.5 for random
scores and 1.0 for a perfectly risk-ordered score.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd

BASE = Path(r'E:\脑膜瘤npj')
IDENT = BASE / 'identity_preprocessing_20261019'
sys.path.insert(0, str(IDENT))
import donor_model_utilities as legacy
sys.path.insert(0, str(BASE / 'scripts_v4'))
from run_calibration_uno import comparable_pairs, uno_c, make_G

inputs = {
    'baylor_published': (IDENT / 'results/donor_aware/input_baylor_published_overlap_guard.csv',
                         ['grade2', 'gtr', 'age', 'male'], ['menG_B', 'menG_C']),
    'tubingen_published': (IDENT / 'results/donor_aware/input_tubingen_published_overlap_guard.csv',
                           ['grade', 'gtr', 'age', 'male'], ['tor_prolif', 'tor_hypermet']),
    'ucsf_published': (BASE / 'results_v4/input_ucsf_published_overlap_guard.csv',
                       ['grade', 'gtr', 'male'], ['hm', 'ie']),
    'hku_published': (IDENT / 'results/donor_aware/input_hku_published_overlap_guard.csv',
                      ['grade', 'gtr', 'age', 'male'], ['published_immune', 'published_hyper']),
}


def flat_G(x):
    return np.ones_like(np.atleast_1d(np.asarray(x, float)))


print(f'{"comparison":22s} {"archived C":>11s} {"local C (w=1)":>14s} {"diff":>11s}  verdict')
allok = True
for name, (path, clin, mol) in inputs.items():
    d = pd.read_csv(path)
    m = legacy.fit(d, clin)
    lp = m.predict_log_partial_hazard(d).to_numpy()
    leg = legacy.concord(d, lp)[0]
    loc = uno_c(d['T'].to_numpy(), d['E'].to_numpy(), lp, flat_G)
    same = abs(leg - loc) < 1e-10
    allok &= same
    print(f'{name:22s} {leg:11.6f} {loc:14.6f} {leg - loc:11.2e}  {"MATCH" if same else "MISMATCH"}')

# null and perfect-score checks
rng = np.random.default_rng(0)
d = pd.read_csv(inputs['ucsf_published'][0])
t, e = d['T'].to_numpy(), d['E'].to_numpy()
null = [uno_c(t, e, rng.normal(size=len(t)), flat_G) for _ in range(200)]
perfect = uno_c(t, e, -t, flat_G)          # higher risk = earlier event
print(f'\nrandom scores  : mean {np.mean(null):.4f} (expect 0.500), 95% range '
      f'{np.percentile(null, 2.5):.4f}-{np.percentile(null, 97.5):.4f}')
print(f'perfect score  : {perfect:.6f} (expect 1.000)')
null_ok = abs(np.mean(null) - 0.5) < 0.01 and np.percentile(null, 2.5) < 0.5 < np.percentile(null, 97.5)
# a perfectly time-ordered score gives 1.0 except for subjects censored at exactly
# the same time as an event, whose tied risk scores contribute 0.5
perfect_ok = perfect >= 0.999
allok &= null_ok and perfect_ok
print('\nVALIDATION', 'PASSED' if allok else 'FAILED')
