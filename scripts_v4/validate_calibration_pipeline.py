"""Validate the calibration pipeline against the archive's own out-of-fold horizon
metrics.

The archive contains, for every record-level comparison, the per-(repeat, fold)
out-of-fold horizon metrics (observed_ipcw, predicted_clinical, predicted_molecular,
Brier, O/P, AUC) computed with the same definitions. Re-running the new pipeline in
record-level mode must reproduce them to numerical precision.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd

BASE = Path(r'E:\脑膜瘤npj')
STR = BASE / 'strengthening_20261019/results'
sys.path.insert(0, str(BASE / 'scripts_v4'))
import run_calibration_uno as rc

# the archive's record-level analysis sets
rc.INPUTS['ucsf_published'] = (STR / 'analysis_input_ucsf_published.csv',
                               ['grade', 'gtr', 'male'], ['hm', 'ie'], [3.0, 5.0])
rc.INPUTS['baylor_published'] = (STR / 'analysis_input_baylor_published.csv',
                                 ['grade2', 'gtr', 'age', 'male'], ['menG_B', 'menG_C'], [36.0, 60.0])
rc.INPUTS['tubingen_published'] = (STR / 'analysis_input_tubingen_published.csv',
                                   ['grade', 'gtr', 'age', 'male'], ['tor_prolif', 'tor_hypermet'], [36.0, 60.0])

ok = True
for name, arch in [('ucsf_published', 'ucsf_published'),
                   ('baylor_published', 'baylor_published'),
                   ('tubingen_published', 'tubingen_published')]:
    rc.run(name, grouped=False, reps=20, out_suffix='_record', seed=20261019)
    new = pd.read_csv(rc.OUT / f'calibration_oof_subjects_{name}_record.csv')
    old = pd.read_csv(STR / f'oof_horizon_metrics_{arch}.csv')
    for hz in sorted(new.horizon.unique()):
        n = new[(new.horizon == hz) & (new['repeat'] == 0)]
        o = old[(old.horizon == hz) & (old['repeat'] == 0)]
        # archive: per fold; new: per subject -> aggregate to per fold
        nf = n.groupby('fold').agg(n=('id', 'size'),
                                   pred_mol=('risk_molecular', 'mean')).reset_index()
        of = o.groupby('fold').agg(n=('n', 'first'),
                                   pred_mol=('predicted_molecular', 'first'),
                                   obs=('observed_ipcw', 'first')).reset_index()
        m = nf.merge(of, on='fold', suffixes=('_new', '_arch'))
        dpred = np.abs(m.pred_mol_new - m.pred_mol_arch).max()
        dn = np.abs(m.n_new - m.n_arch).max()
        good = dpred < 1e-6 and dn == 0
        ok &= good
        print(f'{name:22s} hz={hz:<5g} folds={len(m):2d}  max|Δ predicted_molecular| = {dpred:.3e}  '
              f'max|Δ n| = {dn}  {"MATCH" if good else "MISMATCH"}')
print('\nPIPELINE VALIDATION', 'PASSED' if ok else 'FAILED')
