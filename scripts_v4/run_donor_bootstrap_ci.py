"""Donor-level bootstrap confidence interval for the primary cross-validated dC.

For each bootstrap resample of donor clusters, the ENTIRE validation procedure is
re-executed: the complete repeated cross-validation is re-run inside the resample,
and the resulting dC is recorded. The 2.5th-97.5th percentile of that distribution
is therefore a resampling interval for the primary estimand itself, not for the
apparent in-sample estimator (which is what the archived apparent bootstrap
interval describes).

Each resample is evaluated twice, with donor-grouped folds and with record-level
folds on the identical resampled records, so the donor-grouping counterfactual can
also be given an interval.

Usage:  python run_donor_bootstrap_ci.py <cohort> <layer>
"""
from pathlib import Path
import sys, json, time
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold

BASE = Path(r'E:\脑膜瘤npj')
IDENT = BASE / 'identity_preprocessing_20261019'
sys.path.insert(0, str(IDENT))
import donor_model_utilities as legacy

OUT = BASE / 'results_v4'
OUT.mkdir(exist_ok=True)
SEED = 2026101902
B = 1000          # bootstrap resamples of donor clusters
R = 5             # cross-validation repetitions inside each resample
CLIN = {'ucsf': ['grade', 'gtr', 'male']}
MOL = {'published': ['hm', 'ie']}

INPUTS = {
    'baylor_published': (IDENT / 'results/donor_aware/input_baylor_published_overlap_guard.csv',
                         ['grade2', 'gtr', 'age', 'male'], ['menG_B', 'menG_C']),
    'tubingen_published': (IDENT / 'results/donor_aware/input_tubingen_published_overlap_guard.csv',
                           ['grade', 'gtr', 'age', 'male'], ['tor_prolif', 'tor_hypermet']),
    'ucsf_published': (BASE / 'results_v4/input_ucsf_published_overlap_guard.csv',
                       ['grade', 'gtr', 'male'], ['hm', 'ie']),
    'hku_published': (IDENT / 'results/donor_aware/input_hku_published_overlap_guard.csv',
                      ['grade', 'gtr', 'age', 'male'], ['published_immune', 'published_hyper']),
}


def cv_delta(d, clin, mol, grouped, seed):
    """Mean dC over R repetitions of 5-fold CV on the supplied records."""
    tot = np.zeros(3)
    for rep in range(R):
        if grouped:
            cv = StratifiedGroupKFold(5, shuffle=True, random_state=seed + rep)
            split = cv.split(d, d.E, d.donor)
        else:
            cv = StratifiedKFold(5, shuffle=True, random_state=seed + rep)
            split = cv.split(d, d.E)
        for tr, te in split:
            a, b = d.iloc[tr], d.iloc[te]
            if b.E.sum() == 0 or a.E.sum() == 0:
                continue
            m0 = legacy.fit(a, clin)
            m1 = legacy.fit(a, clin + mol)
            x = m0.predict_log_partial_hazard(b).to_numpy()
            y = m1.predict_log_partial_hazard(b).to_numpy()
            _, c, n = legacy.concord(b, x)
            _, cc, nn = legacy.concord(b, y)
            if n == 0 or n != nn:
                continue
            tot += np.array([c, cc, n])
    if tot[2] == 0:
        return np.nan
    return float((tot[1] - tot[0]) / tot[2])


def run(name):
    path, clin, mol = INPUTS[name]
    d = pd.read_csv(path)
    started = time.time()
    rng = np.random.default_rng(SEED)
    units = d.donor.unique()
    by = {u: d[d.donor == u] for u in units}
    rows = []
    for b in range(B):
        drawn = rng.choice(units, len(units), replace=True)
        db = pd.concat([by[u] for u in drawn], ignore_index=True)
        try:
            dg = cv_delta(db, clin, mol, True, SEED)
            dr = cv_delta(db, clin, mol, False, SEED)
            rows.append({'bootstrap': b, 'n_records': len(db), 'n_events': int(db.E.sum()),
                         'delta_grouped': dg, 'delta_recordlevel': dr})
        except Exception as exc:
            rows.append({'bootstrap': b, 'n_records': len(db), 'n_events': int(db.E.sum()),
                         'delta_grouped': np.nan, 'delta_recordlevel': np.nan,
                         'error': str(exc)[:120]})
        if (b + 1) % 100 == 0:
            print(f'{name}: {b + 1}/{B}  elapsed {time.time() - started:.0f}s', flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / f'donor_bootstrap_ci_{name}.csv', index=False)
    g = df.delta_grouped.dropna()
    r = df.delta_recordlevel.dropna()
    summ = {'comparison': name, 'bootstrap_attempted': B, 'bootstrap_success': int(g.size),
            'records': len(d), 'events': int(d.E.sum()), 'donor_clusters': len(units),
            'delta_grouped_mean': float(g.mean()),
            'delta_grouped_ci_low': float(g.quantile(.025)),
            'delta_grouped_ci_high': float(g.quantile(.975)),
            'delta_grouped_frac_le_zero': float((g <= 0).mean()),
            'delta_recordlevel_mean': float(r.mean()),
            'delta_recordlevel_ci_low': float(r.quantile(.025)),
            'delta_recordlevel_ci_high': float(r.quantile(.975)),
            'seconds': time.time() - started}
    pd.DataFrame([summ]).to_csv(OUT / f'donor_bootstrap_ci_summary_{name}.csv', index=False)
    print('DONE', json.dumps(summ), flush=True)


if __name__ == '__main__':
    run(sys.argv[1])
