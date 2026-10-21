"""Donor-grouped cross-validation for the UCSF/Toronto published-group comparison.

Version 4 of the manuscript reports UCSF at record level because donor grouping was
not applied. This script closes that gap using the identity-proven chain only:

  MOESM19 row index (verified by exact field equality against natgen205_survival.csv)
    -> MID   (ucsf_old_positional_mapping_rejected.csv, M19_record_index)
    -> GSM   (same file; MID<->GSM is ID-based via series title + filelist)
    -> SNP donor cluster (SNP_donor_clusters_current.csv, donor_overlap_guard_key)

The rejected positional bridge (fingerprint_col_index <-> M20_row_index) is NOT used.

Cross-validation, concordance and optimism code is imported unchanged from
donor_model_utilities.py, so the UCSF estimate is produced by the same routine as
the Baylor, Tuebingen and Hong Kong donor-aware analyses.
"""
from pathlib import Path
import sys, json, time
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold
from lifelines import CoxPHFitter

BASE = Path(__file__).resolve().parent.parent
IDENT = BASE / 'identity_preprocessing_20261019'
sys.path.insert(0, str(IDENT))
import donor_model_utilities as legacy

OUT = BASE / 'results_v4'
OUT.mkdir(exist_ok=True)
SEED = 2026101902
CLIN = ['grade', 'gtr', 'male']
MOL = ['hm', 'ie']


def build_input():
    a = pd.read_csv(BASE / 'strengthening_20261019/results/analysis_input_ucsf_published.csv')
    a['idx'] = a['id'].str.extract(r'row_(\d+)').astype(int)

    # verify the row correspondence before relying on it
    n = pd.read_csv(BASE / 'results/natgen205_survival.csv')
    m = a.merge(n, left_on='idx', right_index=True, suffixes=('_a', '_n'))
    assert len(m) == len(a)
    for c in ['T', 'E', 'grade', 'gtr', 'male', 'hm', 'ie']:
        assert np.allclose(m[c + '_a'].astype(float), m[c + '_n'].astype(float), equal_nan=True), c

    chain = pd.read_csv(IDENT / 'results/ucsf_old_positional_mapping_rejected.csv')
    chain = chain[['MID', 'GSM', 'M19_record_index']].rename(columns={'M19_record_index': 'idx'})

    # Two source defects must be handled explicitly, not silently:
    #   (i)  MOESM19 row 128 is claimed by two MIDs (M104, M162), both Immune-enriched;
    #   (ii) two complete-covariate rows (MOESM19 idx 56, 127) carry no MID at all.
    # Both are excluded, and the exclusion ledger is written out.
    ambiguous = chain[chain.duplicated('idx', keep=False)]
    dup_idx = sorted(ambiguous.idx.unique()) if len(ambiguous) else []
    unlinked = sorted(set(a.idx) - set(chain.idx))
    chain = chain.drop_duplicates('idx', keep=False)
    ledger = pd.concat([
        ambiguous.assign(reason='ambiguous: two MIDs claim this MOESM19 row'),
        pd.DataFrame({'MID': [np.nan] * len(unlinked), 'GSM': [np.nan] * len(unlinked),
                      'idx': unlinked,
                      'reason': 'no identifier-proven MID for this MOESM19 row'})
    ], ignore_index=True)
    ledger.to_csv(OUT / 'ucsf_identity_exclusion_ledger.csv', index=False)

    snp = pd.read_csv(IDENT / 'results/SNP_donor_clusters_current.csv')
    snp = snp[snp.cohort == 'ucsf'][['GSM', 'donor_overlap_guard_key', 'donor_strict_key',
                                      'callable_SNPs']]

    d = a.merge(chain, on='idx', how='left').merge(snp, on='GSM', how='left')
    d = d.rename(columns={'donor_overlap_guard_key': 'donor'})
    d = d[['id', 'MID', 'GSM', 'T', 'E', 'donor', 'donor_strict_key',
           'callable_SNPs'] + CLIN + MOL].dropna().reset_index(drop=True)
    print('excluded MOESM19 rows (ambiguous):', dup_idx, '| no MID:', unlinked,
          '| analysis records retained:', len(d), flush=True)
    return d


def analyse():
    started = time.time()
    d = build_input()
    clin, mol = CLIN, MOL
    d.to_csv(OUT / 'input_ucsf_published_overlap_guard.csv', index=False)

    c0 = legacy.fit(d, clin)
    c1 = legacy.fit(d, clin + mol)
    app0 = legacy.concord(d, c0.predict_log_partial_hazard(d))[0]
    app1 = legacy.concord(d, c1.predict_log_partial_hazard(d))[0]

    repeats, predictions, failures = [], [], []
    for rep in range(20):
        cv = StratifiedGroupKFold(5, shuffle=True, random_state=SEED + rep)
        tot = np.zeros(3)
        for fold, (tr, te) in enumerate(cv.split(d, d.E, d.donor)):
            a, b = d.iloc[tr], d.iloc[te]
            assert not set(a.donor) & set(b.donor)
            try:
                m0 = legacy.fit(a, clin)
                m1 = legacy.fit(a, clin + mol)
                x = m0.predict_log_partial_hazard(b).to_numpy()
                y = m1.predict_log_partial_hazard(b).to_numpy()
                _, c, n = legacy.concord(b, x)
                _, cc, nn = legacy.concord(b, y)
                assert n == nn
                tot += np.array([c, cc, n])
                for (_, r), s0, s1 in zip(b.iterrows(), x, y):
                    predictions.append({'repeat': rep, 'fold': fold, 'id': r.id,
                                        'donor': r.donor, 'LP_clinical': s0,
                                        'LP_molecular': s1})
            except Exception as exc:
                failures.append({'stage': 'CV', 'repeat': rep, 'fold': fold, 'error': str(exc)})
        repeats.append({'repeat': rep, 'C_clinical': tot[0] / tot[2],
                        'C_molecular': tot[1] / tot[2],
                        'deltaC': (tot[1] - tot[0]) / tot[2]})

    boot = []
    rng = np.random.default_rng(SEED)
    units = d.donor.unique()
    by = {x: d[d.donor == x] for x in units}
    nboot = 1000
    for k in range(nboot):
        drawn = rng.choice(units, len(units), replace=True)
        db = pd.concat([by[x] for x in drawn], ignore_index=True)
        try:
            a = legacy.fit(db, clin)
            b = legacy.fit(db, clin + mol)
            t0 = legacy.concord(db, a.predict_log_partial_hazard(db))[0]
            t1 = legacy.concord(db, b.predict_log_partial_hazard(db))[0]
            e0 = legacy.concord(d, a.predict_log_partial_hazard(d))[0]
            e1 = legacy.concord(d, b.predict_log_partial_hazard(d))[0]
            boot.append({'bootstrap': k, 'n_records': len(db), 'n_events': int(db.E.sum()),
                         'delta_app': t1 - t0, 'optimism_clinical': t0 - e0,
                         'optimism_molecular': t1 - e1,
                         'optimism_delta': (t1 - t0) - (e1 - e0)})
        except Exception as exc:
            failures.append({'stage': 'bootstrap', 'bootstrap': k, 'error': str(exc)})

    rp = pd.DataFrame(repeats)
    bp = pd.DataFrame(boot)
    for obj, suffix in [(rp, 'cv'), (pd.DataFrame(predictions), 'oof'),
                        (bp, 'bootstrap'), (pd.DataFrame(failures), 'failures')]:
        obj.to_csv(OUT / (suffix + '_ucsf_published_overlap_guard.csv'), index=False)

    r = {'comparison': 'ucsf_published_overlap_guard', 'records': len(d),
         'events': int(d.E.sum()), 'donor_clusters': d.donor.nunique(),
         'C_app_clinical': app0, 'C_app_molecular': app1, 'delta_app': app1 - app0,
         'C_cv_clinical': rp.C_clinical.mean(), 'C_cv_molecular': rp.C_molecular.mean(),
         'delta_cv': rp.deltaC.mean(),
         'delta_cv_split_p025': rp.deltaC.quantile(.025),
         'delta_cv_split_p975': rp.deltaC.quantile(.975),
         'positive_repeats': int((rp.deltaC > 0).sum()),
         'bootstrap_success': len(bp), 'bootstrap_attempted': nboot,
         'n_failures': len(failures)}
    if len(bp):
        r.update(delta_app_CI_low=bp.delta_app.quantile(.025),
                 delta_app_CI_high=bp.delta_app.quantile(.975),
                 C_corrected_clinical=app0 - bp.optimism_clinical.mean(),
                 C_corrected_molecular=app1 - bp.optimism_molecular.mean(),
                 delta_corrected=(app1 - app0) - bp.optimism_delta.mean())
    try:
        m = CoxPHFitter().fit(d[['T', 'E', 'donor'] + clin + mol], 'T', 'E',
                              cluster_col='donor', robust=True)
        m.summary.to_csv(OUT / 'cluster_robust_cox_ucsf_published_overlap_guard.csv')
    except Exception as exc:
        r['cluster_robust_cox_error'] = str(exc)

    pd.DataFrame([r]).to_csv(OUT / 'summary_ucsf_published_overlap_guard.csv', index=False)
    print('DONE', json.dumps(r, indent=2), 'seconds', time.time() - started, flush=True)


if __name__ == '__main__':
    analyse()