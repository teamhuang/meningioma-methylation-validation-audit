"""Cross-cohort external validation of a frozen risk model on a harmonised endpoint.

This is the level-C analysis that versions 1-4 of the manuscript state it does not
provide. It becomes possible only after (a) a common endpoint definition is fixed
in advance, (b) a common predictor set is fixed in advance, and (c) the molecular
predictor is supplied by ONE frozen implementation fitted elsewhere - here the
author's published 403-CpG KNN classifier applied unchanged to raw IDAT.

Harmonisation (fixed before any model was fitted)
--------------------------------------------------
  endpoint   : time to tumour recurrence or progression, in MONTHS. Tuebingen is
               already in months; UCSF and Hong Kong are converted from years (x12).
               The source endpoints are not identical (recurrence-free survival,
               progression, local freedom from recurrence); they are harmonised to
               a common "recurrence or progression" event, and the residual
               semantic mismatch is stated as a limitation, not hidden.
  predictors : WHO grade (ordinal) and extent of resection (binary) are shared and
               verified across all three cohorts. The frozen classifier's class-2
               and class-3 indicators are shared by construction, because the
               classifier is one frozen implementation applied to all three.
  age        : available in Tuebingen and Hong Kong but absent in UCSF, so it is
               excluded from the primary specification to keep the predictor set
               identical, and reported as a secondary specification.
  sex        : EXCLUDED from the primary specification. The UCSF clinical table
               carries a binary sex column whose direction the archive itself
               records as uninterpreted without source coding confirmation, and
               which does not agree with the GEO sex metadata for the same
               specimens. Including it in a frozen model produces an apparently
               inverted external discrimination, which is reported as a
               harmonisation control rather than as a property of the classifier.

  model      : Cox with ridge penalty 0.05, fitted ONCE on the development cohort
               and then applied with frozen coefficients to the external cohort.
               No refitting, no recalibration, no re-standardisation.
"""
from pathlib import Path
import sys
import time
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
import statsmodels.api as sm

BASE = Path(r'E:\脑膜瘤npj')
IDENT = BASE / 'identity_preprocessing_20261019/results'
DA = IDENT / 'donor_aware'
STR = BASE / 'strengthening_20261019/results'
OUT = BASE / 'results_v4'
OUT.mkdir(exist_ok=True)

PRIMARY = ['grade', 'gtr', 'author_2', 'author_3']
SECONDARY = ['grade', 'gtr', 'male', 'author_2', 'author_3']
PENALTY = 0.05
HORIZON_M = 60.0


def build_tubingen():
    d = pd.read_csv(DA / 'input_tubingen_author_KNN_raw_overlap_guard.csv')
    d['T_months'] = d['T'].astype(float)
    return d


def build_hongkong():
    d = pd.read_csv(DA / 'input_hku_author_KNN_raw_overlap_guard.csv')
    d['T_months'] = d['T'].astype(float) * 12.0
    return d


def build_ucsf():
    """UCSF clinical rows -> molecular identifiers -> raw-IDAT classifier output."""
    clin = pd.read_csv(STR / 'analysis_input_ucsf_published.csv')
    clin['idx'] = clin['id'].str.extract(r'row_(\d+)').astype(int)
    chain = pd.read_csv(IDENT / 'ucsf_old_positional_mapping_rejected.csv')
    chain = chain[['MID', 'GSM', 'M19_record_index']].rename(columns={'M19_record_index': 'idx'})
    chain = chain.drop_duplicates('idx', keep=False)          # drop the ambiguous row
    pred = pd.read_csv(IDENT / 'ucsf_raw_author_ID_linked_predictions.csv')
    pred = pred[~pred.tied_top.astype(bool)][['GSM', 'predicted_class']]
    snp = pd.read_csv(IDENT / 'SNP_donor_clusters_current.csv')
    snp = snp[snp.cohort == 'ucsf'][['GSM', 'donor_overlap_guard_key']]
    d = (clin.merge(chain, on='idx', how='inner')
              .merge(pred, on='GSM', how='inner')
              .merge(snp, on='GSM', how='inner'))
    d['T_months'] = d['T'].astype(float) * 12.0
    d['author_2'] = np.where(d.predicted_class.notna(), d.predicted_class.eq(2).astype(float), np.nan)
    d['author_3'] = np.where(d.predicted_class.notna(), d.predicted_class.eq(3).astype(float), np.nan)
    d = d.rename(columns={'donor_overlap_guard_key': 'donor'})
    return d[['id', 'GSM', 'T', 'E', 'donor', 'T_months', 'grade', 'gtr', 'male',
              'author_2', 'author_3']]


COHORTS = {'Tubingen': build_tubingen(), 'UCSF': build_ucsf(), 'Hong Kong': build_hongkong()}


def km_censor(t, e):
    kmf = KaplanMeierFitter().fit(np.asarray(t, float), np.asarray(e, float))
    sf = kmf.survival_function_.iloc[:, 0]
    x = sf.index.to_numpy(); y = sf.to_numpy()

    def G(v, left=False):
        v = np.atleast_1d(np.asarray(v, float))
        idx = np.searchsorted(x, v, side='left' if left else 'right') - 1
        out = np.where(idx < 0, 1.0, y[np.maximum(idx, 0)])
        return out if out.size > 1 else float(out[0])
    return G


def concord(t, e, score):
    t = np.asarray(t, float); e = np.asarray(e, float); score = np.asarray(score, float)
    I, J, C = [], [], []
    for i in np.where(e == 1)[0]:
        later = np.where(t > t[i])[0]
        same_c = np.where((t == t[i]) & (e == 0))[0]
        for j in np.concatenate([later, same_c]):
            I.append(i); J.append(j)
            C.append(1.0 if score[i] > score[j] else (0.5 if score[i] == score[j] else 0.0))
    if not I:
        return np.nan
    return float(np.sum(C) / len(C))


def fit_frozen(dev, cols):
    cols = [c for c in cols if c in dev.columns]
    d = dev[['T_months', 'E'] + cols].dropna()
    m = CoxPHFitter(penalizer=PENALTY).fit(d, 'T_months', 'E')
    sf = m.baseline_survival_.iloc[:, 0]
    times = sf.index.to_numpy(); surv = sf.to_numpy()
    keep = times <= HORIZON_M
    s0 = float(np.interp(HORIZON_M, times[keep], surv[keep])) if keep.sum() else float(surv[0])
    return {'coef': m.params_.to_dict(), 'active': cols, 'S0_60': s0,
            'dev_n': len(d), 'dev_events': int(d['E'].sum()),
            'dev_coefs': {k: round(float(v), 3) for k, v in m.params_.items()}}


def evaluate(frozen, test):
    active = [c for c in frozen['active'] if c in test.columns and test[c].nunique() > 1]
    lp = np.zeros(len(test))
    for c in active:
        lp += frozen['coef'][c] * test[c].to_numpy()
    t = test['T_months'].to_numpy(float)
    e = test['E'].to_numpy(float)
    C = concord(t, e, lp)
    risk = 1 - np.power(np.clip(frozen['S0_60'], 1e-8, 1.0), np.exp(lp))
    G = km_censor(t, e)
    ev = (t <= HORIZON_M) & (e == 1)
    ct = t > HORIZON_M
    w = np.zeros(len(test))
    ge = np.atleast_1d(G(t[ev], left=True))
    gh = float(np.atleast_1d(G([HORIZON_M]))[0])
    supported = gh > 0 and not np.any(ge <= 0)
    if supported:
        w[ev] = 1 / ge
        w[ct] = 1 / gh
        observed = float(np.sum(w * ev) / len(test))
        predicted = float(np.mean(risk))
        brier = float(np.sum(w * (ev - risk) ** 2) / len(test))
        op = observed / predicted if predicted > 0 else np.nan
        citl = observed - predicted
        pk = np.clip(risk, 1e-6, 1 - 1e-6)
        X = sm.add_constant(np.log(pk / (1 - pk)))
        try:
            fit = sm.GLM(ev.astype(float), X, family=sm.families.Binomial(), freq_weights=w).fit()
            slope, intercept = float(fit.params[1]), float(fit.params[0])
        except Exception:
            slope, intercept = np.nan, np.nan
    else:
        observed = predicted = brier = op = citl = slope = intercept = np.nan
    return {'n': len(test), 'events': int(e.sum()), 'c_index': C,
            'observed_60m': observed, 'predicted_60m': predicted,
            'observed_over_predicted': op, 'calibration_in_the_large': citl,
            'calibration_slope': slope, 'calibration_intercept': intercept,
            'brier_60m': brier, 'horizon_supported': bool(supported),
            'predictors_used': ','.join(active)}


def run_spec(spec_name, cols, rows):
    for dev_name, dev in COHORTS.items():
        frozen = fit_frozen(dev, cols)
        for test_name, test in COHORTS.items():
            if test_name == dev_name:
                continue
            need = [c for c in frozen['active'] if c in test.columns]
            t2 = test.dropna(subset=need + ['T_months', 'E']).copy()
            r = evaluate(frozen, t2)
            r.update({'specification': spec_name, 'development': dev_name, 'test': test_name,
                      'dev_n': frozen['dev_n'], 'dev_events': frozen['dev_events'],
                      'dev_coefficients': str(frozen['dev_coefs'])})
            rows.append(r)


def main():
    started = time.time()
    rows = []
    run_spec('primary (no sex, no age)', PRIMARY, rows)
    run_spec('control (sex included, unverified direction)', SECONDARY, rows)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / 'cross_cohort_external_validation.csv', index=False)

    pd.set_option('display.width', 250)
    for spec in df.specification.unique():
        s = df[df.specification == spec]
        print(f'\n=== {spec} ===')
        print(s[['development', 'test', 'n', 'events', 'c_index', 'observed_over_predicted',
                 'calibration_in_the_large', 'calibration_slope', 'brier_60m',
                 'predictors_used']].to_string(index=False))
        print(f'  C-index range              : {s.c_index.min():.3f} to {s.c_index.max():.3f}')
        print(f'  O/P range                   : {s.observed_over_predicted.min():.3f} to {s.observed_over_predicted.max():.3f}')
        print(f'  pairs with O/P in [0.8,1.25]: {int(((s.observed_over_predicted >= 0.8) & (s.observed_over_predicted <= 1.25)).sum())} of {len(s)}')
        print(f'  pairs with C > 0.5          : {int((s.c_index > 0.5).sum())} of {len(s)}')
    print(f'\nseconds {time.time() - started:.0f}')


if __name__ == '__main__':
    main()
