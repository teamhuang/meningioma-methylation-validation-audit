"""Per-subject out-of-fold horizon predictions, Uno's C, time-dependent AUC and
calibration diagnostics.

The horizon weights, risk computation, observed-risk estimator and Brier score
follow the archived definitions in strengthening_20261019/run_strengthening.py
exactly (censoring Kaplan-Meier estimated in the TRAINING fold; event-first
handling of event/censor ties; risk from predict_survival_function at the horizon).
The pipeline is validated by reproducing the archive's own out-of-fold horizon
metrics for a record-level comparison before it is used for the donor-grouped
analyses.

Outputs per comparison:
  uno_c_{name}.csv                      Uno's C and time-dependent AUC per repetition
  calibration_oof_subjects_{name}.csv   per-subject out-of-fold horizon predictions
  calibration_summary_{name}.csv        calibration-in-the-large, slope, intercept
  calibration_curve_{name}_{model}_{h}.csv  predicted-quintile calibration curve
"""
from pathlib import Path
import sys, json, time
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold
from lifelines import KaplanMeierFitter
import statsmodels.api as sm

BASE = Path(r'E:\脑膜瘤npj')
IDENT = BASE / 'identity_preprocessing_20261019'
sys.path.insert(0, str(IDENT))
import donor_model_utilities as legacy

OUT = BASE / 'results_v4'
OUT.mkdir(exist_ok=True)
SEED = 2026101902
REPS = 20

INPUTS = {
    'baylor_published': (IDENT / 'results/donor_aware/input_baylor_published_overlap_guard.csv',
                         ['grade2', 'gtr', 'age', 'male'], ['menG_B', 'menG_C'], [36.0, 60.0]),
    'tubingen_published': (IDENT / 'results/donor_aware/input_tubingen_published_overlap_guard.csv',
                           ['grade', 'gtr', 'age', 'male'], ['tor_prolif', 'tor_hypermet'], [36.0, 60.0]),
    'ucsf_published': (BASE / 'results_v4/input_ucsf_published_overlap_guard.csv',
                       ['grade', 'gtr', 'male'], ['hm', 'ie'], [3.0, 5.0]),
    'hku_published': (IDENT / 'results/donor_aware/input_hku_published_overlap_guard.csv',
                      ['grade', 'gtr', 'age', 'male'], ['published_immune', 'published_hyper'], [3.0, 5.0]),
}


def G_at(train, times, left=False):
    """KM censoring survival with event-first handling of event/censor ties."""
    km = KaplanMeierFitter().fit(train['T'], 1 - train['E'])
    x = km.survival_function_.index.to_numpy()
    y = km.survival_function_.iloc[:, 0].to_numpy()
    idx = np.searchsorted(x, np.asarray(times, float),
                          side='left' if left else 'right') - 1
    return np.where(idx < 0, 1.0, y[np.maximum(idx, 0)])


def weights(train, test, h):
    """Archived horizon weights: event-by-horizon and still-at-risk subjects."""
    ev = (test['T'].to_numpy() <= h) & (test['E'].to_numpy() == 1)
    ct = test['T'].to_numpy() > h
    w = np.zeros(len(test))
    ge = G_at(train, test.loc[ev, 'T'].to_numpy(), True)
    gh = float(G_at(train, [h])[0])
    if h > train['T'].max() or gh <= 0 or np.any(ge <= 0):
        return None, ev, ct, gh
    w[ev] = 1 / ge
    w[ct] = 1 / gh
    return w, ev, ct, gh


def comparable_pairs(t, e, risk):
    """Comparable-pair structure matching lifelines' _concordance_summary_statistics."""
    t = np.asarray(t, float)
    e = np.asarray(e, float)
    risk = np.asarray(risk, float)
    I, J, C = [], [], []
    for i in np.where(e == 1)[0]:
        later = np.where(t > t[i])[0]
        same_censored = np.where((t == t[i]) & (e == 0))[0]
        for j in np.concatenate([later, same_censored]):
            I.append(i); J.append(j)
            C.append(1.0 if risk[i] > risk[j] else (0.5 if risk[i] == risk[j] else 0.0))
    return np.array(I, int), np.array(J, int), np.array(C, float)


def uno_c(t, e, score, gvals):
    """Uno's C: Harrell's comparable pairs with IPCW weight 1/G(t_i-) on the event subject.

    gvals must be the censoring survival evaluated at each subject's event time
    using the left-continuous convention, matching the archived weights().
    """
    I, J, C = comparable_pairs(t, e, score)
    if len(I) == 0:
        return np.nan
    w = 1.0 / np.clip(np.asarray(gvals, float)[I], 1e-6, None)
    return float(np.sum(w * C) / np.sum(w))


def time_dependent_auc(t, e, risk, horizon, gvals):
    """Cumulative (Horvitz) time-dependent AUC at a horizon."""
    t = np.asarray(t, float)
    e = np.asarray(e, float)
    risk = np.asarray(risk, float)
    gvals = np.asarray(gvals, float)
    num = den = 0.0
    for i in np.where((e == 1) & (t <= horizon))[0]:
        j = np.where(t > t[i])[0]
        if len(j) == 0:
            continue
        w = 1.0 / max(gvals[i], 1e-6)
        con = np.where(risk[i] > risk[j], 1.0, np.where(risk[i] == risk[j], 0.5, 0.0))
        num += w * con.sum()
        den += w * len(j)
    return num / den if den > 0 else np.nan


def run(name, grouped=True, reps=REPS, out_suffix='', seed=SEED):
    path, clin, mol, horizons = INPUTS[name]
    d = pd.read_csv(path).reset_index(drop=True)
    started = time.time()
    pred_rows, uno_rows = [], []
    for rep in range(reps):
        if grouped:
            cv = StratifiedGroupKFold(5, shuffle=True, random_state=seed + rep)
            split = cv.split(d, d['E'], d['donor'])
        else:
            cv = StratifiedKFold(5, shuffle=True, random_state=seed + rep)
            split = cv.split(d, d['E'])
        oof0 = np.full(len(d), np.nan)
        oof1 = np.full(len(d), np.nan)
        for fold, (tr, te) in enumerate(split):
            a, b = d.iloc[tr], d.iloc[te]
            if a['E'].sum() == 0 or b['E'].sum() == 0:
                continue
            try:
                m0 = legacy.fit(a, clin)
                m1 = legacy.fit(a, clin + mol)
            except Exception:
                continue
            oof0[te] = m0.predict_log_partial_hazard(b).to_numpy()
            oof1[te] = m1.predict_log_partial_hazard(b).to_numpy()
            for hz in horizons:
                try:
                    r0 = 1 - m0.predict_survival_function(b, times=[hz]).iloc[0].to_numpy()
                    r1 = 1 - m1.predict_survival_function(b, times=[hz]).iloc[0].to_numpy()
                except Exception:
                    continue
                w, ev, ct, gh = weights(a, b, hz)
                for (_, r), p0, p1 in zip(b.iterrows(), r0, r1):
                    pred_rows.append({'repeat': rep, 'fold': fold, 'id': r['id'], 'horizon': hz,
                                      'T': r['T'], 'E': r['E'],
                                      'risk_clinical': float(p0), 'risk_molecular': float(p1)})
        ok = ~np.isnan(oof0)
        t_ok = d.loc[ok, 'T'].to_numpy()
        e_ok = d.loc[ok, 'E'].to_numpy()
        # censoring survival at each event time, left-continuous (archived convention)
        gvals = G_at(d[ok], t_ok, left=True)
        row = {'repeat': rep, 'n_scored': int(ok.sum()),
               'uno_c_clinical': uno_c(t_ok, e_ok, oof0[ok], gvals),
               'uno_c_molecular': uno_c(t_ok, e_ok, oof1[ok], gvals)}
        for hz in horizons:
            for lab, sc in [('clinical', oof0), ('molecular', oof1)]:
                row[f'auc_{lab}_{hz:g}'] = time_dependent_auc(
                    t_ok, e_ok, sc[ok], hz, gvals)
        uno_rows.append(row)
    pr = pd.DataFrame(pred_rows)
    pr.to_csv(OUT / f'calibration_oof_subjects_{name}{out_suffix}.csv', index=False)

    summ = []
    for hz in horizons:
        sub = pr[pr.horizon == hz]
        for rep in sorted(sub['repeat'].unique())[:1]:
            s = sub[sub['repeat'] == rep].reset_index(drop=True)
            for model in ['clinical', 'molecular']:
                rk = s[f'risk_{model}'].to_numpy()
                t = s['T'].to_numpy(float)
                e = s['E'].to_numpy(float)
                ev = (t <= hz) & (e == 1)
                ct = t > hz
                w = np.zeros(len(s))
                ge = G_at(s, t[ev], True)
                gh = float(G_at(s, [hz])[0])
                if gh <= 0 or np.any(ge <= 0):
                    continue
                w[ev] = 1 / ge
                w[ct] = 1 / gh
                observed = float(np.sum(w * ev) / len(s))
                predicted = float(np.mean(rk))
                citl = observed - predicted
                pk = np.clip(rk, 1e-6, 1 - 1e-6)
                lp = np.log(pk / (1 - pk))
                X = sm.add_constant(lp)
                try:
                    fit = sm.GLM(ev.astype(float), X, family=sm.families.Binomial(),
                                 freq_weights=w).fit()
                    slope, intercept = float(fit.params[1]), float(fit.params[0])
                except Exception:
                    slope, intercept = np.nan, np.nan
                q = pd.qcut(pd.Series(rk).rank(method='first'), 5, labels=False).to_numpy()
                curve = []
                for g in range(5):
                    m = q == g
                    if m.sum() == 0:
                        continue
                    curve.append({'horizon': hz, 'model': model, 'group': g + 1, 'n': int(m.sum()),
                                  'mean_predicted': float(rk[m].mean()),
                                  'observed_ipcw': float(np.sum(w[m] * ev[m]) / max(np.sum(w[m]), 1e-9))})
                pd.DataFrame(curve).to_csv(
                    OUT / f'calibration_curve_{name}{out_suffix}_{model}_{hz:g}.csv', index=False)
                summ.append({'comparison': name, 'grouped': grouped, 'horizon': hz, 'model': model,
                             'n': len(s), 'events': int(s['E'].sum()),
                             'observed_ipcw': observed, 'mean_predicted': predicted,
                             'calibration_in_the_large': citl, 'calibration_slope': slope,
                             'calibration_intercept': intercept,
                             'observed_over_predicted': observed / predicted if predicted > 0 else np.nan})
    pd.DataFrame(summ).to_csv(OUT / f'calibration_summary_{name}{out_suffix}.csv', index=False)
    u = pd.DataFrame(uno_rows)
    u.to_csv(OUT / f'uno_c_{name}{out_suffix}.csv', index=False)
    print('DONE', name, out_suffix or '_grouped', json.dumps({
        'uno_c_clinical': float(u.uno_c_clinical.mean()),
        'uno_c_molecular': float(u.uno_c_molecular.mean()),
        'uno_delta': float((u.uno_c_molecular - u.uno_c_clinical).mean()),
        'seconds': time.time() - started}), flush=True)


if __name__ == '__main__':
    name = sys.argv[1]
    grouped = (len(sys.argv) < 3) or (sys.argv[2] != 'record')
    run(name, grouped=grouped, out_suffix='' if grouped else '_record')
