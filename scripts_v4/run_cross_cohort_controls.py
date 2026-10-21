"""Controls for the cross-cohort external validation.

A C-index below 0.5 in an external cohort is a strong claim and has two very
different explanations: (a) the model genuinely does not transport, or (b) the
harmonisation itself is wrong - for example the covariates are coded in opposite
directions in different source tables. This script distinguishes them.

Controls run:
  1. apparent discrimination of the development model in the development cohort
     (self-consistency; must be in the usual range)
  2. the same frozen transfer with the MOLECULAR predictor removed, i.e. a clinical
     -only model (if this also collapses, the collapse is a harmonisation artefact)
  3. the signed coefficient of every shared covariate in each cohort
  4. the event rate and follow-up distribution in each cohort
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

BASE = Path(r'E:\脑膜瘤npj')
sys.path.insert(0, str(BASE / 'scripts_v4'))
import run_cross_cohort_external_validation as xv

rows = []
for name, d in xv.COHORTS.items():
    cols = xv.CLIN + xv.MOL
    sub = d[['T_months', 'E'] + cols].dropna()
    m = CoxPHFitter(penalizer=0.05).fit(sub, 'T_months', 'E')
    lp = m.predict_log_partial_hazard(sub).to_numpy()
    app = xv.concord(sub['T_months'].to_numpy(), sub['E'].to_numpy(), lp)
    # clinical-only model
    csub = d[['T_months', 'E'] + xv.CLIN].dropna()
    mc = CoxPHFitter(penalizer=0.05).fit(csub, 'T_months', 'E')
    lpc = mc.predict_log_partial_hazard(csub).to_numpy()
    appc = xv.concord(csub['T_months'].to_numpy(), csub['E'].to_numpy(), lpc)
    rows.append({
        'cohort': name, 'n': len(sub), 'events': int(sub['E'].sum()),
        'event_rate': round(float(sub['E'].mean()), 3),
        'followup_median_months': round(float(sub['T_months'].median()), 1),
        'followup_max_months': round(float(sub['T_months'].max()), 1),
        'apparent_C_molecular': round(app, 3),
        'apparent_C_clinical_only': round(appc, 3),
        'events_by_grade': ';'.join(
            f"g{int(g)}={int(v)}" for g, v in sub.groupby('grade')['E'].agg(['sum', 'count']).values),
        'gtr_mean': round(float(sub['gtr'].mean()), 3),
        'gtr_event_rate': round(float(sub.loc[sub['gtr'] == 1, 'E'].mean()), 3) if (sub['gtr'] == 1).any() else np.nan,
        'gtr0_event_rate': round(float(sub.loc[sub['gtr'] == 0, 'E'].mean()), 3) if (sub['gtr'] == 0).any() else np.nan,
        'male_mean': round(float(sub['male'].mean()), 3),
        'male_event_rate': round(float(sub.loc[sub['male'] == 1, 'E'].mean()), 3) if (sub['male'] == 1).any() else np.nan,
        'female_event_rate': round(float(sub.loc[sub['male'] == 0, 'E'].mean()), 3) if (sub['male'] == 0).any() else np.nan,
    })
    for c in cols:
        rows[-1][f'coef_{c}'] = round(float(m.params_.get(c, np.nan)), 3)
    for c in xv.CLIN:
        rows[-1][f'clin_coef_{c}'] = round(float(mc.params_.get(c, np.nan)), 3)

df = pd.DataFrame(rows)
pd.set_option('display.width', 250)
print(df[['cohort', 'n', 'events', 'event_rate', 'followup_median_months',
          'apparent_C_molecular', 'apparent_C_clinical_only']].to_string(index=False))
print()
print(df[['cohort', 'gtr_mean', 'gtr_event_rate', 'gtr0_event_rate',
          'male_mean', 'male_event_rate', 'female_event_rate',
          'events_by_grade']].to_string(index=False))
print()
print(df[['cohort', 'coef_grade', 'coef_gtr', 'coef_male', 'coef_author_2', 'coef_author_3',
          'clin_coef_grade', 'clin_coef_gtr', 'clin_coef_male']].to_string(index=False))
df.to_csv(BASE / 'results_v4/cross_cohort_harmonisation_controls.csv', index=False)

# frozen clinical-only transfer, for comparison with the molecular model
print('\nFROZEN CLINICAL-ONLY TRANSFER (no molecular predictor)')
out = []
for dev_name, dev in xv.COHORTS.items():
    frozen = xv.fit_frozen(dev)
    frozen_clin = {k: v for k, v in frozen.items()}
    d = dev[['T_months', 'E'] + xv.CLIN].dropna()
    m = CoxPHFitter(penalizer=0.05).fit(d, 'T_months', 'E')
    sf = m.baseline_survival_.iloc[:, 0]
    times = sf.index.to_numpy(); surv = sf.to_numpy()
    keep = times <= 60.0
    frozen_clin['coef'] = m.params_.to_dict()
    frozen_clin['active'] = xv.CLIN
    frozen_clin['S0_60'] = float(np.interp(60.0, times[keep], surv[keep]))
    for test_name, test in xv.COHORTS.items():
        if test_name == dev_name:
            continue
        need = [c for c in xv.CLIN if c in test.columns]
        r = xv.evaluate(frozen_clin, test.dropna(subset=need + ['T_months', 'E']).copy(),
                        f'{dev_name} -> {test_name}')
        out.append({'development': dev_name, 'test': test_name, 'c_index': round(r['c_index'], 3),
                    'observed_over_predicted': round(r['observed_over_predicted'], 3)})
print(pd.DataFrame(out).to_string(index=False))
pd.DataFrame(out).to_csv(BASE / 'results_v4/cross_cohort_clinical_only_transfer.csv', index=False)
