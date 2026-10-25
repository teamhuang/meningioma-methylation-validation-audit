"""NEW analysis (reviewer plan §6.3C / §11): does being REJECTED by the frozen
classifier's 0.75 confidence threshold associate with a different actual recurrence
risk? Exposure = rejected (tied vote OR confidence < 0.75) vs retained. Unit = donor
(specimens of one donor must share exposure; ambiguous donors dropped and counted).
Adjusts for WHO grade; cluster-robust by donor. Reports event counts + KM/log-rank +
Cox HR with CI. Association only, NOT causal; and rejection is entangled with the
hypermitotic class, which is stated. Writes results_v4/rejection_vs_outcome_v7.json."""
import json, os
import numpy as np, pandas as pd
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

ROOT = r'E:\脑膜瘤npj'
ID = os.path.join(ROOT, 'identity_preprocessing_20261019/results')
DA = os.path.join(ID, 'donor_aware')
V4 = os.path.join(ROOT, 'results_v4')
SEED = 20261021

COH = {
 'hku':      (os.path.join(ID,'hku_raw_author_ID_linked_predictions.csv'),
              os.path.join(DA,'input_hku_published_overlap_guard.csv')),
 'tubingen': (os.path.join(ID,'tubingen_raw_author_ID_linked_predictions.csv'),
              os.path.join(DA,'input_tubingen_published_overlap_guard.csv')),
 'ucsf':     (os.path.join(ID,'ucsf_raw_author_ID_linked_predictions.csv'),
              os.path.join(V4,'input_ucsf_published_overlap_guard.csv')),
}

def run(cohort, predp, inp):
    pr = pd.read_csv(predp)[['GSM','confidence','tied_top','Group']].copy()
    pr['tied_top'] = pr.tied_top.astype(bool)
    pr['rejected'] = pr.tied_top | (pr.confidence < 0.75)
    d = pd.read_csv(inp).merge(pr[['GSM','rejected','confidence','tied_top','Group']], on='GSM', how='inner')
    d = d.dropna(subset=['rejected'])
    # donor-level: keep donors whose specimens all share the same rejected status
    g = d.groupby('donor')['rejected'].nunique()
    unambig = g[g == 1].index
    dd = d[d.donor.isin(unambig)].drop_duplicates('donor').copy()
    dd['rejected'] = dd.rejected.astype(int)
    n_ambig = d.donor.nunique() - len(unambig)
    out = {'cohort': cohort, 'records_joined': len(d), 'donors': int(dd.donor.nunique()),
           'donors_dropped_ambiguous': int(n_ambig),
           'n_rejected': int(dd.rejected.sum()), 'n_retained': int((1-dd.rejected).sum()),
           'events_rejected': int(dd[dd.rejected==1].E.sum()), 'events_retained': int(dd[dd.rejected==0].E.sum())}
    er = dd[dd.rejected==1]; ee = dd[dd.rejected==0]
    out['event_rate_rejected'] = round(float(er.E.mean()),4) if len(er) else None
    out['event_rate_retained'] = round(float(ee.E.mean()),4) if len(ee) else None
    # simple descriptive 2x2 test on observed events (recurrence yes/no) by rejection
    if len(er) and len(ee):
        try:
            from scipy.stats import fisher_exact
            tab = [[int(er.E.sum()), int((1-er.E).sum())],
                   [int(ee.E.sum()), int((1-ee.E).sum())]]
            orr, pf = fisher_exact(tab)
            out['event_odds_ratio_rejected_vs_retained'] = round(float(orr),3)
            out['event_fisher_p'] = float(f"{pf:.3g}")
        except Exception as e:
            out['event_test_error'] = str(e)[:120]
    # grade-adjusted Cox, one record per donor
    cols = ['T','E','rejected','grade']
    cc = dd[cols].dropna().reset_index(drop=True)
    if cc.E.sum()>=10 and cc.rejected.nunique()==2:
        try:
            cph = CoxPHFitter(penalizer=0.0)
            cph.fit(cc, 'T','E', robust=True)
            s = cph.summary.loc['rejected']
            out['cox_HR_rejected_vs_retained'] = round(float(np.exp(s['coef'])),3)
            out['cox_HR_ci'] = [round(float(np.exp(s['coef']-1.96*s['se(coef)'])),3),
                                round(float(np.exp(s['coef']+1.96*s['se(coef)'])),3)]
            out['cox_p'] = float(f"{s['p']:.3g}")
            out['cox_adjusted_for'] = 'WHO grade (ordinal); one record per donor'
        except Exception as e:
            out['cox_error'] = str(e)[:120]
    else:
        out['cox_note'] = 'events too few for a stable grade-adjusted Cox; descriptive only'
    # class entanglement: rejected share by published group
    ent = d.groupby('Group')['rejected'].agg(['mean','size']).round(3)
    out['rejected_fraction_by_group'] = {k: float(v) for k,v in ent['mean'].to_dict().items()}
    return out

res = {}
for coh,(pr,inp) in COH.items():
    try:
        res[coh] = run(coh, pr, inp)
    except Exception as e:
        res[coh] = {'error': str(e)[:200]}

path = os.path.join(V4,'rejection_vs_outcome_v7.json')
tmp = path+'.tmp'
with open(tmp,'w',encoding='utf-8') as fh: json.dump(res, fh, indent=2, ensure_ascii=False); fh.flush(); os.fsync(fh.fileno())
os.replace(tmp,path)
print(json.dumps(res, indent=2, ensure_ascii=False))
