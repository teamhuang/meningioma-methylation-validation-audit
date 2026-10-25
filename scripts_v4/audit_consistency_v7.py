"""audit_consistency_v7.py — PASS/FAIL gate for manuscript_npjDM_v7.md.
Reverse-checks the new rejection-vs-outcome analysis, the 37/20/17/11 reconciliation,
the -39.1 pp bootstrap, title length, retired phrasings, and declaration completeness."""
from pathlib import Path
import re, json, sys
import pandas as pd
ROOT = Path(r'E:\脑膜瘤npj')
M = (ROOT/'manuscript_npjDM_v7.md').read_text(encoding='utf-8')
V6 = ROOT/'npjDM_submission_candidate_v6'
S = ROOT/'strengthening_20261019/results'
fails, checks = [], 0
def chk(l, c, d=''):
    global checks; checks += 1
    if not c: fails.append(f'{l}: {d}')

# title
title = M.splitlines()[0].lstrip('#').replace('\ufeff','').strip()
chk('title <=15 words', len(title.split()) <= 15, f'{len(title.split())}: {title}')
chk('title no punctuation', not re.search(r'[:,;.?!]', title), title)

# retired phrasings
for bad in ['genuinely non-training', '11 of its 37', 'on its 11 non-training profiles']:
    chk(f'retired absent: {bad!r}', bad not in M)
# corrected reconciliation present
chk('11-of-17 reconciliation present', '11 of its 17 non-training profiles' in M)

# reverse-check 37/20/17/11 from the archived audit
kn = pd.read_csv(S/'official_knn_audit_summary.csv').set_index('cohort').loc['release37']
chk('release37 n_scored==37', int(kn.n_scored)==37, str(kn.n_scored))
chk('release37 n_in_model_training==20', int(kn.n_in_model_training)==20, str(kn.n_in_model_training))
chk('non-training = 37-20 = 17', int(kn.n_scored)-int(kn.n_in_model_training)==17)
chk('n_nontraining_known_nontied_labels==11', int(kn.n_nontraining_known_nontied_labels)==11, str(kn.n_nontraining_known_nontied_labels))

# rejection-vs-outcome JSON reverse-check
d = json.load(open(ROOT/'results_v4/rejection_vs_outcome_v7.json', encoding='utf-8'))
for coh, hr, lo, hi, p in [('hku',1.266,0.678,2.363,0.459),('ucsf',1.043,0.563,1.933,0.894),('tubingen',0.62,0.299,1.284,0.198)]:
    r=d[coh]
    chk(f'{coh} HR matches JSON', abs(r['cox_HR_rejected_vs_retained']-hr)<1e-3, str(r['cox_HR_rejected_vs_retained']))
    chk(f'{coh} CI matches JSON', abs(r['cox_HR_ci'][0]-lo)<1e-3 and abs(r['cox_HR_ci'][1]-hi)<1e-3, str(r['cox_HR_ci']))
# the HRs as printed in text (rounded to 2dp)
for txt in ['1.27 (95% CI 0.68–2.36; p = 0.46)', '1.04 (0.56–1.93; p = 0.89)', '0.62 (0.30–1.28; p = 0.20)']:
    chk(f'text HR string present: {txt}', txt in M)
chk('rejected fractions in text', '48.2%' in M and '6.7%' in M and '34.4%' in M and '9.7%' in M)
chk('Supplementary Data 6 referenced', 'Supplementary Data 6' in M)
chk('Supp Data 6 file exists', (V6/'Supplementary_Data_6_rejection_outcome.csv').exists())

# -39.1 pp bootstrap (Supp Data 5) reverse-check
sd5 = pd.read_csv(V6/'Supplementary_Data_5_class_retention_bootstrap.csv').set_index('cohort')
chk('HK difference ~ -39.1 pp', abs(sd5.loc['hku','difference_hypermitotic_minus_merlin_intact']+0.3908)<1e-3, str(sd5.loc['hku','difference_hypermitotic_minus_merlin_intact']))
chk('text carries -39.1 and CI', '39.1 percentage-point' in M and '27.9–50.3' in M)

# declarations present
for sec in ['## Ethics approval','Data availability','Code availability','Competing interests','Author contributions','Generative AI assistance']:
    chk(f'declaration present: {sec}', sec in M)

# bundle figures present
for f in ['Figure1.pdf','Figure4.pdf','Figure6.pdf']:
    chk(f'bundle figure {f}', (V6/f).exists())

print(f'checks run: {checks}')
if fails:
    print(f'FAILURES ({len(fails)}):')
    for f in fails: print('  -', f)
    sys.exit(1)
print('ALL CHECKS PASSED')
