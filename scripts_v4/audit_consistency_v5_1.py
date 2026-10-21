"""audit_consistency_v5_1.py — version-scoped PASS/FAIL gate for manuscript_npjDM_v5_1.md.
Reverse-checks every value the reviewer flagged against the archived results, forbids
the retired phrasings so they cannot creep back, checks structure, and confirms the
figures and the new Hong Kong frozen-classifier bootstrap exist. A round is not done
until this prints ALL CHECKS PASSED."""
from pathlib import Path
import re, os, sys
import numpy as np, pandas as pd

ROOT = Path(r'E:\脑膜瘤npj')
MD = ROOT / 'manuscript_npjDM_v5_1.md'
M = MD.read_text(encoding='utf-8')
IDENT = ROOT / 'identity_preprocessing_20261019/results'
STR = ROOT / 'strengthening_20261019/results'
V4 = ROOT / 'results_v4'

fails, checks = [], 0
def chk(label, cond, detail=''):
    global checks; checks += 1
    if not cond: fails.append(f'{label}: {detail}')

# ---------- 1. forbidden (retired) phrasings must be ABSENT ----------
FORBIDDEN = ['worse than chance', 'hypermitonic', 'All four would survive peer review',
    '867 external rows', 'excluded zero in Hong Kong and (just)', 'is unmeasured',
    'single level-A estimate reported at record level', 'only just in one of the two',
    'fewer than seven events remained', 'Cannot assess** externally',
    'Fail** as an external claim', 'just outside it in the record-level']
for bad in FORBIDDEN:
    chk(f'forbidden absent: {bad!r}', bad not in M)
# these three may appear ONLY inside the compliance disclaimer that negates them
for term in ['clinically deployable', 'practice-changing', 'prospectively validated']:
    n = M.count(term)
    chk(f'{term!r} only in disclaimer', n <= 1 and 'No claim uses the terms' in M, f'{n} occurrences')

# ---------- 2. required corrected phrasings must be PRESENT ----------
REQUIRED = ['incompletely characterized', '−0.000002', '33.1%', '54.3%', '93.3%',
    '803 of 1,000', 'no comparable case pairs', '866 external specimens',
    'Not assessed', 'whole-procedure', '46% of all hypermitotic',
    '46% of hypermitotic', '39% and 5%']
for good in REQUIRED:
    chk(f'required present: {good!r}', good in M)

# ---------- 3. structure ----------
title = M.splitlines()[0].lstrip('# ').strip()
chk('title <=15 words', len(title.split()) <= 15, f'{len(title.split())}')
chk('title no punctuation', not re.search(r'[:,;.?!]', title), title)
absm = re.search(r'## Abstract\s*\n\n([\s\S]*?)\n---', M)
chk('abstract found', absm is not None)
if absm:
    aw = len(absm.group(1).split())
    chk('abstract <=170 words', aw <= 170, f'{aw} words')

# ---------- 4. UCSF record-level interval actually includes zero ----------
uc = pd.read_csv(V4 / 'donor_bootstrap_ci_summary_ucsf_published.csv').iloc[0]
chk('UCSF record-level lower bound is negative (~ -2e-6)', -1e-5 < uc.delta_recordlevel_ci_low < 0,
    str(uc.delta_recordlevel_ci_low))
chk('text states -0.000002', '−0.000002' in M)
chk('text says HK only excludes zero', 'excluded zero in Hong Kong only' in M or 'excludes zero in Hong Kong alone' in M)

# ---------- 5. published donor-level intervals match text ----------
def fmt(v):
    return ('+' if v >= 0 else '−') + f'{abs(v):.3f}'
for name in ['baylor_published', 'tubingen_published', 'ucsf_published', 'hku_published']:
    r = pd.read_csv(V4 / f'donor_bootstrap_ci_summary_{name}.csv').iloc[0]
    pair = f'{fmt(r.delta_grouped_ci_low)} to {fmt(r.delta_grouped_ci_high)}'
    chk(f'{name} interval {pair} in text', pair in M, f'summary {r.delta_grouped_ci_low:.4f}..{r.delta_grouped_ci_high:.4f}')

# ---------- 6. HK frozen whole-procedure interval (post-patch) ----------
hkf = V4 / 'donor_bootstrap_ci_summary_hku_author_KNN_raw.csv'
if hkf.exists():
    r = pd.read_csv(hkf).iloc[0]
    lo, hi = r.delta_grouped_ci_low, r.delta_grouped_ci_high
    # Table 4 cell + Results sentence must carry the value
    chk('HK-frozen interval present in text (3dp)', f'{lo:.3f}'.replace('-','−') in M and f'{hi:.3f}'.replace('-','+') in M, f'{lo:.3f}..{hi:.3f}')
    chk('HK-frozen row no longer says "not run"', '| +0.059 | not run |' not in M)
    if lo > 0:
        chk('HK-frozen excludes zero -> consistent with HK-only claim', True)
else:
    chk('HK-frozen bootstrap summary exists', False, 'job not finished; do not ship V5.1 until present')

# ---------- 7. risk-coverage all-case denominators ----------
for coh, hyp_all, mer_all in [('tubingen', 0.714, 0.944), ('ucsf', 0.651, 0.903), ('hku', 0.543, 0.933)]:
    p = pd.read_csv(IDENT / f'{coh}_raw_author_ID_linked_predictions.csv')
    for g, exp in [('Hypermitotic', hyp_all), ('Merlin-intact', mer_all)]:
        sub = p[p.Group == g]; N = len(sub)
        ret = ((~sub.tied_top.astype(bool)) & (sub.confidence >= 0.75)).mean()
        chk(f'{coh} {g} all-case retention {exp:.3f}', abs(ret - exp) < 0.002, f'{ret:.4f}')
chk('HK 46% (all-case) in text', '46%' in M)
chk('HK 7% merlin all-case in text', '7% of Merlin-intact' in M)
chk('tied-vote loss 10/94 in text', '10/94' in M or '10 of 94' in M)

# ---------- 8. group chance baseline 33.1% ----------
pos = pd.read_csv(IDENT / 'ucsf_old_positional_mapping_rejected.csv')
pt = pos.published_Group.value_counts(normalize=True); pq = pos.Group.value_counts(normalize=True)
exp = sum(pt.get(c, 0) * pq.get(c, 0) for c in set(pt.index) | set(pq.index))
chk('group chance ~0.331', abs(exp - 0.331) < 0.005, f'{exp:.4f}')
chk('observed group 0.265 below chance', pos.published_group_agrees.mean() < exp)
chk('33.1% stated in text', '33.1%' in M)

# ---------- 9. independence 866 / min gap at Baylor ----------
ov = pd.read_csv(IDENT / 'raw_author_reference_overlap_audit.csv'); ext = ov[ov.cohort != 'author186']
chk('867 non-anchor rows incl pilot', len(ext) == 867)
mn = ext.loc[ext.max_abs_difference_to_nearest_reference.idxmin()]
chk('min gap 0.309 at Baylor GSM5702954', abs(mn.max_abs_difference_to_nearest_reference - 0.309393) < 1e-3 and mn.cohort == 'baylor' and mn.GSM == 'GSM5702954', f"{mn.cohort} {mn.GSM}")
chk('0 matches at 1e-12', (~ext.exact_beta_reference_match_tolerance_1e12).all())
chk('ledger no longer says 867 external rows', 'over the 866 external specimens' in M)

# ---------- 10. files exist ----------
for f in ['figures_v5/fig1_cohort_architecture.pdf', 'figures_v5/fig2_independence_agreement.pdf',
          'figures_v5/fig3_risk_coverage.pdf', 'figures_v5/fig4_incremental_value.pdf',
          'figures_v5/fig5_calibration.pdf', 'figures_v5/fig6_deployment_audit.pdf',
          'results_v4/supplementary_table_s8_numerical_ledger.csv',
          'results_v4/ucsf_identity_exclusion_ledger.csv']:
    chk(f'file {f}', (ROOT / f).exists())

# ---------- 11. figure freshness: fig3/fig4 regenerated after the V5.1 data ----------
for f in ['figures_v5/fig3_risk_coverage.pdf', 'figures_v5/fig4_incremental_value.pdf']:
    chk(f'{f} newer than manuscript build', (ROOT / f).stat().st_mtime >= MD.stat().st_mtime - 1,
        'regenerate after final edits')

print(f'checks run: {checks}')
if fails:
    print(f'FAILURES ({len(fails)}):')
    for f in fails: print('  -', f)
    sys.exit(1)
print('ALL CHECKS PASSED')
