import os
p = r'E:\脑膜瘤npj\manuscript_npjDM_v5_1.md'
s = open(p, encoding='utf-8').read()

reps = [
 ("| Hong Kong, frozen classifier (raw IDAT) | 332/77 | 308 | 0.695 → 0.754 | +0.059 | not run | 0.025 to 0.116 | +0.059 |",
  "| Hong Kong, frozen classifier (raw IDAT) | 332/77 | 308 | 0.695 → 0.754 | +0.059 | +0.008 to +0.132 (1,000/1,000) | 0.025 to 0.116 | +0.059 |"),
 ("whose lower bound is −0.000002 and therefore includes zero. The Baylor resampling succeeded",
  "whose lower bound is −0.000002 and therefore includes zero. The whole-procedure interval for the Hong Kong frozen-classifier route (+0.059) likewise excluded zero (+0.008 to +0.132; 1,000/1,000 draws), so both Hong Kong analyses — the published labels and the deployable frozen classifier — are distinguishable from zero. The Baylor resampling succeeded"),
 ("| Hong Kong frozen KNN +0.059, CI 0.025 to 0.116, corrected +0.059 | `identity_preprocessing_20261019/results/donor_aware/summary_hku_author_KNN_raw_overlap_guard.csv` | full row | yes |",
  "| Hong Kong frozen KNN +0.059, CI 0.025 to 0.116, corrected +0.059 | `identity_preprocessing_20261019/results/donor_aware/summary_hku_author_KNN_raw_overlap_guard.csv` | full row | yes |\n| Hong Kong frozen-KNN whole-procedure donor-level interval +0.008 to +0.132 (1,000/1,000 draws) | `results_v4/donor_bootstrap_ci_summary_hku_author_KNN_raw.csv` | full row | yes |"),
]
for old, new in reps:
    c = s.count(old)
    assert c == 1, f'count={c} for {old[:60]!r}'
    s = s.replace(old, new)

tmp = p + '.tmp'
with open(tmp, 'w', encoding='utf-8') as fh:
    fh.write(s); fh.flush(); os.fsync(fh.fileno())
os.replace(tmp, p)
print('patched OK; not-run gone:', '| +0.059 | not run |' not in s,
      '| HK-frozen interval in text:', '+0.008 to +0.132' in s,
      '| ledger row added:', 'donor_bootstrap_ci_summary_hku_author_KNN_raw.csv' in s)
