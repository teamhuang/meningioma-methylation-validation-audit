"""Final consistency check of manuscript_npjDM_v5.md against the archived results."""
from pathlib import Path
import re
import pandas as pd
import numpy as np

ROOT = Path(r'E:\脑膜瘤npj')
M = (ROOT / 'manuscript_npjDM_v5.md').read_text(encoding='utf-8')
IDENT = ROOT / 'identity_preprocessing_20261019/results'
STR = ROOT / 'strengthening_20261019/results'
V4 = ROOT / 'results_v4'

fails, checks = [], 0


def chk(label, cond, detail=''):
    global checks
    checks += 1
    if not cond:
        fails.append(f'{label}: {detail}')


# --- format compliance
title = M.splitlines()[0].lstrip('# ').strip()
chk('title <=15 words', len(title.split()) <= 15, f'{len(title.split())} words: {title}')
chk('title no punctuation', not re.search(r'[:,;.?!]', title), title)
abs_ = re.search(r'## Abstract\s*\n\n([\s\S]*?)\n---', M).group(1)
chk('abstract <=150 words', len(abs_.split()) <= 150, f'{len(abs_.split())} words')
chk('abstract no subheadings', not re.search(r'\*\*[A-Z][a-z]+:?\*\*', abs_))
chk('refs <=60', len(re.findall(r'(?m)^\d+\. ', M)) <= 60)

# --- key numbers vs sources
uc = pd.read_csv(V4 / 'summary_ucsf_published_overlap_guard.csv').iloc[0]
ul = pd.read_csv(V4 / 'summary_ucsf_published_recordlevel_same198.csv').iloc[0]
chk('UCSF dC +0.027 in text', '+0.027' in M, str(uc.delta_cv))
chk('UCSF corrected +0.030', '+0.030' in M, str(uc.delta_corrected))
chk('UCSF 198/76', '198/76' in M, f'{uc.records}/{uc.events}')
chk('UCSF 187 clusters', '187' in M)
chk('UCSF record-level +0.024', '+0.024' in M, str(ul.delta_cv))
chk('UCSF CI 0.011-0.076', '0.011 to 0.076' in M, f'{uc.delta_app_CI_low:.3f}-{uc.delta_app_CI_high:.3f}')

for name, f, key, val in [
    ('Baylor pub', 'summary_baylor_published_overlap_guard.csv', 'delta_cv', 0.0916656079618887),
    ('Tübingen pub', 'summary_tubingen_published_overlap_guard.csv', 'delta_cv', 0.030600862095361413),
    ('HK pub', 'summary_hku_published_overlap_guard.csv', 'delta_cv', 0.06314230341426245),
    ('HK knn', 'summary_hku_author_KNN_raw_overlap_guard.csv', 'delta_cv', 0.0589923980079469),
]:
    r = pd.read_csv(IDENT / 'donor_aware' / f).iloc[0]
    chk(f'{name} dC', abs(r[key] - val) < 1e-9)

# --- independence
ov = pd.read_csv(IDENT / 'raw_author_reference_overlap_audit.csv')
ext = ov[ov.cohort != 'author186']
chk('867 external rows', len(ext) == 867, str(len(ext)))
chk('min gap 0.309 in text', '0.309' in M)
chk('min gap correct', abs(ext.max_abs_difference_to_nearest_reference.min() - 0.309393) < 1e-3,
    str(ext.max_abs_difference_to_nearest_reference.min()))
chk('866 external in text', '866 external' in M)
chk('no overlap at 1e-2', (ext.max_abs_difference_to_nearest_reference < 1e-2).sum() == 0)

# --- HK donor clusters
snp = pd.read_csv(IDENT / 'SNP_donor_clusters_current.csv')
hk = snp[snp.cohort == 'hku']
sizes = hk.groupby('donor_overlap_guard_key').size()
mult = sizes[sizes > 1]
chk('HK 329 clusters', hk.donor_overlap_guard_key.nunique() == 329, str(hk.donor_overlap_guard_key.nunique()))
chk('HK 23 multi clusters', len(mult) == 23, str(len(mult)))
chk('HK 59 specimens', mult.sum() == 59, str(mult.sum()))
chk('HK 329 in text', '329 donor clusters' in M)
chk('HK 59 in text', '59 specimens' in M)

# --- positional bridge
pos = pd.read_csv(IDENT / 'ucsf_old_positional_mapping_rejected.csv')
chk('positional 0.265', abs((pos.published_group_agrees == True).mean() - 0.265) < 0.001)
chk('sex 0.535', abs((pos.sex_agrees == True).mean() - 0.535) < 0.001)
chk('grade 0.410', abs((pos.grade_agrees == True).mean() - 0.410) < 0.001)
chk('grade 41.0% in text', '41.0%' in M)

# --- agreement
ag = pd.read_csv(IDENT / 'raw_author_supervised_label_metrics_with_CI.csv')
for coh, acc, hi in [('tubingen', 0.915, 0.958), ('ucsf', 0.872, 0.917), ('hku', 0.888, 0.919)]:
    r = ag[(ag.cohort == coh) & (ag.subset == 'all_assigned')].iloc[0]
    chk(f'{coh} accuracy', abs(r.accuracy_fixed_training_semantics - acc) < 0.001, str(r.accuracy_fixed_training_semantics))
    chk(f'{coh} CI high', abs(r.accuracy_CI_high - hi) < 0.001, str(r.accuracy_CI_high))

# --- risk coverage retention
for coh, hyp, mer in [('tubingen', 0.781, 0.985), ('ucsf', 0.672, 0.915), ('hku', 0.607, 0.949)]:
    p = pd.read_csv(IDENT / f'{coh}_raw_author_ID_linked_predictions.csv')
    u = p[~p.tied_top.astype(bool)]
    h = (u[u.Group == 'Hypermitotic'].confidence >= 0.75).mean()
    mm = (u[u.Group == 'Merlin-intact'].confidence >= 0.75).mean()
    chk(f'{coh} hyper retention', abs(h - hyp) < 0.002, str(h))
    chk(f'{coh} merlin retention', abs(mm - mer) < 0.002, str(mm))
chk('HK 39% in text', '39%' in M)
chk('HK 5% in text', '5% of Merlin-intact' in M)

# --- association family
at = pd.read_csv(STR / 'unpenalized_association_tests.csv')
chk('16-test family stated', 'sixteen' in M.lower())
chk('six survive stated', 'Six of sixteen' in M)
chk('4 archived Holm(14) <0.05', (at.P_Holm_all_tested < 0.05).sum() == 4,
    str((at.P_Holm_all_tested < 0.05).sum()))
chk('4 archived + 2 HK = 6 stated', 'six' in M.lower())

# --- no stale claims (explanatory mentions of the corrected values are allowed)
for bad, why in [('867 external specimens are', '867 replaced by 866'),
                 ('334 donor clusters', 'HK clusters corrected to 329'),
                 ('across four independent institutions', 'only three cohorts have labels'),
                 ('The UCSF row is reported at record level', 'UCSF now donor-grouped'),
                 ('UCSF, published groups (**record level**)', 'table updated')]:
    chk(f'no stale: {bad}', bad not in M, why)
chk('867 appears only as an explicit correction',
    M.count('867') == 2 and 'the 867-row figure previously quoted included the pilot specimen' in M,
    f'{M.count("867")} occurrences')

# --- files exist
for f in ['figures_v5/fig1_cohort_architecture.pdf', 'figures_v5/fig2_independence_agreement.pdf',
          'figures_v5/fig3_risk_coverage.pdf', 'figures_v5/fig4_incremental_value.pdf',
          'figures_v5/fig5_calibration.pdf', 'figures_v5/fig6_deployment_audit.pdf',
          'results_v4/supplementary_table_s8_numerical_ledger.csv',
          'results_v4/summary_ucsf_published_overlap_guard.csv',
          'results_v4/ucsf_identity_exclusion_ledger.csv']:
    chk(f'file {f}', (ROOT / f).exists())

print(f'checks run: {checks}')
if fails:
    print(f'FAILURES ({len(fails)}):')
    for f in fails:
        print('  -', f)
else:
    print('ALL CHECKS PASSED')
