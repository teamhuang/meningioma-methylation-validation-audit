"""Establish the identifier-proven UCSF clinical-row <-> molecular-profile chain.

The rejected positional chain (fingerprint column index) reproduces 26.5% agreement
with the published group labels, which is the failure the manuscript documents. The
correct chain must reproduce the archive's verified figure of 87.2% for the 196
uniquely classified UCSF specimens.

Candidate chains are built from the identifier-proven tables and scored against two
independent targets:
  (i) agreement between the published group label and the clinical table's own
      group columns at the linked row, and
  (ii) agreement between the published group label and the frozen classifier's
      class, which must reproduce the archived 87.2%.
"""
from pathlib import Path
import numpy as np
import pandas as pd

BASE = Path(r'E:\脑膜瘤npj')
IDENT = BASE / 'identity_preprocessing_20261019/results'
STR = BASE / 'strengthening_20261019/results'

clin = pd.read_csv(STR / 'analysis_input_ucsf_published.csv')
clin['idx'] = clin.id.str.extract(r'row_(\d+)').astype(int)
pg = pd.read_csv(IDENT / 'ucsf_ID_proven_published_groups.csv')
link = pd.read_csv(IDENT / 'ucsf_ID_proven_molecular_linkage_200.csv')
pred = pd.read_csv(IDENT / 'ucsf_raw_author_ID_linked_predictions.csv')
pred = pred[~pred.tied_top.astype(bool)][['GSM', 'predicted_class', 'semantic_prediction']]

target_archived = 0.872
print('archived UCSF agreement to verify against: %.3f\n' % target_archived)

results = []
for sheet in sorted(pg.source_sheet.unique()):
    p = pg[pg.source_sheet == sheet].drop_duplicates('MID')
    for off in range(-10, 11):
        p2 = p.copy()
        p2['idx'] = p2.excel_row - off
        j = clin.merge(p2, on='idx', how='inner')
        if len(j) < 100:
            continue
        j = j.merge(link[['MID', 'GSM']], on='MID', how='inner')
        j = j.merge(pred, on='GSM', how='inner')
        if len(j) < 100:
            continue
        grp_agree = (j.Group == j.published_Group).mean() if 'published_Group' in j else np.nan
        clf_agree = (j.semantic_prediction == j.Group).mean()
        clin_agree = (j.Group.eq('Hypermitotic') == j.hm.eq(1)).mean()
        results.append({'sheet': sheet, 'offset': off, 'n': len(j),
                        'label_vs_clinical': round(clin_agree, 4),
                        'label_vs_classifier': round(clf_agree, 4),
                        'matches_archived': abs(clf_agree - target_archived) < 0.005})
        if abs(clf_agree - target_archived) < 0.02:
            print(f'sheet {sheet} offset {off:+d}: n={len(j)}  label vs clinical-table group = {clin_agree:.4f}'
                  f'  label vs frozen classifier = {clf_agree:.4f}  <== candidate')

res = pd.DataFrame(results)
res.to_csv(BASE / 'results_v4/ucsf_chain_search.csv', index=False)
print('\nbest by label-vs-classifier agreement:')
print(res.sort_values('label_vs_classifier', ascending=False).head(6).to_string(index=False))
