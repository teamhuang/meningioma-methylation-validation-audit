"""Build Supplementary Table S8: a numerical ledger mapping every value reported in
manuscript_npjDM_v5.md to its source file, and the per-specimen sample accounting of
Table 1. Output is a markdown table plus a machine-readable CSV.
"""
from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path(r'E:\脑膜瘤npj')
IDENT = ROOT / 'identity_preprocessing_20261019/results'
STR = ROOT / 'strengthening_20261019/results'
V4 = ROOT / 'results_v4'
OUT = ROOT / 'results_v4'

# value -> (source file, locator)
LEDGER = [
    ('866 external specimens; 870 genotyped', 'identity_preprocessing_20261019/results/raw_author_reference_overlap_audit.csv',
     '870 rows; cohort counts baylor 110, tubingen 150 (149 clinical + 1 pilot), ucsf 200, hku 365, longitudinal 42, author186 3'),
    ('minimum gap to nearest training reference = 0.309',
     'identity_preprocessing_20261019/results/raw_author_reference_overlap_audit.csv',
     'min(max_abs_difference_to_nearest_reference) over the 866 external specimens; the raw file holds 867 non-anchor rows including one pilot specimen, and the 0.309 minimum is attained at a Baylor specimen (GSM5702954), so it is unchanged whether or not the pilot is included'),
    ('0 matches at every threshold 1e-12 to 0.30',
     'identity_preprocessing_20261019/results/raw_author_reference_overlap_audit.csv',
     'exact_beta_reference_match_tolerance_1e12 == False for all external rows'),
    ('121/186 and 20/37 training-reference matches',
     'strengthening_20261019/results/official_knn_audit_summary.csv',
     'rows release186 and release37, column n_in_model_training'),
    ('Hong Kong 329 donor clusters; 23 multi-specimen clusters covering 59 specimens',
     'identity_preprocessing_20261019/results/SNP_donor_clusters_current.csv',
     'group counts of donor_overlap_guard_key for cohort hku'),
    ('Baylor 108 clusters (2 pairs); Tübingen 143 (6 pairs); UCSF 189 (11 pairs)',
     'identity_preprocessing_20261019/results/SNP_donor_clusters_current.csv',
     'group counts of donor_overlap_guard_key by cohort'),
    ('positional bridge agreement 0.265; sex 0.535; grade 0.410',
     'identity_preprocessing_20261019/results/ucsf_old_positional_mapping_rejected.csv',
     'published_group_agrees, sex_agrees, grade_agrees over 200 rows'),
    ('UCSF identity exclusion ledger (2 MIDs claim one row; 2 rows with no MID)',
     'results_v4/ucsf_identity_exclusion_ledger.csv', 'full file'),
    ('Tübingen 91.5% (86.4-95.8) and 96.0% (92.1-99.2) at 0.75',
     'identity_preprocessing_20261019/results/raw_author_supervised_label_metrics_with_CI.csv',
     'rows tubingen/all_assigned and tubingen/highconfidence_075'),
    ('UCSF 87.2% (82.2-91.7) and 93.0% (89.0-96.8) at 0.75',
     'identity_preprocessing_20261019/results/raw_author_supervised_label_metrics_with_CI.csv',
     'rows ucsf/all_assigned and ucsf/highconfidence_075'),
    ('Hong Kong 88.8% (85.1-91.9) and 93.6% (90.6-96.2) at 0.75',
     'identity_preprocessing_20261019/results/raw_author_supervised_label_metrics_with_CI.csv',
     'rows hku/all_assigned and hku/highconfidence_075'),
    ('Hong Kong identity-concordant subset 88.6% / 93.5%',
     'identity_preprocessing_20261019/results/raw_author_supervised_label_metrics_with_CI.csv',
     'rows hku/identity_concordant_*'),
    ('per-class recall and precision; confusion matrices',
     'identity_preprocessing_20261019/results/{tubingen,ucsf,hku}_raw_author_confusion_all_assigned.csv',
     'all rows; accuracy and ARI also in raw_author_supervised_label_metrics_with_CI.csv'),
    ('risk-coverage curves; class-specific retention at 0.75 and 0.95',
     'identity_preprocessing_20261019/results/{tubingen,ucsf,hku}_raw_author_ID_linked_predictions.csv',
     'confidence, tied_top, Group, semantic_prediction per specimen'),
    ('Baylor donor-grouped dC +0.092, CI -0.004 to 0.170, corrected +0.072',
     'identity_preprocessing_20261019/results/donor_aware/summary_baylor_published_overlap_guard.csv', 'full row'),
    ('Baylor surrogate +0.025, corrected +0.033',
     'identity_preprocessing_20261019/results/donor_aware/summary_baylor_centroid_surrogate_overlap_guard.csv', 'full row'),
    ('Baylor frozen KNN +0.016, corrected +0.018',
     'identity_preprocessing_20261019/results/donor_aware/summary_baylor_author_KNN_raw_overlap_guard.csv', 'full row'),
    ('Tuebingen published +0.031, corrected +0.046',
     'identity_preprocessing_20261019/results/donor_aware/summary_tubingen_published_overlap_guard.csv', 'full row'),
    ('Tuebingen surrogate +0.025, corrected +0.024',
     'identity_preprocessing_20261019/results/donor_aware/summary_tubingen_centroid_surrogate_overlap_guard.csv', 'full row'),
    ('Tuebingen frozen KNN +0.041, corrected +0.041',
     'identity_preprocessing_20261019/results/donor_aware/summary_tubingen_author_KNN_raw_overlap_guard.csv', 'full row'),
    ('UCSF donor-grouped dC +0.027, CI 0.011 to 0.076, corrected +0.030 (198 records, 187 clusters)',
     'results_v4/summary_ucsf_published_overlap_guard.csv', 'full row'),
    ('UCSF record-level on the identical 198 records: dC +0.024, 19/20 positive',
     'results_v4/summary_ucsf_published_recordlevel_same198.csv', 'full row'),
    ('Hong Kong published +0.063, CI 0.020 to 0.110, corrected +0.056',
     'identity_preprocessing_20261019/results/donor_aware/summary_hku_published_overlap_guard.csv', 'full row'),
    ('Hong Kong frozen KNN +0.059, CI 0.025 to 0.116, corrected +0.059',
     'identity_preprocessing_20261019/results/donor_aware/summary_hku_author_KNN_raw_overlap_guard.csv', 'full row'),
    ('Hong Kong frozen-KNN whole-procedure donor-level interval +0.008 to +0.132 (1,000/1,000 draws)',
     'results_v4/donor_bootstrap_ci_summary_hku_author_KNN_raw.csv', 'full row'),
    ('risk-coverage retention over all true cases (HK 54.3% hyper / 93.3% Merlin at 0.75)',
     'identity_preprocessing_20261019/results/hku_raw_author_ID_linked_predictions.csv',
     'retention = not tied_top and confidence >= 0.75, over all specimens of each published Group'),
    ('record-level dC for the four counterfactual comparisons (+0.073, +0.029, +0.023, +0.011)',
     'strengthening_20261019/results/validation_summary.csv', 'column cv_delta_mean'),
    ('sixteen-test family: raw P and Holm-adjusted P',
     'strengthening_20261019/results/unpenalized_association_tests.csv',
     'column P_Holm_all_tested for the 14 archived tests; the two Hong Kong rows added with the same Holm procedure'),
    ('Baylor HR 4.99 (1.24-20.05); Tuebingen 3.53 (1.46-8.53); UCSF 3.73 (1.79-7.79)',
     'strengthening_20261019/results/cohort_specific_unpenalized_HR.csv', 'rows baylor_published, tubingen_published, ucsf_published'),
    ('Hong Kong HR 5.29 (2.64-10.59), P 2.6e-06, cluster-robust',
     'identity_preprocessing_20261019/results/donor_aware/cluster_robust_cox_hku_published_overlap_guard.csv',
     'row covariate=published_hyper'),
    ('apparent 5-year O/P 1.105, 1.039, 0.962',
     'strengthening_20261019/results/apparent_horizon_{baylor_published,tubingen_transferred,ucsf_published}.csv',
     'column OP_molecular at the 60-month or 5-year horizon'),
    ('recalibration slopes 1.208 (0.605-1.811), 1.133 (0.503-1.762), 1.149 (0.822-1.476)',
     'strengthening_20261019/results/apparent_calibration_slopes.csv', 'rows with model=molecular'),
    ('cross-fitted slope means 0.796, 0.480, 1.021',
     'strengthening_20261019/results/oof_slope_summary.csv', 'rows with model=molecular'),
    ('five-year Brier 0.264->0.260, 0.248->0.242, 0.199->0.183, 0.248->0.250, 0.268->0.287',
     'strengthening_20261019/results/oof_risk_summary.csv', 'brier_clinical_mean and brier_molecular_mean at the horizon'),
    ('Baylor 5-year: 11 of 99 still under observation; G(60) = 0.13',
     'strengthening_20261019/results/apparent_horizon_baylor_published.csv', 'columns at_risk, G_horizon, n'),
    ('corrected robust panel: Baylor 679 CpGs, 3/110, ARI 0.808, best-map 0.910; UCSF 1,002 CpGs, 13/200, ARI 0.485, best-map 0.795; multicentre 1,999, 31/994',
     'strengthening_20261019/results/corrected_P3_transfer_summary.csv', 'rows baylor_correctedP3, ucsf_discovery200_correctedP3, multicentre_correctedP3'),
    ('erroneous panel UCSF best-map 0.860 (counterfactual only)',
     'analysis_report.md', 'third-round audit entry; superseded value retained only as the counterfactual in Table 6'),
    ('longitudinal: 17 clusters, 16 concordant, 15 all-class-3',
     'strengthening_20261019/results/longitudinal_baseline_correction.csv', 'full file'),
    ('18 source-labelled groups vs 16 SNP donor clusters',
     'identity_preprocessing_20261019/results/longitudinal_identity_resolution_status.json',
     'source_labelled_groups_recovered and SNP clustering'),
]

rows = []


def expand(src):
    """Expand a {a,b,c} pattern into concrete paths and report whether all exist."""
    if '{' in src:
        head, rest = src.split('{', 1)
        body, tail = rest.split('}', 1)
        parts = [head + p + tail for p in body.split(',')]
    else:
        parts = [src]
    ok = all((ROOT / p).exists() for p in parts)
    return parts, ('yes' if ok else 'no')


for value, src, loc in LEDGER:
    parts, present = expand(src)
    shown = src if '{' not in src else parts[0].split('/')[-1] + ' (+ siblings)'
    rows.append({'value_or_finding': value, 'source_file': shown,
                 'locator': loc, 'present': present,
                 'resolved_paths': '; '.join(parts)})
df = pd.DataFrame(rows)
df.to_csv(OUT / 'supplementary_table_s8_numerical_ledger.csv', index=False)

md = ['### Supplementary Table S8. Numerical ledger: every reported value mapped to its source file', '',
      '| Reported value or finding | Source file | Locator | File present |', '| --- | --- | --- | --- |']
for _, r in df.iterrows():
    md.append(f"| {r.value_or_finding} | `{r.source_file}` | {r.locator} | {r.present} |")
md += ['', f'*{len(df)} entries. All paths are relative to the repository root. '
          'The per-specimen sample accounting behind Table 1 is archived in '
          '`identity_preprocessing_20261019/results/raw_author_reference_overlap_audit.csv` '
          '(870 rows) together with the UCSF exclusion ledger in `results_v4/ucsf_identity_exclusion_ledger.csv`.*']
(OUT / 'supplementary_table_s8_numerical_ledger.md').write_text('\n'.join(md), encoding='utf-8')
print('\n'.join(md[:8]))
print('...')
print(f'ledger rows: {len(df)}; missing files: {(df.present == "no").sum()}')
