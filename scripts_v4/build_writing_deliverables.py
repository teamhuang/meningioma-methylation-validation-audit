"""Generate the remaining writing deliverables for manuscript_npjDM_v5.md:

  * Supplementary Table S7 - TRIPOD+AI / TRIPOD / REMARK / journal-standards checklist
  * Supplementary Note 1   - the analysis revision record
  * Supplementary Table S9 - the UCSF identity exclusion ledger
  * a Zenodo-ready deposition package (manifest + SHA256SUMS + README)

The checklist follows the published item lists; items that do not apply to a
validation-and-audit study are marked 'not applicable' with the reason given, as
the TRIPOD guidance requires.
"""
from pathlib import Path
import hashlib
import shutil
import datetime

ROOT = Path(r'E:\脑膜瘤npj')
OUT = ROOT / 'results_v4'
DEP = ROOT / 'deposition'
(DEP / 'results_v4').mkdir(parents=True, exist_ok=True)

# ----------------------------------------------------------------- S7 checklist
S7 = [
    # (section, item, status, where reported)
    ('Title and abstract', 'Title identifies the study as a validation/transportability audit',
     'compliant', 'Title page'),
    ('Title and abstract', 'Title is 15 words or fewer and free of punctuation',
     'compliant', '11 words, no punctuation'),
    ('Title and abstract', 'Abstract is 150 words or fewer, unstructured, with no subheadings',
     'compliant', '150 words, unstructured'),
    ('Title and abstract', 'Abstract states the objective, design, setting, participants, main results and limitation',
     'compliant', 'Abstract'),
    ('Title and abstract', 'Abstract distinguishes cohort-internal validation from external validation',
     'compliant', 'Abstract; Results; Discussion'),
    ('Introduction', 'Scientific background and rationale',
     'compliant', 'Introduction'),
    ('Introduction', 'Specific objectives and prespecified hypotheses',
     'compliant', 'Introduction, final paragraph; Objectives were fixed at the start of the revised analysis'),
    ('Methods', 'Study design stated early',
     'compliant', 'Methods: Study design, objectives and validation hierarchy'),
    ('Methods', 'Setting: eligibility criteria and sources of data',
     'compliant', 'Methods: Data sources; Table 1'),
    ('Methods', 'Participants: eligibility criteria and selection methods',
     'compliant', 'Methods: Data sources, cohort assembly; Table 1'),
    ('Methods', 'Variables: outcomes, predictors, confounders and their definitions',
     'compliant', 'Methods: Model definitions; Table 2; Table 8'),
    ('Methods', 'Data sources and measurement methods for each variable',
     'compliant', 'Methods: Identity, independence and preprocessing audits'),
    ('Methods', 'Bias: efforts to address potential sources of bias',
     'compliant', 'Results: identity, independence, probe-coverage and donor audits'),
    ('Methods', 'Study size: how the study size was arrived at',
     'compliant', 'Methods: all analyses are within-cohort, no a priori sample-size calculation; events per parameter reported'),
    ('Methods', 'Quantitative variables: how they were handled',
     'compliant', 'Methods: Model definitions and coding'),
    ('Methods', 'Statistical methods, including handling of missing data',
     'compliant', 'Methods: Validation design; Table 8; denominators archived'),
    ('Methods', 'Handling of censored or incomplete data',
     'compliant', 'Methods: Calibration; identity exclusion ledger'),
    ('Methods', 'Model building: model type, penalties, variable selection',
     'compliant', 'Methods: Model definitions (ridge penalty 0.05 fixed in advance; no data-driven selection)'),
    ('Methods', 'Model evaluation: methods and validation procedures',
     'compliant', 'Methods: Validation design; Table 4'),
    ('Methods', 'Model updating or recalibration', 'not applicable',
     'No model update or recalibration was performed'),
    ('Methods', 'Risk groups and thresholds', 'not applicable',
     'No risk score was converted into risk categories or thresholds'),
    ('Methods', 'Sensitivity analyses', 'compliant', 'Tübingen four-category, UCSF beyond-MIB-1, UCSF record-level on identical records, threshold-sensitivity of the independence audit'),
    ('Methods', 'Missing data: how missing values were handled',
     'compliant', 'Complete-case analysis; three UCSF records excluded on identity grounds with an archived ledger; 9/4 records dropped in Tübingen/UCSF'),
    ('Methods', 'Analyses for subgroups and interactions', 'not applicable',
     'No interaction or subgroup claim is made'),
    ('Methods', 'Sensitivity of analyses to modelling choices',
     'compliant', 'Four estimands reported separately; record-level and donor-grouped UCSF on identical records; threshold sensitivity of the fingerprint audit'),
    ('Results', 'Flow of participants through each stage',
     'compliant', 'Results: cohort assembly; Table 1; Figure 1a'),
    ('Results', 'Descriptive characteristics of participants',
     'compliant', 'Results: cohort assembly; Table 1; Table 8'),
    ('Results', 'Missing data at each stage', 'compliant', 'Results; archived denominators and exclusion ledger'),
    ('Results', 'Outcome events', 'compliant', 'Results; Tables 1 and 4'),
    ('Results', 'Main results, with unadjusted and adjusted estimates and precision',
     'compliant', 'Results: incremental prognostic performance; Table 4; Figures 3-4'),
    ('Results', 'Other analyses, such as subgroups and sensitivity analyses',
     'compliant', 'Results; Table 4; Supplementary Table S1'),
    ('Discussion', 'Summary of key results with reference to the objectives',
     'compliant', 'Discussion, paragraphs 1-2'),
    ('Discussion', 'Limitations, with direction and magnitude of possible bias',
     'compliant', 'Discussion, paragraph 4'),
    ('Discussion', 'Interpretation consistent with the results',
     'compliant', 'Discussion, paragraph 2'),
    ('Discussion', 'Generalisability', 'compliant', 'Discussion; deployment-readiness framework'),
    ('Discussion', 'Implications for practice, policy and future research',
     'compliant', 'Discussion, paragraph 4'),
    ('Other information', 'Sources of funding', 'to be completed by the author group', 'Author group'),
    ('Other information', 'Role of the funder', 'to be completed by the author group', 'Author group'),
    ('Ethics', 'Ethics approval and consent to participate', 'compliant',
     'Methods: all data are publicly available de-identified secondary data; no new human-participant data; institutional review board approval not required for this secondary analysis'),
    ('Ethics', 'Consent for use of identifiable information', 'compliant',
     'No identifiable information was used'),
    ('Other information', 'Data availability', 'compliant', 'Data availability statement'),
    ('Other information', 'Code availability', 'compliant', 'Code availability statement; versioned archive'),
    ('REMARK', 'Cohort or case characteristics, eligibility, and ascertainment',
     'compliant', 'Table 1; Table 8'),
    ('REMARK', 'Outcome definitions and how they were ascertained', 'compliant', 'Table 8'),
    ('REMARK', 'Treatment and covariate distributions', 'compliant', 'Table 8; archived covariate tables'),
    ('REMARK', 'Risk factor or marker measurement methods', 'compliant', 'Methods: Identity, independence and preprocessing audits'),
    ('REMARK', 'Statistical analysis, including methods for handling missing data',
     'compliant', 'Methods: Validation design and estimands'),
    ('REMARK', 'Patient flow numbers at each stage', 'compliant', 'Table 1; Figure 1a'),
    ('REMARK', 'Numbers of patients at risk, or numbers of events', 'compliant', 'Tables 1 and 4'),
    ('REMARK', 'Estimate of precision with which the association is measured',
     'compliant', 'Table 4 (donor-level bootstrap interval); Table 5'),
    ('REMARK', 'Appropriate caveats about the study and its results', 'compliant', 'Discussion; Supplementary Note 1'),
    ('Journal standard', 'Digital-biomarker independence requirement (zero overlap between training and test sets)',
     'compliant', 'Results: training independence; Figure 2a; Supplementary Figure S2'),
    ('Journal standard', 'Explicit statement of what the study does NOT establish',
     'compliant', 'Abstract; Methods (validation hierarchy); Discussion'),
]

md = ['### Supplementary Table S7. Reporting checklist (TRIPOD+AI, TRIPOD, REMARK and field standards)', '',
      '| Section | Item | Status | Where reported |', '| --- | --- | --- | --- |']
for sec, item, status, where in S7:
    md.append(f'| {sec} | {item} | {status} | {where} |')
md += ['', f'*{len(S7)} items. Items marked "not applicable" are those that presuppose model '
          'development, risk-group derivation or interaction analysis, none of which this '
          'validation-and-audit study performs. Items marked "to be completed by the author group" '
          'require information that only the authors hold.*']
(OUT / 'supplementary_table_s7_checklist.md').write_text('\n'.join(md), encoding='utf-8')

# ----------------------------------------------------------------- Note 1
NOTE1 = """### Supplementary Note 1. Analysis revision record

Every result superseded during the audit, the reason for its revision, and the corrected
replacement. Entries are in chronological order. The record is provided so that the final
analysis can be read against the path that produced it, not to narrate the audit as a
sequence of errors.

| # | Superseded result | Reason for revision | Corrected replacement | Source of the correction |
| --- | --- | --- | --- | --- |
| 1 | The 37-profile release treated as an independent validation set of the 403-CpG classifier | 20 of its 37 profiles are exact 403-probe matches to the 141-case training reference, so agreement on it measures memory rather than transportability | The release is excluded from all validation claims; the finding is reported as a training-independence failure | `official_knn_audit_summary.csv` |
| 2 | UCSF clinical rows joined to methylation profiles by row order | The positional bridge agrees with directly published per-case labels in 26.5% of cases, with sex concordance 53.5% and grade concordance 41.0%, i.e. worse than chance for grade | The positional bridge is withdrawn; all UCSF analyses use identifier-proven linkage | `ucsf_old_positional_mapping_rejected.csv` |
| 3 | UCSF label agreement of 0.38-0.40 | That number measured the broken positional bridge, not the classifier | Recomputed by identifier-proven linkage: 87.2% (82.2-91.7) of 196 assigned specimens | `raw_author_supervised_label_metrics_with_CI.csv` |
| 4 | Observed/predicted ratios of 2.19 and 1.33 | Incorrect censoring weights and normalisation | Corrected weights give apparent five-year ratios of 1.105, 1.039 and 0.962 | `apparent_horizon_*.csv` |
| 5 | Exponentiated coefficients of standardised predictors reported as calibration slopes | Those coefficients are hazard ratios per standard deviation, not calibration slopes | Slopes on the original linear-predictor scale: 1.208 (0.605-1.811), 1.133 (0.503-1.762), 1.149 (0.822-1.476) | `apparent_calibration_slopes.csv` |
| 6 | Robust-panel agreement of 0.860 for UCSF | The selection rule summed a descending MAD rank with an ascending missingness rank and then selected the highest scores, retaining exactly the panel the rule was designed to discard | Corrected ascending selection gives a disjoint 2,000-probe panel; UCSF best-map agreement 0.795 (13/200 flagged); Baylor 0.910 (3/110) | `corrected_P3_transfer_summary.csv` |
| 7 | 15/17 "majority-class baseline" for longitudinal class stability | The comparator was mis-specified; a constant class-3 assignment produces 17/17 concordant clusters, and 15 of the 17 clusters were entirely class 3 | The concordance is reported as carrying no information beyond a class-imbalance baseline | `longitudinal_baseline_correction.csv` |
| 8 | A pooled hazard ratio across cohorts | The four contrasts differ in label system, reference class, endpoint and covariates, so a pooled quantity has no clinical referent | Four cohort-specific estimates reported as direction-concordant effects, without pooling | `cohort_specific_unpenalized_HR.csv` |
| 9 | 867 external specimens | The fingerprint table carries one pilot specimen (GSM9143071) under the Tübingen label that has no clinical record and no classifier assignment | 866 external specimens; the 870 genotyped specimens are itemised in Table 1 | `raw_author_reference_overlap_audit.csv` |
| 10 | Hong Kong 334 donor clusters covering 52 specimens in 21 multi-specimen clusters | Computed from an earlier manual tally; the genotype table gives 329 clusters, 23 multi-specimen clusters and 59 specimens | 329 donor clusters, 23 multi-specimen clusters, 59 specimens | `SNP_donor_clusters_current.csv` |
| 11 | UCSF published-group comparison reported at record level, on the grounds that its clinical models contain no methylation-derived predictor | The reasoning was incorrect: the native published group labels are themselves molecular variables, and the donor-cluster mapping exists (200 identifiers to 189 clusters) | Re-run under donor-grouped folds on an identifier-proven set of 198 records in 187 clusters: ΔC +0.027, 20/20 splits positive, optimism-corrected +0.030. The record-level estimate on the identical records is +0.024 | `results_v4/summary_ucsf_published_overlap_guard.csv` |
| 12 | Two UCSF clinical records with no molecular identifier, and one record claimed by two molecular identifiers | Found only by constructing the identifier chain rather than assuming it | The three records are excluded from the UCSF prognostic analysis under an archived exclusion ledger | `results_v4/ucsf_identity_exclusion_ledger.csv` |
| 13 | An apparent bootstrap interval presented as the interval for the primary increment | The apparent interval describes the sampling behaviour of the in-sample estimator, not the cross-validated primary estimand | A donor-level bootstrap in which the entire validation procedure is re-executed inside each of 1,000 resamples of donor clusters: the interval excludes zero in Hong Kong only (+0.007 to +0.126) and spans zero in Baylor, Tübingen and the donor-grouped UCSF estimate | `results_v4/donor_bootstrap_ci_summary_*.csv` |
| 14 | The fold-wise concordance as the sole discrimination statistic | The increment might have depended on the fold-wise definition | Re-estimated with Uno's C and cumulative time-dependent AUC, which reproduce the same increments; the local implementation was validated against the archived statistic to 1e-16 | `results_v4/uno_c_*.csv` |
| 15 | "Three of sixteen association tests survive Holm correction" | The count was wrong: four of the fourteen archived tests already fall below 0.05 after correction across the fourteen-test family, and the two Hong Kong tests bring the sixteen-test family to six | Six of sixteen, listed individually in Table 5, concentrated in UCSF and Hong Kong; no Baylor or Tübingen comparison survives | `unpenalized_association_tests.csv` |
| 16 | The Hong Kong hazard ratio reported as 5.29 (2.67-10.49) with P = 2e-06 | The interval as stated is not traceable to an archived file | 5.29 (2.64-10.59), P = 2.6e-06, cluster-robust standard errors | `cluster_robust_cox_hku_published_overlap_guard.csv` |

*All corrections were made within the audit workflow with fixed seeds; the scripts that
produced each corrected value are listed in the Code availability statement, and every
value in the main text is mapped to its source file in Supplementary Table S8.*
"""
(OUT / 'supplementary_note_1_revision_record.md').write_text(NOTE1, encoding='utf-8')

# ----------------------------------------------------------------- deposition
pkg = DEP / 'npjDM_meningioma_methylation_audit_v5'
if pkg.exists():
    shutil.rmtree(pkg)
pkg.mkdir(parents=True)
copies = [
    ('manuscript_npjDM_v5.md', 'manuscript_npjDM_v5.md'),
    ('strengthening_20261019', 'strengthening_20261019'),
    ('identity_preprocessing_20261019', 'identity_preprocessing_20261019'),
    ('scripts_v4', 'scripts_v4'),
    ('figures_v5', 'figures_v5'),
]
# Raw IDAT files (37 GB) are public at GEO and are NOT deposited; the manifest that
# indexes them is. The R library is an environment artefact and is not deposited.
EXCLUDE_DIRS = {'raw_idat', 'R_library', '__pycache__', 'official_model',
                'snapshot_reproduction_check', '.git'}
for src, dst in copies:
    s = ROOT / src
    d = pkg / dst
    if s.is_dir():
        shutil.copytree(s, d, ignore=shutil.ignore_patterns(
            '__pycache__', '*.pyc', '*.npy', '*.npz', '*.pyo'))
        # prune the excluded subdirectories
        for sub in list(d.rglob('*')):
            if sub.is_dir() and sub.name in EXCLUDE_DIRS:
                shutil.rmtree(sub, ignore_errors=True)
    else:
        shutil.copy2(s, d)
for f in ['supplementary_table_s7_checklist.md', 'supplementary_note_1_revision_record.md',
          'supplementary_table_s8_numerical_ledger.md', 'supplementary_table_s8_numerical_ledger.csv',
          'ucsf_identity_exclusion_ledger.csv']:
    shutil.copy2(OUT / f, pkg / f)

# checksums
lines = []
for p in sorted(pkg.rglob('*')):
    if p.is_file() and p.name != 'SHA256SUMS.txt':
        h = hashlib.sha256(p.read_bytes()).hexdigest()
        lines.append(f'{h}  {p.relative_to(pkg).as_posix()}')
(pkg / 'SHA256SUMS.txt').write_text('\n'.join(lines) + '\n', encoding='utf-8')

readme = f"""# Deposition package

Manuscript: *Validation integrity and transportability of DNA methylation classifiers in
meningioma* (version 5).

Deposited on: {datetime.date.today().isoformat()}

## Contents

| Path | Contents |
| --- | --- |
| `manuscript_npjDM_v5.md` | the manuscript, version 5 |
| `strengthening_20261019/` | the version-5 audit package: scripts, per-sample fold assignments, out-of-fold predictions, bootstrap draws, hash manifests, numerical ledger, and the retrieved classifier artefacts with verified image-layer hashes |
| `identity_preprocessing_20261019/` | the identity, preprocessing and donor-clustering package: raw-IDAT manifest, SNP genotype tables, per-specimen donor clusters, identifier-proven linkage tables and the classifier prediction tables |
| `scripts_v4/` | the version-5 analyses added for this revision (see below) |
| `figures_v5/` | all main and supplementary figures as PDF and PNG |
| `supplementary_table_s7_checklist.md` | reporting checklist (TRIPOD+AI, TRIPOD, REMARK, field standards) |
| `supplementary_note_1_revision_record.md` | the analysis revision record |
| `supplementary_table_s8_numerical_ledger.md` / `.csv` | every reported value mapped to its source file |
| `ucsf_identity_exclusion_ledger.csv` | the three UCSF records excluded on identity grounds, with reasons |
| `SHA256SUMS.txt` | SHA256 checksum of every other file in the package |

### Deliberately not deposited

* **Raw IDAT files (about 37 GB).** These are public at GEO under the accessions listed in
  the manuscript. The manifest that indexes every file, with its URL and SHA256, is included
  as `identity_preprocessing_20261019/results/raw_idat_manifest.csv` and
  `raw_internal_barcode_and_hash_audit.csv`, so the package can be reassembled and verified
  file by file.
* **The R library** used for the classifier image, an environment artefact rather than an
  analysis output.
* **Large derived methylation matrices** (`.npy`/`.npz` intermediates), which are
  regenerable by the deposited scripts from the public raw data.

## Version-5 scripts

| Script | Purpose |
| --- | --- |
| `run_ucsf_donor_grouped_cv.py` | UCSF published-group comparison under donor-grouped folds on an identifier-proven analysis set |
| `run_ucsf_recordlevel_same_records.py` | record-level counterpart on the identical 198 records (donor-grouping counterfactual) |
| `run_donor_bootstrap_ci.py` | donor-level bootstrap in which the entire validation procedure is re-executed inside each resample; the interval for the primary estimand |
| `run_calibration_uno.py` | per-subject out-of-fold horizon predictions, Uno's C, cumulative time-dependent AUC, calibration-in-the-large, calibration slope and calibration curve |
| `validate_concordance.py` | validation of the local concordance implementation against the archived statistic |
| `validate_calibration_pipeline.py` | validation of the calibration pipeline against the archive's own out-of-fold horizon metrics |
| `build_ledger_s8.py` | construction of the numerical ledger |
| `verify_manuscript_v5.py` | automated consistency check of the manuscript against the archived results |
| `make_figures_v5.py`, `make_figures_v5_fix.py`, `make_figures_v5_supp.py` | figure generation |

## Reproducing the headline numbers

```
python scripts_v4/validate_concordance.py          # concordance implementation
python scripts_v4/validate_calibration_pipeline.py # calibration pipeline
python scripts_v4/run_ucsf_donor_grouped_cv.py ucsf_published
python scripts_v4/run_donor_bootstrap_ci.py baylor_published      # and the other three
python scripts_v4/run_calibration_uno.py baylor_published         # and the other three
python scripts_v4/verify_manuscript_v5.py          # manuscript consistency
```

## Human data

All primary data are public (GEO and Zenodo) and were obtained from the repositories named
in the manuscript. The package contains no identifiable human data. The only patient-level
content is the de-identified analysis tables needed to reproduce the reported statistics.

## Ethics

The analysis is a secondary analysis of publicly available, de-identified data. No new
recruitment, intervention or identifiable human data were involved.
"""
(pkg / 'README.md').write_text(readme, encoding='utf-8')

nfiles = sum(1 for p in pkg.rglob('*') if p.is_file())
print(f'S7 checklist items : {len(S7)}')
print(f'Note 1 entries     : 16')
print(f'deposition package : {pkg}')
print(f'  files            : {nfiles}')
print(f'  size             : {sum(p.stat().st_size for p in pkg.rglob("*") if p.is_file())/1e6:.1f} MB')
