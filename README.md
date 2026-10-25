# Auditing the transportability of a meningioma DNA methylation classifier across independent cohorts

Analysis **code** for a deployment-oriented audit of a frozen 403-CpG meningioma
methylation classifier across four external cohorts (Baylor, Tübingen, UCSF, Hong Kong).
The study asks whether conventional external validation of a molecular digital biomarker
can produce results that look credible but are not, and quantifies what a validation
pipeline would have reported had the identity, training-independence, probe-selection and
donor-structure checks been skipped.

This repository is **code only**. It contains no raw data, no derived result tables and no
figures, so it stays small and free of any patient-level content.

## Authors

Jianhuang Huang¹ # \*, Yangyang Lin² #, Qixiu Wang³ #, Jianning Chen¹ \*

¹ Department of Neurosurgery, Affiliated Hospital of Putian University, Putian, Fujian, China
² Department of Neurosurgery, Chinese People's Armed Police Force Fujian Provincial Corps Hospital, Fuzhou, Fujian, China
³ Department of Cerebral Diseases Rehabilitation III, Affiliated Hospital of Liaoning University of Traditional Chinese Medicine, Shenyang, Liaoning, China

\# Co-first authors (Jianhuang Huang, Yangyang Lin, Qixiu Wang contributed equally).
\* Corresponding authors: Jianhuang Huang (teamhuang@sina.com); Jianning Chen (cjn2006@sohu.com).
Correspondence address: Department of Neurosurgery, Affiliated Hospital of Putian University, No. 399 Xihai Road, Licheng District, Putian, Fujian 351100, China.

## Repository layout

```
identity_preprocessing_20261019/   SNP-based donor recovery, numeric-fingerprint overlap
                                   audit, identifier-proven linkage, donor-aware models
strengthening_20261019/            frozen-classifier retrieval/reproduction, robust-panel
                                   transfer, association tests, counterfactual inputs
scripts/                           earlier cohort-assembly and analysis scripts
scripts_v4/                        version-4/5 re-analysis, donor-level bootstrap of the
                                   whole cross-validation procedure, figure generation,
                                   the numerical ledger and the consistency-audit gate
```

Scripts use absolute paths rooted at the project directory and a pinned environment; they
are provided for transparency and re-execution by reviewers, not as an installable package.

## Data availability

All inputs are public, de-identified, secondary data. No new human-participant data were
generated.

- GEO: `GSE189672`, `GSE189521` (Baylor); `GSE304094` (Tübingen); `GSE183647`
  (UCSF/Toronto and Hong Kong); `GSE292325`, `GSE292326` (classifier training and nominal
  held-out releases); `GSE304096` (longitudinal ancillary series); `GSE270371`
  (multicentre, no public patient-level outcomes).
- Zenodo: Hong Kong numbered clinical-outcome table (`11501567`).
- Frozen classifier: the public 403-CpG, k = 10 KNN model with a 141-case reference
  (Bayley and colleagues; Choudhury and colleagues), retrieved from a versioned image with
  verified layer hashes.

## Reproducing the headline numbers

```
python identity_preprocessing_20261019/audit_raw_reference_overlap.py      # training independence
python identity_preprocessing_20261019/audit_snp_fingerprints.py           # donor structure
python scripts_v4/run_ucsf_donor_grouped_cv.py ucsf_published              # donor-grouped CV
python scripts_v4/run_donor_bootstrap_ci.py baylor_published               # whole-procedure CI
python scripts_v4/run_donor_bootstrap_ci_hk_frozen.py                      # HK frozen-classifier CI
python scripts_v4/audit_consistency_v5_1.py                                # numeric consistency gate
```

All computations use fixed seeds (`20261019`–`20261038`). Every reported value is mapped to
its source file in the manuscript's numerical ledger (Supplementary Table S8).

## Scope and honest limits

This is a retrospective methodological audit of publicly available data, not a completed
clinical external validation and not a deployable clinical tool. It provides cohort-internal
predictive modelling and cross-cohort label transportability; it provides no external
validation of a frozen risk model for absolute risk.

## Citation

If you use this code, please cite:

> Huang J, Lin Y, Wang Q, Chen J. *Identity, training independence and selective
> prediction shape the transportability of a meningioma DNA methylation classifier.*
> npj Digital Medicine. [AUTHOR-SPECIFIC: add manuscript DOI / Zenodo archive DOI on
> acceptance.]

A machine-readable `CITATION.cff` is included. The versioned Zenodo DOI is not yet minted
and is intentionally left as a placeholder rather than invented.

## License

[AUTHOR-SPECIFIC: choose and add a LICENSE file. No license file is included yet.]
