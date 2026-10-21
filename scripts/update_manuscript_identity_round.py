# -*- coding: utf-8 -*-
"""Update manuscript abstract/Results/Methods with the identity-preprocessing round."""
import io

p = r"E:\脑膜瘤npj\manuscript_npjDM.md"
s = io.open(p, encoding="utf-8").read()

# ---- abstract ----
a1o = ("We evaluated molecular risk stratification across three outcome-bearing cohorts "
       "(110 tumours from 108 patients; 149; 201) and a 994-sample distribution-level release")
a1n = ("We evaluated molecular risk stratification across four outcome-bearing cohorts "
       "(110 tumours from 108 patients; 149; 201; and a newly linked Hong Kong series of 365 profiles "
       "with identifier-proven clinical outcomes) and a 994-sample distribution-level release")
assert a1o in s, "a1"
s = s.replace(a1o, a1n)

a2o = ("In repeated, event-stratified cross-validation with identical complete cases and folds, "
       "the per-case published labels increased cross-validated C-index beyond clinical covariates in all three cohorts "
       "(Baylor 0.634\u21920.707; T\u00fcbingen 0.541\u21920.564; UCSF 0.710\u21920.741), "
       "with optimism-corrected \u0394C of +0.072, +0.045 and +0.033; "
       "the UCSF group effect remained significant after Holm correction across all fourteen evaluated comparisons "
       "(P = 0.009), and cohort-specific hazard ratios were 4.99 (95% CI 1.24\u201320.05), 3.53 (1.46\u20138.53) "
       "and 3.73 (1.79\u20137.79).")
a2n = ("Under donor-grouped cross-validation that keeps duplicate-donor specimens within the same fold, "
       "the per-case published labels increased cross-validated C-index beyond clinical covariates in all four cohorts "
       "(Baylor 0.603\u21920.695; T\u00fcbingen 0.532\u21920.563; UCSF 0.710\u21920.741; Hong Kong 0.681\u21920.744), "
       "with mean \u0394C of +0.092, +0.031, +0.032 and +0.063 \u2014 entirely positive across all twenty repeated splits "
       "in Baylor, T\u00fcbingen (author KNN) and Hong Kong \u2014 and optimism-corrected \u0394C of +0.072, +0.045, +0.033 "
       "and +0.056; an unpenalized group effect survived Holm correction across all fourteen evaluated comparisons only "
       "in UCSF (P = 0.009; Hong Kong raw P < 10\u207b\u2076), and cohort-specific hazard ratios for the high-risk contrasts "
       "were 4.99 (95% CI 1.24\u201320.05), 3.53 (1.46\u20138.53), 3.73 (1.79\u20137.79) and 5.29 (2.67\u201310.49).")
assert a2o in s, "a2"
s = s.replace(a2o, a2n)

a3o = ("Direct per-case labels from the source publications agreed with the author's KNN classifier in 91.5% of "
       "uniquely classified T\u00fcbingen samples and 87.2% in UCSF.")
a3n = ("Direct per-case labels from the source publications agreed with the author's KNN classifier in 88.8% "
       "(95% CI 85.1\u201391.9) of 349 assigned Hong Kong samples, 91.5% (86.4\u201395.8) of 141 uniquely classified "
       "T\u00fcbingen samples and 87.2% (82.2\u201391.7) of 196 assigned UCSF samples, rising to 93.6%, 96.0% and 93.0% "
       "at the author's 0.75 confidence threshold; an earlier conflicting UCSF estimate arose from a positional clinical "
       "linkage that the direct labels invalidated and is withdrawn.")
assert a3o in s, "a3"
s = s.replace(a3o, a3n)

# ---- Results: cohort section, add HKU + donor-duplicate findings ----
c1o = ("The UCSF clinical table contained 205 labelled records, of which 201 with 76 local-recurrence events had "
       "complete model covariates; its endpoint was local freedom from recurrence, recorded in years [3]. "
       "The 200 methylation profiles with discovery fingerprints were analysed separately from these clinical models. "
       "The 994-sample multicentre methylation release contributed distribution and quality-control information only "
       "(Fig. 1; Table 1).")
c1n = ("The UCSF clinical table contained 205 labelled records, of which 201 with 76 local-recurrence events had "
       "complete model covariates; its endpoint was local freedom from recurrence, recorded in years [3]. "
       "The 200 methylation profiles with discovery fingerprints were analysed separately from these clinical models. "
       "A fourth outcome-bearing series was linked during this audit: the 365 Hong Kong profiles from the same "
       "GSE183647 release carry case-number identifiers that connect them, by exact case-number normalisation with "
       "independently checked demographics, to a numbered clinical-outcome table published by a follow-up study [8]; "
       "302 records shared by both clinical sources had identical event times and status, 7 records with demographic "
       "conflicts and 2 with fewer than 35 callable SNP probes were excluded, leaving 356 eligible records (81 local-"
       "recurrence events in the primary complete set) (Fig. 1; Table 1). The 994-sample multicentre methylation release "
       "contributed distribution and quality-control information only.")
assert c1o in s, "c1"
s = s.replace(c1o, c1n)

c2o = ("SNP-fingerprint checks identified six pairs of specimens with highly concordant donor profiles, "
       "so the 149 records do not represent 149 independently confirmed patients.")
c2n = ("SNP-fingerprint checks identified six pairs of specimens with highly concordant donor profiles "
       "(twelve specimens with archived GSM and barcode identifiers), so the 149 records represent 143 confirmed "
       "independent donors; all analyses were therefore rerun with donor-grouped folds.")
assert c2o in s, "c2"
s = s.replace(c2o, c2n)

c3o = ("The Baylor series comprised 110 primary meningioma specimens from 108 patients [1], "
       "rather than 110 independent patients. The public clinical table does not establish patient clusters, "
       "so analyses of this series were performed at the tumour-record level and remain exploratory.")
c3n = ("The Baylor series comprised 110 primary meningioma specimens from 108 patients [1], rather than 110 "
       "independent patients. The public clinical table does not label patient clusters, but chip-embedded SNP "
       "fingerprints recovered the correspondence directly: one specimen pair shares a strict donor fingerprint and a "
       "second pair merges at the pre-registered overlap-guard threshold, together accounting for the 110-from-108 "
       "structure; all analyses were therefore rerun with donor-grouped folds under both thresholds.")
assert c3o in s, "c3"
s = s.replace(c3o, c3n)

# ---- Results: donor-grouped CV section (insert before calibration section) ----
sec = """
### Donor-grouped validation removes the duplicate-donor leakage concern

Because ten specimens across the outcome-bearing cohorts share a donor with another specimen (one Baylor pair at the
strict threshold plus one guard-level pair, six T\u00fcbingen pairs, and eleven UCSF fingerprint pairs), every comparison
was re-estimated with StratifiedGroupKFold on donor clusters, so that duplicate-donor specimens always fall in the same
fold. All increments remained positive (Table 2, donor-grouped columns): the Baylor published-label \u0394C rose to +0.092
(0.603\u21920.695; all twenty repeated splits positive), the T\u00fcbingen published-label \u0394C was +0.031, the T\u00fcbingen
author-KNN \u0394C was +0.041 (0.552\u21920.593; all twenty splits positive), the Hong Kong published-label \u0394C was +0.063
(0.681\u21920.744; all twenty splits positive) and the Hong Kong author-KNN \u0394C was +0.059 (0.695\u21920.754; all twenty
splits positive); the T\u00fcbingen surrogate \u0394C was +0.025 with nineteen of twenty splits positive. Cluster-robust Cox
models and donor-cluster bootstraps gave concordant directions. Duplicate-donor leakage therefore does not explain the
increments; if anything, donor grouping strengthened them.

"""
anchor = "### Correct censoring weights revise the calibration findings"
assert anchor in s, "sec"
s = s.replace(anchor, sec.lstrip("\n") + anchor)

# ---- Discussion opening: add HKU + donor-grouped ----
d1o = ("Molecular labels carried prognostic information beyond clinical variables, with the clearest reproducible "
       "increment in the UCSF clinical release: the published group labels increased cross-validated discrimination, "
       "reduced the cross-validated Brier score, and their association with recurrence survived correction across all "
       "fourteen evaluated comparisons.")
d1n = ("Molecular labels carried prognostic information beyond clinical variables, with the clearest reproducible "
       "increments in the UCSF and Hong Kong releases: the published group labels increased donor-grouped "
       "cross-validated discrimination in every outcome-bearing cohort \u2014 entirely positive across repeated splits in "
       "Baylor, T\u00fcbingen (author KNN) and Hong Kong \u2014 the association survived correction across all fourteen "
       "evaluated comparisons in UCSF and was uncorrected-significant in Hong Kong, and the molecular layer reduced the "
       "cross-validated Brier score in both.")
assert d1o in s, "d1"
s = s.replace(d1o, d1n)

io.open(p, "w", encoding="utf-8").write(s)
print("manuscript updated: abstract + cohorts + donor-grouped section + discussion")
