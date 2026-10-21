from pathlib import Path
import re

source = Path(__file__).resolve().parent.parent / 'manuscript_npjDM_v5_2_revised.md'
text = source.read_text(encoding='utf-8-sig')

text = text.replace(
    '# Validation integrity and transportability of DNA methylation classifiers in meningioma',
    '# Identity, training independence and selective prediction shape the transportability of a meningioma DNA methylation classifier',
    1,
)
author_block = '''

Jianhuang Huang¹#, Yangyang Lin²#, Qixiu Wang³#, Jianning Chen¹*

¹ Department of Neurosurgery, Affiliated Hospital of Putian University, Putian, Fujian, China
² Department of Neurosurgery, Chinese People's Armed Police Force Fujian Provincial Corps Hospital, Fuzhou, Fujian, China
³ Department of Cerebral Diseases Rehabilitation III, Affiliated Hospital of Liaoning University of Traditional Chinese Medicine, Shenyang, Liaoning, China

# These authors contributed equally.
* Corresponding authors: Jianhuang Huang, MD, Department of Neurosurgery, Affiliated Hospital of Putian University, No. 399 Xihai Road, Licheng District, Putian, Fujian 351100, China. Email: teamhuang@sina.com; and Jianning Chen, same department and address. Email: cjn2006@sohu.com.
'''
if 'Jianhuang Huang¹#' not in text:
    text = text.replace('\n\n---\n\n## Abstract', author_block + '\n\n---\n\n## Abstract', 1)
text = re.sub(
    r'\n\*Manuscript version 5\.1\..*?The complete analysis-revision record is Supplementary Note 1\.\*\n',
    '\n*Version 5.2. Quantitative claims are drawn from the archived identity, preprocessing and validation analyses cited in the Supplementary Information. Cohort-specific estimates are interpreted as internal validation unless explicitly identified as frozen-classifier label agreement.*\n',
    text, count=1, flags=re.S,
)
replacements = {
'''DNA methylation classifiers improve meningioma risk stratification, but cross-institutional reliability remains incompletely characterized. We audited a frozen 403-CpG classifier and published labels across four cohorts, examining identity, training overlap, probe coverage and performance. Numerical fingerprinting found no overlap between 866 external specimens and the 141-case training reference, whereas the nominally held-out release shared 20 of 37 training profiles. The classifier agreed with published group assignments in 91.5% (Tübingen), 87.2% (UCSF) and 88.8% (Hong Kong) of assigned specimens, but positional joining agreed in only 26.5%, below the 33.1% chance rate. At the 0.75 confidence threshold agreement rose to 93-96% while 46% of all hypermitotic versus 7% of Merlin-intact Hong Kong cases were discarded. Molecular labels raised cohort-internal discrimination, corroborated by Uno's C and time-dependent AUC, but whole-procedure donor-level bootstrap intervals excluded zero only in Hong Kong. A five-dimension audit surfaced errors that conventional validation would leave undetected.''':
'''Cross-institutional evaluation of molecular classifiers depends on more than discrimination: test specimens must be independent of training data, molecular and clinical records must be correctly linked, and selective predictions must retain clinically relevant cases. We audited a frozen 403-CpG meningioma classifier across four outcome-bearing cohorts. Fingerprinting found no exact overlap between 866 external specimens and the 141-profile training reference, but 20 of 37 profiles in a nominally held-out release matched training profiles. Agreement with published molecular groups was 91.5% in Tübingen, 87.2% in UCSF and 88.8% in Hong Kong; a rejected positional UCSF join agreed in 26.5%. At a confidence threshold of 0.75, agreement rose to 93–96%, while only 54.3% of all Hong Kong hypermitotic cases were retained, compared with 93.3% of Merlin-intact cases. Donor-grouped cross-validation showed positive point estimates of incremental discrimination in all four cohorts, but whole-procedure donor-bootstrap intervals excluded zero only in Hong Kong. These findings show how identity, independence and coverage checks change the interpretation of classifier transportability.''',
'''Gross-total resection and WHO grade remain the backbone of postoperative meningioma management, yet recurrence risk within grade spans a wide range and informs decisions on adjuvant radiotherapy and on surveillance intensity. DNA methylation profiling has produced molecular classifications with documented prognostic value: the meningioma classes of Bayley and colleagues, the Merlin-intact, immune-enriched and hypermitotic groups of Choudhury and colleagues, the Toronto nomenclature of Wang and colleagues, and recurrence-oriented classifiers developed in single-centre series. That methylation carries prognostic information beyond histology is not in dispute.''':
'''Postoperative meningioma risk assessment still relies heavily on histological grade and extent of resection, although recurrence risk varies within grade. DNA methylation classifiers have identified biologically and prognostically informative groups, including the meningioma classes of Bayley and colleagues, the Merlin-intact, immune-enriched and hypermitotic groups of Choudhury and colleagues, the Toronto nomenclature of Wang and colleagues, and recurrence-oriented classifiers developed in single-centre series. For these classifiers to support decisions beyond their development cohorts, their implementation and evidence for transport to new specimens must be independently verifiable.''',
'''What is uncertain is a more operational question. A published molecular classifier is a digital biomarker, and is expected to satisfy the conditions that separate a research variable from a transferable clinical instrument: an implementation another laboratory can execute unchanged; test specimens demonstrably independent of those used to fit it; a verifiable link between a molecular profile and the patient whose outcome is recorded; analytical performance that survives differences in probe coverage and preprocessing; discrimination, calibration and absolute-risk accuracy reproducible in new patients; and net benefit in a defined clinical decision. These are properties of a validation chain, not of a correlation. The field's own guidance for digital biomarker manuscripts requires that a separate external test set be used and that zero overlap between training and test data be demonstrated [1]; such demonstrations are rarely reported.''':
'''For a digital biomarker, an apparent validation result is interpretable only if several links in the evidence chain hold: test specimens are independent of the training reference, sample identity is established before clinical outcomes are joined, preprocessing and feature selection are reproducible, and performance is reported with uncertainty and coverage. Discrimination alone does not establish calibrated absolute risk or clinical utility. Reporting guidance for healthcare AI emphasizes independent evaluation and transparent validation methods [1], but these links are not always auditable in published classifier evaluations.''',
'''Three failure modes make this a live problem. A release described as held out can share individuals with the model's training reference, so agreement estimated on it measures memory rather than transportability. Molecular profiles and clinical rows can be joined by position rather than identifier, and where an independent per-case label table exists such a join can be shown to be wrong. Probe-selection and preprocessing choices made in one cohort can silently determine the transfer result in another, so that a single sort-direction error can change a reported agreement materially.''':
'''Several failure modes can change the apparent result: a nominally held-out release may contain training profiles; a positional join may attach a molecular profile to the wrong clinical record; and a feature-selection or preprocessing error may shift transfer metrics. Duplicate specimens from one donor can also cross validation folds, making internal performance appear more stable than it is.''',
'''We report a deployment-oriented evaluation, not a demonstration of clinical utility. The primary objective was to establish, by numerical fingerprinting, whether the four external cohorts are genuinely independent of the frozen classifier's training reference, and to measure how far the frozen classifier reproduces independently published per-case group labels in those cohorts. The secondary objective was to quantify incremental prognostic information under one statistical framework, with donor-grouped validation. A third objective was to quantify what a conventional validation pipeline would have reported had these checks not been performed. We deliberately do not claim external validation of a frozen clinical risk equation, because the four cohorts share neither an endpoint definition nor a covariate set that would permit one.''':
'''We evaluated training-reference overlap, specimen identity, cross-nomenclature agreement and selective prediction for a frozen classifier, then estimated cohort-internal prognostic increments using donor-grouped validation. We also quantified how omitted checks would have changed reported results. Because endpoints and covariate sets differ across cohorts, we do not treat these analyses as external validation of a frozen absolute-risk model or as evidence of clinical utility.''',
'''The four outcome-bearing cohorts are legitimate independent test sets for this frozen classifier, which is the condition the field's digital-biomarker guidance requires [1].''':
'''The fingerprint audit supports independence of these external specimens from the 141-profile training reference; it does not by itself establish clinical validity or generalizability.''',
'''Taken together, these intervals and the point estimates support the following position: the direction of the discrimination increment is positive in all four cohorts, its magnitude is between +0.027 and +0.092, and it is statistically distinguishable from zero in Hong Kong on a donor-level resampling of the whole procedure.''':
'''Across the four primary published-label comparisons, donor-grouped cross-validation produced positive point estimates of incremental discrimination (+0.027 to +0.092). The whole-procedure donor-level interval excluded zero only in Hong Kong; intervals in Baylor, Tübingen and UCSF included zero. These results are compatible with incremental prognostic information, but do not establish a consistent non-zero gain across cohorts.''',
'''The point of an audit is only established by what it changes. We therefore computed, for each failure mode, the result a conventional validation pipeline would have produced had the corresponding check not been performed (Table 6). These are not hypothetical quantities; each is already in the record.''':
'''To show the practical effect of these checks, we compared audited results with the values that would have been reported if each corresponding check had been omitted (Table 6). These counterfactual comparisons quantify the impact of specific analysis choices; they do not estimate how often such failures occur in the broader literature.''',
'''Four of these five errors would have been statistically invisible. None produces a failed test, an implausible estimate or a warning sign; each produces a number that looks reportable. That is the argument for treating identity, independence, probe selection and donor structure as validation steps rather than as quality control, and it is the empirical basis for the framework below.''':
'''The comparisons illustrate that standard performance summaries may not reveal identity, overlap, feature-selection or donor-structure errors. We therefore recommend treating these checks as part of validation design and reporting. The audit criteria used here are study-specific and have not been validated as a general readiness score.''',
'''Molecular labels carried prognostic information beyond clinical variables in every outcome-bearing cohort, and the frozen classifier reproduced independently published group labels in three institutions on specimens with no training-reference overlap. The direction of the prognostic information is stable under every discrimination measure we applied: donor-grouped repeated cross-validation, Uno's C and cumulative time-dependent AUC all give positive increments of +0.024 to +0.092. The precision is more limited than an apparent bootstrap interval would suggest.''':
'''The frozen classifier showed cross-nomenclature agreement with published group labels in three cohorts, and fingerprinting found no exact overlap between these specimens and the training reference. In the separate cohort-internal prognostic analyses, point estimates were positive across the primary comparisons and were directionally similar under donor-grouped cross-validation, Uno's C and time-dependent AUC. Precision was limited: when the full validation procedure was repeated in donor-level bootstrap samples, the interval for the primary increment excluded zero only in Hong Kong.''',
'''Four of five routes''':
'''Three of five evaluated routes''',
'''The complete audit package — scripts, per-sample fold assignments, out-of-fold predictions, bootstrap draws, hash manifests, the numerical ledger and the retrieved classifier artefacts with verified image-layer hashes — accompanies the submission as a frozen, versioned snapshot made available to reviewers, and will be minted with a public archive DOI on acceptance. The authors declare no competing interests. CRediT roles are to be finalised by the author group before submission.''':
'''**Data availability.** The source datasets are available from GEO under the accession numbers listed above and from the linked public source repositories. The derived analysis files and complete audit package must be deposited in a repository and made available to editors and reviewers before submission. **[AUTHOR ACTION: insert repository name, persistent DOI or accession, version, and any third-party data restrictions; confirm that the deposited package contains the scripts, fold assignments, predictions, bootstrap draws, hash manifests and numerical ledger cited here.]**

**Code availability.** **[AUTHOR ACTION: provide the DOI or stable repository link and version for the custom code central to the reported results, and state the terms of access. The journal requires a Code Availability section; a future promise to mint a DOI after acceptance is not sufficient for peer review.]**

**Ethics approval.** **[AUTHOR CONFIRMATION REQUIRED: state the institutional review board/ethics determination or applicable exemption for this secondary analysis, with institution and reference number if applicable. Do not infer exemption solely from public, de-identified source data.]**

**Competing interests.** **[AUTHOR CONFIRMATION REQUIRED: confirm the declaration for all authors and submit the journal's signed statement.]**

**Author contributions.** **[AUTHOR ACTION: complete CRediT roles and confirm all authors have approved the submitted version.]**''',
}
for old, new in replacements.items():
    if old not in text:
        print('NOT FOUND:', old[:100])
    text = text.replace(old, new)

# Correct legacy selective-prediction denominator in the counterfactual table.
text = text.replace(
    '| Selective prediction | Coverage and class-specific retention alongside confidence thresholds | Agreement 93.6% at the 0.75 threshold in Hong Kong | Same agreement, with coverage 0.773 and hypermitotic retention 0.607 | Hides that 39% of the highest-risk class is discarded |',
    '| Selective prediction | Coverage and class-specific retention alongside confidence thresholds | Agreement 93.6% at the 0.75 threshold in Hong Kong | Same agreement, with coverage 0.773 and hypermitotic retention 0.543 when all true class cases are the denominator | Hides that 46% of hypermitotic cases are discarded |',
)
text = text.replace('**Figure 6.** *(figures_v5/fig6_deployment_audit.pdf)* **Deployment-readiness audit and counterfactual evidence.** (a) The five dimensions with audit-defined criteria and three-level verdicts. (b) The five counterfactuals of Table 6 plotted as pairs, showing the reportable number a conventional pipeline would have produced against the audited value.',
                  '**Figure 6.** *(figures_v5/fig6_deployment_audit.pdf)* **Validation audit and counterfactual evidence.** (a) Five audit dimensions with study-defined criteria and three-level verdicts. (b) Six counterfactual comparisons across five audit domains, showing the value a conventional pipeline would have reported against the audited result.')
text = text.replace('None of these five errors produces a failed test, an implausible estimate or a statistical warning.',
                  'These omitted checks did not produce a failed test, an implausible estimate or a statistical warning in the analyses examined here.')
text = text.replace('### Table 5. Global association tests, sixteen-test family',
                  '### Table 5. Exploratory global association tests, sixteen-test family')
text = text.replace('### Table 6. Counterfactual analysis: what conventional validation would have reported',
                  '### Table 6. Counterfactual analysis: impact of omitted validation checks')
text = text.replace('The point of an audit is only established by what it changes.',
                  'To show the practical effect of these checks,')
text = text.replace('The point of an audit is only established by what it changes. We therefore computed, for each failure mode, the result a conventional validation pipeline would have produced had the corresponding check not been performed (Table 6). These are not hypothetical quantities; each is already in the record.',
                  'To show the practical effect of these checks, we compared audited results with the values that would have been reported if each corresponding check had been omitted (Table 6). These counterfactual comparisons quantify the impact of specific analysis choices; they do not estimate how often such failures occur in the broader literature.')

# Avoid overstating comparison to chance as a formal null test.
text = text.replace('below the 33.1% expected by chance from the observed label marginals',
                    'compared with 33.1% expected from the observed label marginals')
text = text.replace('lies below its 33.1% marginal chance rate',
                    'is below the 33.1% agreement expected from the observed label marginals')

# Remove version-history prose from the manuscript and keep the abstract concise.
text = re.sub(r'\n\*Version 5\.2\..*?\*\n', '\n', text, count=1)
text = re.sub(
    r'(## Abstract\n\n).*?(\n---)',
    r'''\1External evaluation of molecular classifiers requires independent test specimens, verified sample identity and reporting of prediction coverage. We audited a frozen 403-CpG meningioma classifier and prognostic analyses across four outcome-bearing cohorts. Fingerprinting found no exact overlap between 866 external specimens and the 141-profile training reference, whereas 20 of 37 profiles in a nominally held-out release matched training. Agreement with published molecular groups was 91.5% in Tübingen, 87.2% in UCSF and 88.8% in Hong Kong. A rejected positional UCSF linkage agreed in 26.5%. At a confidence threshold of 0.75, agreement reached 93–96%, but retention in Hong Kong was 54.3% for all hypermitotic cases and 93.3% for Merlin-intact cases. Donor-grouped cross-validation yielded positive point estimates of incremental discrimination in four cohorts, but whole-procedure donor-bootstrap intervals excluded zero only in Hong Kong. These results show that overlap, identity and coverage checks materially qualify apparent classifier transportability.\2''',
    text, count=1, flags=re.S,
)

# Replace the long discussion with a tighter, evidence-bounded interpretation.
discussion = '''## Discussion

This audit shows that claims about classifier transportability depend on the integrity of the validation chain. None of the 866 external specimens matched the 141-profile training reference, while 20 of 37 profiles in a nominally held-out release did. In three cohorts, the frozen classifier agreed with independently published group labels in 87.2–91.5% of assigned specimens. In UCSF, an identifier-based linkage was essential: the tested positional join agreed with published labels in 26.5% of cases. Confidence filtering increased selective agreement to 93–96%, but retained only 54.3% of all Hong Kong hypermitotic cases, compared with 93.3% of Merlin-intact cases. In cohort-internal prognostic analyses, donor-grouped point estimates were positive, but the whole-procedure donor-level interval excluded zero only in Hong Kong.

These findings extend prior work establishing prognostic associations for meningioma methylation groups by measuring a different property: whether a frozen implementation and its validation evidence can be traced across cohorts. The distinction matters. Molecular labels were added to clinical models whose coefficients were refitted within each cohort, and endpoint definitions, follow-up units and covariate availability differed. The resulting ΔC estimates are internal validation, not evidence that a fixed risk equation transports to a new hospital. The data therefore support cautious evidence of prognostic signal, not external validation of absolute risk or clinical benefit.

Selective prediction also requires a coverage denominator that includes specimens without an assigned class. Reporting only accuracy among confidently classified samples can conceal preferential rejection of a clinically important group. Here, tied votes and low confidence both reduced coverage, with hypermitotic cases disproportionately affected. This pattern is a reason to report class-specific retention alongside selective accuracy; it does not by itself establish how a clinical decision based on the classifier would affect patient outcomes.

Several limitations constrain interpretation. Baylor per-case labels could not be linked to methylation profiles at the required coverage, so that cohort contributed to the training-independence audit but not classifier-label agreement. The prognostic analyses were retrospective, cohort-specific and not preregistered; the donor-level resampling intervals were wide, and the Baylor bootstrap was frequently non-estimable because comparable pairs disappeared within resampled folds. Endpoints and covariates could not be harmonized, adjuvant radiotherapy was unavailable in most public clinical tables, and event counts limited calibration assessment. We did not test a frozen absolute-risk model in an independent institution or evaluate a prospectively defined treatment decision. The audit framework itself is study-defined and has not been validated as a general readiness score.

Future studies should use a prespecified endpoint and covariate set, establish specimen-to-patient linkage before outcomes are analyzed, document training–test independence, keep donor clusters intact during validation, and report discrimination, calibration and class-specific coverage. Clinical utility should be evaluated only after a decision and its threshold have been defined. Until such evidence is available, the classifier's cross-cohort label agreement should not be interpreted as proof of clinical readiness.

'''
text = re.sub(r'## Discussion\n.*?(?=\n---\n\n## Methods)', discussion, text, count=1, flags=re.S)

# Nature Portfolio requires disclosure when an LLM assists with manuscript preparation.
anchor = '### Reporting standards, ethics, data and code\n'
if '**Use of AI assistance.**' not in text:
    text = text.replace(anchor, anchor + '''
**Use of AI assistance.** OpenAI Codex was used to assist with manuscript language and structure. It was not used to generate data or perform analyses. The authors are responsible for reviewing all content and for the accuracy of the manuscript.

''', 1)

# Make the journal-facing statements explicit, separate sections and retain author-only facts as placeholders.
old_statements = '''### Reporting standards, ethics, data and code

**Use of AI assistance.** OpenAI Codex was used to assist with manuscript language and structure. It was not used to generate data or perform analyses. The authors are responsible for reviewing all content and for the accuracy of the manuscript.


Analyses were designed against TRIPOD+AI [2], the original TRIPOD statement [3], REMARK [4] and the field's guidance for digital biomarker manuscripts [1], with extraction of cohort inputs against CHARMS [5]; items specific to model development are marked not applicable with reasons (Supplementary Table S7). No claim uses the terms clinically deployable, practice-changing, prospectively validated or generalizable.

All analysed data are publicly available, de-identified secondary data from GEO and a public Zenodo deposit containing no directly identifying information; no new recruitment, intervention or identifiable human data occurred, so institutional review board approval was not required. **Data availability.** The source datasets are available from GEO under the accession numbers listed above and from the linked public source repositories. The derived analysis files and complete audit package must be deposited in a repository and made available to editors and reviewers before submission. **[AUTHOR ACTION: insert repository name, persistent DOI or accession, version, and any third-party data restrictions; confirm that the deposited package contains the scripts, fold assignments, predictions, bootstrap draws, hash manifests and numerical ledger cited here.]**

**Code availability.** **[AUTHOR ACTION: provide the DOI or stable repository link and version for the custom code central to the reported results, and state the terms of access. The journal requires a Code Availability section; a future promise to mint a DOI after acceptance is not sufficient for peer review.]**

**Ethics approval.** **[AUTHOR CONFIRMATION REQUIRED: state the institutional review board/ethics determination or applicable exemption for this secondary analysis, with institution and reference number if applicable. Do not infer exemption solely from public, de-identified source data.]**

**Competing interests.** **[AUTHOR CONFIRMATION REQUIRED: confirm the declaration for all authors and submit the journal's signed statement.]**

**Author contributions.** **[AUTHOR ACTION: complete CRediT roles and confirm all authors have approved the submitted version.]**
'''
new_statements = '''### Reporting standards

Analyses were designed against TRIPOD+AI [2], the original TRIPOD statement [3], REMARK [4] and the field's guidance for digital biomarker manuscripts [1], with extraction of cohort inputs against CHARMS [5]. Items specific to model development are marked not applicable with reasons in Supplementary Table S7.

---

## Ethics approval

**[AUTHOR CONFIRMATION REQUIRED: state the institutional review board or ethics determination for this secondary analysis, including the institution and reference number, or the applicable exemption. Public, de-identified source data alone do not establish the local ethics determination.]**

## Data availability

The source datasets are available from GEO under the accession numbers listed above and from the linked public source repositories. **[AUTHOR ACTION: provide the repository name, persistent DOI or accession, version, and any third-party data restrictions for the derived analysis files and audit package. Confirm that the deposit includes the scripts, fold assignments, predictions, bootstrap draws, hash manifests and numerical ledger needed to reproduce the reported findings, and make it available to editors and reviewers at submission.]**

## Code availability

**[AUTHOR ACTION: provide the DOI or stable repository link and version for custom code central to the results, and state its access terms. A future promise to mint a DOI after acceptance is not sufficient for peer review.]**

## Competing interests

**[AUTHOR CONFIRMATION REQUIRED: confirm the declaration for all authors and submit the journal's signed statement.]**

## Author contributions

**[AUTHOR ACTION: complete the CRediT roles and confirm that all authors have approved the submitted version.]**
'''
if old_statements not in text:
    print('Could not find original declarations block for section separation')
else:
    text = text.replace(old_statements, new_statements, 1)

# Reconcile internal inconsistencies and avoid conclusions not supported by this design.
text = text.replace('The abstract', 'The abstract')
text = text.replace('Because a bridge that reproduces the published groups in only 26.5% of cases — compared with 33.1% expected from the observed label marginals — and sex at coin-flip level is not an imperfect bridge but a wrong one,',
                    'Because this bridge reproduced the published groups in 26.5% of cases, compared with 33.1% expected from the observed label marginals, and sex concordance was 53.5%,')
text = text.replace('Group agreement 26.5% is below the 33.1% agreement expected from the observed label marginals (sex at coin-flip 53.5%); a mis-linkage, not noise',
                    'Group agreement 26.5% compared with 33.1% expected from observed label marginals, with sex concordance of 53.5%; consistent with a mis-linkage')
text = text.replace('A confidence threshold adopted for convenience would therefore degrade performance exactly where it matters clinically, and this is a property of the classifier and its confidence estimate, not of any single cohort.',
                    'This selective-accuracy gain therefore comes with lower coverage in the hypermitotic group and should be interpreted with class-specific retention, not agreement alone.')
text = text.replace('All analysed data are publicly available, de-identified secondary data from GEO and a public Zenodo deposit containing no directly identifying information; no new recruitment, intervention or identifiable human data occurred, so institutional review board approval was not required. ', '')
text = text.replace('**Figure 6.** *(figures_v5/fig6_deployment_audit.pdf)* Deployment-readiness audit and counterfactual evidence.** (a) The five dimensions with audit-defined criteria and three-level verdicts. (b) The five counterfactuals of Table 6 plotted as pairs, showing the reportable number a conventional pipeline would have produced against the audited value.',
                    '**Figure 6.** *(figures_v5/fig6_deployment_audit.pdf)* **Validation audit and counterfactual evidence.** (a) Five audit dimensions with study-defined criteria and three-level verdicts. (b) Six counterfactual comparisons across five audit domains, showing the value a conventional pipeline would have reported against the audited result.')
text = text.replace('### A five-dimension deployment-readiness audit', '### A five-dimension validation audit')
text = text.replace('### Table 7. Deployment-readiness framework: audit-defined criteria and three-level verdicts', '### Table 7. Five-dimension validation audit: study-defined criteria and verdicts')
text = text.replace('Group agreement 26.5% lies below its 33.1% marginal chance rate (sex at coin-flip 53.5%); a mis-linkage, not noise',
                    'Group agreement 26.5% compared with 33.1% expected from the observed label marginals (sex concordance 53.5%); consistent with a mis-linkage')
text = text.replace('The four external cohorts are legitimate independent test sets for this frozen classifier',
                    'These external specimens are independent of the training reference at the fingerprint level')

source.write_text(text, encoding='utf-8-sig')
print('Updated', source, 'characters=', len(text), 'words=', len(text.split()))
