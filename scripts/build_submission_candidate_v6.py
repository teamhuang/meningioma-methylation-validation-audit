from pathlib import Path
import re, shutil

ROOT = Path(r'E:\脑膜瘤npj')
OUT = ROOT / 'npjDM_submission_candidate_v6'
OUT.mkdir(exist_ok=True)
src = (ROOT / 'manuscript_npjDM_v5_2_revised.md').read_text(encoding='utf-8')
old = (ROOT / 'manuscript_npjDM_v5_1.md').read_text(encoding='utf-8')

# Author block supplied by the corresponding author in this conversation.
src = src.replace(
    'Jianhuang Huang¹#*, Yangyang Lin²#, Qixiu Wang³#, Jianning Chen¹*',
    'Jianhuang Huang¹*#, Yangyang Lin²#, Qixiu Wang³#, Jianning Chen¹*'
)
src = src.replace(
    '¹ Department of Neurosurgery, Affiliated Hospital of Putian University, Putian, Fujian, China\n'
    "² Department of Neurosurgery, Chinese People's Armed Police Force Fujian Provincial Corps Hospital, Fuzhou, Fujian, China\n"
    '³ Department of Cerebral Diseases Rehabilitation III, Affiliated Hospital of Liaoning University of Traditional Chinese Medicine, Shenyang, Liaoning, China\n\n'
    '# Jianhuang Huang, Yangyang Lin and Qixiu Wang contributed equally.\n'
    '* Corresponding authors: Jianhuang Huang, MD, Department of Neurosurgery, Affiliated Hospital of Putian University, No. 399 Xihai Road, Licheng District, Putian, Fujian 351100, China. Email: teamhuang@sina.com; and Jianning Chen, same department and address. Email: cjn2006@sohu.com.',
    '¹ Department of Neurosurgery, Affiliated Hospital of Putian University, Putian, Fujian, China\n'
    "² Department of Neurosurgery, Chinese People's Armed Police Force Fujian Provincial Corps Hospital, Fuzhou, Fujian, China\n"
    '³ Department of Cerebral Diseases Rehabilitation III, Affiliated Hospital of Liaoning University of Traditional Chinese Medicine, Shenyang, Liaoning, China\n\n'
    '# Jianhuang Huang, Yangyang Lin and Qixiu Wang contributed equally.\n'
    '* Correspondence: Jianhuang Huang, MD, and Jianning Chen, Department of Neurosurgery, Affiliated Hospital of Putian University, No. 399 Xihai Road, Licheng District, Putian, Fujian 351100, China. Emails: teamhuang@sina.com (J.H.); cjn2006@sohu.com (J.C.).'
)

# Remove references to unavailable materials, then cite the files included in the package.
src = src.replace(
    'Four distinct quantities must not be conflated, because they answer different questions (Table 4 footnote). The donor-grouped mean cross-validated ΔC is the primary effect estimate. The split-to-split spread across the twenty repetitions describes algorithmic stability under re-partitioning of the same records and is not a sampling interval; it is reported in Supplementary Table S1 together with the count of positive repetitions. The apparent bootstrap percentile interval describes the sampling behaviour of an in-sample ΔC estimator. The optimism-corrected ΔC is a distinct internal-validation estimate.',
    'The donor-grouped mean cross-validated ΔC is the primary effect estimate. Split-to-split variation across the twenty repetitions describes algorithmic stability on these records, not sampling uncertainty. The apparent bootstrap percentile interval describes an in-sample estimator, whereas the optimism-corrected ΔC is a separate internal-validation estimate; neither substitutes for the whole-procedure donor bootstrap interval reported below and in Table 4.'
)
src = src.replace('Five quantities are reported and are not interchangeable. (i) The donor-grouped mean cross-validated ΔC is the primary effect estimate. (ii) The split-to-split spread across repetitions describes algorithmic stability on fixed records and is not a sampling interval; with the positive-repeat counts it is reported in Supplementary Table S1. (iii)', 'Four quantities are reported and are not interchangeable. (i) The donor-grouped mean cross-validated ΔC is the primary effect estimate. (ii) The apparent paired bootstrap percentile interval over 1,000 donor-cluster resamples describes the sampling behaviour of an in-sample estimator. (iii)')
src = src.replace(' (iv) The optimism-corrected ΔC is a distinct internal-validation estimate, computed as the difference between each bootstrap sample and the original data evaluated with the bootstrap-fitted model. (v) The donor-level bootstrap interval', ' (iv) The donor-level bootstrap interval')
src = src.replace('A fifth quantity is the interval that actually corresponds to the primary estimand.', 'The whole-procedure donor-level bootstrap interval is the interval that corresponds to the primary estimand.')
src = src.replace('Supplementary Table S7', 'Supplementary Table S1')
src = src.replace('Supplementary Figure S4', 'Supplementary Figure S2')
src = src.replace('Supplementary Figure S6', 'Supplementary Figure S4')
src = src.replace('Supplementary Table S3', 'Supplementary Data 1')
src = src.replace('Supplementary Data 1. Longitudinal specimen-to-cluster mapping', 'Supplementary Data 2. Longitudinal specimen-to-cluster mapping')
src = src.replace('Supplementary Data 2. Identifier-proven UCSF', 'Supplementary Data 3. Identifier-proven UCSF')
src = src.replace('Supplementary Data 3. Corrected robust-panel', 'Supplementary Data 1. Corrected robust-panel')
src = src.replace('Supplementary Data 3)', 'Supplementary Data 1)')
src = src.replace('in Supplementary Figure S4 and Supplementary Data 1.', 'in Supplementary Figure S4 and Supplementary Data 2.')
src = src.replace('provided as Supplementary Data 2.', 'provided as Supplementary Data 3.')

# Reconcile all main-text table citations with the available tables (old Tables 6-9 become 5-8).
for old_n, marker in [(6, 'TABLE_FIVE'), (7, 'TABLE_SIX'), (8, 'TABLE_SEVEN'), (9, 'TABLE_EIGHT')]:
    src = re.sub(rf'\bTable {old_n}\b', marker, src)
src = src.replace('TABLE_FIVE', 'Table 5').replace('TABLE_SIX', 'Table 6').replace('TABLE_SEVEN', 'Table 7').replace('TABLE_EIGHT', 'Table 8')

# Add figure callouts at first substantive mention, in numerical order.
replacements = [
    ('assembled (Table 1).', 'assembled (Table 1; Fig. 1a,c).'),
    ('the result is robust to the matching threshold.', 'the result is robust to the matching threshold (Fig. 2a; Supplementary Fig. 1).'),
    ('and compared with the per-case group labels published by the cohort authors and linked by verified identifiers.', 'and compared with the per-case group labels published by the cohort authors and linked by verified identifiers (Fig. 2b,c).'),
    ('and it discards patients non-randomly (Table 3).', 'and it discards patients non-randomly (Table 3; Fig. 3).'),
    ('Under donor grouping, every methylation-derived route increased mean cross-validated discrimination (Table 4).', 'Under donor grouping, every methylation-derived route increased mean cross-validated discrimination (Table 4; Fig. 4).'),
    ('Discrimination and absolute-risk accuracy dissociated.', 'Discrimination and absolute-risk accuracy dissociated (Fig. 5).'),
    ('(Table 5). These counterfactual comparisons', '(Table 5; Fig. 6). These counterfactual comparisons'),
    ('were assembled (Table 1; Fig. 1a,c).', 'were assembled (Table 1; Fig. 1a,c).'),
    ('marks definitions that could not be verified from public data. Analyses were not pooled across cohorts.', 'marks definitions that could not be verified from public data (Table 7). Analyses were not pooled across cohorts.'),
]
for a,b in replacements:
    src = src.replace(a,b)

# Correct the current results wording and provide file-specific supplementary citations.
src = src.replace('These comparisons quantify the impact of specific analysis choices;', 'These comparisons quantify the impact of specific analysis choices;')
src = src.replace('The comparisons illustrate that standard performance summaries may not reveal identity, overlap, feature-selection or donor-structure errors.', 'The corrected robust-panel probe list and per-target shared-probe counts are provided as Supplementary Data 3. The comparisons illustrate that standard performance summaries may not reveal identity, overlap, feature-selection or donor-structure errors.')
src = src.replace('Decision-curve quantities were recomputed with corrected weights and are reported as exploratory relative comparisons in Supplementary Figure S2;', 'Decision-curve quantities were recomputed with corrected weights and are reported as exploratory relative comparisons in Supplementary Figure S2;')
src = src.replace('in Supplementary Results. In brief,', 'in Supplementary Figure S4 and Supplementary Data 1. In brief,')
src = src.replace('the same identifier-proven route linked the 365 Hong Kong profiles', 'the identifier-proven UCSF profile linkage is provided as Supplementary Data 2. The same identifier-proven route linked the 365 Hong Kong profiles')

# Insert missing author-facing statements and an LLM-use disclosure required by the journal.
src = src.replace('## Competing interests\n\n**[AUTHOR CONFIRMATION REQUIRED: confirm the declaration for all authors and submit the journal\'s signed statement.]**',
                  '## Acknowledgements\n\n**[AUTHOR CONFIRMATION REQUIRED: state all funding sources and grant numbers, or confirm that the study received no funding; describe the role of each funder. Add any other brief acknowledgements or state that none are applicable.]**\n\n## Competing interests\n\n**[AUTHOR CONFIRMATION REQUIRED: provide the declaration for all four authors and complete the journal\'s required competing-interests form.]**')
src = src.replace('### Reporting standards\n\n', '### Reporting standards\n\n')
anchor = 'Analyses were designed against TRIPOD+AI [6], the original TRIPOD statement [7], REMARK [8] and the field\'s guidance for digital biomarker manuscripts [5], with extraction of cohort inputs against CHARMS [9]. Items specific to model development are marked not applicable with reasons in Supplementary Table S1.'
src = src.replace(anchor, 'Analyses were designed against TRIPOD+AI [6], the original TRIPOD statement [7], REMARK [8] and the field\'s guidance for digital biomarker manuscripts [5], with extraction of cohort inputs against CHARMS [9]. Model-development items are not applicable because no new prediction model was fitted.')
src = src.replace(anchor, anchor + '\n\n### Use of generative AI\n\nOpenAI Codex was used to assist with English-language editing and manuscript organization. It was not used to generate or alter source data, perform statistical analyses or select results. **[AUTHOR CONFIRMATION REQUIRED: confirm that all authors reviewed and approved the final text and remain responsible for its content.]**')

# Add data files to the SI legends and renumber legacy figure assets as a contiguous, supported set.
src = re.sub(r'\*\*Supplementary Figure S1\.\*\*.*?(?=\*\*Supplementary Table)',
'''**Supplementary Figure S1.** Training-reference fingerprint audit: per-cohort distances to the nearest reference profile and identifier checks.\n**Supplementary Figure S2.** Exploratory decision curves with corrected censoring weights, presented as relative comparisons without a clinical threshold.\n**Supplementary Figure S3.** Robust-panel transfer: the corrected panel, shared-probe counts and rejection flags.\n**Supplementary Figure S4.** Longitudinal ancillary series: per-cluster labels against the constant-class comparator, paired Δβ distribution and donor-fingerprint conflict flags.\n\n**Supplementary Data 1.** Corrected robust-panel profile-level results and probe-level metadata.\n**Supplementary Data 2.** Longitudinal specimen-to-cluster mapping, donor-fingerprint conflict flags and the constant-class baseline.\n**Supplementary Data 3.** Identifier-proven UCSF molecular-profile linkage and cohort metadata used in the identity audit.\n''', src, flags=re.S)

# References and legends in v5.2 end the article; restore the actual display tables from v5.1,
# omitting the legacy unpenalized association Table 5 and internal numerical ledger S8.
old_lines = old.splitlines()
start = next(i for i,l in enumerate(old_lines) if l.strip() == '## Tables')
end = next(i for i,l in enumerate(old_lines[start+1:], start+1) if l.strip() == '## Figure legends')
block = old_lines[start+1:end]
sections = {}
cur = None
for line in block:
    m = re.match(r'### Table (\d+)\.', line)
    if m:
        cur = int(m.group(1)); sections[cur] = [line]
    elif cur is not None:
        sections[cur].append(line)
tables = []
for n in (1,2,3,4,6,7,8,9):
    lines = sections[n]
    if n == 6:
        lines = [re.sub(r'^### Table 6\.', '### Table 5.', x) for x in lines]
        lines = [x.replace('five errors', 'six comparisons across five audit domains') for x in lines]
    elif n > 6:
        lines = [re.sub(rf'^### Table {n}\.', f'### Table {n-1}.', x) for x in lines]
    tables.extend(lines)
    tables.append('')
table_block = '\n'.join(['## Tables', ''] + tables).strip()
table_block = table_block.replace('with the positive-repeat counts it is reported in Supplementary Table S1, not here;', 'with split-to-split stability summarized in the Results rather than in a separate supplementary table;')
table_block = re.sub(
    r'Twenty repetitions of five-fold cross-validation with StratifiedGroupKFold.*?(?=\n\n### Table 5)',
    'Twenty repetitions of five-fold cross-validation with StratifiedGroupKFold on SNP donor clusters (runtime-asserted); identical complete cases and folds were used for both models, and concordance was computed from comparable pairs within each test fold. **Four estimands are reported and are not interchangeable:** (i) the donor-grouped mean cross-validated ΔC is the primary effect estimate; (ii) split-to-split variation describes algorithmic stability on these fixed records and is not a sampling interval; (iii) the apparent paired bootstrap interval describes the sampling behaviour of an in-sample ΔC estimator; and (iv) the optimism-corrected ΔC is a distinct internal-validation estimate. The donor-level bootstrap interval in the table is the only interval corresponding to the primary estimand: the complete repeated cross-validation was re-executed in each of 1,000 donor-cluster resamples. These percentile intervals are not assumed to have exactly nominal 95% coverage because resamples with collapsed Cox models or comparable-pair sets were counted as failures. Baylor succeeded in 803/1,000 draws; 196 failures arose from loss of comparable pairs and one from Cox non-convergence. The Baylor interval is therefore an exploratory summary of a partially degenerate distribution; Tübingen succeeded in 996/1,000 draws. Whole-procedure intervals were computed for the four published-label routes and the Hong Kong frozen-classifier route only. The UCSF donor-grouped and record-level rows use the same 198 records and are not distinguishable in donor-level resampling. The Baylor frozen-classifier routes differ by one specimen because the raw-IDAT and processed-beta routes exclude different tied-vote specimens.',
    table_block,
    flags=re.S,
    count=1
)

# Journal-facing availability and author declarations remain after Methods; display tables
# are placed at the end of the manuscript, after references and figure legends.

# In a final draft, correct former Table 8/9 mentions in Methods and table captions are contiguous.
src = src.replace('Table 8 cross-walks', 'Table 7 cross-walks')
src = src.replace('Table 9 summarizes', 'Table 8 summarizes')
src = src.replace('None of these five errors produces', 'None of these six counterfactual comparisons produces')
src = src.replace('and marks definitions that could not be verified from public data. Analyses were not pooled across cohorts.', 'and marks definitions that could not be verified from public data (Table 7). Analyses were not pooled across cohorts.')
src = src.replace('the erroneous robust panel would have reported a UCSF best-map agreement of 0.860 where the corrected panel gives 0.795, so the coding error would have flattered the method by 6.5 points.', 'the erroneous robust panel would have reported a UCSF best-map agreement of 0.860 where the corrected panel gives 0.795 (Table 8; Supplementary Fig. 3; Supplementary Data 1), so the coding error would have flattered the method by 6.5 points.')
src = src.replace('Three evidential levels are distinguished throughout and never substituted for one another:', 'Three evidential levels are distinguished throughout and never substituted for one another (Fig. 1b):')
src = src.replace('Table 7 cross-walks endpoint definition, time unit, radiotherapy status, covariate set, array type and methylation input for each cohort, and marks definitions that could not be verified from public data (Table 7).', 'Table 7 cross-walks cohort-specific endpoint, follow-up, radiotherapy, covariate, array and methylation-input definitions, including items that could not be verified from public data.')
src = src.replace('**Figure 4. Incremental prognostic performance and its uncertainty.** (a) Donor-grouped mean cross-validated ΔC for all nine comparisons with apparent bootstrap intervals, in blue where the interval excludes zero. (b) The four estimands for the four comparisons in which all are available, displayed on one ΔC axis and labelled by value rank, so that stability, sampling and optimism are not read as significance. (c) Cohort-specific unpenalized hazard ratios with 95% CIs; no pooled estimate.', '**Figure 4. Incremental prognostic performance and its uncertainty.** (a) Donor-grouped mean cross-validated ΔC for all nine comparisons, with apparent bootstrap intervals shown separately from the primary estimand. (b) The four estimands for the four comparisons in which all are available, displayed on one ΔC axis and labelled by value rank. (c) Cohort-specific unpenalized hazard ratios with 95% CIs; no pooled estimate.')

# Repair the estimand paragraph after the earlier prose reduction and retain its uncertainty caveat.
src = re.sub(
    r'Four quantities are reported and are not interchangeable\..*?no resample was discarded silently\.',
    'Four quantities are reported and are not interchangeable. (i) The donor-grouped mean cross-validated ΔC is the primary effect estimate. (ii) The apparent paired bootstrap percentile interval over 1,000 donor-cluster resamples describes the sampling behaviour of an in-sample ΔC estimator. (iii) The optimism-corrected ΔC is a separate internal-validation estimate. (iv) The donor-level bootstrap interval corresponds to the primary estimand: in each of 1,000 donor-cluster resamples the complete repeated cross-validation was re-executed and the 2.5th–97.5th percentiles are reported. These percentile intervals are not assumed to have exactly nominal 95% coverage because event counts vary and collapsed resamples are excluded rather than imputed. The Baylor interval is particularly exploratory (803/1,000 usable draws; 196 failures from loss of comparable pairs and one Cox non-convergence), whereas Tübingen had 996/1,000 usable draws. All failures and success rates are reported; no resample was discarded silently.',
    src,
    flags=re.S,
    count=1
)

legacy_supp = '''**Supplementary Table S1.** All fourteen originally evaluated comparisons at record level, with the split-to-split 2.5th–97.5th percentiles and positive-repeat counts, explicitly labelled as algorithmic stability rather than sampling intervals.
**Supplementary Table S2.** Within-grade increments in the external cohorts.
**Supplementary Data 1.** Longitudinal specimen-to-cluster mapping with donor-fingerprint conflict flags and the constant-class baseline correction.
**Supplementary Table S5.** UCSF/Toronto per-sample identity chain with the positional assumption disclosed and its rejection documented.
**Supplementary Table S6.** Corrected robust-panel probe list and per-target shared-probe counts.
**Supplementary Table S1.** TRIPOD+AI, TRIPOD, REMARK and field-standards reporting checklist, with items not applicable marked and reasoned.'''
src = src.replace(legacy_supp, '')
for stale in (
    '**Supplementary Table S1.** All fourteen originally evaluated comparisons at record level, with the split-to-split 2.5th–97.5th percentiles and positive-repeat counts, explicitly labelled as algorithmic stability rather than sampling intervals.',
    '**Supplementary Table S2.** Within-grade increments in the external cohorts.',
    '**Supplementary Data 1.** Longitudinal specimen-to-cluster mapping with donor-fingerprint conflict flags and the constant-class baseline correction.',
    '**Supplementary Table S5.** UCSF/Toronto per-sample identity chain with the positional assumption disclosed and its rejection documented.',
    '**Supplementary Table S6.** Corrected robust-panel probe list and per-target shared-probe counts.',
    '**Supplementary Table S1.** TRIPOD+AI, TRIPOD, REMARK and field-standards reporting checklist, with items not applicable marked and reasoned.'
):
    src = src.replace(stale, '')

# Keep all figure legends in the manuscript, with repaired title syntax and packaged relative paths.
for i in range(1,7):
    src = src.replace(f'figures_v5/fig{i}_', f'Figure{i}.')
src = src.replace('**Figure 2.** *(Figure2.independence_agreement.pdf)* Frozen-classifier reproducibility, independence and cross-cohort agreement.**', '**Figure 2. Frozen-classifier reproducibility, independence and cross-cohort agreement.**')
src = src.replace('**Figure 3.** *(Figure3.risk_coverage.pdf)* Risk–coverage behaviour of the confidence threshold.**', '**Figure 3. Risk–coverage behaviour of the confidence threshold.**')
src = src.replace('**Figure 4.** *(Figure4.incremental_value.pdf)* Incremental prognostic performance and its uncertainty.**', '**Figure 4. Incremental prognostic performance and its uncertainty.**')
src = src.replace('**Figure 5.** *(Figure5.calibration.pdf)* Calibration, absolute-risk error and follow-up sufficiency.**', '**Figure 5. Calibration, absolute-risk error and follow-up sufficiency.**')
src = src.replace('**Figure 6.** *(Figure6.deployment_audit.pdf)* **Validation audit and counterfactual evidence.**', '**Figure 6. Validation audit and counterfactual evidence.**')
src = src.replace('**Figure 1.** *(Figure1.cohort_architecture.pdf)* **Cohort assembly, sample accounting and validation architecture.**', '**Figure 1. Cohort assembly, sample accounting and validation architecture.**')

# Make the title-page marks explicit and remove stale path-like figure references in captions.
src = src.replace('## Figure legends\n', '## Figure legends\n')
src = src.replace('\n---\n\n## Abstract', '\n\n## Abstract', 1)
src = src.rstrip() + '\n\n' + table_block + '\n'

# Add normalized figure files and the CSV/GZ source data as a coherent package.
for i in range(1,7):
    matches = list((ROOT / 'figures_v5').glob(f'fig{i}_*.pdf'))
    if len(matches) != 1: raise RuntimeError(f'Expected one Figure {i} PDF; got {matches}')
    shutil.copy2(matches[0], OUT / f'Figure{i}.pdf')
    png = matches[0].with_suffix('.png')
    if png.exists(): shutil.copy2(png, OUT / f'Figure{i}.png')
supp_map = {
    1: 'figS2_fingerprint_audit.pdf', 2: 'figS4_decision_curves.pdf',
    3: 'figS5_robust_panel.pdf', 4: 'figS6_longitudinal.pdf'
}
for new, oldname in supp_map.items():
    shutil.copy2(ROOT / 'figures_v5' / oldname, OUT / f'Supplementary_Figure_S{new}.pdf')
    png = ROOT / 'figures_v5' / oldname.replace('.pdf','.png')
    if png.exists(): shutil.copy2(png, OUT / f'Supplementary_Figure_S{new}.png')
shutil.copy2(ROOT / 'strengthening_20261019' / 'results' / 'official_panel_ucsf_correctedP3.csv.gz', OUT / 'Supplementary_Data_1_corrected_panel_results.csv.gz')
shutil.copy2(ROOT / 'results' / 'supptable_s3_longitudinal_mapping.csv', OUT / 'Supplementary_Data_2_longitudinal_mapping.csv')
shutil.copy2(ROOT / 'identity_preprocessing_20261019' / 'results' / 'ucsf_ID_proven_molecular_linkage_200.csv', OUT / 'Supplementary_Data_3_ucsf_identity_linkage.csv')
(OUT / 'manuscript_npjDM_submission_candidate_v6.md').write_text(src, encoding='utf-8')
print(f'Wrote {OUT / "manuscript_npjDM_submission_candidate_v6.md"}')
print(f'Restored tables: {[n for n in (1,2,3,4,6,7,8,9)]}')
print(f'Manuscript characters: {len(src):,}; lines: {len(src.splitlines()):,}')
