from pathlib import Path
import re

path = Path(__file__).resolve().parent.parent / 'manuscript_npjDM_v5_2_revised.md'
s = path.read_text(encoding='utf-8-sig')

s = s.replace('Jianhuang Huang¹#, Yangyang Lin²#, Qixiu Wang³#, Jianning Chen¹*',
              'Jianhuang Huang¹#*, Yangyang Lin²#, Qixiu Wang³#, Jianning Chen¹*')
s = s.replace('# These authors contributed equally.', '# Jianhuang Huang, Yangyang Lin and Qixiu Wang contributed equally.')

# Add primary citations in order of first appearance and cite the classifier under review.
s = s.replace('and recurrence-oriented classifiers developed in single-centre series.',
              'and recurrence-oriented classifiers developed in single-centre series [1–4].')
s = s.replace('Reporting guidance for healthcare AI emphasizes independent evaluation and transparent validation methods [1]',
              'Reporting guidance for healthcare AI emphasizes independent evaluation and transparent validation methods [5]')
s = s.replace('The frozen classifier is the public 403-CpG, k = 10 KNN model with a 141-case balanced reference,',
              'The frozen classifier is the public 403-CpG, k = 10 KNN model with a 141-case balanced reference described by Duckett et al. [4],')
s = s.replace('TRIPOD+AI [2], the original TRIPOD statement [3], REMARK [4] and the field\'s guidance for digital biomarker manuscripts [1], with extraction of cohort inputs against CHARMS [5]',
              'TRIPOD+AI [6], the original TRIPOD statement [7], REMARK [8] and the field\'s guidance for digital biomarker manuscripts [5], with extraction of cohort inputs against CHARMS [9]')

# Remove an exploratory association block that still uses the pre-identity-resolution UCSF denominator.
s = re.sub(r'Cohort-specific unpenalized Cox models gave hazard ratios.*?These are exploratory association tests on the estimation records and do not test out-of-sample ΔC\.\n\n', '', s, count=1, flags=re.S)
s = re.sub(r'### Association testing and diagnostics\n.*?(?=### Calibration, absolute-risk error and decision analysis)', '', s, count=1, flags=re.S)
s = re.sub(r'### Table 5\. Exploratory global association tests, sixteen-test family\n.*?(?=### Table 6\.)', '', s, count=1, flags=re.S)
s = s.replace('Across sixteen unpenalized likelihood-ratio tests, six fall below 0.05 after Holm correction (Table 5);',
              'Across sixteen unpenalized likelihood-ratio tests, six fell below 0.05 after Holm correction in the earlier record-level audit;')

# Keep publication-facing references limited to cited works and order them by first citation.
refs = '''## References

1. Bayley, J. C. T. et al. Multiple approaches converge on three biological subtypes of meningioma and extract new insights from published studies. *Sci. Adv.* **8**, eabm6247 (2022). https://doi.org/10.1126/sciadv.abm6247
2. Choudhury, A. et al. Meningioma DNA methylation groups identify biological drivers and therapeutic vulnerabilities. *Nat. Genet.* **54**, 649–659 (2022). https://doi.org/10.1038/s41588-022-01061-8
3. Wang, J. Z. et al. Molecular classification to refine surgical and radiotherapeutic decision-making in meningioma. *Nat. Med.* **30**, 3173–3183 (2024). https://doi.org/10.1038/s41591-024-03167-4
4. Duckett, D. et al. Development, validation, and utility of a clinically applicable methylation classifier for recurrence risk prediction in meningiomas. *Acta Neuropathol. Commun.* **13**, 259 (2025). https://doi.org/10.1186/s40478-025-02184-9
5. Kakarmath, S., Esteva, A., Arnaout, R. et al. Best practices for authors of healthcare-related artificial intelligence manuscripts. *npj Digit. Med.* **3**, 47 (2020). https://doi.org/10.1038/s41746-020-00336-w
6. Collins, G. S., Moons, K. G. M., Dhiman, P. et al. TRIPOD+AI statement: updated guidance for reporting clinical prediction models that use regression or machine learning methods. *BMJ* **385**, e078378 (2024).
7. Collins, G. S., Reitsma, J. B., Altman, D. G. & Moons, K. G. M. Transparent reporting of a multivariable prediction model for individual prognosis or diagnosis (TRIPOD): the TRIPOD statement. *BMJ* **350**, g7594 (2015). https://doi.org/10.1136/bmj.g7594
8. McShane, L. M., Altman, D. G., Sauerbrei, W., Taube, S. E., Gion, M. & Clark, G. M. Reporting recommendations for tumor marker prognostic studies (REMARK). *J. Natl. Cancer Inst.* **97**, 1180–1184 (2005). https://doi.org/10.1093/jnci/dji237
9. Moons, K. G. M., de Groot, J. G. H., Bouwmeester, W. et al. Critical appraisal and data extraction for systematic reviews of prediction modelling studies: the CHARMS checklist. *PLoS Med.* **11**, e1001744 (2014). https://doi.org/10.1371/journal.pmed.1001744

'''
s = re.sub(r'## References\n.*?(?=\n---\n\n## Figure legends)', refs, s, count=1, flags=re.S)

# Internal numerical file paths and historical analysis notes are not publication-facing supplementary content.
s = re.sub(r'### Supplementary Table S8\. Numerical ledger: every reported value mapped to its source file\n.*$', '', s, count=1, flags=re.S)
s = s.replace('**Supplementary Table S4.** Unpenalized association tests with raw P, Holm-adjusted P and df for the five-core, fourteen- and sixteen-comparison families.\n', '')
s = s.replace('**Figure 6.** *(figures_v5/fig6_deployment_audit.pdf)* Deployment-readiness audit and counterfactual evidence.** (a) The five dimensions with audit-defined criteria and three-level verdicts. (b) The five counterfactuals of Table 6 plotted as pairs, showing the reportable number a conventional pipeline would have produced against the audited value.',
              '**Figure 6.** *(figures_v5/fig6_deployment_audit.pdf)* **Validation audit and counterfactual evidence.** (a) Five audit dimensions with study-defined criteria and three-level verdicts. (b) Six counterfactual comparisons across five audit domains, showing the value a conventional pipeline would have reported against the audited result.')
s = s.replace('Deployment-readiness', 'Validation')
s = s.replace('preconditions for deployment', 'checks for interpreting classifier validation')
s = s.replace('confidence filtering therefore raises reported agreement while preferentially discarding the highest-risk class',
              'confidence filtering therefore raises reported agreement while preferentially discarding hypermitotic cases')
s = s.replace('will therefore degrade performance exactly where it matters clinically',
              'may reduce coverage in a potentially higher-risk subgroup')

path.write_text(s, encoding='utf-8-sig')
print('Finalized', path, 'words=', len(s.split()))
