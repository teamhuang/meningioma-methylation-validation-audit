from pathlib import Path
p=Path(r'E:\脑膜瘤npj\npjDM_submission_candidate_v6\manuscript_npjDM_submission_candidate_v6.md')
s=p.read_text(encoding='utf-8')
s=s.replace('reported in Supplementary Figure S4 and Supplementary Data 1.', 'reported in Supplementary Figure S4 and Supplementary Data 2.')
s=s.replace('provided as Supplementary Data 3. The comparisons illustrate', 'provided in Supplementary Data 1. The comparisons illustrate')
s=s.replace('net-benefit curves largely indistinguishable over the examined range | **Cannot assess** — exploratory quantities reported in Supplementary Figure S4 only', 'net-benefit curves largely indistinguishable over the examined range | **Cannot assess** — exploratory quantities reported in Supplementary Figure S2 only')
s=s.replace('Sample accounting is itemised in the Supplementary ledger, separating the 866 external specimens, three reference anchors and one pilot specimen.', 'Sample accounting separating the 866 external specimens, three reference anchors and one pilot specimen is reported in Table 1 and Fig. 1a.')
old="For UCSF, the identifier-proven chain runs from the clinical table's source row, through a molecular identifier, to a specimen and its genotype, and the chain was verified to be field-exact against the source clinical table before use. Records that the chain could not attribute to a single specimen — one claimed by two molecular identifiers and two carrying none — were excluded from the prognostic analysis and listed in an archived exclusion ledger."
new="For UCSF, the identifier-proven chain runs from the clinical table's source row, through a molecular identifier, to a specimen and its genotype, and the chain was verified to be field-exact against the source clinical table before use. The profile-level linkage is provided in Supplementary Data 3. Records that the chain could not attribute to a single specimen — one claimed by two molecular identifiers and two carrying none — were excluded from the prognostic analysis and listed in an archived exclusion ledger."
s=s.replace(old,new)
p.write_text(s,encoding='utf-8')
print('Finalized text cross-references.')
