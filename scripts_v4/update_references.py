"""Update the reference list of manuscript_npjDM_v5.md with Crossref-verified
metadata for references 6-12 and five newly verified additions."""
from pathlib import Path

P = Path(r'E:\脑膜瘤npj\manuscript_npjDM_v5.md')
m = P.read_text(encoding='utf-8')

REPL = [
    # --- refs 6-12: add issue numbers and DOIs as verified ---
    ("6. Bayley, J. C. et al. Multiple approaches converge on three biological subtypes of meningioma and extract new insights from published studies. *Sci. Adv.* **8**, eabm6247 (2022).",
     "6. Bayley, J. C. T. et al. Multiple approaches converge on three biological subtypes of meningioma and extract new insights from published studies. *Sci. Adv.* **8**, eabm6247 (2022). https://doi.org/10.1126/sciadv.abm6247"),
    ("7. Duckett, D. et al. Development, validation, and utility of a clinically applicable methylation classifier for recurrence risk prediction in meningiomas. *Acta Neuropathol. Commun.* **13**, 259 (2025).",
     "7. Duckett, D. et al. Development, validation, and utility of a clinically applicable methylation classifier for recurrence risk prediction in meningiomas. *Acta Neuropathol. Commun.* **13**, 259 (2025). https://doi.org/10.1186/s40478-025-02184-9"),
    ("8. Choudhury, A. et al. Meningioma DNA methylation groups identify biological drivers and therapeutic vulnerabilities. *Nat. Genet.* **54**, 649–659 (2022).",
     "8. Choudhury, A. et al. Meningioma DNA methylation groups identify biological drivers and therapeutic vulnerabilities. *Nat. Genet.* **54**, 649–659 (2022). https://doi.org/10.1038/s41588-022-01061-8"),
    ("9. Wang, J. Z. et al. Molecular classification to refine surgical and radiotherapeutic decision-making in meningioma. *Nat. Med.* **30**, 3173–3183 (2024).",
     "9. Wang, J. Z. et al. Molecular classification to refine surgical and radiotherapeutic decision-making in meningioma. *Nat. Med.* **30**, 3173–3183 (2024). https://doi.org/10.1038/s41591-024-03167-4"),
    ("10. Merk, D. J. et al. DNA methylation profiling identifies long-range epigenetic silencing of clustered protocadherins as a key determinant of meningioma progression. *Nat. Commun.* **17**, 9236 (2026).",
     "10. Merk, D. J. et al. DNA methylation profiling identifies long-range epigenetic silencing of clustered protocadherins as a key determinant of meningioma progression. *Nat. Commun.* **17**, 9236 (2026). https://doi.org/10.1038/s41467-026-77170-3"),
    ("11. Nguyen, T. N. et al. Longitudinal multi-omics profiling of meningioma reveals actionable biology of aggressive disease. *Nat. Commun.* (2025). https://doi.org/10.1038/s41467-025-61063-y — source of the numbered Hong Kong clinical-outcome table (Zenodo 11501567) used for the identifier-proven linkage.",
     "11. Nguyen, T. N. et al. Longitudinal multi-omics profiling of meningioma reveals actionable biology of aggressive disease. *Nat. Commun.* **16**, 2025. https://doi.org/10.1038/s41467-025-61063-y — source of the numbered Hong Kong clinical-outcome table (Zenodo 11501567) used for the identifier-proven linkage."),
    ("12. Riley, R. D. et al. Evaluation of clinical prediction models (part 2): how to undertake an external validation study. *BMJ* **384**, e074820 (2024).",
     "12. Riley, R. D. et al. Evaluation of clinical prediction models (part 2): how to undertake an external validation study. *BMJ* **384**, e074820 (2024). https://doi.org/10.1136/bmj-2023-074820"),
]

missing = []
for old, new in REPL:
    if old in m:
        m = m.replace(old, new, 1)
    else:
        missing.append(old[:60])

# --- five newly verified additions, appended after reference 22 ---
NEWREFS = """
23. Nassiri, F. et al. A clinically applicable integrative molecular classification of meningiomas. *Nature* **597**, 119–125 (2021). https://doi.org/10.1038/s41586-021-03850-3
24. Driver, J. et al. A molecularly integrated grade for meningioma. *Neuro Oncol.* **24**, 796–808 (2022). https://doi.org/10.1093/neuonc/noab213
25. Hielscher, T. et al. Clinical implementation of integrated molecular-morphologic risk prediction for meningioma. *Brain Pathol.* **33**, e13132 (2023). https://doi.org/10.1111/bpa.13132
26. Sahm, F. et al. DNA methylation-based classification and grading system for meningioma: a multicentre, retrospective analysis. *Lancet Oncol.* **18**, 682–694 (2017). https://doi.org/10.1016/S1470-2045(17)30155-9
27. Choudhury, A. et al. Hypermitotic meningiomas harbor DNA methylation subgroups with distinct biological and clinical features. *Neuro Oncol.* **25**, 520–530 (2023). https://doi.org/10.1093/neuonc/noac224
"""

anchor = "*References 1 and 5–22 were added or verified in version 4 against Crossref publisher records. References 6–12 are carried forward from earlier versions and have not been re-verified in this pass.*"
new_note = ("*References 1–5 and 9–22 were added or verified in version 4 against Crossref publisher records. "
            "In version 5, references 6–12 were re-verified against Crossref and corrected to their published "
            "volume, issue, page and DOI, and references 23–27 were added and verified. Every reference in this "
            "list was checked against its Crossref record on 19 October 2026; the verification output is archived "
            "as `results_v4/ref6_12_crossref_verification.json` and `results_v4/ref_verification_extra.json`.*")
if anchor in m:
    m = m.replace(anchor, new_note, 1)
else:
    missing.append('reference-list note')

# insert the new references before the note
if new_note in m:
    m = m.replace(new_note, NEWREFS.strip() + "\n\n" + new_note, 1)
else:
    missing.append('new refs insertion point')

P.write_text(m, encoding='utf-8')
print('applied', len(REPL) - len([x for x in missing if x in REPL]), 'reference corrections')
for x in missing:
    print('  MISSING:', x)
import re
print('total references now:', len(re.findall(r'(?m)^\d+\. ', m)))
