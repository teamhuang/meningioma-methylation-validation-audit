"""Replace the V5 abstract with the version reflecting the donor-level bootstrap,
Uno's C and out-of-fold calibration results, and check the npj DM word limit."""
from pathlib import Path
import re

P = Path(r'E:\脑膜瘤npj\manuscript_npjDM_v5.md')
m = P.read_text(encoding='utf-8')

NEW = (
    "DNA methylation classifiers improve meningioma risk stratification, but "
    "cross-institutional reliability is unmeasured. We audited a frozen 403-CpG classifier "
    "and published labels across four cohorts, examining sample identity, "
    "training overlap, probe coverage and performance. Numerical fingerprinting found no "
    "overlap between 866 external specimens and its 141-case training reference, whereas "
    "the nominally held-out release contained 20 of 37 training profiles. The frozen "
    "classifier agreed with independently published group assignments in 91.5% of "
    "Tübingen, 87.2% of UCSF and 88.8% of Hong Kong assigned specimens, but agreement "
    "fell to 26.5% under positional joining and to 41.0% for grade. Raising agreement to "
    "93–96% at the 0.75 confidence threshold discarded 39% of hypermitotic versus 5% of "
    "Merlin-intact Hong Kong cases. Molecular labels increased cohort-internal "
    "discrimination, corroborated by Uno's C and time-dependent AUC, but donor-level "
    "bootstrap intervals excluded zero only in Hong Kong; absolute-risk gains were "
    "inconsistent. A five-dimension audit identified errors conventional validation would "
    "not detect."
)

pat = re.compile(r'(## Abstract\s*\n\n)([\s\S]*?)(\n---)', re.M)
assert pat.search(m), 'abstract block not found'
m = pat.sub(lambda g: g.group(1) + NEW + g.group(3), m, count=1)
P.write_text(m, encoding='utf-8')

n = len(NEW.split())
title_words = len(m.splitlines()[0].lstrip('# ').split())
print(f'abstract words : {n}  (npj DM limit 150) -> {"OK" if n <= 150 else "TOO LONG"}')
print(f'title words    : {title_words}  (npj DM limit 15) -> {"OK" if title_words <= 15 else "TOO LONG"}')
print(f'subheadings    : {"**" in NEW}')
