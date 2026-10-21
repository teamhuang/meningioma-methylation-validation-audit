"""Assemble a code-only GitHub payload from the meningioma npjDM project.
Copies every .py/.R analysis script, preserving relative structure, EXCLUDING raw
data, results, figures, manuscripts, the deposition copies, caches and run logs.
Nothing outside the new payload folder is modified."""
import os, shutil
from pathlib import Path

ROOT = Path(r'E:\脑膜瘤npj')
OUT = ROOT / 'github_repo'
EXCLUDE_DIR_PARTS = {'.git', 'data', 'deposition', '__pycache__', '.openscience',
                     '.wisp', 'node_modules', 'runs', '.venv', 'venv'}
# third-party / vendored / duplicated code that is NOT this study's own analysis code
EXCLUDE_SUBSTR = ('R_library', 'snapshot_reproduction_check', 'official_model',
                  'meningiomas_Analysis', 'minfi')

if OUT.exists():
    shutil.rmtree(OUT)   # OUT is a freshly-generated payload dir we own
OUT.mkdir(parents=True)

copied = []
total = 0
for p in ROOT.rglob('*'):
    if not p.is_file():
        continue
    rel = p.relative_to(ROOT)
    if any(part in EXCLUDE_DIR_PARTS for part in rel.parts):
        continue
    if any(s in str(rel) for s in EXCLUDE_SUBSTR):
        continue
    if p.suffix.lower() not in ('.py', '.r'):
        continue
    dst = OUT / rel
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(p, dst)
    copied.append(str(rel))
    total += p.stat().st_size

# also capture the classifier model config / image manifest that the code depends on
# (small JSON descriptors, not data) — optional provenance for reviewers
print(f'code files copied: {len(copied)}   total: {total/1e6:.2f} MB')
for c in sorted(copied):
    print('  ', c)
