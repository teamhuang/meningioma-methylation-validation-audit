"""Build the V5.1 Zenodo-ready deposition snapshot from the V5 one, refreshing every
changed artefact and regenerating SHA256SUMS. Deterministic; never deletes the V5 tree."""
import os, shutil, hashlib, subprocess
from pathlib import Path

ROOT = Path(r'E:\脑膜瘤npj')
SRC = ROOT / 'deposition/npjDM_meningioma_methylation_audit_v5'
DST = ROOT / 'deposition/npjDM_meningioma_methylation_audit_v5_1'

def to_recycle_bin(path: Path):
    """Never permanently delete under the user tree; send to Windows Recycle Bin."""
    ps = ("Add-Type -AssemblyName Microsoft.VisualBasic; "
          "[Microsoft.VisualBasic.FileIO.FileSystem]::DeleteFile("
          f"'{path}','OnlyErrorDialogs','SendToRecycleBin')")
    subprocess.run(['powershell.exe', '-NoProfile', '-Command', ps], check=False)

shutil.copytree(SRC, DST, dirs_exist_ok=True)

# 1. manuscript: retire v5 copy (Recycle Bin), add v5.1 md + docx
old_ms = DST / 'manuscript_npjDM_v5.md'
if old_ms.exists(): to_recycle_bin(old_ms)
shutil.copy2(ROOT / 'manuscript_npjDM_v5_1.md', DST / 'manuscript_npjDM_v5_1.md')
shutil.copy2(ROOT / 'manuscript_npjDM_v5_1.docx', DST / 'manuscript_npjDM_v5_1.docx')

# 2. figures: refresh the regenerated ones (fig3, fig4) + keep the rest
for f in ['fig3_risk_coverage', 'fig4_incremental_value']:
    for ext in ('pdf', 'png'):
        src = ROOT / f'figures_v5/{f}.{ext}'
        if src.exists():
            shutil.copy2(src, DST / f'figures_v5/{f}.{ext}')

# 3. supplementary text artefacts
shutil.copy2(ROOT / 'results_v4/supplementary_note_1_revision_record.md',
             DST / 'supplementary_note_1_revision_record.md')
for ext in ('csv', 'md'):
    shutil.copy2(ROOT / f'results_v4/supplementary_table_s8_numerical_ledger.{ext}',
                 DST / f'supplementary_table_s8_numerical_ledger.{ext}')

# 4. add the new HK-frozen bootstrap outputs (results_v4 subdir in the package)
(rv4 := DST / 'results_v4').mkdir(exist_ok=True)
for name in ['donor_bootstrap_ci_hku_author_KNN_raw.csv',
             'donor_bootstrap_ci_summary_hku_author_KNN_raw.csv']:
    src = ROOT / f'results_v4/{name}'
    if src.exists(): shutil.copy2(src, rv4 / name)

# 5. add the V5.1 scripts that produced the new artefacts
for name in ['run_donor_bootstrap_ci_hk_frozen.py', 'make_figures_v5_1_fig3.py',
             'make_figures_v5_1_fig4.py', 'build_v5_1.py', 'patch_v5_1_hk_frozen.py',
             'audit_consistency_v5_1.py']:
    src = ROOT / f'scripts_v4/{name}'
    if src.exists(): shutil.copy2(src, DST / f'scripts_v4/{name}')

# 6. README version bump + V5.1 note
rd = DST / 'README.md'
if rd.exists():
    t = rd.read_text(encoding='utf-8')
    t = t.replace('Version-5 scripts', 'Version-5.1 scripts')
    note = ('\n\n## Version 5.1 (submission cleanup)\n\n'
            'V5.1 removes the residual UCSF record-level/level-A contradiction, retracts the '
            'claim that the UCSF record-level interval excluded zero (lower bound -0.000002, so it '
            'includes zero; Hong Kong only), reports risk-coverage over all true class cases with '
            'tie- vs confidence-rejection separated, adds the whole-procedure donor-level interval '
            'for the Hong Kong frozen-classifier route '
            '(+0.008 to +0.132, 1000/1000 draws; run_donor_bootstrap_ci_hk_frozen.py), corrects the '
            'Baylor bootstrap-failure mechanism (196/197 = loss of comparable pairs), and re-words '
            'two counterfactual claims. The audit gate is scripts_v4/audit_consistency_v5_1.py.\n')
    if 'Version 5.1 (submission cleanup)' not in t:
        t += note
    rd.write_text(t, encoding='utf-8')

# 7. regenerate SHA256SUMS over every file in the package (relative paths, forward slashes)
def sha256(p):
    h = hashlib.sha256()
    with open(p, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 20), b''): h.update(b)
    return h.hexdigest()

lines = []
for p in sorted(DST.rglob('*')):
    if p.is_file() and p.name != 'SHA256SUMS.txt':
        rel = p.relative_to(DST).as_posix()
        lines.append(f'{sha256(p)}  {rel}')
(DST / 'SHA256SUMS.txt').write_text('\n'.join(lines) + '\n', encoding='utf-8')

nfiles = sum(1 for p in DST.rglob('*') if p.is_file())
print(f'V5.1 package: {DST}')
print(f'  files: {nfiles}   sha256 lines: {len(lines)}')
print(f'  manuscript present: {(DST/"manuscript_npjDM_v5_1.md").exists()} '
      f'docx: {(DST/"manuscript_npjDM_v5_1.docx").exists()} '
      f'hk_frozen boot: {(rv4/"donor_bootstrap_ci_summary_hku_author_KNN_raw.csv").exists()}')
