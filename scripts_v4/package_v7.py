"""Build the v7 submission bundle from v6, swapping in the v7 manuscript + docx and
adding Supplementary Data 6. Regenerates SHA256SUMS for the bundle."""
import os, shutil, hashlib, subprocess
from pathlib import Path
ROOT = Path(r'E:\脑膜瘤npj')
SRC = ROOT / 'npjDM_submission_candidate_v6'
DST = ROOT / 'npjDM_submission_candidate_v7'

def rb(p):
    subprocess.run(['powershell.exe','-NoProfile','-Command',
      "Add-Type -AssemblyName Microsoft.VisualBasic; [Microsoft.VisualBasic.FileIO.FileSystem]::DeleteFile("
      f"'{p}','OnlyErrorDialogs','SendToRecycleBin')"], check=False)

shutil.copytree(SRC, DST, dirs_exist_ok=True)
# retire old v6 manuscript/docx inside the v7 bundle
for old in [DST/'manuscript_npjDM_submission_candidate_v6.md', DST/'manuscript_npjDM_submission_candidate_v6.docx']:
    if old.exists(): rb(old)
# bring in v7 md + its docx (docx built separately first)
shutil.copy2(ROOT/'manuscript_npjDM_v7.md', DST/'manuscript_npjDM_submission_candidate_v7.md')
v7docx = ROOT/'manuscript_npjDM_v7.docx'
if v7docx.exists(): shutil.copy2(v7docx, DST/'manuscript_npjDM_submission_candidate_v7.docx')
# ensure Supp Data 6 present (build_v7 wrote it into v6 folder; copy forward)
sd6 = ROOT/'npjDM_submission_candidate_v6/Supplementary_Data_6_rejection_outcome.csv'
if sd6.exists(): shutil.copy2(sd6, DST/'Supplementary_Data_6_rejection_outcome.csv')

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda: f.read(1<<20), b''): h.update(b)
    return h.hexdigest()
lines=[f'{sha(p)}  {p.relative_to(DST).as_posix()}' for p in sorted(DST.rglob('*')) if p.is_file() and p.name!='SHA256SUMS.txt']
(DST/'SHA256SUMS.txt').write_text('\n'.join(lines)+'\n', encoding='utf-8')
print('v7 bundle:', DST, '| files:', sum(1 for p in DST.rglob("*") if p.is_file()))
print('has v7 md:', (DST/'manuscript_npjDM_submission_candidate_v7.md').exists(),
      '| has v7 docx:', (DST/'manuscript_npjDM_submission_candidate_v7.docx').exists(),
      '| has SuppData6:', (DST/'Supplementary_Data_6_rejection_outcome.csv').exists())
