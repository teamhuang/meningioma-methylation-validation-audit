from pathlib import Path
import shutil
root=Path(r'E:\脑膜瘤npj')
src=root/'figures_v5'
out=root/'npjDM_submission_candidate_v6'
for i in range(1,7):
    p=next(src.glob(f'fig{i}_*.pdf'))
    shutil.copy2(p,out/f'Figure{i}.pdf')
    shutil.copy2(p.with_suffix('.png'),out/f'Figure{i}.png')
for new,old in {1:'figS2_fingerprint_audit',2:'figS4_decision_curves',3:'figS5_robust_panel',4:'figS6_longitudinal'}.items():
    shutil.copy2(src/f'{old}.pdf',out/f'Supplementary_Figure_S{new}.pdf')
    shutil.copy2(src/f'{old}.png',out/f'Supplementary_Figure_S{new}.png')
print('Copied regenerated figures into the submission package.')
