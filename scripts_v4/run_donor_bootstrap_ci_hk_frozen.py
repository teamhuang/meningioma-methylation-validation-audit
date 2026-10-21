"""V5.1 addition: donor-level bootstrap of the WHOLE repeated-CV procedure for the
Hong Kong FROZEN-KNN route (the single most clinically relevant / deployable-classifier
comparison, Table 4 delta +0.059). Reuses run_donor_bootstrap_ci.py verbatim:
same SEED (identical donor resampling stream as the HK-published route), same
StratifiedGroupKFold folds, same Cox fit/concord estimators, same molecular
covariates (frozen predictions author_2/author_3 vs Merlin-intact) already present
in the input file (frozen => not refit per fold). Originals are NOT modified."""
from pathlib import Path
import sys
BASE = Path(r'E:\脑膜瘤npj')
IDENT = BASE / 'identity_preprocessing_20261019'
sys.path.insert(0, str(BASE / 'scripts_v4'))
import run_donor_bootstrap_ci as mod

# inject the frozen-KNN input: clinical + fixed frozen-class prediction covariates
mod.INPUTS['hku_author_KNN_raw'] = (
    IDENT / 'results/donor_aware/input_hku_author_KNN_raw_overlap_guard.csv',
    ['grade', 'gtr', 'age', 'male'],      # clinical
    ['author_2', 'author_3'],             # frozen-classifier predicted-class (vs Merlin-intact)
)
print('Launching donor-level bootstrap for hku_author_KNN_raw ...', flush=True)
mod.run('hku_author_KNN_raw')
