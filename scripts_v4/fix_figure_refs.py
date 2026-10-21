"""Repair the figure-legend file references in manuscript_npjDM_v5.md."""
from pathlib import Path

p = Path(r'E:\脑膜瘤npj\manuscript_npjDM_v5.md')
m = p.read_text(encoding='utf-8')

paths = {
    1: 'fig1_cohort_architecture',
    2: 'fig2_independence_agreement',
    3: 'fig3_risk_coverage',
    4: 'fig4_incremental_value',
    5: 'fig5_calibration',
    6: 'fig6_deployment_audit',
}
fixed = 0
for n, fn in paths.items():
    bad = '**Figure {}.** *($(System.Collections.Hashtable[Figure {}]))*'.format(n, n)
    good = '**Figure {}.** *(figures_v5/{}.pdf)*'.format(n, fn)
    if bad in m:
        m = m.replace(bad, good)
        fixed += 1
p.write_text(m, encoding='utf-8')
print('repaired', fixed, 'figure references')
