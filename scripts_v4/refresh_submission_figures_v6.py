import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import make_figures_v5 as f1f2
import make_figures_v5_1_fig3 as f3
import make_figures_v5_1_fig4 as f4
import make_figures_v5_fix as f6
import make_figures_v5_supp as supp

f1f2.figure1()
f1f2.figure2()
f3.figure3()
f4.figure4()
f6.figure6()
supp.figure5()
supp.figure_s2()
supp.figure_s4()
supp.figure_s5()
supp.figure_s6()
print('Submission figures refreshed.')
