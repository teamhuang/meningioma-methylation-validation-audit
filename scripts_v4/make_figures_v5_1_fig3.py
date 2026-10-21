"""V5.1 figure updates: make Figure 3 reflect the corrected DUAL-DENOMINATOR
risk-coverage message (retention over ALL true class cases, with tie-rejection shown
separately from confidence-rejection), instead of the old uniquely-assigned-only bars.
Reuses the v5 house style / save() / load_agreement() from make_figures_v5.py.
Figure 1 and 2 graphics need no change (their plotted numbers are unchanged); the
Figure 1 fix was a legend-text fix already applied in the manuscript.
"""
from pathlib import Path
import warnings, sys
import numpy as np
sys.path.insert(0, r'E:\脑膜瘤npj\scripts_v4')
import make_figures_v5 as m   # applies style, exposes save(), load_agreement(), consts

SHORT, COH, CLASSES, NCOL = m.SHORT, m.COH, m.CLASSES, m.NCOL
BLUE, AQUA, YELLOW, VIOLET, RED, GREY, INK = m.BLUE, m.AQUA, m.YELLOW, m.VIOLET, m.RED, m.GREY, m.INK
SEQ = m.SEQ
import matplotlib.pyplot as plt


def figure3():
    ag = m.load_agreement()
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 3.8), gridspec_kw={'wspace': 0.42})

    # (a) risk-coverage: coverage over all eligible profiles, selective accuracy over assigned
    ax = axes[0]
    for i, c in enumerate(COH):
        a = ag[c][~ag[c].tied]
        order = np.argsort(a.confidence.to_numpy())
        ok = a.ok.to_numpy()[order].astype(float)
        cum = np.cumsum(ok) / np.arange(1, len(ok) + 1)
        cov = np.arange(1, len(ok) + 1) / NCOL[c]
        ax.plot(cov, cum, lw=2, color=[BLUE, AQUA, VIOLET][i], label=SHORT[c])
    ax.axvline(0.773, color=GREY, lw=0.9, ls=':')
    ax.text(0.785, 0.56, '0.75 confidence\nthreshold', fontsize=7.5, color=GREY, linespacing=1.3)
    ax.set_xlabel('coverage (fraction of all eligible profiles)', fontsize=8.5)
    ax.set_ylabel('selective accuracy', fontsize=8.5)
    ax.set_ylim(0.5, 1.02)
    ax.set_title('(a) Risk–coverage of the classifier confidence', loc='left', fontsize=9.5)
    ax.legend(fontsize=8, frameon=False, loc='lower right')
    ax.grid(axis='y', lw=0.4, alpha=0.3); ax.set_axisbelow(True)

    # (b) retention at 0.75 computed over ALL true cases (tied = discarded), assigned overlaid
    ax = axes[1]
    xs = np.arange(3); w = 0.26
    for i, g in enumerate(CLASSES):
        all_ret, assigned_ret = [], []
        for c in COH:
            a = ag[c]
            sub_all = a[a.Group == g]                       # all true cases incl tied
            retained = ((~sub_all.tied) & (sub_all.confidence >= 0.75)).mean()
            sub_as = sub_all[~sub_all.tied]                 # uniquely assigned only
            ret_as = (sub_as.confidence >= 0.75).mean()
            all_ret.append(retained); assigned_ret.append(ret_as)
        b = ax.bar(xs + (i - 1) * w, all_ret, width=w, color=[RED, YELLOW, BLUE][i], label=g)
        for xi, v, va in zip(xs, all_ret, assigned_ret):
            ax.text(xi + (i - 1) * w, v + 0.02, f'{v:.2f}', ha='center', fontsize=7, color=INK)
            ax.plot([xi + (i - 1) * w], [va], marker='_', ms=12, mew=1.8,
                    color=INK, alpha=0.55)  # tick = uniquely-assigned-only value
    ax.set_xticks(xs); ax.set_xticklabels([SHORT[c] for c in COH], fontsize=8.5)
    ax.set_ylabel('retention at confidence 0.75', fontsize=8.5)
    ax.set_ylim(0, 1.16)
    ax.set_title('(b) Bars = all true cases; tick = assigned only', loc='left', fontsize=9.0)
    ax.legend(fontsize=7.5, frameon=False, ncol=1, loc='upper center',
              bbox_to_anchor=(0.5, 1.03), handlelength=1.1)
    ax.grid(axis='y', lw=0.4, alpha=0.3); ax.set_axisbelow(True)

    # (c) discarded share per class over ALL cases, split tie-rejection vs confidence-rejection
    # NOTE: bars are drawn at y[i] (not i) so each cohort's data sits on its own label row.
    ax = axes[2]
    y = np.arange(3)[::-1]
    for i, c in enumerate(COH):
        a = ag[c]; pos = y[i]
        for g, sign in [('Hypermitotic', +1), ('Merlin-intact', -1)]:
            sub = a[a.Group == g]; N = len(sub)
            tie = sub.tied.astype(bool).sum()                       # lost before threshold
            conf_fail = ((~sub.tied) & (sub.confidence < 0.75)).sum()  # lost by threshold
            tie_f, conf_f = tie / N, conf_fail / N
            base = tie_f if sign > 0 else -tie_f
            span = conf_f if sign > 0 else -conf_f
            ax.barh(pos, base, color=GREY, height=0.5)              # tie loss (grey)
            ax.barh(pos, span, left=base, color=RED if sign > 0 else SEQ[1], height=0.5)
            tot = base + span
            # short label at the outer end; tie/conf split is in the legend
            ax.text(tot + 0.015 * sign, pos, f'{abs(tot):.0%}', va='center',
                    ha='left' if sign > 0 else 'right', fontsize=7.5, color=INK)
    ax.axvline(0, color=INK, lw=1)
    ax.set_yticks(y); ax.set_yticklabels([SHORT[c] for c in COH], fontsize=8.5)
    ax.set_xlim(-0.55, 0.62)
    ax.set_xlabel('Merlin-intact discarded  ←  |  →  Hypermitotic discarded', fontsize=7.0)
    ax.set_title('(c) Tie- + confidence-rejection, all cases', loc='left', fontsize=9.0)
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=GREY),
                       plt.Rectangle((0, 0), 1, 1, color=RED)],
              labels=['tied / unassigned', 'confidence < 0.75'],
              fontsize=7, frameon=False, loc='lower left', handlelength=1.0)
    ax.grid(axis='x', lw=0.4, alpha=0.3); ax.set_axisbelow(True)
    ax.tick_params(axis='y', pad=2)

    m.save(fig, 'fig3_risk_coverage')


if __name__ == '__main__':
    warnings.filterwarnings('ignore')
    figure3()
    print('FIG3_DONE', flush=True)
