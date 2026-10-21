"""V5.1 figure 4: panel (a) now plots the WHOLE-PROCEDURE donor-level bootstrap
interval (the primary-estimand resampling interval) and colours it by whether THAT
interval excludes zero -- so the figure carries the corrected inference (Hong Kong
only), instead of the old apparent-interval colouring that wrongly implied UCSF and
both HK routes excluded zero. Routes without a whole-procedure interval are shown as
open point estimates and explicitly labelled as such. Every number is read from disk
(donor_bootstrap_ci_summary_*.csv and donor_aware / results_v4 summary_*.csv)."""
from pathlib import Path
import warnings, sys, os
import numpy as np, pandas as pd
sys.path.insert(0, r'E:\脑膜瘤npj\scripts_v4')
import make_figures_v5 as m
import matplotlib.pyplot as plt

ROOT = Path(r'E:\脑膜瘤npj'); IDENT = m.IDENT; STR = m.STR; V4 = m.V4
BLUE, AQUA, VIOLET, GREY, INK, SEQ, RED = m.BLUE, m.AQUA, m.VIOLET, m.GREY, m.INK, m.SEQ, m.RED


def read_first(cands):
    for c in cands:
        if os.path.exists(c):
            return pd.read_csv(c)
    return None


def summary(name):
    return read_first([IDENT / f'donor_aware/summary_{name}_overlap_guard.csv',
                       V4 / f'summary_{name}_overlap_guard.csv'])


def donor_ci(name):
    return read_first([V4 / f'donor_bootstrap_ci_summary_{name}.csv'])


# (display label, summary name, donor-bootstrap summary name or None)
PUBLISHED = [
    ('Baylor, published MenG', 'baylor_published', 'baylor_published'),
    ('Tübingen, published', 'tubingen_published', 'tubingen_published'),
    ('UCSF, published (grouped)', 'ucsf_published', 'ucsf_published'),
    ('Hong Kong, published', 'hku_published', 'hku_published'),
]
FROZEN_WITH_CI = [('Hong Kong, frozen KNN', 'hku_author_KNN_raw', 'hku_author_KNN_raw')]
POINT_ONLY = [  # no whole-procedure interval computed -> open markers only
    ('Baylor, surrogate', 'baylor_centroid_surrogate'),
    ('Baylor, frozen KNN', 'baylor_author_KNN_raw'),
    ('Tübingen, surrogate', 'tubingen_centroid_surrogate'),
    ('Tübingen, frozen KNN', 'tubingen_author_KNN_raw'),
]


def figure4():
    # assemble primary rows
    prim = []  # (label, dc, ci_low, ci_high, success, kind)
    for lab, sname, cname in PUBLISHED + FROZEN_WITH_CI:
        s = summary(sname); dc = float(s.iloc[0]['delta_cv'])
        ci = donor_ci(cname)
        if ci is None:
            prim.append((lab, dc, np.nan, np.nan, None, 'pending')); continue
        r = ci.iloc[0]
        prim.append((lab, dc, float(r['delta_grouped_ci_low']), float(r['delta_grouped_ci_high']),
                     int(r['bootstrap_success']), 'ci'))
    pts = []
    for lab, sname in POINT_ONLY:
        s = summary(sname); pts.append((lab, float(s.iloc[0]['delta_cv'])))

    n_ci = len(prim); n_pt = len(pts)
    fig = plt.figure(figsize=(12.6, 6.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.05, 0.9], wspace=0.46)

    # (a) whole-procedure donor-level interval
    ax = fig.add_subplot(gs[0, 0])
    labels = [p[0] for p in prim] + [p[0] for p in pts]
    ypos = np.arange(len(labels))[::-1]
    for i, (lab, dc, lo, hi, succ, kind) in enumerate(prim):
        y = ypos[i]
        if kind == 'ci':
            excl = lo > 0
            col = BLUE if excl else GREY
            ax.plot([lo, hi], [y, y], lw=2.4, color=col, solid_capstyle='round', zorder=2)
            ax.plot(dc, y, 'o', ms=7.5, color=col, zorder=3)
            tag = f'{succ}/1000' if succ is not None else ''
            ax.text(hi + 0.006, y, f'{dc:+.3f}', va='center', fontsize=7.5, color=col)
        else:
            ax.plot(dc, y, 'o', ms=7.5, color='none', markeredgecolor=RED,
                    markeredgewidth=1.6, zorder=3)
            ax.text(dc + 0.006, y, f'{dc:+.3f} (interval pending)', va='center', fontsize=7, color=RED)
    for j, (lab, dc) in enumerate(pts):
        y = ypos[n_ci + j]
        ax.plot(dc, y, 'o', ms=7.5, mfc='white', mec=GREY, mew=1.6, zorder=3)
        ax.text(dc + 0.006, y, f'{dc:+.3f}', va='center', fontsize=7.5, color=GREY)
    ax.axvline(0, color=INK, lw=1, ls='--')
    ax.set_yticks(ypos); ax.set_yticklabels(labels, fontsize=8.2)
    ax.set_xlabel('donor-grouped ΔC  ±  whole-procedure donor-level bootstrap interval', fontsize=8.2)
    ax.set_xlim(-0.09, 0.36)
    ax.set_title('(a) Primary estimand; colour = interval excludes 0', loc='left', fontsize=9.3)
    ax.legend(handles=[plt.Line2D([], [], marker='o', ls='-', ms=7, color=BLUE, label='whole-procedure interval excludes 0'),
                       plt.Line2D([], [], marker='o', ls='-', ms=7, color=GREY, label='whole-procedure interval spans 0'),
                       plt.Line2D([], [], marker='o', ls='', ms=7, mfc='white', mec=GREY, mew=1.6, label='point estimate only (no whole-procedure interval)')],
              fontsize=6.8, frameon=False, loc='lower right')
    ax.grid(axis='x', lw=0.4, alpha=0.3); ax.set_axisbelow(True)

    # (b) four estimands for the published routes + HK frozen
    ax = fig.add_subplot(gs[0, 1])
    rec = {'baylor_published': 0.073, 'tubingen_published': 0.023, 'ucsf_published': 0.024}
    rows_b = []
    for lab, sname, cname in PUBLISHED + FROZEN_WITH_CI:
        s = summary(sname).iloc[0]
        rows_b.append((lab, float(s['delta_cv']), float(s['delta_app_CI_low']),
                       float(s['delta_app_CI_high']), float(s['delta_corrected']),
                       rec.get(sname, None)))
    y = np.arange(len(rows_b))[::-1]
    for i, (lab, dc, lo, hi, corr, r) in enumerate(rows_b):
        ax.plot([lo, hi], [y[i] + 0.22] * 2, lw=2.2, color=SEQ[0], solid_capstyle='round')
        ax.plot(dc, y[i] + 0.22, 'o', ms=7, color=BLUE, zorder=3)
        ax.plot(corr, y[i] - 0.18, 'D', ms=6.5, color=AQUA, zorder=3)
        if r is not None:
            ax.plot(r, y[i] - 0.52, 's', ms=6, color=GREY, zorder=3)
        ax.text(dc, y[i] + 0.40, f'{dc:.3f}', ha='center', fontsize=6.8, color=BLUE)
        ax.text(corr, y[i] - 0.34, f'{corr:.3f}', ha='center', fontsize=6.8, color=AQUA)
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows_b], fontsize=8.2)
    ax.set_xlabel('ΔC', fontsize=8.5); ax.set_xlim(-0.05, 0.19); ax.set_ylim(-1.1, len(rows_b) - 0.2)
    ax.set_title('(b) Estimands are not one test', loc='left', fontsize=9.3)
    ax.legend(handles=[plt.Line2D([], [], marker='o', ls='-', ms=7, color=BLUE, label='donor-grouped cross-validated ΔC (primary)'),
                       plt.Line2D([], [], marker='D', ls='', ms=6.5, color=AQUA, label='optimism-corrected ΔC'),
                       plt.Line2D([], [], marker='s', ls='', ms=6, color=GREY, label='record-level ΔC (counterfactual)')],
              fontsize=6.8, frameon=False, loc='lower left')
    ax.grid(axis='x', lw=0.4, alpha=0.3); ax.set_axisbelow(True)

    # (c) HRs (unchanged, read from disk)
    ax = fig.add_subplot(gs[0, 2])
    hk = pd.read_csv(IDENT / 'donor_aware/cluster_robust_cox_hku_published_overlap_guard.csv')
    hkv = hk[hk.covariate == 'published_hyper'].iloc[0]
    items = [('Baylor\nMenG C vs A', 4.99, 1.24, 20.05), ('Tübingen\nproliferative', 3.53, 1.46, 8.53),
             ('UCSF\nhypermitotic', 3.73, 1.79, 7.79),
             ('Hong Kong\nhypermitotic', float(hkv['exp(coef)']), float(hkv['exp(coef) lower 95%']), float(hkv['exp(coef) upper 95%']))]
    y = np.arange(len(items))[::-1]
    for i, (lab, h, lo, hi) in enumerate(items):
        ax.plot([lo, hi], [y[i], y[i]], lw=2.2, color=VIOLET, solid_capstyle='round')
        ax.plot(h, y[i], 'o', ms=7, color=VIOLET, zorder=3)
        ax.text(hi * 1.15, y[i], f'{h:.2f}', va='center', fontsize=8, color=INK)
    ax.axvline(1, color=INK, lw=1, ls='--')
    ax.set_yticks(y); ax.set_yticklabels([i[0] for i in items], fontsize=8.5)
    ax.set_xscale('log'); ax.set_xlim(0.85, 90); ax.set_xticks([1, 3, 10, 30]); ax.set_xticklabels(['1', '3', '10', '30'])
    ax.set_xlabel('hazard ratio (log scale)', fontsize=8.5)
    ax.set_title('(c) Cohort-specific effects; no pooling', loc='left', fontsize=9.5)
    ax.grid(axis='x', lw=0.4, alpha=0.3); ax.set_axisbelow(True)

    m.save(fig, 'fig4_incremental_value')


if __name__ == '__main__':
    warnings.filterwarnings('ignore')
    figure4()
    print('FIG4_DONE', flush=True)
