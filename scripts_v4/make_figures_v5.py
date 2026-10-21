"""Publication figures for manuscript_npjDM_v5.md (revision 2).

Reads only from the archived audit packages and results_v4; produces PDF + PNG at
two-column final size using the Open Science style. Layout rules enforced here:
one y-axis per panel, legends outside the data area, no text below 7 pt at final
size, no overlapping labels, and no two incompatible scales on one axis.
"""
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from sklearn.metrics import adjusted_rand_score

STYLE = Path(r'C:\Users\Admin\AppData\Roaming\com.ai4s.workbench\runtime\xdg-config\opencode\skills\publication-figures\openscience.mplstyle')
plt.style.use(str(STYLE) if STYLE.exists() else 'default')

ROOT = Path(r'E:\脑膜瘤npj')
IDENT = ROOT / 'identity_preprocessing_20261019/results'
STR = ROOT / 'strengthening_20261019/results'
V4 = ROOT / 'results_v4'
FIG = ROOT / 'figures_v5'
FIG.mkdir(exist_ok=True)

BLUE, AQUA, YELLOW, GREEN, VIOLET, RED = '#2a78d6', '#1baf7a', '#eda100', '#008300', '#4a3aa7', '#e34948'
GREY = '#8a8f98'
INK = '#1b1d21'
SEQ = ['#cde2fb', '#9ec5f4', '#6da7ec', '#3987e5', '#256abf', '#184f95', '#104281']

COH = ['tubingen', 'ucsf', 'hku']
SHORT = {'tubingen': 'Tübingen', 'ucsf': 'UCSF', 'hku': 'Hong Kong',
         'baylor': 'Baylor', 'hku2': 'Hong Kong'}
CLASSES = ['Hypermitotic', 'Immune-enriched', 'Merlin-intact']
NCOL = {'tubingen': 149, 'ucsf': 200, 'hku': 365}


def save(fig, name):
    for ext in ('pdf', 'png'):
        fig.savefig(FIG / f'{name}.{ext}', bbox_inches='tight')
    plt.close(fig)
    print('saved', name, flush=True)


def load_agreement():
    out = {}
    for c in COH:
        p = pd.read_csv(IDENT / f'{c}_raw_author_ID_linked_predictions.csv')
        p['ok'] = (p.semantic_prediction == p.Group)
        p['tied'] = p.tied_top.astype(bool)
        out[c] = p
    return out


# ================================================================ figure 1
def figure1():
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 3.9),
                             gridspec_kw={'wspace': 0.34})

    ax = axes[0]
    snp = pd.read_csv(IDENT / 'SNP_donor_clusters_current.csv')
    rows = [('baylor', 'Baylor'), ('tubingen', 'Tübingen'), ('ucsf', 'UCSF'), ('hku', 'Hong Kong')]
    y = np.arange(len(rows))[::-1]
    for i, (c, lab) in enumerate(rows):
        g = snp[snp.cohort == c]
        sizes = g.groupby('donor_overlap_guard_key').size()
        mult = sizes[sizes > 1]
        singles = len(g) - mult.sum()
        ax.barh(y[i], singles, color=SEQ[1], height=0.55)
        left = singles
        for size, n in mult.value_counts().sort_index().items():
            c_ = BLUE if size == 2 else (AQUA if size == 3 else VIOLET)
            ax.barh(y[i], n * size, left=left, color=c_, height=0.55)
            left += n * size
    ax.set_yticks(y); ax.set_yticklabels([r[1] for r in rows], fontsize=9)
    ax.set_xlabel('specimens', fontsize=9)
    ax.set_xlim(0, 430)
    ax.set_title('(a) Donor clusters recovered by SNP fingerprint', loc='left', fontsize=9.5)
    handles = [plt.Rectangle((0, 0), 1, 1, color=SEQ[1]),
               plt.Rectangle((0, 0), 1, 1, color=BLUE),
               plt.Rectangle((0, 0), 1, 1, color=AQUA),
               plt.Rectangle((0, 0), 1, 1, color=VIOLET)]
    ax.legend(handles, ['single specimen', 'pair', 'triple', '4–6 members'],
              fontsize=7.0, frameon=False, loc='upper right',
              bbox_to_anchor=(1.0, 0.90), ncol=2, handlelength=1.0)
    ax.grid(axis='x', lw=0.4, alpha=0.3); ax.set_axisbelow(True)
    ax.annotate('Hong Kong: 23 multi-specimen clusters\ncovering 59 of 365 specimens',
                xy=(365, y[3]), xytext=(185, 0.48), fontsize=7.0, color=INK,
                arrowprops=dict(arrowstyle='->', color=GREY, lw=1.0))

    ax = axes[1]
    ax.axis('off')
    levels = [('Level C', 'External validation of a\nfrozen risk model', RED, 'not performed'),
               ('Level B', 'Frozen classifier,\ncross-cohort labels', AQUA, '3 cohorts'),
               ('Level A', 'Cohort-internal\npredictive modelling', BLUE, '4 cohorts')]
    for i, (tag, body, c, note) in enumerate(levels):
        yb = 0.66 - i * 0.30
        ax.add_patch(Rectangle((0.02, yb), 0.66, 0.24, facecolor=c, alpha=0.14,
                               edgecolor=c, linewidth=1.4))
        ax.text(0.05, yb + 0.185, tag, va='top', fontsize=8.5, color=c, fontweight='bold')
        ax.text(0.05, yb + 0.115, body, va='top', fontsize=7.8, color=INK, linespacing=1.35)
        ax.text(0.70, yb + 0.12, note, va='center', fontsize=7.0, color=c, fontweight='bold')
    ax.annotate('', xy=(1.12, 0.12), xytext=(1.12, 0.80),
                arrowprops=dict(arrowstyle='-|>', color=GREY, linewidth=1.2))
    ax.set_xlim(0, 1.20); ax.set_ylim(0, 1)
    ax.set_title('(b) Evidence hierarchy actually provided', loc='left', fontsize=9.5)

    ax = axes[2]
    labels = ['Baylor', 'Tübingen', 'UCSF', 'Hong\nKong', 'longitud.', 'ref.\nanchor', 'pilot']
    vals = [110, 149, 200, 365, 42, 3, 1]
    cols = [BLUE] * 5 + [GREY, YELLOW]
    b = ax.bar(labels, vals, color=cols, width=0.68)
    for r, v in zip(b, vals):
        ax.text(r.get_x() + r.get_width() / 2, v + 9, str(v), ha='center', fontsize=8, color=INK)
    ax.set_ylabel('specimens', fontsize=9)
    ax.set_ylim(0, 440)
    ax.set_title('(c) Genotyping inventory (870 specimens)', loc='left', fontsize=9.5)
    ax.tick_params(axis='x', labelsize=7.5, rotation=28)
    for tick in ax.get_xticklabels():
        tick.set_ha('right')
    ax.grid(axis='y', lw=0.4, alpha=0.3); ax.set_axisbelow(True)
    ax.text(0.03, 0.90, '866 external cohort specimens', transform=ax.transAxes,
            fontsize=8.5, color=BLUE, fontweight='bold')

    save(fig, 'fig1_cohort_architecture')


# ================================================================ figure 2
def figure2():
    ag = load_agreement()
    fig = plt.figure(figsize=(12.6, 7.4))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.35], hspace=0.62, wspace=0.40)

    ax = fig.add_subplot(gs[0, 0])
    ov = pd.read_csv(IDENT / 'raw_author_reference_overlap_audit.csv')
    d = ov[(ov.cohort != 'author186') & (ov.GSM != 'GSM9143071')].max_abs_difference_to_nearest_reference.to_numpy()
    ax.hist(d, bins=42, color=SEQ[2], edgecolor='white', linewidth=0.4)
    ax.axvline(0.309, color=RED, lw=1.8)
    ax.annotate('minimum gap to the nearest\ntraining-reference profile = 0.309',
                xy=(0.312, 8), xytext=(0.45, 46), fontsize=7.5, color=RED,
                arrowprops=dict(arrowstyle='->', color=RED, lw=1.1))
    ax.set_xlabel('max |Δβ| to nearest training reference', fontsize=8.5)
    ax.set_ylabel('external specimens', fontsize=8.5)
    ax.set_title('(a) No overlap at any threshold ≤ 0.30', loc='left', fontsize=9.5)
    ax.grid(axis='y', lw=0.4, alpha=0.3); ax.set_axisbelow(True)

    ax = fig.add_subplot(gs[0, 1])
    kn = pd.read_csv(STR / 'official_knn_audit_summary.csv').set_index('cohort')
    vals = [kn.loc['release186', 'n_in_model_training'], kn.loc['release37', 'n_in_model_training']]
    tot = [186, 37]
    x = np.arange(2)
    ax.bar(x, vals, color=RED, width=0.5, label='in training reference')
    ax.bar(x, [t - v for t, v in zip(tot, vals)], bottom=vals, color=SEQ[1], width=0.5,
           label='not in reference')
    for xi, v, t in zip(x, vals, tot):
        ax.text(xi, t + 5, f'{int(v)}/{t}', ha='center', fontsize=8.5, color=INK, fontweight='bold')
    ax.set_xticks(x); ax.set_xticklabels(['training\nrelease', 'nominal held-out\nrelease'], fontsize=8)
    ax.set_ylabel('profiles', fontsize=8.5)
    ax.set_ylim(0, 215)
    ax.set_title('(b) Held-out release is 54% training data', loc='left', fontsize=9.5)
    ax.legend(fontsize=7.5, frameon=False, loc='upper right')
    ax.grid(axis='y', lw=0.4, alpha=0.3); ax.set_axisbelow(True)

    ax = fig.add_subplot(gs[0, 2])
    vals = []
    for c in COH:
        a = ag[c][~ag[c].tied]
        vals.append(adjusted_rand_score(a.Group, a.semantic_prediction))
    b = ax.bar([SHORT[c] for c in COH], vals, color=[BLUE, AQUA, VIOLET], width=0.55)
    for r, v in zip(b, vals):
        ax.text(r.get_x() + r.get_width() / 2, v + 0.018, f'{v:.3f}', ha='center',
                fontsize=8.5, color=INK)
    ax.set_ylabel('adjusted Rand index', fontsize=8.5)
    ax.set_ylim(0, 1.0)
    ax.tick_params(axis='x', labelsize=8.5)
    ax.set_title('(c) Cross-nomenclature agreement', loc='left', fontsize=9.5)
    ax.grid(axis='y', lw=0.4, alpha=0.3); ax.set_axisbelow(True)

    for k, c in enumerate(COH):
        ax = fig.add_subplot(gs[1, k])
        cm = pd.read_csv(IDENT / f'{c}_raw_author_confusion_all_assigned.csv', index_col=0)
        cm = cm.loc[CLASSES, CLASSES]
        v = cm.to_numpy()
        ax.imshow(v, cmap=plt.get_cmap('Blues'), vmin=0, vmax=v.max() * 0.75)
        short = ['Hyper-\nmitotic', 'Immune-\nenriched', 'Merlin-\nintact']
        ax.set_xticks(range(3)); ax.set_xticklabels(short, fontsize=7.5)
        ax.set_yticks(range(3)); ax.set_yticklabels(short, fontsize=7.5)
        ax.set_xlabel('predicted class', fontsize=8)
        if k == 0:
            ax.set_ylabel('published group', fontsize=8)
        for i in range(3):
            for j in range(3):
                ax.text(j, i, f'{v[i, j]}', ha='center', va='center', fontsize=9,
                        color='white' if v[i, j] > v.max() * 0.4 else INK)
        ax.set_title(f'({chr(100 + k)}) {SHORT[c]}  (n={int(v.sum())}, accuracy {np.trace(v) / v.sum():.3f})',
                     loc='left', fontsize=9.5)
        for s in ax.spines.values():
            s.set_visible(False)
    fig.text(0.005, 0.005,
             'Residual error is structured: the aggressive and the indolent extremes are each confused with intermediate disease, not with one another.',
             fontsize=8, color=GREY)

    save(fig, 'fig2_independence_agreement')


# ================================================================ figure 3
def figure3():
    ag = load_agreement()
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 3.8),
                             gridspec_kw={'wspace': 0.42})

    ax = axes[0]
    for i, c in enumerate(COH):
        a = ag[c][~ag[c].tied]
        order = np.argsort(a.confidence.to_numpy())
        conf = a.confidence.to_numpy()[order]
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

    ax = axes[1]
    xs = np.arange(3); widths = 0.26
    for i, g in enumerate(CLASSES):
        vals = []
        for c in COH:
            a = ag[c][~ag[c].tied]
            sub = a[a.Group == g]
            vals.append((sub.confidence >= 0.75).mean())
        b = ax.bar(xs + (i - 1) * widths, vals, width=widths,
                   color=[RED, YELLOW, BLUE][i], label=g)
        for r, v in zip(b, vals):
            ax.text(r.get_x() + r.get_width() / 2, v + 0.018, f'{v:.2f}', ha='center',
                    fontsize=7.5, color=INK)
    ax.set_xticks(xs); ax.set_xticklabels([SHORT[c] for c in COH], fontsize=8.5)
    ax.set_ylabel('retention at confidence 0.75', fontsize=8.5)
    ax.set_ylim(0, 1.16)
    ax.set_title('(b) Filtering drops the highest-risk class', loc='left', fontsize=9.5)
    ax.legend(fontsize=7.5, frameon=False, ncol=1, loc='upper center',
              bbox_to_anchor=(0.5, 1.02), handlelength=1.1)
    ax.grid(axis='y', lw=0.4, alpha=0.3); ax.set_axisbelow(True)

    ax = axes[2]
    y = np.arange(3)[::-1]
    for i, c in enumerate(COH):
        a = ag[c][~ag[c].tied]
        hyp = 1 - (a[a.Group == 'Hypermitotic'].confidence >= 0.75).mean()
        mer = 1 - (a[a.Group == 'Merlin-intact'].confidence >= 0.75).mean()
        ax.barh(i, hyp, color=RED, height=0.5)
        ax.barh(i, -mer, color=SEQ[1], height=0.5)
        ax.text(hyp + 0.02, i, f'{hyp:.0%} hypermitotic', va='center', fontsize=8, color=INK)
        ax.text(-mer - 0.02, i, f'{mer:.0%} Merlin-intact', va='center', ha='right',
                fontsize=8, color=INK)
    ax.axvline(0, color=INK, lw=1)
    ax.set_yticks(y); ax.set_yticklabels([SHORT[c] for c in COH], fontsize=8.5)
    ax.set_xlim(-0.42, 0.85)
    ax.set_xlabel('fraction discarded by the 0.75 threshold', fontsize=8.5)
    ax.set_title('(c) Asymmetric class loss', loc='left', fontsize=9.5)
    ax.grid(axis='x', lw=0.4, alpha=0.3); ax.set_axisbelow(True)
    ax.tick_params(axis='y', pad=2)

    save(fig, 'fig3_risk_coverage')


# ================================================================ figure 4
def figure4():
    rows = [
        ('Baylor, published MenG', 0.092, -0.004, 0.170, 0.072, 0.073),
        ('Baylor, surrogate', 0.025, -0.009, 0.132, 0.033, 0.029),
        ('Baylor, frozen KNN', 0.016, -0.016, 0.116, 0.018, None),
        ('Tübingen, published Toronto', 0.031, -0.004, 0.107, 0.046, 0.023),
        ('Tübingen, surrogate', 0.025, -0.014, 0.097, 0.024, 0.011),
        ('Tübingen, frozen KNN', 0.041, -0.006, 0.118, 0.041, None),
        ('UCSF, published groups', 0.027, 0.011, 0.076, 0.030, 0.024),
        ('Hong Kong, published groups', 0.063, 0.020, 0.110, 0.056, None),
        ('Hong Kong, frozen KNN', 0.059, 0.025, 0.116, 0.059, None),
    ]
    fig = plt.figure(figsize=(12.6, 6.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.05, 0.9], wspace=0.46)

    # (a)
    ax = fig.add_subplot(gs[0, 0])
    y = np.arange(len(rows))[::-1]
    for i, (lab, dc, lo, hi, _, _) in enumerate(rows):
        col = BLUE if lo > 0 else GREY
        ax.plot([lo, hi], [y[i], y[i]], lw=2.2, color=col, solid_capstyle='round')
        ax.plot(dc, y[i], 'o', ms=7, color=col, zorder=3)
    ax.axvline(0, color=INK, lw=1, ls='--')
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=8.5)
    ax.set_xlabel('donor-grouped ΔC (apparent bootstrap interval)', fontsize=8.5)
    ax.set_xlim(-0.045, 0.185)
    ax.set_title('(a) Primary effect estimate', loc='left', fontsize=9.5)
    ax.legend(handles=[plt.Line2D([], [], marker='o', ls='-', ms=7, color=BLUE,
                                  label='interval excludes zero'),
                       plt.Line2D([], [], marker='o', ls='-', ms=7, color=GREY,
                                  label='interval includes zero')],
              fontsize=7.5, frameon=False, loc='lower left', bbox_to_anchor=(0.0, -0.01))
    ax.grid(axis='x', lw=0.4, alpha=0.3); ax.set_axisbelow(True)

    # (b) three estimands on one comparable axis
    ax = fig.add_subplot(gs[0, 1])
    sub = [r for r in rows if r[5] is not None]
    y = np.arange(len(sub))[::-1]
    for i, (lab, dc, lo, hi, corr, rec) in enumerate(sub):
        ax.plot([lo, hi], [y[i] + 0.22] * 2, lw=2.2, color=SEQ[0], solid_capstyle='round')
        ax.plot(dc, y[i] + 0.22, 'o', ms=7, color=BLUE, zorder=3)
        ax.plot(corr, y[i] - 0.18, 'D', ms=6.5, color=AQUA, zorder=3)
        ax.plot(rec, y[i] - 0.52, 's', ms=6, color=GREY, zorder=3)
        ax.text(dc, y[i] + 0.40, f'{dc:.3f}', ha='center', fontsize=7, color=BLUE)
        ax.text(corr, y[i] - 0.34, f'{corr:.3f}', ha='center', fontsize=7, color=AQUA)
        ax.text(rec, y[i] - 0.66, f'{rec:.3f}', ha='center', fontsize=7, color=GREY)
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in sub], fontsize=8.5)
    ax.set_xlabel('ΔC', fontsize=8.5)
    ax.set_xlim(-0.05, 0.19); ax.set_ylim(-1.1, len(sub) - 0.3)
    ax.set_title('(b) Four estimands are not one test', loc='left', fontsize=9.5)
    ax.legend(handles=[plt.Line2D([], [], marker='o', ls='-', ms=7, color=BLUE,
                                  label='donor-grouped cross-validated ΔC (primary)'),
                       plt.Line2D([], [], marker='D', ls='', ms=6.5, color=AQUA,
                                  label='optimism-corrected ΔC'),
                       plt.Line2D([], [], marker='s', ls='', ms=6, color=GREY,
                                  label='record-level ΔC')],
              fontsize=7.2, frameon=False, loc='lower left')
    ax.grid(axis='x', lw=0.4, alpha=0.3); ax.set_axisbelow(True)

    # (c)
    ax = fig.add_subplot(gs[0, 2])
    hk = pd.read_csv(IDENT / 'donor_aware/cluster_robust_cox_hku_published_overlap_guard.csv')
    hkv = hk[hk.covariate == 'published_hyper'].iloc[0]
    items = [('Baylor\nMenG C vs A', 4.99, 1.24, 20.05),
             ('Tübingen\nproliferative', 3.53, 1.46, 8.53),
             ('UCSF\nhypermitotic', 3.73, 1.79, 7.79),
             ('Hong Kong\nhypermitotic', float(hkv['exp(coef)']),
              float(hkv['exp(coef) lower 95%']), float(hkv['exp(coef) upper 95%']))]
    y = np.arange(len(items))[::-1]
    for i, (lab, h, lo, hi) in enumerate(items):
        ax.plot([lo, hi], [y[i], y[i]], lw=2.2, color=VIOLET, solid_capstyle='round')
        ax.plot(h, y[i], 'o', ms=7, color=VIOLET, zorder=3)
        ax.text(hi * 1.15, y[i], f'{h:.2f}', va='center', fontsize=8, color=INK)
    ax.axvline(1, color=INK, lw=1, ls='--')
    ax.set_yticks(y); ax.set_yticklabels([i[0] for i in items], fontsize=8.5)
    ax.set_xscale('log'); ax.set_xlim(0.85, 90)
    ax.set_xticks([1, 3, 10, 30]); ax.set_xticklabels(['1', '3', '10', '30'])
    ax.set_xlabel('hazard ratio (log scale)', fontsize=8.5)
    ax.set_title('(c) Cohort-specific effects; no pooling', loc='left', fontsize=9.5)
    ax.grid(axis='x', lw=0.4, alpha=0.3); ax.set_axisbelow(True)

    save(fig, 'fig4_incremental_value')


# ================================================================ figure 5
def figure5():
    risk = pd.read_csv(STR / 'oof_risk_summary.csv')
    slope = pd.read_csv(STR / 'apparent_calibration_slopes.csv')
    oof = pd.read_csv(STR / 'oof_slope_summary.csv')
    key = [('baylor_published', 'Baylor'), ('tubingen_transferred', 'Tübingen'),
           ('ucsf_published', 'UCSF')]
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 3.6),
                             gridspec_kw={'wspace': 0.46})

    ax = axes[0]
    y = np.arange(len(key))
    for i, (k, lab) in enumerate(key):
        a = slope[(slope.comparison == k) & (slope.model == 'molecular')].iloc[0]
        o = oof[(oof.comparison == k) & (oof.model == 'molecular')].iloc[0]
        ax.plot([a['CI_low'], a['CI_high']], [y[i] + 0.17] * 2, lw=2.2, color=BLUE,
                solid_capstyle='round')
        ax.plot(a['slope'], y[i] + 0.17, 'o', ms=7, color=BLUE, zorder=3)
        ax.plot([o['min'], o['max']], [y[i] - 0.17] * 2, lw=2.2, color=AQUA,
                solid_capstyle='round')
        ax.plot(o['mean'], y[i] - 0.17, 'D', ms=6.5, color=AQUA, zorder=3)
        ax.text(a['CI_high'] + 0.12, y[i] + 0.17, f"{a['slope']:.2f}", fontsize=7.5,
                color=BLUE, va='center')
        ax.text(o['max'] + 0.12, y[i] - 0.17, f"{o['mean']:.2f}", fontsize=7.5,
                color=AQUA, va='center')
    ax.axvline(1, color=INK, lw=1, ls='--')
    ax.set_yticks(y); ax.set_yticklabels([k[1] for k in key], fontsize=8.5)
    ax.set_xlabel('recalibration slope', fontsize=8.5)
    ax.set_xlim(0, 3.6)
    ax.set_title('(a) Apparent vs cross-fitted slopes', loc='left', fontsize=9.5)
    ax.legend(handles=[plt.Line2D([], [], marker='o', ls='-', ms=7, color=BLUE,
                                  label='apparent (95% CI)'),
                       plt.Line2D([], [], marker='D', ls='-', ms=6.5, color=AQUA,
                                  label='cross-fitted (range)')],
              fontsize=7.5, frameon=False, loc='lower right')
    ax.grid(axis='x', lw=0.4, alpha=0.3); ax.set_axisbelow(True)

    ax = axes[1]
    ops = []
    for k, lab in key:
        r = risk[risk.comparison == k]
        hz = 5.0 if k == 'ucsf_published' else 60.0
        ops.append(float(r[r.horizon == hz].OP_molecular_mean.iloc[0]))
    b = ax.bar([k[1] for k in key], ops, color=[BLUE, AQUA, VIOLET], width=0.55)
    ax.axhline(1, color=INK, lw=1, ls='--')
    for r, v in zip(b, ops):
        ax.text(r.get_x() + r.get_width() / 2, v + 0.012, f'{v:.3f}', ha='center',
                fontsize=8.5, color=INK)
    ax.set_ylabel('observed / predicted (5 years)', fontsize=8.5)
    ax.set_ylim(0.9, 1.19)
    ax.tick_params(axis='x', labelsize=8.5)
    ax.set_title('(b) Aggregate calibration near 1', loc='left', fontsize=9.5)
    ax.grid(axis='y', lw=0.4, alpha=0.3); ax.set_axisbelow(True)

    ax = axes[2]
    xs = np.arange(len(key)); w = 0.34
    for j, (col, lab, c) in enumerate([('brier_clinical_mean', 'clinical only', GREY),
                                       ('brier_molecular_mean', 'clinical + molecular', BLUE)]):
        vals = []
        for k, _ in key:
            r = risk[risk.comparison == k]
            hz = 5.0 if k == 'ucsf_published' else 60.0
            vals.append(float(r[r.horizon == hz][col].iloc[0]))
        bb = ax.bar(xs + (j - 0.5) * w, vals, width=w, color=c, label=lab)
        for rr, v in zip(bb, vals):
            ax.text(rr.get_x() + rr.get_width() / 2, v + 0.005, f'{v:.3f}', ha='center',
                    fontsize=7, color=INK)
    ax.set_xticks(xs); ax.set_xticklabels([k[1] for k in key], fontsize=8.5)
    ax.set_ylabel('cross-validated 5-year Brier score', fontsize=8.5)
    ax.set_ylim(0, 0.33)
    ax.legend(fontsize=7.5, frameon=False, loc='upper right')
    ax.set_title('(c) Absolute-risk error moves both ways', loc='left', fontsize=9.5)
    ax.grid(axis='y', lw=0.4, alpha=0.3); ax.set_axisbelow(True)

    save(fig, 'fig5_calibration')


# ================================================================ figure 6
def figure6():
    dims = ['Identity\nintegrity', 'Training\nindependence', 'Analytical\ntransportability',
            'Clinical predictive\nvalidity', 'Decision\nutility']
    PASS, FAIL, NA = 1, 0, 2
    sym = {PASS: '✓', FAIL: '✗', NA: 'n/a'}
    col = {PASS: AQUA, FAIL: RED, NA: GREY}
    grid = [
        [PASS, PASS, PASS, PASS, NA],
        [PASS, PASS, PASS, PASS, NA],
        [PASS, PASS, PASS, PASS, NA],
        [PASS, PASS, PASS, PASS, NA],
        [PASS, PASS, PASS, FAIL, NA],
    ]
    labs = ['Baylor', 'Tübingen', 'UCSF', 'Hong Kong', 'overall']
    notes = ['label linkage not assessable', '', '', '', 'no level-C analysis']

    fig = plt.figure(figsize=(12.6, 4.3))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.95, 1.15], wspace=0.30)

    ax = fig.add_subplot(gs[0, 0])
    for i in range(len(labs)):
        for j in range(5):
            s = grid[i][j]
            ax.add_patch(Rectangle((j + 0.04, len(labs) - i - 0.96), 0.92, 0.88,
                                   facecolor=col[s], alpha=0.18, edgecolor=col[s], lw=1.2))
            ax.text(j + 0.5, len(labs) - i - 0.5, sym[s], ha='center', va='center',
                    fontsize=11, color=col[s], fontweight='bold')
    ax.set_xlim(0, 5); ax.set_ylim(0, len(labs))
    ax.set_xticks(np.arange(5) + 0.5); ax.set_xticklabels(dims, fontsize=7.5)
    ax.set_yticks(np.arange(len(labs))[::-1] + 0.5); ax.set_yticklabels(labs, fontsize=8.5)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    ax.set_title('(a) Deployment-readiness verdicts (audit-defined criteria)',
                 loc='left', fontsize=9.5)
    ax.legend(handles=[plt.Line2D([], [], marker='o', ls='', ms=9, color=AQUA, label='pass'),
                       plt.Line2D([], [], marker='o', ls='', ms=9, color=RED, label='fail'),
                       plt.Line2D([], [], marker='o', ls='', ms=9, color=GREY, label='cannot assess')],
              fontsize=7.5, frameon=False, ncol=3, loc='upper center',
              bbox_to_anchor=(0.5, -0.02), handlelength=0.8)

    # (b) text panel: values are on different scales, so no shared axis is drawn
    ax = fig.add_subplot(gs[0, 1])
    ax.axis('off')
    items = [
        ('Held-out release treated as independent', 'best-map 0.857', 'release excluded (54% training)'),
        ('Positional instead of identifier join', 'agreement 0.265', 'agreement 0.872'),
        ('Erroneous probe-selection rule', 'best-map 0.860', 'best-map 0.795'),
        ('Record-level instead of donor grouping', 'ΔC 0.011', 'ΔC 0.025'),
        ('Confidence threshold without coverage', 'agreement 0.936', 'hypermitotic retention 0.607'),
    ]
    ax.text(0.02, 0.99, 'check omitted  →  value reported', fontsize=8, color=RED,
            fontweight='bold', va='top')
    ax.text(0.52, 0.99, 'audit performed  →  value obtained', fontsize=8, color=BLUE,
            fontweight='bold', va='top')
    for i, (name, naive, audited) in enumerate(items):
        y = 0.86 - i * 0.165
        ax.text(0.02, y, name, fontsize=8.2, color=INK, va='center')
        ax.text(0.02, y - 0.062, naive, fontsize=8.2, color=RED, va='center')
        ax.text(0.52, y - 0.062, audited, fontsize=8.2, color=BLUE, va='center')
        ax.plot([0.46], [y - 0.062], marker='>', ms=6, color=GREY)
        if i < len(items) - 1:
            ax.plot([0.02, 0.98], [y - 0.115] * 2, color='#e6e8eb', lw=0.8)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.02)
    ax.set_title('(b) What a conventional validation pipeline would have reported',
                 loc='left', fontsize=9.5)
    ax.text(0.02, 0.015, 'Rows are on different scales; no shared axis is drawn. '
                         'None of the five errors produces a failed test or an implausible value.',
            fontsize=7.5, color=GREY, va='bottom')

    save(fig, 'fig6_deployment_audit')


if __name__ == '__main__':
    warnings.filterwarnings('ignore')
    figure1(); figure2(); figure3(); figure4(); figure5(); figure6()
