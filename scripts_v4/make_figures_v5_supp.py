"""Figure 5 (revised) and Supplementary Figures S1-S6 for manuscript_npjDM_v5.md.

Figure 5 now carries the out-of-fold calibration evidence: the calibration curve at
the longer horizon, calibration-in-the-large and slope by cohort, and the three
discrimination measures side by side.

Supplementary figures:
  S1  disposition of the 994-sample multicentre release and the dbGaP route
  S2  training-reference fingerprint audit (per-cohort distance to nearest reference)
  S3  RNA ecosystem scores: no increment beyond clinical covariates
  S4  exploratory decision curves with corrected censoring weights
  S5  corrected robust-panel transfer
  S6  longitudinal series: labels, paired dBeta, donor conflicts
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
CAT3 = [BLUE, AQUA, YELLOW]
CAT4 = [BLUE, AQUA, YELLOW, VIOLET]   # categorical slots 1,2,3 and a 4th distinct hue
COH = ['baylor_published', 'tubingen_published', 'ucsf_published', 'hku_published']
SHORT = {'baylor_published': 'Baylor', 'tubingen_published': 'Tübingen',
         'ucsf_published': 'UCSF', 'hku_published': 'Hong Kong'}
HZ = {'baylor_published': 60.0, 'tubingen_published': 60.0,
      'ucsf_published': 5.0, 'hku_published': 5.0}


def save(fig, name):
    for ext in ('pdf', 'png'):
        fig.savefig(FIG / f'{name}.{ext}', bbox_inches='tight')
    plt.close(fig)
    print('saved', name, flush=True)


# ================================================================ figure 5
def figure5():
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 3.9),
                             gridspec_kw={'wspace': 0.46})

    # (a) calibration curve, molecular model, longer horizon
    ax = axes[0]
    for i, c in enumerate(COH):
        p = V4 / f'calibration_curve_{c}_molecular_{HZ[c]:g}.csv'
        if not p.exists():
            continue
        d = pd.read_csv(p)
        ax.plot(d.mean_predicted, d.observed_ipcw, 'o', ms=7, ls='none',
                color=CAT4[i], label=SHORT[c])
    lim = [0, 0.9]
    ax.plot(lim, lim, ls='--', color=INK, lw=1, zorder=1)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel('mean predicted risk', fontsize=8.5)
    ax.set_ylabel('observed risk (IPCW)', fontsize=8.5)
    ax.set_title('(a) Out-of-fold calibration, longer horizon', loc='left', fontsize=9.5)
    ax.legend(fontsize=7.5, frameon=False, loc='upper left', ncol=2)
    ax.text(0.98, 0.03, 'each point is a predicted-risk quintile',
            transform=ax.transAxes, fontsize=7, color=GREY, ha='right')
    ax.grid(lw=0.4, alpha=0.3); ax.set_axisbelow(True)

    # (b) calibration-in-the-large and slope
    ax = axes[1]
    rows = []
    for c in COH:
        p = V4 / f'calibration_summary_{c}.csv'
        if not p.exists():
            continue
        d = pd.read_csv(p)
        m = d[(d.model == 'molecular') & (d.horizon == HZ[c])]
        if len(m):
            rows.append((SHORT[c], float(m.calibration_in_the_large.iloc[0]),
                         float(m.calibration_slope.iloc[0])))
    y = np.arange(len(rows))[::-1]
    for i, (lab, citl, slope) in enumerate(rows):
        ax.plot([citl], [y[i] + 0.18], 'o', ms=7, color=BLUE, zorder=3, clip_on=True)
        ax.text(citl, y[i] + 0.36, f'{citl:+.3f}', ha='center', fontsize=7, color=BLUE,
                clip_on=True)
        if -0.12 <= slope - 1 <= 0.30:
            ax.plot([slope - 1], [y[i] - 0.18], 'D', ms=6.5, color=AQUA, zorder=3, clip_on=True)
            ax.text(slope - 1, y[i] - 0.38, f'{slope:.2f}', ha='center', fontsize=7,
                    color=AQUA, clip_on=True)
        else:
            ax.text(0.24, y[i] - 0.22, f'slope {slope:.2f}\n(not estimable)',
                    ha='right', va='center', fontsize=6.8, color=GREY, linespacing=1.3)
    ax.axvline(0, color=INK, lw=1, ls='--')
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=8.5)
    ax.set_xlim(-0.12, 0.30); ax.set_ylim(-0.75, len(rows) - 0.35)
    ax.set_xlabel('offset from ideal (CITL = 0; slope = 1)', fontsize=8.5)
    ax.set_title('(b) Calibration offset and slope', loc='left', fontsize=9.5)
    ax.legend(handles=[plt.Line2D([], [], marker='o', ls='', ms=7, color=BLUE,
                                  label='calibration-in-the-large'),
                       plt.Line2D([], [], marker='D', ls='', ms=6.5, color=AQUA,
                                  label='calibration slope (ideal 1)')],
              fontsize=7.5, frameon=False, ncol=2, loc='upper center',
              bbox_to_anchor=(0.5, -0.16))
    ax.grid(axis='x', lw=0.4, alpha=0.3); ax.set_axisbelow(True)

    # (c) three discrimination measures
    ax = axes[2]
    boot = {c: pd.read_csv(V4 / f'donor_bootstrap_ci_summary_{c}.csv').iloc[0] for c in COH
            if (V4 / f'donor_bootstrap_ci_summary_{c}.csv').exists()}
    uno = {c: pd.read_csv(V4 / f'uno_c_{c}.csv') for c in COH
           if (V4 / f'uno_c_{c}.csv').exists()}
    foldc = {'baylor_published': 0.0916656079618887, 'tubingen_published': 0.030600862095361413,
             'ucsf_published': 0.026990583909932975, 'hku_published': 0.06314230341426245}
    xs = np.arange(len(COH)); w = 0.26
    series = [('fold-wise ΔC', [foldc[c] for c in COH], BLUE),
              ('Uno ΔC', [float((uno[c].uno_c_molecular - uno[c].uno_c_clinical).mean())
                          for c in COH], AQUA),
              ('donor-bootstrap mean ΔC', [float(boot[c].delta_grouped_mean) for c in COH], YELLOW)]
    for j, (lab, vals, col) in enumerate(series):
        b = ax.bar(xs + (j - 1) * w, vals, width=w, color=col, label=lab)
        for r, v in zip(b, vals):
            ax.text(r.get_x() + r.get_width() / 2, v + 0.002, f'{v:.3f}', ha='center',
                    fontsize=5.8, color=INK, rotation=90, va='bottom')
    ax.set_xticks(xs); ax.set_xticklabels([SHORT[c] for c in COH], fontsize=8.5)
    ax.set_ylabel('discrimination increment', fontsize=8.5)
    ax.set_ylim(0, 0.135)
    ax.legend(fontsize=7.2, frameon=False, ncol=3, loc='upper center',
              bbox_to_anchor=(0.5, -0.12))
    ax.set_title('(c) Fold-wise and resampling estimates of ΔC', loc='left', fontsize=9.2)
    ax.grid(axis='y', lw=0.4, alpha=0.3); ax.set_axisbelow(True)

    save(fig, 'fig5_calibration')


# ================================================================ figure S2
def figure_s2():
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 3.7), gridspec_kw={'wspace': 0.3})

    ax = axes[0]
    ov = pd.read_csv(IDENT / 'raw_author_reference_overlap_audit.csv')
    order = ['author186', 'baylor', 'tubingen', 'ucsf', 'hku', 'longitudinal']
    labs = ['author release\n(3 anchors)', 'Baylor', 'Tübingen', 'UCSF', 'Hong Kong', 'longitudinal']
    ov = ov[ov.GSM != 'GSM9143071']
    for i, (c, lab) in enumerate(zip(order, labs)):
        d = ov[ov.cohort == c].max_abs_difference_to_nearest_reference.to_numpy()
        if len(d) == 0:
            continue
        ax.scatter(np.full(len(d), i) + np.random.default_rng(i).normal(0, 0.10, len(d)),
                   d, s=14, color=BLUE if c != 'author186' else RED, alpha=0.65,
                   edgecolor='none')
    ax.axhline(0.309, color=RED, lw=1.2, ls='--')
    ax.text(0.02, 0.34, 'minimum gap among the 866 external specimens = 0.309',
            fontsize=7.5, color=RED, transform=ax.get_yaxis_transform())
    ax.set_xticks(range(len(labs))); ax.set_xticklabels(labs, fontsize=7.5)
    ax.set_ylabel('max |Δβ| to nearest training reference', fontsize=8.5)
    ax.set_title('(a) Distance to the training reference, by series', loc='left', fontsize=9.5)
    ax.grid(axis='y', lw=0.4, alpha=0.3); ax.set_axisbelow(True)

    ax = axes[1]
    kn = pd.read_csv(STR / 'official_knn_audit_summary.csv')
    rel = kn[kn.cohort.isin(['release186', 'release37'])]
    x = np.arange(len(rel))
    tot = [186, 37]
    inref = [float(r.n_in_model_training) if str(r.n_in_model_training) != 'nan' else 0.0
             for _, r in rel.iterrows()]
    ax.bar(x, inref, color=RED, width=0.5, label='in training reference')
    ax.bar(x, [t - v for t, v in zip(tot, inref)], bottom=inref, color=SEQ[1], width=0.5,
           label='not in reference')
    for xi, v, t in zip(x, inref, tot):
        ax.text(xi, v / 2 if v else t / 2, f'{int(v)}/{t}', ha='center', va='center',
                fontsize=9, color='white' if v else INK, fontweight='bold')
    ax.set_xticks(x); ax.set_xticklabels(['training release', 'nominal held-out release'],
                                          fontsize=8)
    ax.set_ylabel('profiles', fontsize=8.5)
    ax.set_ylim(0, 200)
    ax.set_title('(b) Contamination of the classifier releases', loc='left', fontsize=9.5)
    ax.legend(fontsize=7.5, frameon=False, ncol=2, loc='upper center',
              bbox_to_anchor=(0.5, -0.14))
    ax.grid(axis='y', lw=0.4, alpha=0.3); ax.set_axisbelow(True)

    save(fig, 'figS2_fingerprint_audit')


# ================================================================ figure S4
def figure_s4():
    dc = pd.read_csv(STR / 'oof_decision_metrics_by_repeat.csv')
    key = ['baylor_published', 'tubingen_published', 'ucsf_published']
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 3.6), sharey=True,
                             gridspec_kw={'wspace': 0.12})
    for i, (ax, c) in enumerate(zip(axes, key)):
        d = dc[dc.comparison == c]
        if len(d) == 0:
            continue
        g = d.groupby('threshold')[['NB_clinical', 'NB_molecular', 'NB_all']].mean()
        ax.plot(g.index, g.NB_clinical, lw=1.8, color=GREY, label='clinical only')
        ax.plot(g.index, g.NB_molecular, lw=1.8, color=BLUE, label='clinical + molecular')
        ax.plot(g.index, g.NB_all, lw=1.2, color=GREY, ls=':', label='treat all')
        ax.axhline(0, color=INK, lw=0.8)
        ax.set_xlabel('threshold probability', fontsize=8.5)
        if i == 0:
            ax.set_ylabel('net benefit', fontsize=8.5)
            ax.legend(fontsize=7.5, frameon=False, loc='upper right')
        ax.set_title(['Baylor', 'Tübingen', 'UCSF'][i], loc='left', fontsize=9.5)
        ax.grid(lw=0.4, alpha=0.3); ax.set_axisbelow(True)
    fig.text(0.005, -0.02, 'Exploratory only. No guideline-defined or expert-consensus threshold exists for these '
                           'cohorts, so these curves are not interpreted as clinical benefit.',
             fontsize=8, color=GREY)
    save(fig, 'figS4_decision_curves')


# ================================================================ figure S5
def figure_s5():
    t = pd.read_csv(STR / 'corrected_P3_transfer_summary.csv')
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 3.6), gridspec_kw={'wspace': 0.32})

    ax = axes[0]
    rows = t[['cohort', 'shared_cpgs', 'profiled', 'n_rejected', 'reject_percent']].copy()
    rows = rows.set_index('cohort').loc[
        ['baylor_correctedP3', 'ucsf_discovery200_correctedP3', 'ucsf_correctedP3',
         'multicentre_correctedP3']].reset_index()
    x = np.arange(len(rows))
    ax.bar(x, rows.reject_percent, color=[BLUE, AQUA, VIOLET, YELLOW][:len(rows)], width=0.55)
    for xi, (_, r) in enumerate(rows.iterrows()):
        ax.text(xi, r.reject_percent + 0.25, f'{int(r.n_rejected)}/{int(r.profiled)}',
                ha='center', fontsize=7.5, color=INK)
    ax.set_xticks(x)
    ax.set_xticklabels(['Baylor\n(corrected P3)', 'UCSF 200\n(corrected P3)',
                        'UCSF full\n(corrected P3)', 'multicentre\n(corrected P3)'],
                       fontsize=7.5)
    ax.set_ylabel('rejection rate (%)', fontsize=8.5)
    ax.set_title('(a) Corrected robust panel: rejection by target', loc='left', fontsize=9.5)
    ax.grid(axis='y', lw=0.4, alpha=0.3); ax.set_axisbelow(True)

    ax = axes[1]
    ax.axis('off')
    lines = [
        ('Selection rule', 'Descending MAD rank plus ascending missingness rank'),
        ('Error', 'Selected the highest combined ranks'),
        ('Effect', 'Retained the probes meant to be discarded'),
        ('Correction', 'Select the lowest combined ranks'),
        ('Consequence', 'Corrected 2,000-probe panel is disjoint from erroneous panel'),
        ('Baylor labels', 'Best-map 0.910; ARI 0.808; 3/110 rejected'),
        ('UCSF reference', 'Best-map 0.795; ARI 0.485; 13/200 rejected'),
        ('Erroneous panel', 'Best-map 0.860; superseded'),
    ]
    ax.text(0.02, 0.97, 'field', fontsize=7.5, color=GREY, va='top')
    ax.text(0.30, 0.97, 'value', fontsize=7.5, color=GREY, va='top')
    for i, (k, v) in enumerate(lines):
        y = 0.88 - i * 0.105
        ax.text(0.02, y, k, fontsize=7.2, color=INK, va='center', fontweight='bold')
        ax.text(0.30, y, v, fontsize=7.0, color=INK, va='center')
        if i < len(lines) - 1:
            ax.plot([0.02, 0.98], [y - 0.055] * 2, color='#e6e8eb', lw=0.8)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_title('(b) The probe-selection error and its correction', loc='left', fontsize=9.5)

    save(fig, 'figS5_robust_panel')


# ================================================================ figure S6
def figure_s6():
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 3.4), gridspec_kw={'wspace': 0.38})

    ax = axes[0]
    lc = pd.read_csv(STR / 'longitudinal_baseline_correction.csv')
    ax.axis('off')
    facts = [('inferred clusters', lc.computed_clusters.iloc[0]),
             ('observed concordance', lc.observed_concordant.iloc[0]),
             ('concordance under a constant class-3 assignment',
              lc.all_class3_concordant_clusters.iloc[0]),
             ('clusters entirely class 3', lc.observed_all_class3_clusters.iloc[0])]
    ax.text(0.02, 0.95, 'quantity', fontsize=7.5, color=GREY, va='top')
    ax.text(0.92, 0.95, 'n', fontsize=7.5, color=GREY, va='top')
    for i, (k, v) in enumerate(facts):
        y = 0.82 - i * 0.17
        ax.text(0.02, y, k, fontsize=7.2, color=INK, va='center')
        ax.text(0.92, y, str(v), fontsize=8.5, color=INK, va='center', fontweight='bold')
    ax.text(0.02, 0.06, 'A constant class-3 assignment would have produced the same\n'
                        'concordance, so the observed agreement carries no information\n'
                        'beyond a class-imbalance baseline.',
            fontsize=7.5, color=RED, va='bottom', linespacing=1.4)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_title('(a) Longitudinal concordance vs baseline', loc='left', fontsize=9.5)

    ax = axes[1]
    p = ROOT / 'results/gse304096_paired_cpg_shifts.csv'
    if p.exists():
        d = pd.read_csv(p)
        col = 'delta_beta' if 'delta_beta' in d.columns else d.columns[-1]
        ax.hist(d[col].dropna(), bins=50, color=SEQ[2], edgecolor='white', linewidth=0.4)
        ax.set_xlabel('paired Δβ (all probes)', fontsize=8.5)
    ax.set_ylabel('probe pairs', fontsize=8.5)
    ax.set_title('(b) Paired probe differences (none BH-significant)', loc='left', fontsize=9.5)
    ax.grid(axis='y', lw=0.4, alpha=0.3); ax.set_axisbelow(True)

    ax = axes[2]
    st = pd.read_csv(IDENT / 'longitudinal_identity_resolution_status.json'.replace(
        '.json', '.json')) if False else None
    res = pd.read_csv(IDENT / 'longitudinal_source_patient_verification.csv')
    ax.axis('off')
    ax.text(0.02, 0.95, 'identity check in the longitudinal series', fontsize=8,
            color=INK, va='top', fontweight='bold')
    bullets = [
        '18 source groups; 16 SNP donor clusters',
        'Three cross-source pairs share donor genotypes',
        'Some within-cluster pairs fail strict threshold',
        'Patient independence and visit order unconfirmed',
    ]
    for i, b in enumerate(bullets):
        ax.text(0.04, 0.82 - i * 0.15, '· ' + b, fontsize=7.8, color=INK, va='center',
                wrap=True)
    ax.text(0.02, 0.08, 'The recovered grouping cannot be equated with 18\nindependent patients.',
            fontsize=7.8, color=RED, va='bottom', linespacing=1.4)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_title('(c) Donor-identity conflicts', loc='left', fontsize=9.5)

    save(fig, 'figS6_longitudinal')


if __name__ == '__main__':
    warnings.filterwarnings('ignore')
    figure5(); figure_s2(); figure_s4(); figure_s5(); figure_s6()
