import os, glob, pandas as pd, numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def P(*a): return os.path.join(ROOT, *a)

print("=== ROOT:", ROOT)

# ---- 1. overlap audit: 866 vs 867 external rows, min gap 0.309 ----
f = P("identity_preprocessing_20261019/results/raw_author_reference_overlap_audit.csv")
if os.path.exists(f):
    df = pd.read_csv(f)
    print("\n[overlap audit] cols:", list(df.columns))
    if 'cohort' in df.columns:
        print("[overlap audit] cohort value counts:\n", df['cohort'].value_counts())
    # identify external vs reference-anchor/pilot
    for c in df.columns:
        if 'external' in c.lower() or 'reference' in c.lower() or 'match' in c.lower() or 'abs' in c.lower():
            print("   col sample:", c, "->", df[c].head(3).tolist())
    print("   nrows:", len(df))

# ---- 2. donor-level bootstrap summary files (the published intervals) ----
print("\n=== donor_bootstrap_ci_summary files in results_v4 ===")
for fp in sorted(glob.glob(P("results_v4/donor_bootstrap_ci_summary_*.csv"))):
    print("\n----", os.path.basename(fp))
    print(pd.read_csv(fp).to_string(index=False))

# ---- 3. do HK frozen-KNN / Tubingen frozen / surrogate bootstrap DRAWS exist? ----
print("\n=== bootstrap_*.csv in donor_aware: check dC columns + count ===")
da = P("identity_preprocessing_20261019/results/donor_aware")
for name in ["bootstrap_hku_author_KNN_raw_overlap_guard.csv",
             "bootstrap_tubingen_author_KNN_raw_overlap_guard.csv",
             "bootstrap_tubingen_centroid_surrogate_overlap_guard.csv",
             "bootstrap_baylor_centroid_surrogate_overlap_guard.csv",
             "bootstrap_baylor_author_KNN_raw_overlap_guard.csv"]:
    fp = os.path.join(da, name)
    if os.path.exists(fp):
        b = pd.read_csv(fp)
        print(f"\n{name}: nrows={len(b)} cols={list(b.columns)}")
        # print any dC/delta column head
        for c in b.columns:
            if 'd' in c.lower() and ('c' in c.lower() or 'delta' in c.lower() or 'concord' in c.lower()):
                print("   ", c, "n_null=", b[c].isna().sum(), "quantiles=",
                      np.nanpercentile(b[c].dropna(), [2.5,50,97.5]).round(5).tolist() if b[c].notna().any() else "NA")
    else:
        print(f"\n{name}: MISSING")

# ---- 4. UCSF record-level precise lower quantile (the -0.000) ----
for fp in [P("results_v4/donor_bootstrap_ci_summary_ucsf_published.csv"),
           P("results_v4/bootstrap_ucsf_published_overlap_guard.csv"),
           P("results_v4/donor_bootstrap_ci_ucsf_published.csv")]:
    if os.path.exists(fp):
        print("\n[ucsf boot]", os.path.basename(fp))
        d = pd.read_csv(fp)
        print(d.head(30).to_string(index=False))

# ---- 5. positional mapping -> Cohen kappa for grade ----
fp = P("identity_preprocessing_20261019/results/ucsf_old_positional_mapping_rejected.csv")
if os.path.exists(fp):
    d = pd.read_csv(fp)
    print("\n[positional mapping] cols:", list(d.columns), "n=", len(d))
    for c in d.columns:
        if 'grade' in c.lower():
            print("   grade col:", c, "->", d[c].head(10).tolist())
