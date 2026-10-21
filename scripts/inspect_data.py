import gzip, io, os, sys
import pandas as pd
import numpy as np

DATA = r"E:\脑膜瘤npj\data"
META = os.path.join(DATA, "meta")
os.makedirs(META, exist_ok=True)

def sec(t): print("\n" + "=" * 70 + "\n" + t + "\n" + "=" * 70)

# ---------- 1. Clinical xlsx ----------
sec("GSE189672_Clinical_data.xlsx")
xl = pd.ExcelFile(os.path.join(DATA, "GSE189672_Clinical_data.xlsx"))
print("sheets:", xl.sheet_names)
df = pd.read_excel(xl, sheet_name=xl.sheet_names[0], header=None)
print("raw shape:", df.shape)
print(df.head(8).to_string())

sec("GSE189521_Clinical_data.xlsx")
xl2 = pd.ExcelFile(os.path.join(DATA, "GSE189521_Clinical_data.xlsx"))
print("sheets:", xl2.sheet_names)
df2 = pd.read_excel(xl2, sheet_name=xl2.sheet_names[0], header=None)
print("raw shape:", df2.shape)
print(df2.head(8).to_string())

# ---------- 2. Methylation sample sheet ----------
sec("GSE189521_Methylation_sample_sheet.csv.gz")
ss = pd.read_csv(os.path.join(DATA, "GSE189521_Methylation_sample_sheet.csv.gz"))
print("shape:", ss.shape)
print("columns:", list(ss.columns))
print(ss.head(5).to_string())

# ---------- 3. RNA raw counts ----------
sec("GSE189672_RNAseq_raw_counts.csv.gz")
with gzip.open(os.path.join(DATA, "GSE189672_RNAseq_raw_counts.csv.gz"), "rt") as f:
    head = [next(f) for _ in range(3)]
print("header cols:", head[0].rstrip().split("\t")[:12], "...")
print("n_header_cols:", len(head[0].rstrip().split("\t")))
print("line2:", head[1][:120])

# ---------- 4. Series matrices -> sample maps ----------
sec("Series matrix sample maps")
for gse in ["GSE292325", "GSE292326", "GSE304096", "GSE270371", "GSE189672"]:
    p = os.path.join(DATA, gse + "_series_matrix.txt.gz")
    if not os.path.exists(p):
        print(gse, ": MISSING"); continue
    geo_acc, titles, chars, srcs = [], [], [], []
    with gzip.open(p, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"): break
            if line.startswith("!Sample_geo_accession"):
                geo_acc = [x.strip() for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_title"):
                titles = [x.strip() for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                chars.append([x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]])
            elif line.startswith("!Sample_source_name_ch1"):
                srcs = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
    print(f"\n{gse}: n={len(geo_acc)}")
    for i, c in enumerate(chars):
        uniq = sorted(set(c))
        print(f"  CHAR{i} [{len(uniq)} uniq]: " + " | ".join(u[:60] for u in uniq[:12]))
    # build long table
    tab = pd.DataFrame({"GSM": geo_acc, "title": titles})
    for i, c in enumerate(chars):
        tab[f"char{i}"] = c
    if srcs: tab["source"] = srcs
    tab.to_csv(os.path.join(META, f"{gse}_sample_map.csv"), index=False)
    print("  saved:", f"{gse}_sample_map.csv")

# ---------- 5. Processed methylation matrices: structure ----------
sec("Processed methylation matrix structure (streaming)")
for name, path in [("GSE292325", "GSE292325_matrix_processed.txt.gz"),
                   ("GSE292326", "GSE292326_matrix_processed.txt.gz"),
                   ("GSE304096", "GSE304096_longitudinal_processed.txt.gz"),
                   ("GSE270371", "GSE270371_Matrix_processed.txt.gz")]:
    p = os.path.join(DATA, path)
    with gzip.open(p, "rt", encoding="utf-8", errors="replace") as f:
        h = f.readline().rstrip("\n")
        r1 = f.readline().rstrip("\n")
        r2 = f.readline().rstrip("\n")
    hc = h.split("\t")
    print(f"\n{name}: n_cols={len(hc)}; first cols: {hc[:6]}; last col: {hc[-1]}")
    print("  row1 fields:", r1.split("\t")[:6])
    v1 = r1.split("\t")[1:6]
    print("  row1 values (first 5):", v1)
    try:
        vv = np.array([float(x) for x in r1.split("\t")[1:101]])
        print(f"  row1 val min={vv.min():.3f} max={vv.max():.3f} (n={len(vv)})")
    except Exception as e:
        print("  value parse note:", e)
print("\nINSPECT_DONE")
