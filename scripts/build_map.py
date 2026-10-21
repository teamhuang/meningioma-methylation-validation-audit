import gzip, re, os
import pandas as pd

gtf = r"E:\脑膜瘤npj\data\gencode.v38.annotation.gtf.gz"
out = r"E:\脑膜瘤npj\data\gencode_v38_ensg_symbol.csv"
rows = []
with gzip.open(gtf, "rt") as f:
    for line in f:
        if line.startswith("#"): continue
        p = line.split("\t")
        if p[2] != "gene": continue
        m = re.search(r'gene_id "([^"]+)"', p[8])
        n = re.search(r'gene_name "([^"]+)"', p[8])
        t = re.search(r'gene_type "([^"]+)"', p[8])
        if m and n:
            rows.append((m.group(1), n.group(1), t.group(1) if t else ""))
df = pd.DataFrame(rows, columns=["ensg_id", "symbol", "gene_type"])
df.to_csv(out, index=False)
print("genes:", len(df), " protein_coding:", (df["gene_type"] == "protein_coding").sum())
print(df.head(3).to_string())
