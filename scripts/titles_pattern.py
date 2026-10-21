import gzip, re
from collections import Counter

p = r"E:\脑膜瘤npj\data\GSE183647_series_matrix.txt.gz"
pats = Counter()
examples = {}
with gzip.open(p, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("!Sample_title"):
            tit = [x.strip().strip('"') for x in line.rstrip().split("\t")[1:]]
            for t in tit:
                key = re.sub(r"\d+", "N", t)
                pats[key] += 1
                examples.setdefault(key, t)
        if line.startswith("!series_matrix_table_begin"):
            break
for k, v in pats.most_common(10):
    print(v, "|", k, "| e.g.", examples[k])
