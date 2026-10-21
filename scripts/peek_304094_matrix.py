import gzip
import pandas as pd

p = r"E:\脑膜瘤npj\data\GSE304094_discovery_processed.txt.gz"
with gzip.open(p, "rt", encoding="utf-8", errors="replace") as f:
    h = f.readline().rstrip("\n").split("\t")
    r1 = f.readline().rstrip("\n").split("\t")
    r2 = f.readline().rstrip("\n").split("\t")
print("cols:", len(h), "first:", h[:6], "last:", h[-3:])
print("row1 id:", r1[0], "vals:", r1[1:4])
print("row2 id:", r2[0], "vals:", r2[1:4])
vals = [float(x) for x in r1[1:51]]
print("row1 value range: %.3f-%.3f" % (min(vals), max(vals)))
