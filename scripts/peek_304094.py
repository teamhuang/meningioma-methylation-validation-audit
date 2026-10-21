import gzip

p = r"E:\脑膜瘤npj\data\GSE304094_series_matrix.txt.gz"
acc, tit, chars, srcs = [], [], [], []
with gzip.open(p, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("!series_matrix_table_begin"):
            break
        if line.startswith("!Sample_geo_accession"):
            acc = line.rstrip().split("\t")[1:]
        elif line.startswith("!Sample_title"):
            tit = [x.strip().strip('"') for x in line.rstrip().split("\t")[1:]]
        elif line.startswith("!Sample_characteristics_ch1"):
            chars.append([x.strip().strip('"') for x in line.rstrip().split("\t")[1:]])
        elif line.startswith("!Sample_source_name_ch1"):
            srcs = [x.strip().strip('"') for x in line.rstrip().split("\t")[1:]]
print("n =", len(acc))
print("titles head:", tit[:6])
for i, c in enumerate(chars):
    u = sorted(set(c))
    print(f"char{i} uniq[{len(u)}]:", u[:6])
if srcs:
    print("sources head:", srcs[:4])
