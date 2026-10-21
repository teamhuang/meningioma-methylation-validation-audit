import gzip, os
import pandas as pd

# 1. GSE183647 series matrix titles
p = os.path.join(r"E:\脑膜瘤npj\data", "GSE183647_series_matrix.txt.gz")
if not os.path.exists(p):
    import urllib.request
    url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183647/matrix/GSE183647_series_matrix.txt.gz"
    urllib.request.urlretrieve(url, p)
acc, tit = [], []
with gzip.open(p, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("!series_matrix_table_begin"):
            break
        if line.startswith("!Sample_geo_accession"):
            acc = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        elif line.startswith("!Sample_title"):
            tit = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
print("183647 series matrix: n =", len(acc))
print("titles head:", tit[:8])

# 2. GSE183647 RAW filelist (GSM -> barcode)
fl_path = os.path.join(r"E:\脑膜瘤npj\data\suppl", "GSE183647_filelist.txt")
if not os.path.exists(fl_path):
    import urllib.request
    url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183647/suppl/filelist.txt"
    urllib.request.urlretrieve(url, fl_path)
gsm2bar = {}
for line in open(fl_path, encoding="utf-8", errors="replace"):
    if line.startswith("File\t"):
        parts = line.split("\t")[1].strip().split("_")
        if parts[0].startswith("GSM") and len(parts) >= 3:
            gsm2bar[parts[0]] = parts[1] + "_" + parts[2]
print("183647 filelist: GSM->barcode entries =", len(gsm2bar))

# 3. bridge: GSM -> title(M-ID?) x GSM -> barcode
bridged = [(g, t, gsm2bar.get(g)) for g, t in zip(acc, tit)]
sample = bridged[:6]
for s in sample:
    print(s)

# 4. GSE183653 (RNA, M-ID columns): series matrix titles
p3 = os.path.join(r"E:\脑膜瘤npj\data", "GSE183653_series_matrix.txt.gz")
acc3, tit3 = [], []
with gzip.open(p3, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("!series_matrix_table_begin"):
            break
        if line.startswith("!Sample_geo_accession"):
            acc3 = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        elif line.startswith("!Sample_title"):
            tit3 = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
print("183653 titles head:", tit3[:6])
