import os, re, warnings
import pandas as pd
warnings.filterwarnings("ignore")

D = r"E:\脑膜瘤npj\data\suppl\natgen2022"
kw = re.compile(r"recurren|progression|pfs|rfs|follow|month|survival|outcome|event", re.I)
rows_kw = []
for fn in sorted(os.listdir(D)):
    if not fn.endswith(".xlsx"):
        continue
    p = os.path.join(D, fn)
    try:
        xl = pd.ExcelFile(p)
    except Exception as e:
        print(fn, "OPEN_FAIL", str(e)[:60]); continue
    for sh in xl.sheet_names:
        try:
            df = xl.parse(sh, nrows=3)
        except Exception:
            continue
        cols = " ".join(str(c) for c in df.columns)
        if kw.search(cols):
            rows_kw.append((fn, sh, len(df.columns), cols[:160]))
        # also check row count hints (parse without data)
        try:
            full_ncols = len(xl.parse(sh, nrows=0).columns)
        except Exception:
            full_ncols = len(df.columns)
print("=== sheets whose COLUMNS mention outcome keywords ===")
for r in rows_kw:
    print(f"{r[0]} [{r[1]}] ncols={r[2]}: {r[3]}")
print("SCAN_DONE")
