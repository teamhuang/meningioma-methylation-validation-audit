import os, re, warnings
import pandas as pd
warnings.filterwarnings("ignore")

D = r"E:\脑膜瘤npj\data\suppl\natgen2022"
pat = re.compile(r"^M\d+$")
for fn in sorted(os.listdir(D)):
    if not fn.endswith(".xlsx"):
        continue
    p = os.path.join(D, fn)
    try:
        xl = pd.ExcelFile(p)
    except Exception:
        continue
    for sh in xl.sheet_names:
        try:
            df = xl.parse(sh)
        except Exception:
            continue
        if df.shape[0] < 100:
            continue
        # check every column for M-ID-like content
        for c in df.columns:
            vals = df[c].astype(str).str.strip()
            n_m = int(vals.str.match(pat).sum())
            if n_m >= 100:
                kw_cols = [str(x) for x in df.columns if re.search(r"recurren|pfs|rfs|follow|month|survival|LF|grade|Group", str(x), re.I)]
                print(f"{fn} [{sh}] col='{c}' M-ids={n_m} shape={df.shape} clinical_cols={kw_cols[:8]}")
                break
print("SCAN2_DONE")
