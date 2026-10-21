import pandas as pd
import re
import warnings
warnings.filterwarnings("ignore")

pat = re.compile(r"\d{9,12}_R\d{2}C\d{2}")
targets = {
    "suppl": r"E:\脑膜瘤npj\data\suppl\protocadherin_suppldata.xlsx",
    "source": r"E:\脑膜瘤npj\data\suppl\protocadherin_sourcedata.xlsx",
}
for name, p in targets.items():
    xl = pd.ExcelFile(p)
    for sh in xl.sheet_names:
        try:
            df = xl.parse(sh, header=None)
        except Exception:
            continue
        if df.empty:
            continue
        flat = " ".join(str(v) for v in df.values.ravel()[:200000])
        hits = pat.findall(flat)
        if hits:
            uniq = sorted(set(hits))
            print(f"{name} [{sh}]: {len(hits)} sentrix mentions, {len(uniq)} unique; e.g. {uniq[:5]}")
        # find longitudinal-style tables
        first = " ".join(str(v) for v in df.values.ravel()[:60])
        for kw in ["Longitudinal", "longitudinal", "paired", "Pair", "T26_00"]:
            if kw in first:
                print(f"  {name} [{sh}]: keyword '{kw}' in header area")
                break
print("SCAN_DONE")
