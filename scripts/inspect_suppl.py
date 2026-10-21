import os
import pandas as pd

D = r"E:\脑膜瘤npj\data\suppl"
for name in ["protocadherin_suppldata.xlsx", "protocadherin_sourcedata.xlsx"]:
    p = os.path.join(D, name)
    try:
        xl = pd.ExcelFile(p)
        print(f"=== {name}: sheets = {xl.sheet_names}")
        for sh in xl.sheet_names[:12]:
            df = xl.parse(sh, nrows=4, header=None)
            print(f"  [{sh}] shape~{df.shape} | first cells: {[str(x)[:28] for x in df.iloc[0, :6].tolist()]}")
    except Exception as e:
        print(f"=== {name} FAIL: {type(e).__name__} {str(e)[:120]}")
