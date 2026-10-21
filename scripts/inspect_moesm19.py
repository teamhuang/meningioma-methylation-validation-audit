import pandas as pd, os
NG = r"E:\脑膜瘤npj\data\suppl\natgen2022"
c19 = pd.read_excel(os.path.join(NG, "41588_2022_1061_MOESM19_ESM.xlsx"))
c20 = pd.read_excel(os.path.join(NG, "41588_2022_1061_MOESM20_ESM.xlsx"))
m21 = pd.read_excel(os.path.join(NG, "41588_2022_1061_MOESM21_ESM.xlsx"), index_col=0)
print("MOESM19 shape:", c19.shape)
print("MOESM19 cols:", list(c19.columns))
print("MOESM20 shape:", c20.shape)
print("MOESM20 cols:", list(c20.columns))
print("MOESM21 shape:", m21.shape)
print("MOESM21 index head:", [str(x) for x in m21.index[:5]])
print("MOESM21 col head:", [str(x) for x in m21.columns[:5]])
