import pandas as pd, numpy as np
p = r'E:\脑膜瘤npj\data\GSE270371_Matrix_processed.txt.gz'
df = pd.read_csv(p, sep='\t', index_col=0, nrows=5)
print('cols[:10]:', [str(c) for c in df.columns[:10]])
print('index[:5]:', [str(x) for x in df.index[:5]])
print('shape (nrows=5):', df.shape)
full_cols = pd.read_csv(p, sep='\t', index_col=0, nrows=0).columns
print('total cols:', len(full_cols))
int_cols = [c for c in full_cols if not str(c).startswith('Detection Pval')]
print('intensity cols:', len(int_cols), 'first/last:', int_cols[0], int_cols[-1])
dups = pd.Series(int_cols).duplicated().sum()
print('dup intensity col names:', dups)
