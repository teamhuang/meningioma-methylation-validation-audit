import gzip
p = r'E:\脑膜瘤npj\data\GSE270371_Matrix_processed.txt.gz'
targets = {'cg19311470'}
found = None
with gzip.open(p, 'rt') as f:
    hdr = f.readline()
    for i, line in enumerate(f, 1):
        if i > 1500000:
            break
        if line.split('\t', 1)[0] in targets:
            found = i
            break
print('cg19311470 at data row:', found)
