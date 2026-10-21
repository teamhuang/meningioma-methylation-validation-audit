import gzip
p189 = r'E:\脑膜瘤npj\data\GSE189521_series_matrix.txt.gz'
with gzip.open(p189, 'rt', encoding='utf-8', errors='replace') as f:
    for line in f:
        if line.startswith('!Sample_title'):
            t = [x.strip().strip('"') for x in line.rstrip().split('\t')[1:]]
            print('titles[:6]:', t[:6])
            print('n titles:', len(t))
        elif line.startswith('!Sample_characteristics_ch1'):
            c = [x.strip().strip('"') for x in line.rstrip().split('\t')[1:]]
            print('char[:3]:', c[:3])
        elif line.startswith('!Sample_geo_accession'):
            g = [x.strip().strip('"') for x in line.rstrip().split('\t')[1:]]
            print('gsm[:3]:', g[:3])
