from pathlib import Path
from xml.etree import ElementTree as ET
from run_strengthening import OUT
for file in ['PMC12750935.xml','PMC13526022.xml','PMC11313601.xml']:
    root=ET.parse(OUT/file).getroot()
    hits=[]
    words=['random forest','centroid','support vector','local freedom','108 patients','patient identifier','213','217','clinical classifier']
    if file=='PMC12750935.xml':words+=['classifier','training','github']
    for el in root.iter('p'):
        s=''.join(el.itertext())
        if any(x in s.lower() for x in words):hits.append(s)
    print('\nFILE',file)
    for s in hits: print(s[:4500])
    (OUT/(file+'.extract.txt')).write_text('\n\n'.join(hits),encoding='utf-8')
