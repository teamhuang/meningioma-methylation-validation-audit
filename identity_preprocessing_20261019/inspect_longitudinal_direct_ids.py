from pathlib import Path
import json
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results'
data=json.loads((OUT/'protocadherin_source_inventory_fast.json').read_text(encoding='utf-8'))
for x in data:
    if '14d' in x['sheet'] or '15d' in x['sheet'] or x['sheet']=='Fig. 5f':
        print(x['sheet'])
        for row in x['preview']:print(row)
