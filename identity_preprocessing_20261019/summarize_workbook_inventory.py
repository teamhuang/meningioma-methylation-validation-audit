import json
from pathlib import Path
p=Path(__file__).resolve().parent/'results/workbook_inventory.json'
data=json.loads(p.read_text(encoding='utf-8'))
for x in data:
    if '1061_' in x['file']:
        print(x['file'],x['sheet'],repr(x['preview'][:3])[:1600])
