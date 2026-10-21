import json
from pathlib import Path
OUT=Path(__file__).resolve().parent/'results'
p=OUT/'protocadherin_source_inventory.json'
if p.exists():
    data=json.loads(p.read_text(encoding='utf-8'))
    for x in data:
        print(x['sheet'],repr(x['preview'][:4])[:1600])
