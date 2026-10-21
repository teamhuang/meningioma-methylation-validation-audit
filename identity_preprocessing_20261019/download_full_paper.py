from pathlib import Path
import requests
dest=Path(__file__).resolve().parent/'source_data/choudhurynatgenet2022.pdf'
url='https://open-proposals.ucsf.edu/sites/default/files/chv/cw-projects-2025/1709303946/choudhurynatgenet2022.pdf'
r=requests.get(url,timeout=60);r.raise_for_status();assert r.content.startswith(b'%PDF');dest.write_bytes(r.content);print('Full paper downloaded',len(r.content))
