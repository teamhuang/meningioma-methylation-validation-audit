"""Verify reference 11 by DOI and check four additional candidate references that the
Crossref record for Wang 2024 (doi 10.1038/s41591-024-03167-4) cites."""
import json
import time
import urllib.request

DOIS = [
    ('11', '10.1038/s41467-025-61063-y', 'Nguyen, longitudinal multi-omics meningioma'),
    ('new1', '10.1038/s41586-021-03850-3', 'Nassiri, clinically applicable integrative molecular classification of meningiomas'),
    ('new2', '10.1093/neuonc/noab213', 'Driver, a molecularly integrated grade for meningioma'),
    ('new3', '10.1111/bpa.13132', 'Hielscher, clinical implementation of integrated molecular-morphologic risk prediction for meningioma'),
    ('new4', '10.1016/S1470-2045(17)30155-9', 'Sahm, DNA methylation-based classification and grading system for meningioma'),
    ('new5', '10.1093/neuonc/noac224', 'Choudhury, hypermitotic meningiomas harbor DNA methylation subgroups'),
]
HDRS = {'User-Agent': 'npjDM-manuscript-refcheck/1.0 (mailto:corresponding.author@example.org)'}
out = []
for tag, doi, label in DOIS:
    try:
        req = urllib.request.Request(f'https://api.crossref.org/works/{doi}', headers=HDRS)
        with urllib.request.urlopen(req, timeout=40) as r:
            m = json.load(r)['message']
        rec = {'tag': tag, 'doi': doi, 'label': label,
               'title': (m.get('title') or [''])[0],
               'container': (m.get('container-title') or [''])[0],
               'volume': m.get('volume', ''), 'issue': m.get('issue', ''),
               'page': m.get('page', ''),
               'year': m.get('issued', {}).get('date-parts', [[None]])[0][0],
               'authors': [a.get('family', '') for a in (m.get('author') or [])[:3]],
               'status': 'verified'}
    except Exception as e:
        rec = {'tag': tag, 'doi': doi, 'label': label, 'status': 'ERROR', 'error': str(e)}
    out.append(rec)
    print(f"[{tag}] {rec.get('status')} | {rec.get('authors','?')} {rec.get('year','?')} | "
          f"{rec.get('container','?')} {rec.get('volume','')}({rec.get('issue','')}) "
          f"{rec.get('page','')} | {rec.get('doi')}", flush=True)
    time.sleep(3)

with open(r'E:\脑膜瘤npj\results_v4\ref_verification_extra.json', 'w', encoding='utf-8') as f:
    json.dump(out, f, indent=2, ensure_ascii=False)
print('\nsaved results_v4/ref_verification_extra.json')
