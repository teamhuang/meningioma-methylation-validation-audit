"""Verify references 6-12 of manuscript_npjDM_v5.md against Crossref.

Each reference is looked up by bibliographic query; the top hit is compared with the
citation as printed, and the verified metadata (volume, issue, page, DOI) is returned
for correction. Requests are rate-limited.
"""
from pathlib import Path
import json
import time
import urllib.parse
import urllib.request

QUERIES = [
    (6, 'Bayley meningioma three biological subtypes Science Advances 2022',
     'Multiple approaches converge on three biological subtypes of meningioma and extract new insights from published studies'),
    (7, 'Duckett clinically applicable methylation classifier recurrence risk meningioma Acta Neuropathologica Communications 2025',
     'Development, validation, and utility of a clinically applicable methylation classifier for recurrence risk prediction in meningiomas'),
    (8, 'Choudhury meningioma DNA methylation groups biological drivers therapeutic vulnerabilities Nature Genetics 2022',
     'Meningioma DNA methylation groups identify biological drivers and therapeutic vulnerabilities'),
    (9, 'Wang molecular classification surgical radiotherapeutic decision-making meningioma Nature Medicine 2024',
     'Molecular classification to refine surgical and radiotherapeutic decision-making in meningioma'),
    (10, 'Merk DNA methylation long-range epigenetic silencing protocadherins meningioma progression Nature Communications',
     'DNA methylation profiling identifies long-range epigenetic silencing of clustered protocadherins as a key determinant of meningioma progression'),
    (11, 'Nguyen longitudinal multi-omics profiling meningioma actionable biology aggressive disease Nature Communications 2025',
     'Longitudinal multi-omics profiling of meningioma reveals actionable biology of aggressive disease'),
    (12, 'Riley evaluation of clinical prediction models part 2 how to undertake an external validation study BMJ 2024',
     'Evaluation of clinical prediction models (part 2): how to undertake an external validation study'),
]

HDRS = {'User-Agent': 'npjDM-manuscript-refcheck/1.0 (mailto:corresponding.author@example.org)'}
out = []
for num, q, expect in QUERIES:
    url = ('https://api.crossref.org/works?query.bibliographic='
           + urllib.parse.quote(q) +
           '&rows=1&select=title,container-title,volume,issue,page,DOI,issued,author')
    rec = {'ref': num, 'query': q, 'expected_title': expect}
    try:
        req = urllib.request.Request(url, headers=HDRS)
        with urllib.request.urlopen(req, timeout=40) as r:
            data = json.load(r)
        it = data['message']['items'][0]
        title = (it.get('title') or [''])[0]
        rec.update({
            'title': title,
            'container': (it.get('container-title') or [''])[0],
            'volume': it.get('volume', ''),
            'issue': it.get('issue', ''),
            'page': it.get('page', ''),
            'doi': it.get('DOI', ''),
            'year': (it.get('issued', {}).get('date-parts', [[None]])[0][0]),
            'first_author': (it.get('author') or [{}])[0].get('family', ''),
            'title_match': title.strip().lower().replace(' ', '')[:40] ==
                           expect.strip().lower().replace(' ', '')[:40],
        })
    except Exception as e:
        rec['error'] = str(e)
    out.append(rec)
    print(f"[{num}] {rec.get('first_author','?')} {rec.get('year','?')} | "
          f"{rec.get('container','?')} {rec.get('volume','')}({rec.get('issue','')}) "
          f"{rec.get('page','')} | {rec.get('doi','ERR')} | title_match={rec.get('title_match')}",
          flush=True)
    time.sleep(3)

Path(r'E:\脑膜瘤npj\results_v4\ref6_12_crossref_verification.json').write_text(
    json.dumps(out, indent=2, ensure_ascii=False), encoding='utf-8')
print('\nsaved results_v4/ref6_12_crossref_verification.json')
