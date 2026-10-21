from pathlib import Path
import requests,json
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results'
s=requests.Session();s.headers['Accept']='application/vnd.github+json'
for user in ['william-c-chen','jbayley']:
    r=s.get('https://api.github.com/users/'+user+'/repos',params={'per_page':100},timeout=40);r.raise_for_status()
    repos=r.json();(OUT/(user+'_repos.json')).write_text(json.dumps(repos,indent=2),encoding='utf-8')
    print(user,[(x['name'],x['default_branch']) for x in repos],flush=True)
    for repo in repos:
        if not any(x in repo['name'].lower() for x in ['mening','methyl']):continue
        r=s.get('https://api.github.com/repos/'+user+'/'+repo['name']+'/git/trees/'+repo['default_branch'],params={'recursive':1},timeout=40);r.raise_for_status()
        tree=r.json();(OUT/(user+'_'+repo['name']+'_tree.json')).write_text(json.dumps(tree,indent=2),encoding='utf-8')
        print(repo['name'],[(x['path'],x.get('size')) for x in tree.get('tree',[]) if x['type']=='blob'][:150],flush=True)
