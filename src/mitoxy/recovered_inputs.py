"""Public input preparation for recovered original notebooks."""
import csv
import json
import shutil
import requests
from .data import ROOT,sha256

def prepare(refresh_gdc=False, include_gdc=False):
    rows=list(csv.DictReader((ROOT/'config/original_public_source_manifest.csv').open()))
    sources=json.loads((ROOT/'config/public_sources.json').read_text())
    by_sha={s['sha256']:ROOT/'data/raw'/s['path'] for s in sources if s.get('sha256')}
    for r in rows:
        src=by_sha.get(r['sha256']);dest=ROOT/r['relative_path']
        if src and src.exists() and not dest.exists():
            dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dest)
    raw=ROOT/'data/raw/tcga_kirc_public';raw.mkdir(parents=True,exist_ok=True)
    p=raw/'gdc_clinical_query_response.json'
    if include_gdc and (refresh_gdc or not p.exists()):
        query={'filters':json.dumps({'op':'in','content':{'field':'project.project_id','value':['TCGA-KIRC']}}),
          'format':'JSON','size':2000,'expand':'diagnoses,demographic',
          'fields':'submitter_id,case_id,primary_site,disease_type'}
        response=requests.get('https://api.gdc.cancer.gov/cases',params=query,timeout=(30,180))
        response.raise_for_status();obj=response.json()
        hits=obj.get('data',{}).get('hits',[]);pagination=obj.get('data',{}).get('pagination',{})
        if not hits or pagination.get('total',len(hits))!=len(hits):raise ValueError('Empty or truncated GDC clinical response')
        p.write_text(json.dumps(obj,indent=2)+'\n')
        (raw/'gdc_clinical_query_payload.json').write_text(json.dumps(query,indent=2)+'\n')
    for r in rows:
        if r['role']=='clinical API response' and p.exists():
            r.update(size_bytes=p.stat().st_size,sha256=sha256(p))
    dest=ROOT/'data/provenance/paper_public_source_manifest.csv';dest.parent.mkdir(parents=True,exist_ok=True)
    with dest.open('w') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
    return rows
