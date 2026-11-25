"""Public-source downloads, with integrity checks and no synthetic fallbacks."""
from pathlib import Path
import hashlib
import json
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

ROOT = Path(__file__).resolve().parents[2]

def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def download_public(force=False, groups=None):
    manifest = json.loads((ROOT/'config/public_sources.json').read_text())
    if groups is not None:
        groups = set(groups)
        manifest = [s for s in manifest if s.get('group', s['path'].split('/')[0]) in groups]
        if not manifest:
            raise ValueError(f'No sources match groups: {sorted(groups)}')
    manifest_rows = []
    session = requests.Session()
    session.mount('https://', HTTPAdapter(max_retries=Retry(
        total=3, backoff_factor=1, status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=['GET'])))
    for source in manifest:
        dest = ROOT/'data/raw'/source['path']
        dest.parent.mkdir(parents=True, exist_ok=True)
        if source.get('pmcid') and not dest.exists():
            from .supplements import download_attachment
            download_attachment(source['pmcid'],source['attachment'],dest,source['bytes'],source['sha256'])
        if force or not dest.exists():
            tmp = dest.with_suffix(dest.suffix+'.partial')
            try:
                with session.get(source['url'], stream=True, timeout=(20, 180)) as response:
                    response.raise_for_status()
                    if 'text/html' in response.headers.get('Content-Type', ''):
                        raise ValueError(f"HTML received instead of data: {source['path']}")
                    with tmp.open('wb') as f:
                        for block in response.iter_content(1024*1024):
                            f.write(block)
                check_source(tmp, source)
                tmp.replace(dest)
            finally:
                tmp.unlink(missing_ok=True)
        check_source(dest, source)
        manifest_rows.append({**source, 'observed_sha256': sha256(dest), 'observed_bytes': dest.stat().st_size})
    output = ROOT/'results/download_manifest.json'
    previous = json.loads(output.read_text()) if output.exists() else []
    by_path = {s['path']: s for s in previous}
    by_path.update({s['path']: s for s in manifest_rows})
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(list(by_path.values()), indent=2)+'\n')
    return manifest_rows

def check_source(path, source):
    if path.stat().st_size == 0:
        raise ValueError(f'Empty input: {path}')
    with path.open('rb') as f:
        head = f.read(256).lstrip().lower()
    if head.startswith((b'<!doctype html', b'<html')):
        raise ValueError(f'HTML input: {path}')
    if source.get('sha256') and sha256(path) != source['sha256']:
        raise ValueError(f'Public source changed: {path}. Review before updating manifest.')
