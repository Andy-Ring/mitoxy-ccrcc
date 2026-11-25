"""Retrieve public study attachments from the documented article and Europe PMC."""
from pathlib import Path
import io,zipfile,requests
from .data import check_source

def download_attachment(pmcid,member_name,destination,expected_size,expected_sha256):
    destination=Path(destination);spec={'sha256':expected_sha256}
    if destination.exists():
        check_source(destination,spec);return destination
    urls=[f'https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/supplementaryFiles',
          f'https://pmc.ncbi.nlm.nih.gov/articles/instance/{pmcid.removeprefix("PMC")}/bin/{member_name}']
    errors=[]
    for url in urls:
        tmp=destination.with_suffix(destination.suffix+'.partial')
        try:
            response=requests.get(url,timeout=(20,180));response.raise_for_status();payload=response.content
            if url.endswith('supplementaryFiles'):
                with zipfile.ZipFile(io.BytesIO(payload)) as z:
                    names=[n for n in z.namelist() if n.rsplit('/',1)[-1]==member_name]
                    if len(names)!=1:raise ValueError('Attachment not uniquely identified')
                    payload=z.read(names[0])
            destination.parent.mkdir(parents=True,exist_ok=True);tmp.write_bytes(payload)
            check_source(tmp,spec)
            if len(payload)!=expected_size:raise ValueError('Source-size mismatch')
            tmp.replace(destination);return destination
        except Exception as exc:errors.append(f'{url}: {exc}')
        finally:tmp.unlink(missing_ok=True)
    raise RuntimeError('Public attachment could not be retrieved. Supply the checksum-matching public file at '+str(destination)+'\n'+'\n'.join(errors))
