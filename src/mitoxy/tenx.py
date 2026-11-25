"""Read public 10x matrices and spatial coordinates without extracting archives."""
import io
import tarfile
import h5py
import numpy as np
import pandas as pd
from scipy import sparse


def read_counts(path):
    """Return cells × genes CSR, barcodes, symbols and stable gene identifiers."""
    decode = lambda a: [v.decode() if isinstance(v, bytes) else str(v) for v in a]
    with h5py.File(path, 'r') as f:
        if 'matrix' in f:
            g = f['matrix']; names = g['features']['name'][:]; ids = g['features']['id'][:]
            mask = (np.asarray(decode(g['features']['feature_type'][:])) == 'Gene Expression'
                    if 'feature_type' in g['features'] else np.ones(len(names), dtype=bool))
        else:
            candidates = [k for k in f if isinstance(f[k], h5py.Group) and 'barcodes' in f[k]]
            if len(candidates) != 1:
                raise ValueError(f'Ambiguous legacy genome groups: {candidates}')
            g = f[candidates[0]]; names = g['gene_names'][:]; ids = g['genes'][:]
            mask = np.ones(len(names), dtype=bool)
        x = sparse.csc_matrix((g['data'][:], g['indices'][:], g['indptr'][:]), shape=tuple(g['shape'][:])).T.tocsr()
        bars = pd.Index(decode(g['barcodes'][:]), name='barcode')
    if not bars.is_unique or np.any(x.data < 0) or not np.isfinite(x.data).all():
        raise ValueError(f'Invalid counts or duplicated barcodes: {path}')
    if not np.equal(x.data, np.rint(x.data)).all():
        raise ValueError(f'Noninteger UMI counts: {path}')
    return x[:, mask], bars, np.asarray(decode(names))[mask], np.asarray(decode(ids))[mask]


def positions(path):
    with tarfile.open(path, 'r:gz') as tf:
        members = [m for m in tf.getmembers() if m.isfile() and
                   m.name.rsplit('/', 1)[-1] in ['tissue_positions_list.csv', 'tissue_positions.csv']]
        if len(members) != 1:
            raise ValueError(f'Expected one coordinate file in {path}')
        payload = tf.extractfile(members[0]).read()
    if payload.split(b',', 1)[0].strip().lower() == b'barcode':
        p = pd.read_csv(io.BytesIO(payload)).rename(columns={
            'pxl_row_in_fullres':'pxl_row', 'pxl_col_in_fullres':'pxl_col'})
    else:
        p = pd.read_csv(io.BytesIO(payload), header=None, names=[
            'barcode','in_tissue','array_row','array_col','pxl_row','pxl_col'])
    if not p.barcode.is_unique:
        raise ValueError('Duplicated coordinate barcodes')
    return p.set_index('barcode')
