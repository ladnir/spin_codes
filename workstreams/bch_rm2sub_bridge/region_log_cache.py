"""Small write-once discovery cache; never substitutes for outward replay."""
from pathlib import Path
import math
import numpy as np
from flint import arb,ctx
import bridge as base
import tightened_occupancy as tight
import polynomial_regions as poly


def get(tilt):
    assert isinstance(tilt,int) and -120<=tilt<=20
    path=base.HERE/'generated'/f'region_logs_t128_s15_tilt_{tilt}.npy'
    receipt=path.with_suffix('.json')
    if receipt.exists():
        saved=base.read(receipt);assert base.sha(path)==saved['array_sha256']
        for name,digest in saved['local_sha256'].items():assert base.sha(base.HERE/name)==digest
        return np.load(path,allow_pickle=False)
    assert not path.exists()
    ctx.prec=256;t,s,spectrum=base.load_map(tight.NAME)
    lam=(arb(tilt)/10).exp()
    region=poly.regions(t,s,spectrum,tight.kernel_spectrum(),(-lam).exp(),8192)
    values=np.array([[float(v.log()) if v>0 else -math.inf for v in row] for row in region]).reshape(-1,3,3)
    with path.open('xb') as stream:np.save(stream,values,allow_pickle=False)
    paths=[Path(__file__),Path(poly.__file__),Path(tight.__file__),base.HERE/'general_occupancy.py',base.HERE/'activation_bridge.py',
           base.HERE/'inputs/t128_s15_selection.json',base.HERE/'inputs/t128_s15_a_spectrum.json',base.HERE/'inputs/t128_s15_b_kernel_spectrum.json']
    base.write_new(receipt,dict(status='DISCOVERY_LOG_CACHE_ONLY',tilt=tilt,precision_bits=256,array_sha256=base.sha(path),
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in paths}))
    return values
