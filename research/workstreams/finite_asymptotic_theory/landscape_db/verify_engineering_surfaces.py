"""Independent high-precision checks at selected RM and random witnesses."""
import hashlib
import json
from pathlib import Path

import mpmath as mp
import numpy as np

import activation_q1_refresh as transfer
import activation_refresh_native as native
import verify_bch_refresh as reference
from report_engineering_surfaces import read,select

HERE=Path(__file__).resolve().parent


def main():
    rows,payload=read();root=HERE.parents[2]
    for name,digest in payload['source_sha256'].items():
        if hashlib.sha256((root/name).read_bytes()).hexdigest()!=digest:
            raise ValueError(f'study dependency changed: {name}')
    path=HERE/'activation_pilot_v1/maps/t64_s20.json'
    if hashlib.sha256(path.read_bytes()).hexdigest()!=payload['source_sha256'][path.relative_to(root).as_posix()]:
        raise ValueError('inner map changed')
    for path_source in (Path(native.__file__),native.LIBRARY,Path(transfer.__file__)):
        if hashlib.sha256(path_source.read_bytes()).hexdigest()!=payload['source_sha256'][path_source.relative_to(root).as_posix()]:
            raise ValueError('production source changed')
    ac={w:n for w,n in enumerate(json.loads(path.read_text())['a_counts']) if w and n}
    kernel=native.RefreshKernel();checks=[]
    with mp.workdps(90):
        for family,b in (('rm',128),('rm',512),('random_mean',512)):
            row=select(rows,family=family,block_bits=b,step_bits=64,state_bits=20,message_exponent=16)[0]
            lam=float(np.exp(row['dominant_log_tilt']));epochs=row['epochs_per_region']
            expected=reference.coefficients(64,20,ac,lam,epochs,b)
            actual=kernel.coefficients(*transfer.epoch_logs(64,20,ac,np.array([lam])),epochs,b)[0]
            error=max(abs(float(x)-y) for x,y in zip(expected,actual))
            if error>2e-9: raise ArithmeticError('high-precision replay failed')
            check=dict(family=family,block_bits=b,epochs=epochs,coefficients_checked=b+1,maximum_log_error=float(error))
            checks.append(check);print(check,flush=True)
    sources=(Path(__file__),Path(reference.__file__),Path(native.__file__),native.LIBRARY,Path(transfer.__file__),
             HERE/'engineering_surfaces.csv',HERE/'engineering_surfaces.json',path)
    result=dict(decimal_digits=90,checks=checks,study_dependencies_checked=len(payload['source_sha256']),
                source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
                note='Independent positive arithmetic at selected witnesses, not outward certification.')
    (HERE/'engineering_surface_verification.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__': main()
