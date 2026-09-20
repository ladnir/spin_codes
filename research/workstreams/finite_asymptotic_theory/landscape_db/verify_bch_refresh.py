"""Independent positive-arithmetic replay of actual-size four-state Q1 kernels.

The reference uses 90 decimal digits, linear epoch iteration, and ordinary
positive polynomial arithmetic. It is intentionally separate from the batched
log-domain production implementation. This is not outward interval arithmetic.
"""
import csv
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np

import activation_q1_refresh as production

HERE=Path(__file__).resolve().parent


def multiply(a,b):
    return [[sum(a[i][k]*b[k][j] for k in range(4)) for j in range(4)] for i in range(4)]


def add(a,b):return [[x+y for x,y in zip(ar,br)] for ar,br in zip(a,b)]


def vector(v,m):
    return [v[0]*m[0][j]+v[1]*m[1][j]+v[2]*m[2][j]+v[3]*m[3][j] for j in range(4)]


def coefficients(t,s,counts,lam,epochs,block):
    z=mp.exp(-mp.mpf(float(lam)));m=mp.mpf((1<<s)-1);kappa=m/(m-1)
    d=min(counts)
    m0=sum(mp.mpf(n)*z**w for w,n in counts.items())/m
    m1=sum(mp.mpf(n)*(mp.mpf(w)*z**(w-1)+(t-w)*z**(w+1))/t for w,n in counts.items())/m
    zero=[[mp.mpf(0) for _ in range(4)] for _ in range(4)]
    one=[[mp.mpf(0) for _ in range(4)] for _ in range(4)]
    zero[0][0]=1;zero[1][2]=z**d;zero[2][2]=m0;zero[3][2]=min(z**d,kappa*m0)
    one[0][1]=z
    for state,r in ((1,z**(d-1)),(2,m1),(3,min(z**(d-1),kappa*m1))):
        one[state][0]=r/m;one[state][3]=r*(m-1)/m
    rz=[[mp.mpf(i==j) for j in range(4)] for i in range(4)]
    ra=[[mp.mpf(0) for _ in range(4)] for _ in range(4)]
    for _ in range(epochs):
        ra=add(multiply(ra,zero),multiply(rz,one));rz=multiply(rz,zero)
    ra=[[v/epochs for v in row] for row in ra]
    current=[[mp.mpf(1),mp.mpf(0),mp.mpf(0),mp.mpf(0)]]
    for _ in range(block):
        updated=[[mp.mpf(0) for _ in range(4)] for _ in range(len(current)+1)]
        for w,row in enumerate(current):
            a=vector(row,rz);b=vector(row,ra)
            for j in range(4):updated[w][j]+=a[j];updated[w+1][j]+=b[j]
        current=updated
    return [mp.log(sum(row)/math.comb(block,w)) for w,row in enumerate(current)]


def main():
    csv_path=HERE/'bch_growth.csv';study_path=HERE/'bch_growth.json'
    study=json.loads(study_path.read_text())
    digest=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
    if digest(csv_path)!=study['csv_sha256']:
        raise ValueError('growth CSV does not match study receipt')
    map_path=HERE/'activation_pilot_v1/maps/t64_s20.json'
    for path in (map_path,Path(production.__file__)):
        key=path.relative_to(HERE.parents[2]).as_posix()
        if digest(path)!=study['source_sha256'][key]:
            raise ValueError(f'growth source changed: {key}')
    rows=list(csv.DictReader(csv_path.open(newline='')))
    maps=json.loads(map_path.read_text())
    ac={w:n for w,n in enumerate(maps['a_counts']) if w and n};checks=[]
    with mp.workdps(90):
        for b in (8,32,64,128):
            row=next(r for r in rows if int(r['block_bits'])==b and r['step_bits']=='64'
                     and r['state_bits']=='20' and r['message_exponent']=='16')
            lam=float(np.exp(float(row['dominant_log_surprisal'])));epochs=int(row['epochs_per_region'])
            expected=coefficients(64,20,ac,lam,epochs,b)
            actual=production.coefficient_logs(*production.epoch_logs(64,20,ac,np.array([lam])),epochs,b)[0]
            errors=[abs(float(x)-y) for x,y in zip(expected,actual)]
            if max(errors)>2e-9:raise ArithmeticError('high-precision coefficient replay failed')
            check=dict(block_bits=b,epochs=epochs,coefficients_checked=b+1,maximum_log_error=max(errors))
            checks.append(check);print(check,flush=True)
    source=Path(__file__)
    (HERE/'bch_refresh_verification.json').write_text(json.dumps(dict(decimal_digits=90,checks=checks,
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        production_sha256=hashlib.sha256(Path(production.__file__).read_bytes()).hexdigest(),
        study_sha256=digest(study_path),csv_sha256=digest(csv_path),map_sha256=digest(map_path),
        note='Independent positive arithmetic, not outward certification.'),indent=2)+'\n')


if __name__=='__main__':main()
