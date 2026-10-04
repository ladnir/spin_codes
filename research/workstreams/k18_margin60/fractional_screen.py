"""Floating proposal for clipping after potential-route conditioning.

No output of this script is a certificate. Inputs are physical active-label
comparison operators; zero labels are averaged before fractional powers.
"""
from __future__ import annotations

import os
for _name in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[_name] = '1'

import argparse
from fractions import Fraction as F
from math import comb, log
from pathlib import Path
import hashlib
import json
import sys
from time import monotonic

import numpy as np
from scipy.special import logsumexp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'k16_design'))
from packet_rs_length_proposal import scaled_placement
from packet_rs_s20_proposal import log_power


def thin(physical, packet_bits=4):
    physical = np.asarray(physical, dtype=float)
    if (physical.ndim != 3 or physical.shape[1] != physical.shape[2]
            or not np.isfinite(physical).all() or np.any(physical < 0)):
        raise ValueError('finite nonnegative square operator family required')
    labels = (1 << packet_bits)-1
    result = np.zeros_like(physical)
    for j in range(len(result)):
        for k in range(j+1):
            result[j] += float(F(comb(j,k)*labels**k, 1 << (packet_bits*j))) * physical[k]
    return result


def grouped(physical, tilt, alpha, steps=2):
    """Power each chronological tuple before averaging its placement."""
    if steps not in (1,2,4) or not 0 < alpha <= 1 or tilt <= 0:
        raise ValueError('steps1/2/4, positive tilt, alpha in(0,1] required')
    potential = thin(physical)
    windows = len(potential)-1
    a = ((1 + np.exp(-float(tilt)))/2)**4
    potential /= a**np.arange(windows+1)[:,None,None]
    products = potential.copy()
    totals = np.arange(windows+1)
    multiplicities = np.array([comb(windows,j) for j in totals], dtype=np.int64)
    choices = multiplicities.copy()
    size = potential.shape[1]
    for _ in range(1,steps):
        products = (products[:,None] @ potential[None]).reshape(-1,size,size)
        totals = (totals[:,None]+np.arange(windows+1)[None]).reshape(-1)
        multiplicities = (multiplicities[:,None]*choices[None]).reshape(-1)
    result = np.zeros((steps*windows+1,size,size))
    for j in range(len(result)):
        mask = totals == j
        weights = multiplicities[mask].astype(float)/comb(steps*windows,j)
        result[j] = np.sum(weights[:,None,None]*products[mask]**alpha, axis=0)
    if not np.isfinite(result).all() or np.any(result < 0):
        raise FloatingPointError('invalid grouped operator')
    return result, a


def score(physical, tilt, alpha, qs, steps=2):
    qs = tuple(qs)
    if not qs or min(qs)<1 or max(qs)>2048:
        raise ValueError('K18 occupancies in1..2048 required')
    with np.errstate(over='raise', invalid='raise', divide='raise', under='raise'):
        operators, a = grouped(physical, tilt, alpha, steps)
        windows = 16*steps
        regional = scaled_placement(operators, max(qs), epochs=2048//windows, windows=windows)
        log_beta = 256*log(2)-8*log(65535)
        result = {}
        for q in qs:
            moment = log_power(regional.matrices[q],64)+64*regional.log_scales[q]
            bound = log(comb(2048,q))+float(alpha)*(q*log_beta+float(tilt)*52428+64*q*log(a))+moment
            result[str(q)] = -bound/log(2)
    return result


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_physical(path):
    record = json.loads(Path(path).read_text(encoding='utf-8'))
    for name, expected in record['source_sha256'].items():
        if digest(name) != expected:
            raise ArithmeticError(f'pin mismatch: {name}')
    if record['map_record']['map_sha256'] != 'e2467b9d5de8403ccd8263f0f55c73df128a8a2dd4c4cc405bc3a8e057fc3581':
        raise ValueError('exact certified K18 maps required')
    raw = record['matrices']
    physical = np.array([[[float(F(int(n))*F(2)**int(e)) for n,e in row]
        for row in matrix] for matrix in raw])
    if physical.shape[0] != 17:
        raise ValueError('physical sixteen-slot family required')
    return physical, F(record['tilt']), record


def run(local, output, alphas, qs, steps):
    if output.exists():
        raise ValueError('fresh output path required')
    physical, tilt, source = load_physical(local)
    pins = dict(source['source_sha256'])
    pins[str(Path(local).resolve())] = digest(local)
    pins[str(Path(__file__).resolve())] = digest(__file__)
    trials=[]
    for alpha in alphas:
        start=monotonic()
        margins=score(physical,tilt,F(alpha),qs,steps)
        trials.append(dict(alpha=str(F(alpha)), margins=margins, seconds=monotonic()-start))
        print(json.dumps(trials[-1]),flush=True)
    best={str(q):max(t['margins'][str(q)] for t in trials) for q in qs}
    if any(digest(name)!=value for name,value in pins.items()):
        raise ArithmeticError('sources changed during screen')
    result=dict(schema='k18-potential-conditioned-fractional-proposal-1',
        proposal_only=True,whole_code_certificate=False,local_receipt=str(local),
        tilt=str(tilt),steps_per_macro=steps,activity=source['local_selection_activity'],
        selected_occupancies=list(qs),trials=trials,best_margins=best,
        selected_union_margin_display=-float(logsumexp([-v*log(2) for v in best.values()]))/log(2),
        source_sha256=pins)
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--local',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--alphas',nargs='+',default=['1','.95','.9','.8'])
    parser.add_argument('--qs',type=int,nargs='+',default=[74,80,89,95,100,110,115])
    parser.add_argument('--steps',type=int,choices=[1,2,4],default=2)
    args=parser.parse_args()
    run(args.local,args.output,args.alphas,args.qs,args.steps)
