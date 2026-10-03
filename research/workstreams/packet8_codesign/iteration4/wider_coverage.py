"""Replay every occupancy for the wider-outer actual24 floating proposal.

No interpolation between sampled q values is used. The single-active-group
term uses its exact expected outer support counts and ordinary local moments.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import log
from pathlib import Path
from time import monotonic

import numpy as np
from scipy.special import logsumexp
import cached_wider24 as cw

wf = cw.wf


def validate_receipt(receipt,shape,record):
    if (receipt.get('schema') != 'packet8-iteration4-wider-actual24-fractional-proposal-1'
            or receipt['geometry'] != shape or receipt['map_record'] != record
            or receipt['group_steps'] != 4
            or not receipt['source_pins_verified_at_finish']
            or not cw.cg.variant_gate.checked_pins(receipt['source_sha256'])):
        raise ValueError('fresh authenticated same-construction wider24 proposal required')


def run(args):
    if args.output.exists() or not args.receipts:
        raise ValueError('input receipts and fresh output required')
    shape = wf.wider.geometry(symbol_bits=32)
    shape['state_bits'] = 24
    result = dict(schema='packet8-iteration4-wider24-all-q-proposal-1',
        geometry=shape,proposal_only=True,whole_code_certificate=False,
        has_outward_endpoints=False,all_occupancies_checked=False,
        source_pins_verified_at_finish=False,trials=[])
    pins, choices, record = {}, {}, None
    for path in args.receipts:
        saved=json.loads(path.read_text(encoding='utf-8'))
        if record is None:
            record=saved['map_record']
        validate_receipt(saved,shape,record)
        for name,digest in saved['source_sha256'].items():
            if name in pins and pins[name] != digest:
                raise ArithmeticError('conflicting source hashes')
            pins[name]=digest
        pins[str(path.resolve())]=hashlib.sha256(path.read_bytes()).hexdigest()
        for trial in saved['trials']:
            theta,alpha=Fraction(trial['theta']),Fraction(trial['alpha'])
            choices.setdefault((theta,alpha),trial)
    for path in (Path(__file__).resolve(),wf.HERE/'test_wider_coverage.py'):
        pins[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
    result['map_record'],result['source_sha256']=record,pins
    envelope,counts=wf.wider.outer(shape)
    beta=envelope.beta
    beta_log=log(beta.numerator)-log(beta.denominator)
    best=np.full(257,-np.inf)
    best_choice={}
    started=monotonic()
    for (theta,alpha),saved in choices.items():
        macro=np.asarray(saved['powered_macro_operators'])
        if macro.shape != (33,10,10) or np.any(macro<0) or not np.isfinite(macro).all():
            raise ValueError('invalid powered macro operators')
        regional,backend=wf.birth.prior.placement(macro,256,epochs=8,windows=32,force_log=True)
        margins={}
        for q in range(2,257):
            moment=wf.birth.prior.logarithmic.log_power_matrix(regional[q],64)
            value=-wf.exponent(shape,q,theta,alpha,beta_log,moment)/log(2)
            margins[str(q)]=value
            if value>best[q]:
                best[q]=value
                best_choice[str(q)]=dict(theta=str(theta),alpha=str(alpha))
        discrepancy=margins[str(saved['q'])]-saved['margin_bits']
        if abs(discrepancy)>1e-7:
            raise ArithmeticError('full-range replay disagrees with saved point')
        result['trials'].append(dict(theta=str(theta),alpha=str(alpha),backend=backend,
            point_replay_difference_bits=discrepancy,margins=margins))
        print(f'all-q replay theta={theta}, alpha={alpha}; '
              f'current q>=2 minimum {min(best[2:]):.6f}bits',flush=True)

    # Recover authenticated active-label caches, not the thinned/powered family.
    caches={}
    for name in list(pins):
        path=Path(name)
        if path.name.startswith('theta_') and path.suffix=='.json':
            local,saved=cw.cg.load_local(path)
            if saved['map_record'] != record:
                raise ArithmeticError('q1 cache map does not match')
            caches[Fraction(saved['tilt'])]=local
    if not caches:
        raise ValueError('no authenticated active-label caches for q1')
    q1_support=np.zeros(65)
    for theta,local in sorted(caches.items()):
        regional,_=wf.birth.prior.placement(local,1,epochs=32,windows=8,force_log=True)
        logs=wf.birth.prior.q1_support_logs(regional,regions=64,tilt=theta,cutoff=13107)
        q1_support=np.minimum(q1_support,logs)
    count_logs=np.array([log(v.numerator)-log(v.denominator) if v else -np.inf for v in counts])
    q1_terms=log(256)+count_logs+q1_support
    best[1]=-float(logsumexp(q1_terms))/log(2)
    if not np.isfinite(best[1:]).all():
        raise ArithmeticError('uncovered or nonfinite occupancy bound')
    result.update(best_margins={str(q):float(best[q]) for q in range(1,257)},
        best_choices=best_choice,q1_exact_expected_shells=True,
        q1_support_log_uppers_floating=q1_support.tolist(),
        q1_margin_bits=float(best[1]),weakest_q=int(np.argmin(best[1:])+1),
        weakest_q_margin_bits=float(min(best[1:])),
        combined_margin_bits=-float(logsumexp(-best[1:]*log(2)))/log(2),
        all_occupancies_checked=True,elapsed_seconds=monotonic()-started)
    if not cw.cg.variant_gate.checked_pins(pins):
        raise ArithmeticError('source changed during all-q replay')
    result['source_pins_verified_at_finish']=True
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(f'combined floating margin={result["combined_margin_bits"]:.9f}bits; '
          f'q1={best[1]:.9f}; weakest q={result["weakest_q"]}',flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipts',nargs='+',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args())
