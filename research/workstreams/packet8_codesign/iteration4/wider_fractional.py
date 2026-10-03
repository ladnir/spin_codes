"""Revisit the existing wider outer using the fractional routing bound.

Both physical local families and the larger-outer geometry are imported from
frozen sources. No setup law is approximated by editing metadata alone.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import comb, log
from pathlib import Path
import sys
from time import monotonic

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'iteration3'))
import potential_birth_gate as birth
sys.path.insert(0, str(HERE.parent/'wider_outer'))
import wider_screen as wider
import numpy as np


def point(value):
    q, theta, alpha = value.split(':')
    q, theta, alpha = int(q), Fraction(theta), Fraction(alpha)
    if q < 2 or theta <= 0 or not 0 < alpha <= 1:
        raise ValueError('q>=2, positive theta and 0<alpha<=1 required')
    return q, theta, alpha


def exponent(shape, q, theta, alpha, beta_log, moment):
    return (log(comb(shape['outer_groups'],q))
            +float(alpha)*(q*beta_log+float(theta)*shape['cutoff'])+moment)


def sources(name):
    pins = birth.source_pins()
    pins.update(wider.source_pins())
    pins.update(birth.prior.sources(name))
    for p in (Path(__file__).resolve(), HERE/'test_wider_fractional.py'):
        pins[str(p)] = hashlib.sha256(p.read_bytes()).hexdigest()
    return pins


def run(args):
    shape = wider.geometry(symbol_bits=args.symbol_bits)
    points = tuple(map(point,args.points))
    if args.output.exists() or not points or any(q>shape['outer_groups'] for q,_,_ in points):
        raise ValueError('fresh output and valid points required')
    pins, started = sources(args.map), monotonic()
    data, record = birth.maps.prepare(args.map)
    envelope, _ = wider.outer(shape)
    beta = envelope.beta
    beta_log = log(beta.numerator)-log(beta.denominator)
    result = dict(schema='packet8-iteration4-wider-fractional-proposal-1',
        geometry=shape, map_record=record, envelope=envelope.metadata(),
        proposal_only=True, whole_code_certificate=False, has_outward_endpoints=False,
        group_steps=4, source_sha256=pins, source_pins_verified_at_finish=False,
        subset_union_outside_fractional_power=True, trials=[])
    for q,theta,alpha in points:
        weighted, emission, diagnostics = birth.prior.moments(data,np.exp(-float(theta)))
        active = birth.prior.operators(weighted,emission)
        old = birth.marked.potential_operators(active,0.)
        potential,change,algebra = birth.rebase_potential(old)
        tuples = birth.fine.tuple_products(potential,4)
        macro = birth.fine.fine_operators(tuples,alpha)
        epochs = shape['physical_steps_per_region']
        regional, backend = birth.prior.placement(macro,q,epochs=epochs//4,windows=32,force_log=True)
        moment = birth.prior.logarithmic.log_power_matrix(regional[q],shape['regions'])
        margin = -exponent(shape,q,theta,alpha,beta_log,moment)/log(2)
        alpha1 = birth.fine.fine_operators(tuples,1)
        grouped,_ = birth.prior.placement(alpha1,q,epochs=epochs//4,windows=32,force_log=True)
        physical,_ = birth.prior.placement(potential,q,epochs=epochs,windows=8,force_log=True)
        finite = np.isfinite(physical)
        if not np.array_equal(finite,np.isfinite(grouped)):
            raise ArithmeticError('alpha1 structural zeros differ')
        error = float(np.max(np.abs(grouped[finite]-physical[finite]),initial=0))
        if error > 2e-8:
            raise ArithmeticError('alpha1 geometry regression failed')
        result['trials'].append(dict(q=q,theta=str(theta),alpha=str(alpha),margin_bits=margin,
            log_moment=moment,alpha1_max_regional_log_error=error,backend=backend,
            algebra=algebra,local_diagnostics=diagnostics,potential_operators=potential.tolist(),
            powered_macro_operators=macro.tolist()))
        if sources(args.map) != pins:
            raise ArithmeticError('source changed during screen')
        result['elapsed_seconds']=monotonic()-started
        args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
        print(f'symbol{args.symbol_bits} q={q} theta={theta} alpha={alpha}: '
              f'{margin:.9f} bits; alpha1 log error={error:.3g}',flush=True)
    result['source_pins_verified_at_finish']=True
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--symbol-bits',type=int,choices=(16,32),default=32)
    parser.add_argument('--map',choices=('byte_native','byte_native_A_scaled'),default='byte_native')
    parser.add_argument('--points',nargs='+',default=['8:.06:.4','32:.25:.4','60:.475:.35'])
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args())
