"""Combine authenticated actual24 local maps with the existing wider outer.

This is a different construction from both the measured16-state encoder and
the small-outer actual24 proposal. All results are floating proposals only.
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
import wider_fractional as wf
import capped_gate as cg


def run(args):
    points = tuple(map(wf.point, args.points))
    shape = wf.wider.geometry(symbol_bits=32)
    shape['state_bits'] = 24
    if args.output.exists() or not points or any(q > 256 for q,_,_ in points):
        raise ValueError('fresh output and valid wider-outer occupancies required')
    pins = wf.sources('byte_native')
    for path in (Path(__file__).resolve(), Path(cg.__file__).resolve(),
                 Path(cg.variant_gate.__file__).resolve()):
        pins[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    caches = {}
    for theta in sorted({theta for _,theta,_ in points}):
        basename = f'theta_{theta.numerator}_{theta.denominator}.json'
        found = [d/basename for d in args.caches if (d/basename).exists()]
        if not found:
            raise ValueError(f'no authenticated local cache for theta={theta}')
        path = found[0].resolve()
        active, saved = cg.load_local(path)
        if Fraction(saved['tilt']) != theta:
            raise ArithmeticError('cache filename and tilt disagree')
        caches[theta] = active, saved
        pins.update(saved['source_sha256'])
        pins[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    records = [saved['map_record'] for _,saved in caches.values()]
    if any(record != records[0] for record in records):
        raise ValueError('all local caches must describe the same literal map')
    if not cg.variant_gate.checked_pins(pins):
        raise ArithmeticError('source pin check failed')
    envelope, _ = wf.wider.outer(shape)
    beta = envelope.beta
    beta_log = log(beta.numerator)-log(beta.denominator)
    result = dict(schema='packet8-iteration4-wider-actual24-fractional-proposal-1',
        geometry=shape, map_record=records[0], envelope=envelope.metadata(),
        proposal_only=True, whole_code_certificate=False, has_outward_endpoints=False,
        source_sha256=pins, source_pins_verified_at_finish=False, group_steps=4,
        subset_union_outside_fractional_power=True, trials=[])
    started = monotonic()
    for q,theta,alpha in points:
        active,_ = caches[theta]
        old = wf.birth.marked.potential_operators(active,0.)
        potential,_,algebra = wf.birth.rebase_potential(old)
        tuples = wf.birth.fine.tuple_products(potential,4)
        macro = wf.birth.fine.fine_operators(tuples,alpha)
        regional,backend = wf.birth.prior.placement(macro,q,epochs=8,windows=32,force_log=True)
        moment = wf.birth.prior.logarithmic.log_power_matrix(regional[q],64)
        margin = -wf.exponent(shape,q,theta,alpha,beta_log,moment)/log(2)
        alpha1 = wf.birth.fine.fine_operators(tuples,1)
        grouped,_ = wf.birth.prior.placement(alpha1,q,epochs=8,windows=32,force_log=True)
        physical,_ = wf.birth.prior.placement(potential,q,epochs=32,windows=8,force_log=True)
        finite = np.isfinite(physical)
        if not np.array_equal(finite,np.isfinite(grouped)):
            raise ArithmeticError('alpha1 structural zeros disagree')
        error = float(np.max(np.abs(grouped[finite]-physical[finite]),initial=0))
        if error > 2e-8:
            raise ArithmeticError('alpha1 geometry regression failed')
        result['trials'].append(dict(q=q,theta=str(theta),alpha=str(alpha),
            margin_bits=margin,log_moment=moment,backend=backend,algebra=algebra,
            alpha1_max_regional_log_error=error,powered_macro_operators=macro.tolist()))
        if not cg.variant_gate.checked_pins(pins):
            raise ArithmeticError('source changed during screen')
        result['elapsed_seconds'] = monotonic()-started
        args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
        print(f'actual24 symbol32 q={q} theta={theta} alpha={alpha}: '
              f'{margin:.9f} bits; alpha1 error={error:.3g}',flush=True)
    result['source_pins_verified_at_finish'] = True
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--caches',nargs='+',type=Path,default=[
        wf.HERE.parent/'iteration3'/f'cache24_v{i}' for i in (1,2,3)])
    parser.add_argument('--points',nargs='+',default=['8:.06:.4','32:.3:.35','60:.5:.4'])
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args())
