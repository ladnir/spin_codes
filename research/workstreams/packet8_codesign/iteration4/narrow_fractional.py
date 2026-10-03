"""Revisit the literal t32/s16 maps using fractional routing proposals.

The byte route and small outer stay unchanged. Each physical step processes
four bytes, with a fresh nonzero-transitive state update per physical step.
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
import long_grouped_gate as long_gate


def geometry():
    shape=wf.wider.geometry(symbol_bits=16)
    shape.update(physical_t=32,physical_packet_slots=4,
        physical_steps_per_region=128,physical_steps_total=4096)
    return shape


def source_pins():
    pins=wf.sources('byte_native_t32')
    pins.update(long_gate.source_pins())
    for path in (Path(__file__).resolve(),wf.HERE/'test_narrow_fractional.py'):
        pins[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
    return pins


def run(args):
    points=tuple(map(wf.point,args.points))
    if args.output.exists() or any(q>512 for q,_,_ in points):
        raise ValueError('valid points and fresh output required')
    pins,started=source_pins(),monotonic()
    data,record=wf.birth.maps.prepare('byte_native_t32')
    if (data['bits'],data['windows'],data['width']) != (16,4,32):
        raise ArithmeticError('literal t32/s16 map required')
    shape=geometry()
    envelope,_=wf.wider.outer(shape)
    beta=envelope.beta
    beta_log=log(beta.numerator)-log(beta.denominator)
    result=dict(schema='packet8-iteration4-t32-fractional-proposal-1',
        geometry=shape,map_record=record,group_steps=args.group,
        source_sha256=pins,source_pins_verified_at_finish=False,
        proposal_only=True,whole_code_certificate=False,has_outward_endpoints=False,
        subset_union_outside_fractional_power=True,trials=[])
    for q,theta,alpha in points:
        weighted,emission,diag=wf.birth.prior.moments(data,np.exp(-float(theta)))
        active=wf.birth.prior.operators(weighted,emission)
        old=wf.birth.marked.potential_operators(active,0.)
        potential,_,algebra=wf.birth.rebase_potential(old)
        if args.group==4:
            tuples=wf.birth.fine.tuple_products(potential,4)
            macro=wf.birth.fine.fine_operators(tuples,alpha)
            alpha1=wf.birth.fine.fine_operators(tuples,1)
            grouped,_=wf.birth.prior.placement(alpha1,q,epochs=32,windows=16,force_log=True)
            physical,_=wf.birth.prior.placement(potential,q,epochs=128,windows=4,force_log=True)
            finite=np.isfinite(physical)
            if not np.array_equal(finite,np.isfinite(grouped)):
                raise ArithmeticError('alpha1 structural zeros disagree')
            error=float(np.max(np.abs(grouped[finite]-physical[finite]),initial=0))
            if error>2e-8:
                raise ArithmeticError('alpha1 regression failed')
            group_diag=dict(alpha1_max_regional_log_error=error)
        else:
            macro,group_diag=long_gate.long_operators(potential,8,alpha,memory_mib=64)
        regional,backend=wf.birth.prior.placement(macro,q,
            epochs=128//args.group,windows=4*args.group,force_log=True)
        moment=wf.birth.prior.logarithmic.log_power_matrix(regional[q],32)
        margin=-wf.exponent(shape,q,theta,alpha,beta_log,moment)/log(2)
        result['trials'].append(dict(q=q,theta=str(theta),alpha=str(alpha),
            margin_bits=margin,log_moment=moment,backend=backend,algebra=algebra,
            group_diagnostics=group_diag,local_diagnostics=diag,
            powered_macro_operators=macro.tolist()))
        if source_pins()!=pins:
            raise ArithmeticError('source changed during proposal')
        result['elapsed_seconds']=monotonic()-started
        args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
        print(f't32 g{args.group} q={q} theta={theta} alpha={alpha}: {margin:.9f}bits',flush=True)
    result['source_pins_verified_at_finish']=True
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--points',nargs='+',default=['16:.06:.4','64:.25:.4','119:.475:.35'])
    parser.add_argument('--group',type=int,choices=(4,8),default=4)
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args())
