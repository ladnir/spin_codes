"""Probe new scalar means using exactly one freshly checked saved comparison.

No mixture search or coefficient pruning occurs here. Each rational witness
and directed endpoint is saved separately. These points do not cover the
intervals between them and do not certify the whole code.
"""
import argparse
import copy
import hashlib
import json
from fractions import Fraction as Q
from pathlib import Path

from flint import arb,ctx
from shared_relaxed_strategy import build_model,validate_record
from shared_relaxed_dense import save
from shared_relaxed_parallel import retarget_record


def prepare_record(record,source_sha256,distance=None):
    """Copy one claim, optionally change its exact cutoff, and clear results."""
    validate_record(record)
    if len(record['results'])!=1:
        raise ValueError('one unambiguous distance claim required')
    result=(copy.deepcopy(record) if distance is None
            else retarget_record(record,distance,0,source_sha256))
    row=result['results'][0]
    row.pop('cover',None)
    row.pop('rechecked_precision',None)
    row['probes']=[]
    result['max_cells']=0
    result['screen_only']=True
    result['proof_status']='Fresh selected-point bounds only; not a complete certificate.'
    validate_record(result)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--expected-sha256',required=True)
    parser.add_argument('--means',nargs='+',required=True)
    parser.add_argument('--retarget-distance',
        help='Optional exact distance; retain the saved comparison but rebuild all point bounds')
    parser.add_argument('--precision',type=int,default=256)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():parser.error('refusing to overwrite a receipt')
    means=list(map(Q,args.means))
    if args.precision<256 or any(not 0<=x<=1 for x in means):
        parser.error('precision>=256 and unit-interval means required')
    raw=args.source.read_bytes();digest=hashlib.sha256(raw).hexdigest()
    if digest!=args.expected_sha256.lower():parser.error('source hash does not match')
    record=prepare_record(json.loads(raw),digest,args.retarget_distance)
    model=build_model(record,args.precision)
    model.proposal_stop_bits=record['target_bits']+2
    record['source_record']=dict(path=str(args.source),sha256=digest,
                                scope='Exact saved comparison; no coefficient search or pruning.')
    record['precision']=args.precision
    row=record['results'][0]
    for mean in means:
        cell=(mean,mean)
        if model.empty(cell):
            result=dict(mean=str(mean),empty=True)
        else:
            score,witness=model.proposal(cell)
            ctx.prec=args.precision
            upper=model.outward(cell,witness)
            if ctx.prec!=args.precision or not upper>0:
                raise ArithmeticError('positive fresh bound at the requested precision required')
            result=dict(mean=str(mean),proposal=score,witness=witness,
                upper=[int(v) for v in upper.upper().man_exp()],
                log2_upper=str(upper.log()/arb(2).log()))
        row['probes'].append(result)
        if hashlib.sha256(args.source.read_bytes()).hexdigest()!=digest:
            raise RuntimeError('source record changed during fresh probes')
        save(args.output,record)
        print('SAVED COMPARISON POINT','distance',row['distance'],str(mean),
              result.get('log2_upper','empty'),flush=True)


if __name__=='__main__':main()
