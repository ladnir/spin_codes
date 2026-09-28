"""Reconstruct and replay specified complete two-bit occupancy covers.

Unlisted occupancies remain open. A full-code claim requires every q=1..4096
and the final outward aggregate, including the separately verified q=1 class.
"""
import argparse
import json
from fractions import Fraction as Q
from pathlib import Path

import numpy as np
from flint import arb,ctx

import model
from cover import complete_cover,replay_cover,check_partition,DENOMINATOR
from bch_joint_support import authenticated_caps
from support import weighted_cdf_upper,weighted_union_shells
from mixture_cover import resolve_threshold,checked_threshold


ENSEMBLE_PREFIX=[1048576,2097152,4096,256,64,64,128,19]


def saved_updates(record):
    ensemble=record.get('ensemble')
    if (record.get('schema')!='two-bit-support-cover-1' or not isinstance(ensemble,list)
            or len(ensemble)!=9 or ensemble[:8]!=ENSEMBLE_PREFIX
            or any(type(x) is not int for x in ensemble)):
        raise ValueError('unsupported witness ensemble')
    return model.checked_updates(ensemble[-1])


def saved_threshold(record):
    saved_updates(record)
    return checked_threshold(record.get('threshold',model.THRESHOLD))


def validate_record(record):
    """Check saved scope and the exact support partition before costly replay."""
    updates=saved_updates(record);threshold=saved_threshold(record)
    q=record['q'];cutoff=record.get('cutoff',8);degree=record.get('output_degree')
    penalty=Q(record['penalty'])
    if (type(q) is not int or not 2<=q<=model.GROUPS
            or type(cutoff) is not int or not 1<=cutoff<=10
            or (degree is not None and (type(degree) is not int or not 0<=degree<=64))
            or not 0<penalty<=1 or type(record.get('density_fold',False)) is not bool):
        raise ValueError('valid sparse occupancy and local-operator settings required')
    check_partition(record['leaves'],q)
    for leaf in record['leaves']:
        if (Q(leaf['tilt'])<=0 or len(leaf['numerators'])!=q
                or any(type(p) is not int or not 0<p<DENOMINATOR for p in leaf['numerators'])):
            raise ValueError('positive tilt and complete interior support witnesses required')
    return updates,threshold,q,cutoff,degree,penalty


def build_operators(tilts,maximum,cutoff,penalty,degree,updates,*,floating=True):
    data=model.census(min(cutoff,maximum));operators={}
    model.prepare_density(data,tilts)
    for tilt in tilts:
        ops=model.epoch_operators(data,tilt,penalty,maximum=min(maximum,64),output_degree=degree,updates=updates)
        regions=model.region(ops,maximum)
        arrays=np.array([[[float(t[i,j]) for j in range(9)] for i in range(9)] for t in regions]) if floating else None
        operators[tilt]=regions,arrays
        print('Built outward two-bit regions through',maximum,'tilt',tilt,'precision',ctx.prec,flush=True)
    return operators


def replay_saved_covers(records,caps):
    """Reconstruct every saved bound; cached numerical upper bounds are ignored."""
    groups={};seen=set();results={}
    for record in records:
        updates,_,q,cutoff,degree,penalty=validate_record(record)
        if q in seen:raise ValueError('duplicate sparse occupancy')
        seen.add(q)
        groups.setdefault((updates,cutoff,degree,penalty),[]).append(record)
    for (updates,cutoff,degree,penalty),batch in groups.items():
        maximum=max(record['q'] for record in batch)
        tilts=sorted({leaf['tilt'] for record in batch for leaf in record['leaves']},key=Q)
        operators=build_operators(tilts,maximum,cutoff,str(penalty),degree,updates,floating=False)
        cdf=weighted_cdf_upper(caps,1<<128,rows=2,full_weight=1/penalty)
        shells=weighted_union_shells(caps,rows=2,full_weight=1/penalty);shells[0]-=1
        for record in sorted(batch,key=lambda r:r['q']):
            q=record['q']
            bound=replay_cover(operators,cdf,shells,q,record['leaves'],
                               density_fold=record.get('density_fold',False),output_cutoff=saved_threshold(record))
            if not bound>0:raise ArithmeticError('positive replayed sparse bound required')
            results[q]=bound
            print('REPLAYED SPARSE BRIDGE q',q,'margin',-bound.log()/arb(2).log(),flush=True)
    return results


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occupancies',type=int,nargs='+',default=[2,3,4,8,16])
    parser.add_argument('--tilts',nargs='+',default=['.00064','.001','.0016','.0025','.004','.0064','.01'])
    parser.add_argument('--cutoff',type=int,default=8)
    parser.add_argument('--penalty',default='1')
    parser.add_argument('--precision',type=int,default=192)
    parser.add_argument('--updates',type=int,help='Independent transvections per step; default two or the saved value')
    parser.add_argument('--threshold',type=int,help='Inclusive bad-output-weight cutoff; default 209715 or the saved cutoff')
    parser.add_argument('--target-bits',type=int,default=55)
    parser.add_argument('--max-splits',type=int,default=100)
    parser.add_argument('--joint',action='store_true')
    parser.add_argument('--screen-only',action='store_true')
    parser.add_argument('--output-degree',type=int,default=None)
    parser.add_argument('--density-fold',action='store_true',help='Dominate the complete outer support measure pointwise')
    parser.add_argument('--certificate-dir',type=Path,help='Write small checked witness files here; do not commit generated data')
    saved_mode=parser.add_mutually_exclusive_group()
    saved_mode.add_argument('--replay',type=Path,nargs='+',help='Reconstruct and replay saved support-cover witnesses, never trusting their stored upper bounds')
    saved_mode.add_argument('--retarget',type=Path,nargs='+',help='Recompute saved partitions after explicitly changing updates or the output cutoff')
    args=parser.parse_args()
    saved={};old_updates=set();old_thresholds=set()
    if args.replay or args.retarget:
        if args.screen_only:parser.error('saved witnesses require outward replay')
        if args.retarget and args.updates is None and args.threshold is None:
            parser.error('retarget requires explicit updates or output cutoff')
        for path in args.replay or args.retarget:
            record=json.loads(path.read_text())
            try:
                old_updates.add(saved_updates(record));old_thresholds.add(saved_threshold(record))
            except ValueError as error:parser.error(str(error))
            if Q(record['penalty'])!=Q(args.penalty):parser.error('witness penalty must match --penalty')
            q=record['q']
            if q in saved:parser.error('duplicate witness occupancy')
            saved[q]=record
        args.occupancies=sorted(saved)
        args.tilts=sorted({leaf['tilt'] for record in saved.values() for leaf in record['leaves']},key=Q)
        settings={(record.get('cutoff',args.cutoff),record.get('output_degree',args.output_degree)) for record in saved.values()}
        if len(settings)!=1:parser.error('replay witnesses with different local-operator settings separately')
        args.cutoff,args.output_degree=settings.pop()
    if len(old_updates)>1:parser.error('replay witnesses with different update counts separately')
    if len(old_thresholds)>1:parser.error('replay witnesses with different output cutoffs separately')
    try:
        args.updates=model.resolve_updates(args.updates,{'updates':next(iter(old_updates))} if saved else None,retarget=bool(args.retarget))
        args.threshold=resolve_threshold(args.threshold,{'threshold':next(iter(old_thresholds))} if saved else None,
                                         retarget=bool(args.retarget and args.threshold is not None))
    except ValueError as error:parser.error(str(error))
    if (any(not 2<=q<=4096 for q in args.occupancies) or args.precision<128 or not 1<=args.cutoff<=10
            or not 0<Q(args.penalty)<=1 or any(Q(t)<=0 for t in args.tilts)):
        parser.error('valid occupancies, census cutoff, precision, tilt and penalty required')
    print('ISOLATED TWO-BIT COMPLETE SUPPORT COVERS; full-code closure not presumed.',flush=True)
    print('INDEPENDENT UPDATES PER STEP',args.updates,flush=True)
    print('INCLUSIVE OUTPUT WEIGHT CUTOFF',args.threshold,flush=True)
    caps=authenticated_caps();ctx.prec=args.precision
    cdf=weighted_cdf_upper(caps,1<<128,rows=2,full_weight=1/Q(args.penalty))
    shells=weighted_union_shells(caps,rows=2,full_weight=1/Q(args.penalty));shells[0]-=1
    maximum=max(args.occupancies)
    operators=build_operators(args.tilts,maximum,args.cutoff,args.penalty,args.output_degree,args.updates)
    results={}
    for q in dict.fromkeys(args.occupancies):
        record=dict(schema='two-bit-support-cover-1',ensemble=[*ENSEMBLE_PREFIX,args.updates],
                    penalty=args.penalty,cutoff=args.cutoff,output_degree=args.output_degree,threshold=args.threshold)
        if saved:
            bound=replay_cover(operators,cdf,shells,q,saved[q]['leaves'],density_fold=saved[q].get('density_fold',False),output_cutoff=args.threshold)
            if not 0<bound<arb(2)**-args.target_bits:raise ArithmeticError('saved witness misses the replay budget')
            results[q]=bound
            record.update(q=q,leaves=saved[q]['leaves'],density_fold=saved[q].get('density_fold',False),
                          upper_dyadic=[int(v) for v in bound.upper().man_exp()])
            print('REPLAYED COMPLETE TWO-BIT COVER q',q,'margin',-bound.log()/arb(2).log(),flush=True)
        else:
            results[q]=complete_cover(operators,cdf,shells,q,target_bits=args.target_bits,
                                     max_splits=args.max_splits,joint=args.joint,screen_only=args.screen_only,certificate=record,
                                     density_fold=args.density_fold,output_cutoff=args.threshold)
        if results[q] is not None and args.certificate_dir:
            args.certificate_dir.mkdir(parents=True,exist_ok=True)
            destination=args.certificate_dir/f'q{q:04d}.json'
            if any(destination.resolve()==p.resolve() for p in (args.replay or args.retarget or [])):
                parser.error('write retargeted witnesses to a separate directory; preserve the original')
            destination.write_text(json.dumps(record,sort_keys=True,separators=(',',':'))+'\n')
    successes={q:bound for q,bound in results.items() if bound is not None}
    print('Certified occupancies this run:',sorted(successes),'unresolved:',[q for q,v in results.items() if v is None],flush=True)
    if successes:
        total=model.up(sum(successes.values(),arb(0)))
        print('OUTWARD LISTED-OCCUPANCY AGGREGATE',total,'margin',-total.log()/arb(2).log(),flush=True)
    print('No full-code certificate: the complete occupancy range and aggregate remain required.',flush=True)


if __name__=='__main__':main()
