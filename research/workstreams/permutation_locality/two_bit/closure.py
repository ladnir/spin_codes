"""Replay a dense cover and combine it with the documented sparse lemmas.

This driver replays the dense witnesses and any additional sparse bridge.
The earlier q<=400 sparse calculations are conservative bounds from
COVER_STATUS.md; its commands reproduce them independently. No cached JSON
upper bound is trusted here.
"""
import argparse
import json
from fractions import Fraction as Q
from pathlib import Path

from flint import arb,ctx

import model,iid_kernel
from mixture_cover import CoverModel,SCHEMA as BASE_SCHEMA,replay,resolve_threshold,G,PACKETS
from lifted_mixture import LiftedModel,SCHEMA as LIFTED_SCHEMA
from affine_mixture import AffineModel,SCHEMA as AFFINE_SCHEMA
from row_mixture import envelope,pair_components
from bch_joint_support import authenticated_caps
from probe import aq
from full_cover import validate_record,replay_saved_covers

# Inclusive occupancy ranges, with lower bounds on their failure margins.
# All were verified at cutoff 209715, so apply to every smaller cutoff.
SPARSE_LEMMAS=(
    (1,1,'49.11287738250'),(2,32,'94'),(33,64,'564'),
    (65,128,'357'),(129,168,'70'),(169,169,'163'),
    (170,170,'210'),(171,320,'87'),(321,321,'3623'),(322,400,'348'))
MODELS={BASE_SCHEMA:CoverModel,LIFTED_SCHEMA:LiftedModel,AFFINE_SCHEMA:AffineModel}


def sparse_upper(lemmas=SPARSE_LEMMAS):
    previous=0;total=arb(0)
    for lo,hi,bits in lemmas:
        if lo!=previous+1 or hi<lo or hi>400 or Q(bits)<=0:
            raise ValueError('complete consecutive sparse range and positive margins required')
        total=model.up(total+(-aq(bits)*arb(2).log()).exp());previous=hi
    if previous!=400:raise ValueError('sparse lemmas must cover every occupancy through 400')
    return total


def dense_scope(record,bridge=()):
    if record.get('schema') not in MODELS:raise ValueError('unsupported dense cover')
    if model.resolve_updates(None,record)!=2:raise ValueError('documented sparse lemmas require two updates')
    if record.get('unresolved') or record.get('screen_only'):
        raise ValueError('a complete outward dense cover is required')
    threshold=resolve_threshold(None,record)
    if threshold>model.THRESHOLD:raise ValueError('sparse lemmas do not cover this larger distance')
    q=record['minimum_groups']
    if type(q) is not int or not 1<=q<=G:raise ValueError('bounded dense occupancy required')
    covered=[]
    for extra in bridge:
        updates,cutoff,occupancy,_,_,penalty=validate_record(extra)
        if updates!=2 or penalty!=1:raise ValueError('sparse bridge must use the same two-update, unweighted ensemble')
        if cutoff<threshold:raise ValueError('sparse bridge cutoff does not cover the dense bad event')
        covered.append(occupancy)
    if sorted(covered)!=list(range(401,q)):
        raise ValueError('sparse bridge must fill every occupancy gap exactly once')
    return threshold


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dense',type=Path,required=True)
    parser.add_argument('--bridge',type=Path,nargs='+',help='Complete sparse witnesses filling q=401 through the dense threshold minus one')
    parser.add_argument('--precision',type=int,default=256)
    parser.add_argument('--target-bits',type=int,default=40)
    args=parser.parse_args()
    if args.precision<192 or args.target_bits<1:parser.error('at least 192 bits of outward precision and a positive target required')
    record=json.loads(args.dense.read_text());bridge=[json.loads(p.read_text()) for p in args.bridge or []]
    threshold=dense_scope(record,bridge)
    caps=authenticated_caps();ctx.prec=args.precision
    components=pair_components(envelope(caps,1<<record['central_bits'],Q(record['theta'])))
    images,columns,_=model.maps();data=iid_kernel.prepare(images,columns,19)
    instance=MODELS[record['schema']](components,data,record['minimum_groups'],
                                   [1,*map(Q,record['input_tilt'])],threshold=threshold)
    dense=replay(instance,record)
    if not dense>0:raise ArithmeticError('positive dense enclosure required')
    print('DENSE COMPONENT REPLAYED at precision',ctx.prec,'margin',-dense.log()/arb(2).log(),flush=True)
    extra=replay_saved_covers(bridge,caps)
    sparse=sparse_upper();bridge_total=model.up(sum(extra.values(),arb(0)))
    total=model.up(sparse+dense+bridge_total)
    if not 0<total<arb(2)**-args.target_bits:raise ArithmeticError('full-code aggregate misses its requested margin')
    print('SPARSE COMPONENTS: documented prior replay bounds, cutoff',model.THRESHOLD,flush=True)
    if extra:print('SPARSE BRIDGE REPLAYED q',min(extra),'through',max(extra),'margin',-bridge_total.log()/arb(2).log(),flush=True)
    print('FULL OCCUPANCY RANGE 1..',G,'inclusive bad-weight cutoff',threshold,'of',2*PACKETS,flush=True)
    print('AGGREGATE UPPER',total,'margin',-total.log()/arb(2).log(),flush=True)


if __name__=='__main__':main()
