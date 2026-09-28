"""Find the largest certified cutoff supported by a fixed complete partition.

For fixed witnesses, each cell bound is A*exp(lambda*d). Recompute A at
d=0, then vary only the inclusive output-weight cutoff d. The sparse
lemmas still limit this calculation to d<=209715 and two updates.
This is a frontier of the saved witnesses, not of the actual code.
"""
import argparse
import hashlib
import json
from fractions import Fraction as Q
from pathlib import Path

from flint import arb,ctx

import model,iid_kernel,mixture_probe
from closure import MODELS,dense_scope,sparse_upper
from mixture_cover import partition,PACKETS
from row_mixture import envelope,pair_components
from bch_joint_support import authenticated_caps
from probe import aq


def output_tilt(witness):
    while isinstance(witness,dict):
        if 'old' in witness:witness=witness['old']
        elif 'base' in witness:witness=witness['base']
        elif 'plane' in witness:
            witness=[witness['plane']['tilt']];break
        else:witness=witness['parameters']
    value=Q(witness[0])
    if value<=0:raise ValueError('positive output tilt required')
    return value


def upper_at(terms,cutoff):
    if type(cutoff) is not int or cutoff<0:raise ValueError('nonnegative integer cutoff required')
    total=arb(0)
    for lam,coefficient in terms:
        if Q(lam)<=0 or not coefficient>0:raise ValueError('positive coefficients and tilts required')
        total=model.up(total+coefficient*aq(Q(lam)*cutoff).exp())
    return total


def frontier(terms,sparse,maximum,target_bits):
    if type(maximum) is not int or maximum<0 or type(target_bits) is not int or target_bits<1:
        raise ValueError('nonnegative maximum and positive integer margin required')
    if not sparse>=0:raise ValueError('nonnegative sparse bound required')
    target=arb(2)**-target_bits
    if not model.up(sparse+upper_at(terms,0))<target:
        raise ArithmeticError('saved witnesses do not close even at zero cutoff')
    low,high=0,maximum+1
    while high-low>1:
        mid=(low+high)//2
        if model.up(sparse+upper_at(terms,mid))<target:low=mid
        else:high=mid
    dense=upper_at(terms,low)
    return low,dense,model.up(sparse+dense)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dense',type=Path,required=True)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--target-bits',type=int,default=40)
    parser.add_argument('--precision',type=int,default=256)
    args=parser.parse_args()
    if args.precision<192 or args.target_bits<1:parser.error('positive margin and at least 192-bit precision required')
    if args.output and args.output.resolve()==args.dense.resolve():parser.error('preserve the original witness file')
    raw=args.dense.read_bytes();record=json.loads(raw);old_threshold=dense_scope(record)
    caps=authenticated_caps();ctx.prec=args.precision
    components=pair_components(envelope(caps,1<<record['central_bits'],Q(record['theta'])))
    images,columns,_=model.maps();data=iid_kernel.prepare(images,columns,19)
    instance=MODELS[record['schema']](components,data,record['minimum_groups'],
                                    [1,*map(Q,record['input_tilt'])],threshold=0)
    cells=partition(instance,record['leaves'],{});terms=[]
    for index,(path,leaf) in enumerate(record['leaves'].items(),1):
        cell=cells[path]
        if not instance.empty(cell):
            if leaf.get('empty'):raise ValueError('incorrectly pruned cell')
            witness=leaf['witness']
            if isinstance(witness,dict) and 'compositions' in witness:
                for row in instance.composition_rows(instance.integer_cell(cell),witness):
                    coefficient=mixture_probe.outward(components,data,row['counts'],row['parameters'],threshold=0,
                                                      density_anchors=row.get('density_anchors'),capped=row.get('capped',False))
                    terms.append((output_tilt(row['parameters']),coefficient))
            else:
                coefficient=instance.outward(cell,witness)
                terms.append((output_tilt(witness),coefficient))
        if index%200==0:print('RECOMPUTED THRESHOLD-INDEPENDENT COEFFICIENTS',index,'/',len(cells),flush=True)
    cutoff,dense,total=frontier(terms,sparse_upper(),model.THRESHOLD,args.target_bits)
    print('FIXED-WITNESS FRONTIER cutoff',cutoff,'of',2*PACKETS,'fraction',Q(cutoff,2*PACKETS),flush=True)
    print('DENSE margin',-dense.log()/arb(2).log(),'FULL margin',-total.log()/arb(2).log(),flush=True)
    if args.output:
        record.update(threshold=cutoff,updates=2,upper_dyadic=[int(v) for v in dense.upper().man_exp()],
                      threshold_search=dict(source_threshold=old_threshold,target_bits=args.target_bits,
                                            source_sha256=hashlib.sha256(raw).hexdigest()))
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(record,sort_keys=True,separators=(',',':'))+'\n')
    print('Sparse inputs are prior verified lemmas. This is not a maximum-distance or optimal-proof claim.',flush=True)


if __name__=='__main__':main()
