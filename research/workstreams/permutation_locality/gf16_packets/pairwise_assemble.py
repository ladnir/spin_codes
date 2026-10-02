"""Fresh full verification for the fixed pairwise4-gf16-r4 target.

Reject incomplete dense partitions before expensive work. Rebuild the
comparison and inner maps, replay every dense witness, and regenerate
occupancies 1..48. Saved sparse bounds are never accepted as proof inputs.
"""
import argparse
import hashlib
import json
from fractions import Fraction as Q
from pathlib import Path
from types import SimpleNamespace

from flint import arb,ctx
import pairwise_cover as cover
import pairwise_sparse as sparse
from audit_complete import dyadic

SCHEMA='pairwise-gf16-complete-replay-1'
TILTS=('.00016','.00024','.00028','.00032','.0004','.00064','.001','.0016',
       '.0032','.005','.008','.012','.016','.02','.024','.028','.032','.036',
       '.04','.048','.064','.096')


def validate_dense(record):
    cover.validate_record(record)
    if (record['unresolved'] or not record['leaves'] or record.get('screen_only',False)
            or any(not isinstance(row,dict) or not isinstance(row.get('witness'),dict)
                   for row in record['leaves'].values())):
        raise ValueError('complete dense partition with actual witnesses required')
    # Completeness of binary subdivision paths does not depend on the root
    # endpoints. The actual model reconstructs those endpoints again below.
    cover.sc.partition(SimpleNamespace(root=(Q(0),Q(1))),record['leaves'],{})


def sum_sparse(rows,through):
    if (not isinstance(rows,list) or len(rows)!=through
            or any(not isinstance(row,dict) or type(row.get('occupancy')) is not int for row in rows)
            or sorted(row['occupancy'] for row in rows)!=list(range(1,through+1))):
        raise ValueError('exactly one fresh result for every sparse occupancy required')
    total=Q(0)
    for row in rows:
        if row.get('upper') is None:raise ValueError('fresh sparse computation is incomplete')
        upper=dyadic(row['upper'])
        if upper<=0:raise ValueError('strictly positive sparse upper endpoint required')
        total+=upper
    return total


def validate_geometry(model):
    expected=dict(G=2048,REGIONS=256,PACKETS=1<<19,EPOCHS=1<<14,N=1<<21)
    if (any(getattr(cover.sc,key)!=value for key,value in expected.items())
            or model.threshold!=209715 or model.q_min!=49
            or any(model.data.get(key)!=value for key,value in dict(bits=19,windows=32,updates=4).items())):
        raise ValueError('the actual model must have the fixed K=2^20, four-update IMT(128,19) geometry')


def exact_endpoint(value):
    value=Q(value);den=value.denominator
    if value<0 or den&(den-1):raise ValueError('nonnegative exact dyadic value required')
    return [value.numerator,-(den.bit_length()-1)]


def combine(dense_result,sparse_rows,through,bits):
    if type(bits) is not int or bits<40:raise ValueError('at least 40 margin bits required')
    if (dense_result.get('dense_complete') is not True or dense_result.get('unresolved')!=0):
        raise ValueError('fresh dense replay did not close')
    dense_upper=dyadic(dense_result['upper'])
    if dense_upper<=0:raise ValueError('strictly positive dense upper endpoint required')
    sparse_upper=sum_sparse(sparse_rows,through)
    total=dense_upper+sparse_upper
    if total>=Q(2)**-bits:raise ValueError('fresh aggregate does not meet the requested margin')
    return dict(dense_upper=exact_endpoint(dense_upper),sparse_upper=exact_endpoint(sparse_upper),
                total_upper=exact_endpoint(total))


def verify(record,precision=384,bits=40,max_splits=128,sparse_output=None):
    validate_dense(record)
    if (type(precision) is not int or precision<256 or type(bits) is not int or bits<40
            or type(max_splits) is not int or max_splits<0):
        raise ValueError('precision >=256, margin >=40, and nonnegative search budget required')
    p=record['parameters']
    model=cover.build_model(precision,p)
    validate_geometry(model)
    dense_result=cover.replay(model,record,bits)
    if ctx.prec!=precision:raise ArithmeticError('precision changed during dense replay')
    if not dense_result['dense_complete']:raise ValueError('fresh dense replay is insufficient')
    del model  # Release dense search scratch before building sparse operators.
    through=p['minimum_groups']-1
    rows=sparse.run(list(range(1,through+1)),list(TILTS),precision=precision,
        max_splits=max_splits,target_bits=48,threshold=p['threshold'],updates=p['updates'],
        refined_counts=True,coupled_counts=True,output=sparse_output)
    if ctx.prec!=precision:raise ArithmeticError('precision changed during sparse regeneration')
    result=combine(dense_result,rows,through,bits)
    result.update(schema=SCHEMA,ensemble=cover.ENSEMBLE,complete=True,
        precision=precision,requested_bits=bits,message_length=1<<20,output_length=1<<21,
        updates=p['updates'],threshold=p['threshold'],minimum_distance=p['threshold']+1,
        sparse_occupancies=[1,through],dense_occupancies=[p['minimum_groups'],2048],
        dense_leaves=len(record['leaves']),dense_parameters=p,
        sparse_parameters=dict(tilts=list(TILTS),max_splits=max_splits,target_bits=48,
            refined_counts=True,coupled_counts=True,joint_return_through=3,lazy_density_through=6),
        note='Fresh dense replay and fresh sparse regeneration for the specified pairwise ensemble. '
             'No saved sparse upper bound or dense acceptance label was trusted.')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('dense',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--precision',type=int,default=384)
    parser.add_argument('--bits',type=int,default=40)
    parser.add_argument('--max-splits',type=int,default=128)
    args=parser.parse_args()
    prefix=args.output.with_suffix('.sparse.json')
    if args.output.exists() or prefix.exists() or args.dense.resolve() in (args.output.resolve(),prefix.resolve()):
        parser.error('use new output paths; existing proofs are preserved')
    raw=args.dense.read_bytes();record=json.loads(raw)
    result=verify(record,args.precision,args.bits,args.max_splits,prefix)
    result['dense_sha256']=hashlib.sha256(raw).hexdigest()
    result['sparse_receipt_sha256']=hashlib.sha256(prefix.read_bytes()).hexdigest()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    m,e=result['total_upper'];margin=-(arb(m).log()/arb(2).log()+e)
    print('PAIRWISE COMPLETE: distance >10%, fresh aggregate margin',margin,flush=True)


if __name__=='__main__':main()
