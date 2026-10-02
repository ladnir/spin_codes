"""Fresh process-parallel checking of one saved full-domain partition.

Every cell receives a fresh proposal and outward check. A saved rational
witness is only a fallback, rechecked for the actual requested claim.
Failed cells remain unresolved for the adaptive scalar-cover search.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
from contextlib import redirect_stdout
from copy import deepcopy
from fractions import Fraction as Q
import hashlib
import json
from multiprocessing import get_context
import os
from pathlib import Path
import sys
from flint import arb,ctx
# Spawn inherits paths modified by the older proof modules. Resolve this
# workstream's scalar_cover, not another workstream's namesake.
sys.path.insert(0,str(Path(__file__).resolve().parent))
import scalar_cover as sc
import birth_classes


def make_model(record,precision,*,geometry_only=False):
    if (record.get('schema')!=sc.SCHEMA or record.get('kernel')!='birth-classes'
            or type(record.get('updates')) is not int or record['updates'] not in (2,3,4)
            or type(precision) is not int or precision<128
            or not 0<Q(record.get('row_bias','2/5'))<Q(1,2)):
        raise ValueError('valid GF16 birth-class claim and precision required')
    components=sc.actual_components(record.get('central_bits',130),Q(record.get('row_bias','2/5')),
        Q(record.get('central_scale','1')),record.get('row_parity',False))
    data={'windows':32} if geometry_only else birth_classes.actual(record['updates'])
    ctx.prec=precision
    return sc.Model(components,data,record['threshold'],record['minimum_groups'],
        tilt=Q(record.get('base_tilt','1/32')),inner=birth_classes,
        variance_shuffle=record.get('variance_shuffle',False),
        variance_bins=record.get('variance_bins',0),regional_count=record.get('regional_count',False))


def check_cell(model,cell,old,target_bits):
    score,witness=model.proposal(cell)
    if score<-(target_bits+2):
        upper=model.outward(cell,witness)
        if not 0<upper<arb(2)**-target_bits:
            raise ArithmeticError('fresh outward proposal misses target')
    elif old is not None:
        witness=old['witness'];upper=model.outward(cell,witness)
        if not upper>0:raise ArithmeticError('positive fallback bound required')
        if not upper<arb(2)**-target_bits:return dict(proposal=float(score),upper=None)
    else:return dict(proposal=float(score),upper=None)
    return dict(witness=witness,proposal=float(upper.log()/arb(2).log()),
                upper=[int(x) for x in upper.upper().man_exp()])


def initialize_worker(record,precision,target_bits,log_base):
    global _model,_cells,_record,_target,_log,_precision
    _log=open(f'{log_base}-worker-{os.getpid()}.log','w')
    with redirect_stdout(_log):
        _model=make_model(record,precision)
        _cells=sc.partition(_model,record['leaves'],record['unresolved'])
    _record=record;_target=target_bits;_precision=precision


def worker(path):
    if ctx.prec!=_precision:raise ArithmeticError('worker precision changed before checking')
    with redirect_stdout(_log):
        print('CHECK CELL',path,flush=True)
        result=check_cell(_model,_cells[path],_record['leaves'].get(path),_target)
        if ctx.prec!=_precision:raise ArithmeticError('worker precision changed during checking')
        print('CELL RESULT',path,'verified',result['upper'] is not None,flush=True)
    return path,result


def accept_result(model,leaves,unresolved,expected,result,target_bits):
    """Account for exactly one fresh result; accepted labels alone never suffice."""
    path,row=result
    if path!=expected or path not in unresolved or path in leaves:
        raise ArithmeticError('duplicate or mismatched worker cell')
    bound=row['upper']
    if bound is None:
        unresolved[path]=dict(unresolved[path],proposal=row['proposal'])
        return False
    if (not isinstance(bound,(list,tuple)) or len(bound)!=2
            or any(type(v) is not int for v in bound) or bound[0]<=0
            or 'witness' not in row):
        raise ArithmeticError('positive fresh dyadic bound and witness required')
    upper=sc.kernel.up(arb(bound[0])*arb(2)**bound[1])
    if not 0<upper<arb(2)**-target_bits:raise ArithmeticError('worker bound misses target')
    leaves[path]=dict(witness=row['witness'],proposal=row['proposal'])
    del unresolved[path]
    return True


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--resume',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--retarget-updates',type=int,choices=(2,3,4))
    parser.add_argument('--retarget-distance')
    parser.add_argument('--precision',type=int,default=256)
    parser.add_argument('--target-bits',type=int,default=60)
    parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args()
    if (args.resume.resolve()==args.output.resolve() or not 1<=args.workers<=4
            or args.precision<128 or args.target_bits<40):
        parser.error('separate output, one through four workers, and valid precision/target required')
    raw=args.resume.read_bytes();record=json.loads(raw)
    old_updates=record['updates'];old_threshold=record['threshold']
    threshold,updates=sc.retarget_claim(record,args.retarget_updates,args.retarget_distance)
    record=dict(record,threshold=threshold,updates=updates)
    model=make_model(record,args.precision,geometry_only=True)
    cells=sc.partition(model,record['leaves'],record['unresolved'])
    leaves={};unresolved={p:dict(cell=list(map(str,c))) for p,c in cells.items()};seen=set()
    previous=record.get('visited',0)
    if type(previous) is not int or previous<0:raise ValueError('valid previous work count required')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    def save():
        sc.partition(model,leaves,unresolved)
        result=deepcopy(record)
        result.update(leaves=leaves,unresolved=unresolved,visited=previous+len(seen),precision=args.precision,
            retargeted_from=dict(sha256=hashlib.sha256(raw).hexdigest(),updates=old_updates,threshold=old_threshold),
            parallel_workers=args.workers)
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    save()
    print('PARALLEL FRESH COVER',len(cells),'cells; updates',updates,'cutoff',threshold,flush=True)
    # Each process authenticates the inputs and reconstructs its own operators.
    # Only fresh dyadic results return through this invocation's process channel.
    with ProcessPoolExecutor(max_workers=args.workers,mp_context=get_context('spawn'),
            initializer=initialize_worker,
            initargs=(record,args.precision,args.target_bits,str(args.output.with_suffix('')))) as pool:
        tasks={pool.submit(worker,path):path for path in sorted(cells,key=lambda p:(len(p),p))}
        for task in as_completed(tasks):
            path=tasks[task]
            if path in seen:raise ArithmeticError('worker cell accounted twice')
            accept_result(model,leaves,unresolved,path,task.result(),args.target_bits);seen.add(path)
            if len(seen)%20==0:
                save();print('PARALLEL COVER checked',len(seen),'accepted',len(leaves),
                    'unresolved',len(unresolved),'total',len(cells),flush=True)
    if seen!=set(cells):raise ArithmeticError('missing worker result')
    save()
    print('PARALLEL COVER FINISHED accepted',len(leaves),'unresolved',len(unresolved),flush=True)
    print('INCOMPLETE: continue adaptive search.' if unresolved else
          'All cells checked; fresh aggregate replay still required.',flush=True)


if __name__=='__main__':main()
