"""Replay dense GF16 coverage and regenerate matching sparse support covers.

Run receipts are never accepted as sparse proof inputs. The selected
sparse argument is regenerated for every missing occupancy.
"""
import argparse
from contextlib import contextmanager,redirect_stdout
import hashlib
import json
import sys
from pathlib import Path
from fractions import Fraction as Q
from flint import arb,ctx
import scalar_cover as dense
import refresh_kernel
import feedback_refresh
import gf_refresh
import weighted_return
import profile_return
import shape_return
import conditioned_return
import trimmed_return
import rank_return
import birth_refresh
import birth_classes
import sparse_cover


def sparse_chunks(lo,hi,workers):
    """Disjoint contiguous ranges, approximately balanced by squared occupancy."""
    if (type(lo) is not int or type(hi) is not int or not 1<=lo<hi<=2049
            or type(workers) is not int or not 1<=workers<=8):
        raise ValueError('valid sparse range and one through eight workers required')
    parts=min(workers,hi-lo);remaining=sum(q*q for q in range(lo,hi));start=lo;result=[]
    for count in range(parts,1,-1):
        end=start;weight=0
        while end<hi-count+1 and (end==start or count*weight<remaining):
            weight+=end*end;end+=1
        result.append((start,end));start=end;remaining-=weight
    result.append((start,hi))
    return result


def sparse_range(lo,hi,tilts,args):
    return sparse_cover.run(list(range(lo,hi)),args.threshold,tilts,
        args.precision,args.max_splits,getattr(args,'sparse_target_bits',52),inner=args.sparse_inner,
        shape_penalties=args.shape_penalties,shape_weight_tilt=args.shape_weight_tilt,
        exact_feedback=args.exact_feedback,updates=getattr(args,'updates',2),
        joint_return_through=getattr(args,'sparse_joint_return_through',None),
        lazy_density_through=getattr(args,'sparse_lazy_density_through',None),
        analytic_gradient=getattr(args,'sparse_analytic_gradient',False))


def sparse_worker(task):
    """Fresh child computation; no saved operator or bound is an input."""
    lo,hi,tilts,args,log_path=task
    if log_path is None:upper=sparse_range(lo,hi,tilts,args)
    else:
        log_path=Path(log_path);log_path.parent.mkdir(parents=True,exist_ok=True)
        with log_path.open('w') as log,redirect_stdout(log):upper=sparse_range(lo,hi,tilts,args)
    if upper is None:return lo,hi,None
    if not upper>0:raise ArithmeticError('positive fresh sparse bound required')
    return lo,hi,tuple(int(x) for x in upper.upper().man_exp())


@contextmanager
def sparse_spawn_path():
    # Legacy proof modules prepend their own directories. A spawned child
    # must resolve this assemble.py, not another workstream's namesake.
    original=list(sys.path)
    sys.path.insert(0,str(Path(__file__).resolve().parent))
    try:yield
    finally:sys.path[:]=original


def parallel_sparse(ranges,args,workers):
    """Combine fresh process results only after checking their full ranges."""
    from concurrent.futures import ProcessPoolExecutor
    from multiprocessing import get_context
    tasks=[];log_base=getattr(args,'sparse_log',None)
    for lo,hi,tilts in ranges:
        if lo>=hi:continue
        for left,right in sparse_chunks(lo,hi,workers):
            path=(None if log_base is None else Path(log_base).with_name(
                f'{Path(log_base).stem}-q{left}-{right-1}{Path(log_base).suffix}'))
            tasks.append((left,right,tilts,args,path))
    if not tasks:return arb(0)
    total=arb(0);complete=True;seen=0
    # Spawn isolates Arb precision and the numeric libraries on every OS.
    # Each child constructs its own checked operators; only fresh dyadic
    # results return through this invocation's process channel.
    with sparse_spawn_path(),ProcessPoolExecutor(max_workers=workers,mp_context=get_context('spawn')) as pool:
        for task,result in zip(tasks,pool.map(sparse_worker,tasks)):
            lo,hi=task[:2]
            if result[:2]!=(lo,hi):raise ArithmeticError('sparse worker coverage mismatch')
            seen+=1;bound=result[2]
            if bound is None:
                complete=False;print('Sparse worker INCOMPLETE',lo,'through',hi-1,flush=True);continue
            if (len(bound)!=2 or any(type(x) is not int for x in bound) or bound[0]<=0):
                raise ArithmeticError('invalid fresh dyadic sparse bound')
            upper=dense.kernel.up(arb(bound[0])*arb(2)**bound[1]);total=dense.kernel.up(total+upper)
            print('Sparse worker VERIFIED',lo,'through',hi-1,'log2 upper',upper.log()/arb(2).log(),flush=True)
    if seen!=len(tasks):raise ArithmeticError('missing sparse worker result')
    return total if complete else None


def validate(record):
    if (record.get('schema')!=dense.SCHEMA or type(record.get('updates')) is not int or record['updates'] not in (2,3,4)
            or type(record.get('minimum_groups')) is not int or not 1<=record['minimum_groups']<=dense.G
            or type(record.get('threshold')) is not int or not 0<=record['threshold']<dense.N
            or record.get('unresolved')!={} or record.get('screen_only')
            or not isinstance(record.get('leaves'),dict) or not record['leaves']
            or record.get('kernel','two-state') not in ('two-state','refresh','feedback-refresh','gf-refresh','weighted-return','profile-return','shape-return','conditioned-return','trimmed-return','rank-return','birth-refresh','birth-classes')):
        raise ValueError('complete GF16 cover with two through four updates and valid metadata required')
    if Q(record.get('base_tilt','1/32'))<=0:raise ValueError('positive base activity tilt required')
    bits=record.get('central_bits',130)
    if (type(bits) is not int or not 1<=bits<=256 or Q(record.get('central_scale','1'))<=0
            or type(record.get('row_parity',False)) is not bool
            or type(record.get('variance_shuffle',False)) is not bool
            or type(record.get('variance_bins',0)) is not int or not 0<=record.get('variance_bins',0)<=64
            or (record.get('variance_bins',0) and not record.get('variance_shuffle',False))
            or type(record.get('regional_count',False)) is not bool
            or (record.get('regional_count',False) and (not record.get('variance_bins',0) or record.get('kernel')!='birth-classes'))
            or not 0<Q(record.get('row_bias','2/5'))<Q(1,2)):
        raise ValueError('valid row comparison parameters required')


def dense_model(record,precision,*,geometry_only=False):
    validate(record)
    if type(precision) is not int or precision<128:raise ValueError('precision >=128 required')
    inner={'two-state':dense.kernel,'refresh':refresh_kernel,'feedback-refresh':feedback_refresh,
           'gf-refresh':gf_refresh,'weighted-return':weighted_return,'profile-return':profile_return,
           'shape-return':shape_return,'conditioned-return':conditioned_return,
           'trimmed-return':trimmed_return,'rank-return':rank_return,'birth-refresh':birth_refresh,
           'birth-classes':birth_classes}[record.get('kernel','two-state')]
    components=dense.actual_components(record.get('central_bits',130),Q(record.get('row_bias','2/5')),
                                      Q(record.get('central_scale','1')),record.get('row_parity',False))
    # The parent uses only integer/rational geometry; every worker constructs
    # the actual inner maps before evaluating any witness.
    data={'windows':32} if geometry_only else inner.actual(record['updates'])
    ctx.prec=precision
    return dense.Model(components,data,record['threshold'],record['minimum_groups'],
        tilt=Q(record.get('base_tilt','1/32')),inner=inner,
        variance_shuffle=record.get('variance_shuffle',False),
        variance_bins=record.get('variance_bins',0),regional_count=record.get('regional_count',False))


def initialize_dense_worker(record,precision,log_base):
    global _dense_model,_dense_cells,_dense_record,_dense_precision,_dense_log
    import os
    _dense_log=None if log_base is None else open(f'{log_base}-worker-{os.getpid()}.log','w')
    with redirect_stdout(_dense_log if _dense_log is not None else sys.stdout):
        _dense_model=dense_model(record,precision)
        _dense_cells=dense.partition(_dense_model,record['leaves'],record['unresolved'])
    _dense_record=record;_dense_precision=precision


def dense_worker(path):
    if ctx.prec!=_dense_precision:raise ArithmeticError('dense worker precision changed')
    with redirect_stdout(_dense_log if _dense_log is not None else sys.stdout):
        upper=_dense_model.outward(_dense_cells[path],_dense_record['leaves'][path]['witness'])
        if ctx.prec!=_dense_precision or not upper>0:
            raise ArithmeticError('positive fresh dense bound at requested precision required')
        print('DENSE CELL VERIFIED',path,flush=True)
    return path,tuple(int(x) for x in upper.upper().man_exp())


def sum_dense_results(paths,results):
    """Require exactly one fresh positive result per leaf, in fixed order."""
    if len(set(paths))!=len(paths):raise ArithmeticError('duplicate assigned dense leaf')
    total=arb(0);seen=0
    for result in results:
        if seen>=len(paths) or result[0]!=paths[seen]:
            raise ArithmeticError('dense worker coverage mismatch')
        bound=result[1]
        if (not isinstance(bound,(tuple,list)) or len(bound)!=2
                or any(type(x) is not int for x in bound) or bound[0]<=0):
            raise ArithmeticError('invalid fresh dyadic dense bound')
        upper=dense.kernel.up(arb(bound[0])*arb(2)**bound[1])
        total=dense.kernel.up(total+upper);seen+=1
        if seen%20==0:print('DENSE REPLAY checked',seen,'of',len(paths),flush=True)
    if seen!=len(paths):raise ArithmeticError('missing dense worker result')
    return total


def replay_dense(record,precision,workers=1,log_base=None):
    if type(workers) is not int or not 1<=workers<=4:
        raise ValueError('one through four dense workers required')
    model=dense_model(record,precision,geometry_only=workers>1)
    cells=dense.partition(model,record['leaves'],record['unresolved'])
    paths=list(record['leaves'])
    if workers==1:
        def results():
            for path in paths:
                upper=model.outward(cells[path],record['leaves'][path]['witness'])
                if ctx.prec!=precision or not upper>0:
                    raise ArithmeticError('positive fresh dense bound at requested precision required')
                yield path,tuple(int(x) for x in upper.upper().man_exp())
        return sum_dense_results(paths,results())
    from concurrent.futures import ProcessPoolExecutor
    from multiprocessing import get_context
    if log_base is not None:
        log_base=Path(log_base);log_base.parent.mkdir(parents=True,exist_ok=True)
    # Fixed result order preserves the serial order of outward additions.
    with sparse_spawn_path(),ProcessPoolExecutor(max_workers=workers,mp_context=get_context('spawn'),
            initializer=initialize_dense_worker,initargs=(record,precision,log_base)) as pool:
        return sum_dense_results(paths,pool.map(dense_worker,paths))


def regenerate_sparse(minimum_groups,args):
    """Regenerate disjoint sparse ranges; never import sparse receipts."""
    if minimum_groups==1:return arb(0)
    small_through=getattr(args,'sparse_small_through',0)
    small_tilts=getattr(args,'sparse_small_tilts',None)
    target_bits=getattr(args,'sparse_target_bits',52)
    updates=getattr(args,'updates',2)
    workers=getattr(args,'sparse_workers',1)
    if type(workers) is not int or not 1<=workers<=8:raise ValueError('one through eight sparse workers required')
    joint=getattr(args,'sparse_joint_return_through',None)
    density=getattr(args,'sparse_lazy_density_through',None)
    analytic=getattr(args,'sparse_analytic_gradient',False)
    if type(analytic) is not bool:raise ValueError('boolean sparse analytic-gradient option required')
    if (type(small_through) is not int or not 0<=small_through<=2048
            or bool(small_through)!=bool(small_tilts) or type(target_bits) is not int or target_bits<40
            or type(updates) is not int or updates not in (2,3,4)
            or (small_tilts and any(Q(t)<=0 for t in small_tilts))):
        raise ValueError('valid sparse budget and paired small-occupancy range/grid required')
    for value,limit in ((joint,4),(density,32)):
        if value is not None and (args.sparse_inner!='gf-birth-classes' or type(value) is not int or not 0<=value<=limit):
            raise ValueError('sparse local refinements require GF birth classes and valid integer cutoffs')
    first=1;total=arb(0)
    if args.single_group_exact:
        import single_group
        tilts=sorted(set(args.sparse_tilts+(small_tilts or [])),key=Q)
        total=single_group.run(args.threshold,tilts,args.precision,updates=updates)
        first=2
    boundary=min(minimum_groups,max(first,small_through+1))
    ranges=[(first,boundary,small_tilts),(boundary,minimum_groups,args.sparse_tilts)]
    if workers>1:
        remaining=parallel_sparse(ranges,args,workers)
        return None if remaining is None else dense.kernel.up(total+remaining)
    for lo,hi,tilts in ranges:
        if lo>=hi:continue
        remaining=sparse_range(lo,hi,tilts,args)
        if remaining is None:return None
        total=dense.kernel.up(total+remaining)
    return total


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('dense',type=Path)
    parser.add_argument('--precision',type=int,default=256)
    parser.add_argument('--dense-workers',type=int,default=1,
        help='Fresh dense witness replay in one through four processes; defaults to serial')
    parser.add_argument('--dense-log',type=Path,help='Optional prefix for per-process dense replay logs')
    parser.add_argument('--sparse-tilts',nargs='+',default=['.00032','.0032','.032'])
    parser.add_argument('--sparse-target-bits',type=int,default=52,
                        help='Per-occupancy stopping budget; final aggregate must still exceed 40 bits')
    parser.add_argument('--sparse-small-through',type=int,default=0)
    parser.add_argument('--sparse-small-tilts',nargs='+',
                        help='Separate tilt grid through the requested small occupancy; regenerated independently')
    parser.add_argument('--sparse-inner',choices=('shape-uniform','gf-shape-tilt','gf-occupancy','gf-density','gf-density-refined','gf-classes','gf-rank','gf-birth-classes'),default='shape-uniform')
    parser.add_argument('--shape-penalties',nargs='+',default=['1','.75','.5'])
    parser.add_argument('--shape-weight-tilt',default='1')
    parser.add_argument('--exact-feedback',action='store_true')
    parser.add_argument('--sparse-joint-return-through',type=int,choices=(0,1,2,3,4))
    parser.add_argument('--sparse-lazy-density-through',type=int)
    parser.add_argument('--sparse-analytic-gradient',action='store_true')
    parser.add_argument('--sparse-workers',type=int,default=1,
        help='Fresh independent occupancy ranges in one through eight processes; defaults to serial')
    parser.add_argument('--single-group-exact',action='store_true',
                        help='Use exact GF-rank placement at q=1, then regenerate q>=2 separately')
    parser.add_argument('--max-splits',type=int,default=100)
    parser.add_argument('--sparse-log',type=Path)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    if (args.precision<128 or args.max_splits<0 or not 1<=args.sparse_workers<=8
            or not 1<=args.dense_workers<=4):
        parser.error('precision >=128, nonnegative work limit, 1--8 sparse workers, and 1--4 dense workers required')
    raw=args.dense.read_bytes();record=json.loads(raw);validate(record)
    args.threshold=record['threshold']
    args.updates=record['updates']
    if args.updates>2 and (args.sparse_inner not in ('gf-rank','gf-birth-classes') or not args.exact_feedback):
        parser.error('three/four-update assembly requires an exact-feedback GF rank or birth-class sparse kernel')
    dense_upper=replay_dense(record,args.precision,args.dense_workers,args.dense_log)
    print('GF16 dense replay',record['minimum_groups'],'through',dense.G,
          'log2 upper',dense_upper.log()/arb(2).log(),flush=True)
    def regenerate():
        return regenerate_sparse(record['minimum_groups'],args)
    if args.sparse_log:
        args.sparse_log.parent.mkdir(parents=True,exist_ok=True)
        print('Regenerating sparse covers without operator cache; log',args.sparse_log,flush=True)
        with args.sparse_log.open('w') as log,redirect_stdout(log):sparse_upper=regenerate()
    else:sparse_upper=regenerate()
    if sparse_upper is None:raise ArithmeticError('sparse cover incomplete; no whole-code certificate')
    total=dense.kernel.up(sparse_upper+dense_upper)
    if not 0<total<arb(2)**-40:raise ArithmeticError('complete sum misses 40 bits')
    print('COMPLETE GF16 FOUR-BIT CERTIFICATE; ideal setup; updates',args.updates,flush=True)
    print('minimum distance >=',record['threshold']+1,'of',dense.N,flush=True)
    print('bad-setup probability upper',total,flush=True)
    print('margin enclosure',-total.log()/arb(2).log(),flush=True)
    if args.output:
        report=dict(schema='gf16-complete-replay-1',dense_sha256=hashlib.sha256(raw).hexdigest(),
            precision=args.precision,updates=args.updates,threshold=record['threshold'],minimum_distance=record['threshold']+1,
            sparse_inner=args.sparse_inner,
            sparse_tilts=args.sparse_tilts,shape_penalties=args.shape_penalties,
            sparse_target_bits=args.sparse_target_bits,sparse_small_through=args.sparse_small_through,
            sparse_small_tilts=args.sparse_small_tilts,
            shape_weight_tilt=args.shape_weight_tilt,max_splits=args.max_splits,
            exact_feedback=args.exact_feedback,
            sparse_joint_return_through=args.sparse_joint_return_through,
            sparse_lazy_density_through=args.sparse_lazy_density_through,
            sparse_analytic_gradient=args.sparse_analytic_gradient,
            sparse_workers=args.sparse_workers,
            dense_workers=args.dense_workers,
            single_group_exact=args.single_group_exact,
            output_length=dense.N,sparse_through=record['minimum_groups']-1,dense_through=dense.G,
            dense_upper=[int(x) for x in dense_upper.upper().man_exp()],
            sparse_upper=[int(x) for x in sparse_upper.upper().man_exp()],
            total_upper=[int(x) for x in total.upper().man_exp()],
            note='Regenerated sparse operators and support covers, replayed complete dense partition; rerun this program to verify.')
        args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
