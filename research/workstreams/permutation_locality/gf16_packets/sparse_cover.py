"""Full-support sparse covers for the GF16 ensemble.

The direct bridge uses penalty=1 and input-weight tilt=1. The separate
gf-shape-tilt option pays an explicit GF change-of-measure normalizer.
Both use unweighted outer counts; old certificates are never imported.
"""
import argparse
from pathlib import Path
import sys
import json
from types import SimpleNamespace
from fractions import Fraction as Q
from flint import arb,ctx

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'independent_rows'))
import verify as sparse


def build_args(groups,tilts,precision,max_splits,target_bits,directory,updates=2):
    if type(updates) is not int or updates not in (2,3,4):raise ValueError('two through four inner updates required')
    return SimpleNamespace(groups=groups,tilts=tilts,penalties=['1'],weight_tilt='1',precision=precision,
        updates=updates,
        full_feedback=0,window_histogram=0,joint_cancellation=False,column_density=False,feedback_density=0,
        operator_cache=directory,retain_parents=True,joint_witness=True,joint_top=2,
        check_interval=1,
        probe_supports=[],probe_vector=[],max_splits=max_splits,target_bits=target_bits,screen_only=False)


def run(occupancies,threshold,tilts,precision,max_splits,target_bits,directory=None,output=None,inner='shape-uniform',
        shape_penalties=('1','.75','.5'),shape_weight_tilt='1',operator_through=None,exact_feedback=False,updates=2,
        joint_return_through=None,lazy_density_through=None,analytic_gradient=False):
    if (not occupancies or len(set(occupancies))!=len(occupancies) or any(type(q) is not int or not 1<=q<=2048 for q in occupancies)
            or not 0<=threshold<(1<<21) or precision<128 or not tilts or any(Q(t)<=0 for t in tilts)
            or inner not in ('shape-uniform','gf-shape-tilt','gf-occupancy','gf-density','gf-density-refined','gf-classes','gf-rank','gf-birth-classes')):
        raise ValueError('valid occupancies, cutoff, precision, and positive tilts required')
    if not shape_penalties or any(not 0<Q(p)<=1 for p in shape_penalties) or Q(shape_weight_tilt)<=0:
        raise ValueError('positive shape weights required')
    degree=max(occupancies) if operator_through is None else operator_through
    if type(degree) is not int or not max(occupancies)<=degree<=2048:
        raise ValueError('operator degree must cover every requested occupancy')
    if type(exact_feedback) is not bool or (exact_feedback and inner in ('shape-uniform','gf-shape-tilt')):
        raise ValueError('exact feedback requires a GF fixed-occupancy kernel')
    if inner in ('gf-rank','gf-birth-classes') and not exact_feedback:raise ValueError('GF rank/class operators require the exact-feedback option')
    if type(updates) is not int or updates not in (2,3,4) or (updates>2 and inner not in ('gf-rank','gf-birth-classes')):
        raise ValueError('three/four updates require a regenerated GF rank or birth-class kernel')
    if directory and inner in ('gf-occupancy','gf-density','gf-density-refined','gf-classes','gf-rank','gf-birth-classes'):
        raise ValueError('GF occupancy operators do not use the old operator cache')
    for value,limit in ((joint_return_through,4),(lazy_density_through,32)):
        if value is not None and (inner!='gf-birth-classes' or type(value) is not int or not 0<=value<=limit):
            raise ValueError('local refinements require a GF birth-class kernel and valid integer cutoffs')
    if type(analytic_gradient) is not bool:raise ValueError('boolean analytic-gradient option required')
    ctx.prec=precision
    print('GF16 ENSEMBLE: unweighted support covers, inner',inner,'updates',updates,'; no old certificate reuse',flush=True)
    sparse.mixing_test();sparse.folding_test();sparse.geometry_test();sparse.retained_test();sparse.placement_prefix_test()
    counts=sparse.integer_cdf(sparse.weighted_cdf_upper(sparse.authenticated_caps(),1<<128,full_weight=Q(1)))
    args=build_args(degree,tilts,precision,max_splits,target_bits,directory,updates)
    args.exact_feedback=exact_feedback
    args.joint_return_through=joint_return_through;args.lazy_density_through=lazy_density_through
    args.analytic_gradient=analytic_gradient
    if inner=='gf-shape-tilt':
        import shape_tilt
        operators=shape_tilt.build_operators(args,shape_penalties,shape_weight_tilt);terminal=sparse.TAIL_TERMINAL
    elif inner=='gf-birth-classes':
        import numpy as np
        import occupancy_birth_classes
        operators=occupancy_birth_classes.build_operators(args)
        terminal=np.ones(next(iter(operators.values()))[0][0].nrows())
    elif inner in ('gf-occupancy','gf-density','gf-density-refined','gf-classes','gf-rank'):
        import occupancy_kernel
        import density_kernel
        import occupancy_rank
        engine=occupancy_rank if inner=='gf-rank' else occupancy_kernel if inner=='gf-occupancy' else density_kernel
        if inner=='gf-classes':
            import numpy as np
            operators=engine.build_operators(args,refined=True,classes=True)
            terminal=np.ones(next(iter(operators.values()))[0][0].nrows());terminal[2]=0
        else:
            operators=engine.build_operators(args,refined=True) if inner=='gf-density-refined' else engine.build_operators(args)
            terminal=engine.TERMINAL
    else:operators=sparse.build_operators(args);terminal=sparse.TAIL_TERMINAL
    results={};total=arb(0);unresolved=[]
    for q in occupancies:
        args.groups=q
        print('GF16 occupancy',q,'cutoff',threshold,flush=True)
        upper=sparse.cover(args,operators,{p:counts for _,p in operators},terminal,cutoff=threshold)
        if upper is None:unresolved.append(q)
        else:
            total=sparse.up(total+upper)
            results[str(q)]=[int(x) for x in upper.upper().man_exp()]
        if output:
            record=dict(schema='gf16-sparse-run-1',threshold=threshold,precision=precision,tilts=tilts,inner=inner,
                target_bits=target_bits,max_splits=max_splits,requested=occupancies,
                shape_penalties=list(shape_penalties),shape_weight_tilt=shape_weight_tilt,operator_through=degree,
                exact_feedback=exact_feedback,updates=updates,
                joint_return_through=joint_return_through,lazy_density_through=lazy_density_through,
                analytic_gradient=analytic_gradient,
                verified_uppers=results,unresolved=unresolved+[x for x in occupancies if x>q],
                note='Run receipt only; regenerate support covers for independent verification.')
            output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(record,indent=2)+'\n')
    print('GF16 SPARSE RESULT verified',sorted(map(int,results)),'unresolved',unresolved,
          'verified aggregate log2',total.log()/arb(2).log() if total>0 else 'none',flush=True)
    print('This run does not cover other occupancies.',flush=True)
    return None if unresolved else total


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--through',type=int,default=64)
    parser.add_argument('--updates',type=int,choices=(2,3,4),default=2)
    parser.add_argument('--occupancies',nargs='+',type=int)
    parser.add_argument('--threshold',type=int,default=10485)
    parser.add_argument('--tilts',nargs='+',default=['.00032','.0032','.032'])
    parser.add_argument('--precision',type=int,default=192)
    parser.add_argument('--max-splits',type=int,default=100)
    parser.add_argument('--target-bits',type=int,default=52)
    parser.add_argument('--operator-cache')
    parser.add_argument('--operator-through',type=int,help='Build a longer checked placement prefix, for memo reuse')
    parser.add_argument('--inner',choices=('shape-uniform','gf-shape-tilt','gf-occupancy','gf-density','gf-density-refined','gf-classes','gf-rank','gf-birth-classes'),default='shape-uniform')
    parser.add_argument('--shape-penalties',nargs='+',default=['1','.75','.5'])
    parser.add_argument('--shape-weight-tilt',default='1')
    parser.add_argument('--exact-feedback',action='store_true',help='Use exact GF feedback atoms through eight packets')
    parser.add_argument('--joint-return-through',type=int,choices=(0,1,2,3,4),
        help='Regenerate joint output/return counts through this occupancy; requires gf-birth-classes')
    parser.add_argument('--lazy-density-through',type=int,
        help='Retain density through this local occupancy; requires gf-birth-classes')
    parser.add_argument('--analytic-gradient',action='store_true',help='Use analytic gradients for joint witness proposals only')
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    if not 1<=args.through<=2048 or args.max_splits<0 or args.target_bits<40:parser.error('invalid search limits')
    run(sorted(set(args.occupancies or range(1,args.through+1))),args.threshold,args.tilts,args.precision,
        args.max_splits,args.target_bits,args.operator_cache,args.output,args.inner,args.shape_penalties,args.shape_weight_tilt,
        args.operator_through,args.exact_feedback,args.updates,args.joint_return_through,args.lazy_density_through,args.analytic_gradient)


if __name__=='__main__':main()
