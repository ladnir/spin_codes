"""Exact-repaired positive-mixture alternatives for shared GF16 packets.

Normal runs regenerate authenticated BCH premises. The explicitly exploratory
cached-caps mode reuses a receipt instead and requires later regeneration.
Every proposed positive majorant is checked exactly against its shell caps.
Selected scalar means are diagnostics, never a whole-code certificate.
"""
import argparse
import json
from fractions import Fraction as Q
from pathlib import Path
from flint import ctx

import shared_support
import shared_mixture
import positive_prune
import scalar_cover
import birth_classes
from shared_relaxed_alternative import shell_caps


def candidates(caps,step=8,zero_bits=64,methods=None):
    centers=sorted({Q(u,256) for u in range(38,257,step)}|{Q(1)})
    costs=(Q(1,8),Q(3,16),Q(1,4),Q(3,8),Q(1,2))
    supported={'greedy-quarter'}|{prefix+str(cost) for prefix in ('pruned-','pooled-') for cost in costs}
    selected=supported if methods is None else set(methods)
    if not selected or not selected<=supported:raise ValueError('known nonempty mixture method selection required')
    baseline=shared_mixture.envelope(caps,centers,zero_bits,Q(1,4))
    if 'greedy-quarter' in selected:yield 'greedy-quarter',baseline,{}
    for cost in costs:
        if 'pruned-'+str(cost) not in selected:continue
        pruned,info=positive_prune.prune(caps,baseline,cost)
        yield 'pruned-'+str(cost),pruned,info
    if not any(name.startswith('pooled-') for name in selected):return
    majorants=[baseline]+[shared_mixture.envelope(caps,centers,zero_bits,cost)
               for cost in (Q(1,8),Q(3,16),Q(1,2),Q(1))]
    pooled=positive_prune.pool(caps,majorants)
    for cost in costs:
        if 'pooled-'+str(cost) not in selected:continue
        pruned,info=positive_prune.prune(caps,pooled,cost)
        yield 'pooled-'+str(cost),pruned,info


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--distance',default='.07')
    parser.add_argument('--means',nargs='+',default=['.008','.032'])
    parser.add_argument('--minimum-groups',type=int,default=33)
    parser.add_argument('--base-tilt',default='3/16')
    parser.add_argument('--precision',type=int,default=192)
    parser.add_argument('--variance-bins',type=int,default=16)
    parser.add_argument('--moment-shells',action='store_true')
    parser.add_argument('--methods',nargs='+')
    parser.add_argument('--cached-caps',type=Path,help='Exploratory only: reuse a prior shell-cap receipt; later replay must regenerate premises')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if not 0<Q(args.distance)<1 or any(not 0<Q(m)<1 for m in args.means) or args.precision<128:
        parser.error('valid distance, interior means and precision>=128 required')
    if args.output.exists():parser.error('refusing to overwrite an existing receipt')
    if args.cached_caps:
        cached=json.loads(args.cached_caps.read_text())
        caps=list(map(int,cached['shell_caps']))
        if len(caps)!=257 or caps[0]!=0:parser.error('valid nonzero shared-group shell caps required')
    else:
        caps,_=shared_support.shared_counts(True,True)
    if args.moment_shells:
        moments=shell_caps(256,128,30)
        caps=[min(a,b) for a,b in zip(caps,moments)]
    data=birth_classes.actual(2)
    ctx.prec=args.precision
    rows=[]
    for name,mixture,info in candidates(caps,methods=args.methods):
        shared_mixture.verify(caps,mixture)
        print('EXACT SHARED MIXTURE',name,len(mixture),'components',info,flush=True)
        model=scalar_cover.Model(shared_mixture.as_components(mixture),data,int(Q(args.distance)*scalar_cover.N),
             args.minimum_groups,Q(args.base_tilt),inner=birth_classes,variance_shuffle=True,
             variance_bins=args.variance_bins,regional_count=True)
        probes=[]
        row=dict(name=name,info=info,mixture=[dict(mass=str(c),activity=str(p)) for c,p in mixture],probes=probes)
        rows.append(row)
        for mean in map(Q,args.means):
            cell=(mean,mean)
            if model.empty(cell):continue
            proposal,witness=model.proposal(cell)
            upper=model.outward(cell,witness)
            print('PROBE',name,args.distance,str(mean),'proposal',proposal,'outward',upper.log()/scalar_cover.arb(2).log(),flush=True)
            probes.append(dict(mean=str(mean),proposal=proposal,witness=witness,
                               upper=[int(v) for v in upper.upper().man_exp()]))
            args.output.parent.mkdir(parents=True,exist_ok=True)
            args.output.write_text(json.dumps(dict(schema='shared-relaxed-alternative-mixtures-1',
                parameters={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},
                shell_caps=caps,rows=rows,proof_status='exact majorants and selected-mean bounds only; not a full certificate'),indent=2)+'\n')


if __name__=='__main__':main()
