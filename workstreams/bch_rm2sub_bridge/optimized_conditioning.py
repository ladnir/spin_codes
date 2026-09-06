"""Optimize auxiliary probabilities, retaining exact rational conditioning."""
import argparse
import itertools
import math
from fractions import Fraction as F
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from flint import arb,ctx
import bridge as base
import occupation_three as q3


def log_moment(region,ps):
    a,b,c=ps
    mass=((1-a)*(1-b)*(1-c),a*(1-b)*(1-c)+(1-a)*b*(1-c)+(1-a)*(1-b)*c,
          a*b*(1-c)+a*(1-b)*c+(1-a)*b*c,a*b*c)
    matrix=sum(p*r for p,r in zip(mass,region)); scale=0.
    for _ in range(8):
        matrix=matrix@matrix
        maximum=float(matrix.max()); matrix/=maximum
        scale=2*scale+math.log(maximum)
    return math.log(float(matrix[0].sum()))+scale


def screen(probe=False):
    t,s,spectrum=base.load_map(q3.NAME); caps=q3.q2.deterministic_caps()
    weights=[np.array(band) for band in q3.BANDS]
    logs=[np.array([math.log(caps[w])-math.log(math.comb(256,w)) for w in band]) for band in q3.BANDS]
    mid=[(band[0]+band[-1])/512 for band in q3.BANDS]
    boxes=[(0,0,0),(3,4,4),(4,4,4),(6,6,6)] if probe else list(itertools.combinations_with_replacement(range(len(mid)),3))
    best={box:(math.inf,None,None) for box in boxes}
    for tenth in range(-85,-54,5):
        lam=math.exp(tenth/10)
        region=[np.array(m).reshape(3,3) for m in q3.region_matrices(t,s,spectrum,math.exp(-lam),float)]
        for box in boxes:
            def evaluate(shift,ret=False):
                factor=math.exp(shift)
                ps=[mid[i]*factor/(1-mid[i]+mid[i]*factor) for i in box]
                cost=0.
                for i,p in zip(box,ps):
                    if p==1:
                        assert q3.BANDS[i]==(256,)
                        cost+=math.log(caps[256])
                    else:
                        cost+=float(np.logaddexp.reduce(logs[i]-weights[i]*math.log(p)-(256-weights[i])*math.log1p(-p)))
                value=log_moment(region,ps)+cost+209716*lam
                return (value,ps) if ret else value
            opt=minimize_scalar(evaluate,bounds=(-2.,8.),method='bounded',options={'xatol':0.002})
            value,ps=evaluate(float(opt.x),True)
            if value<best[box][0]:best[box]=(value,tenth,ps)
        print('Optimized Q3 tilt',tenth,flush=True)
    logsum=base.lse([value[0]+math.log(q3.multiplicity(box)) for box,value in best.items()])+math.log(math.comb(8192,3))
    print('Optimized conditioning margin',-logsum/math.log(2),flush=True)
    if probe:
        print(best);return
    payload=dict(status='OPTIMIZED_Q3_CONDITIONING_SCREEN_ONLY',margin_bits_diagnostic=-logsum/math.log(2),
        boxes=[dict(bands=box,witness_tenth=v[1],p=[base.encode(F.from_float(p)) for p in v[2]]) for box,v in best.items()],
        source_sha256={p.name:base.sha(p) for p in (Path(__file__),Path(q3.__file__))})
    base.write_new(base.HERE/'generated'/'t128_s15_q3_optimized_screen.json',payload)


def certify():
    import sys
    sys.path.insert(0,str(base.BCH/'code'))
    from audit_bch_q1_full_arb import rational
    path=base.HERE/'generated'/'t128_s15_q3_optimized_screen.json'; discovery=base.read(path)
    for filename,digest in discovery['source_sha256'].items():assert base.sha(base.HERE/filename)==digest
    t,s,spectrum=base.load_map(q3.NAME); caps=q3.q2.deterministic_caps(); ctx.prec=256
    kernel=base.read(base.HERE/'inputs'/f'{q3.NAME}_b_kernel_spectrum.json')['by_total_weight']
    assert all(row['kernel_words']==0 for row in kernel if 1<=row['total_weight']<=3)
    assert sorted(w for band in q3.BANDS for w in band)==list(base.WEIGHTS)
    assert [tuple(r['bands']) for r in discovery['boxes']]==list(itertools.combinations_with_replacement(range(len(q3.BANDS)),3))
    tables={}; total=F(0); rows=[]
    for row in discovery['boxes']:
        box=tuple(row['bands']); tenth=row['witness_tenth']; ps=[base.decode(p) for p in row['p']]
        assert len(ps)==3 and -120<=tenth<=0
        assert all(0<p<=1 and (p<1 or q3.BANDS[i]==(256,)) for i,p in zip(box,ps))
        if tenth not in tables:
            lam=(arb(tenth)/10).exp()
            tables[tenth]=(q3.region_matrices(t,s,spectrum,(-lam).exp(),arb),(209716*lam).exp())
        region,correction=tables[tenth]
        upper=rational((q3.coefficient(region,ps,arb)*correction).upper())
        costs=[q3.condition_cost(q3.BANDS[i],p,caps) for i,p in zip(box,ps)]
        contribution=math.comb(8192,3)*q3.multiplicity(box)*math.prod(costs)*upper
        assert contribution>0
        total+=contribution
        rows.append(dict(bands=box,witness_tenth=tenth,p=row['p'],moment_upper=base.encode(upper),
                         contribution_upper=base.encode(contribution)))
    local,outer=q3.q2.source_paths(); local += [Path(__file__),Path(q3.__file__),path]
    payload=dict(status='OUTWARD_OPTIMIZED_Q3_CONDITIONING_BOUND',configuration=q3.NAME,
        parameters=dict(message_bits=1<<20,output_bits=1<<21,outer_length=256,outer_rows=8192,
                        step_bits=t,state_bits=s,occupation=3,distance_cutoff=209716),
        bands=q3.BANDS,boxes=rows,Q3_upper=base.encode(total),
        margin_bits_diagnostic=math.log2(total.denominator)-math.log2(total.numerator),
        Q3_below_2_to_minus_60=total<F(1,1<<60),all_occupations_certified=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in local},
        outer_sha256={str(p.relative_to(base.BCH)):base.sha(p) for p in outer})
    base.write_new(base.HERE/'generated'/'t128_s15_q3_optimized_outward.json',payload)
    print('OUTWARD optimized Q3 margin',payload['margin_bits_diagnostic'],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=('probe','screen','certify'))
    args=parser.parse_args()
    if args.mode=='certify':certify()
    else:screen(args.mode=='probe')
