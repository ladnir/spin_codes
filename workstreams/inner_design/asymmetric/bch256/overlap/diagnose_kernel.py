"""Separate exact zero-state mass from transfer slack at a fixed witness.

Alternative feedback maps are diagnostic candidates, not replacements or full
distance certificates. The outer scalar upper bound is common to all rows.
"""
import argparse
from fractions import Fraction as F
import itertools
import random
from pathlib import Path
from flint import arb,ctx
import overlap_transfer as overlap
model=overlap.model


class Candidate(overlap.original.Checker):
    def __init__(self,columns):
        super().__init__(20)
        engine=self.engine
        assert len(columns)==128 and len(set(columns))==128
        assert model.independent.search.rank(columns)==19
        maps=model.independent.g.tv.fixed.maps
        spectrum=maps.spectrum(maps.generators(columns,19))
        kernel=maps.dual_spectrum(spectrum,128,19)
        engine.columns=list(columns)
        engine.b_spectrum={w:n for w,n in spectrum.items() if w}
        engine.kernel=[kernel.get(j,0) for j in range(129)]

    def _fixed_moment(self,index,r):
        lam=(arb(index)/40).exp()
        matrix=self.engine.bernoulli(model.number(r),lam)
        moment=model.independent.terminal(matrix,self.engine.n,1<<self.power)
        return model.up(self.cutoff*lam+moment.log())


def weight_three_lower(rho):
    """Convexity at integer weights, within every dual input-weight shell."""
    import math
    assert 0<rho<1
    result=arb(0)
    for j in range(20):
        count=math.comb(19,j)
        # Every weight-three column is odd on exactly this many weight-j states.
        weight_sum=128*sum(math.comb(3,h)*math.comb(16,j-h) for h in (1,3) if 0<=j-h<=16)
        low,remainder=divmod(weight_sum,count)
        result+=(count-remainder)*rho**low+remainder*rho**(low+1)
    return result


def evaluate(columns,witness,q,v):
    checker=Candidate(columns)
    r=model.number(model.base.decode(witness['fixed_r']))
    lam=(arb(witness['tilt'])/40).exp()
    matrix=checker.engine.bernoulli(r,lam)
    upper=checker.bound(q,q,v,v,witness)
    moment=checker._fixed_moment(witness['tilt'],model.base.decode(witness['fixed_r']))
    scalar=upper-moment
    zero=checker.cutoff*lam+(1<<checker.power)*matrix[0].log()
    z=(-lam).exp();g0=1-r+r*z
    kernel_excess=(matrix[0].log()-128*g0.log())/arb(2).log()+19
    return checker,dict(columns=columns,dual_minimum_weight=min(checker.engine.b_spectrum),
        zero_kernel_weight4=checker.engine.kernel[4],margin_bits=float(-upper/arb(2).log()),
        log_upper=model.pack(upper.exp()),kernel_excess_bits_per_epoch=float(kernel_excess),
        moment_over_zero_path_bits=float((moment-zero)/arb(2).log()),
        zero_path_expression_margin_bits=float(-(scalar+zero)/arb(2).log()))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=overlap.HERE/'KERNEL_DIAGNOSIS.json')
    p.add_argument('--verify',action='store_true')
    args=p.parse_args();saved=model.base.read(args.output) if args.verify else None
    if saved:model.authenticate(saved)
    else:assert not args.output.exists()
    ctx.prec=512 if saved else 256
    gap=model.base.read(model.HERE/'DENSE_GAP_POINT.json')
    witness=gap['witness'];q=gap['occupation'];v=F(gap['coordinate'])
    reference=model.Engine(20)
    if saved:
        candidates=[(row['name'],row['columns']) for row in saved['results']]
    else:
        weight5=[sum(1<<i for i in support) for support in itertools.combinations(range(19),5)]
        candidates=[('current_weight3',reference.columns),('balanced_feedback',reference.a_columns)]
        candidates += [(f'weight5_seed{seed}',random.Random(seed).sample(weight5,128)) for seed in (0,1)]
    results=[]
    for i,(name,columns) in enumerate(candidates):
        checker,result=evaluate(columns,witness,q,v)
        result['name']=name
        if saved:assert model.unpack(result['log_upper'])<=model.unpack(saved['results'][i]['log_upper'])
        results.append(result)
        print(ctx.prec,name,'margin',result['margin_bits'],'zero expression',result['zero_path_expression_margin_bits'],flush=True)
    # This floor concerns the fixed-reference proof expression, not the actual
    # probability of a bad code or the optimum over all outer/input witnesses.
    checker=overlap.Checker(20)
    r=model.number(model.base.decode(witness['fixed_r']));lam=(arb(witness['tilt'])/40).exp();z=(-lam).exp()
    g0=1-r+r*z;rho=1-2*r*z/g0
    universal=weight_three_lower(rho).log()/arb(2).log()
    new_moment=checker._fixed_moment(witness['tilt'],model.base.decode(witness['fixed_r']))
    scalar=checker.bound(q,q,v,v,witness)-new_moment
    ideal_zero=checker.cutoff*lam+(1<<checker.power)*(128*g0.log()-19*arb(2).log())
    floor=scalar+ideal_zero+(1<<checker.power)*universal*arb(2).log()
    diagnosis=dict(all_weight3_kernel_excess_lower_bits=float(universal),
                   all_weight3_fixed_expression_margin_ceiling=float(-floor/arb(2).log()),
                   claim='Fixed outer scalar, reference probability, and output tilt only; not a code impossibility claim')
    if saved:
        model.base.write_new(args.output.with_name(args.output.stem+'_replay.json'),
                             dict(status='KERNEL_DIAGNOSIS_512_BIT_REPLAY_PASSED',producer_sha256=model.base.sha(args.output)))
    else:
        sources=model.sources()
        for path in (overlap.HERE/'PACKING_AUDIT.json',model.HERE/'DENSE_GAP_POINT.json'):
            sources[path.relative_to(model.ROOT).as_posix()]=model.base.sha(path)
        model.base.write_new(args.output,dict(status='FIXED_POINT_KERNEL_DIAGNOSIS_NOT_FULL_CERTIFICATE',
            witness=witness,occupation=q,coordinate=str(v),results=results,weight3_floor=diagnosis,
            full_distance_proved=False,source_sha256=sources))
    print(diagnosis,flush=True)


if __name__=='__main__':main()
