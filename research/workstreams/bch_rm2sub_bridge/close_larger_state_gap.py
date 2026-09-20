"""Discover short shared witnesses, certify every depth, and replay without optimization.

The construction is fixed at t64_s20. Failed depths remain in each shard;
only outward bounds below 2^-70 enter the final complete-range receipt.
"""
import argparse
import math
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import numpy as np
from scipy.optimize import minimize_scalar
from flint import arb, ctx
import bridge as base
import larger_state_maps as maps
import polynomial_regions as poly
import occupation_three as groups
import screen_exponential_modes as cap_loader
import scaled_adaptive as scaled
import positive_line_hull as hull
import activation_bridge as q1
from screen_constant_row_split import moment
from certify_larger_state_range import costs_for
from audit_bch_q1_full_arb import rational


def dyadic_power(value):
    assert value > 0
    n,d = value.numerator,value.denominator
    power = n.bit_length()-d.bit_length()
    candidate = F(1 << power) if power >= 0 else F(1,1 << (-power))
    return max(-80, power+int(candidate < value))


def as_fraction(power):
    return F(1 << power) if power >= 0 else F(1,1 << (-power))


def choose(region, q, tilt, bands, caps):
    selected = np.array([[float(v.log()) if v > 0 else -math.inf for v in row]
                         for row in region[:q+1]]).reshape(-1,3,3)
    ps,gammas = [],[]
    for band in bands:
        if tuple(band) == (256,):
            ps.append(1.);gammas.append(0.);continue
        w = np.array(band)
        v = np.array([math.log(caps[x])-math.log(math.comb(256,x)) for x in band])
        def objective(theta):
            p = 1/(1+math.exp(-theta))
            gamma = float(np.max(v-w*math.log(p)-(256-w)*math.log1p(-p)))
            return moment(selected,q,p)+q*gamma
        optimum = minimize_scalar(objective,bounds=(-6,14),method='bounded',options={'xatol':1e-6})
        p = 1/(1+math.exp(-float(optimum.x)))
        ps.append(p)
        gammas.append(float(np.max(v-w*math.log(p)-(256-w)*math.log1p(-p))))
    ps = np.array(ps);roots = np.exp(np.array(gammas)/256)
    keep = hull.indices(roots*(1-ps),roots*ps)
    mantissas,exponents = scaled.initial(region[:q+1],outward=False)
    value = scaled.terminal_log(mantissas,exponents,ps[keep],roots[keep])
    value += base.CUTOFF*math.exp(tilt/10)+math.log(math.comb(8192,q))+q*math.log(len(bands))
    return dict(anchor=q, tilt=tilt, p=[base.encode(F.from_float(float(p))) for p in ps],
                margin_bits_diagnostic=-value/math.log(2))


def evaluate(region, witness, lower, upper, bands, caps):
    ps = [base.decode(v) for v in witness['p']]
    costs = costs_for(bands,ps,caps)
    roots = [((arb(c.numerator)/c.denominator).log()/256).exp().upper() for c in costs]
    probs = [arb(p.numerator)/p.denominator for p in ps]
    left = np.nextafter(np.array([float((r*(1-p)).upper()) for r,p in zip(roots,probs)]),np.inf)
    right = np.nextafter(np.array([float((r*p).upper()) for r,p in zip(roots,probs)]),np.inf)
    keep = hull.indices(left,right)
    mantissas,exponents = scaled.initial(region[:upper+1])
    lam = (arb(witness['tilt'])/10).exp();correction = (base.CUTOFF*lam).exp()
    result = []
    for q,matrix,exponent in scaled.matrices(mantissas,exponents,left[keep],right[keep]):
        if q < lower:continue
        value = tuple(arb(float(v)) for v in matrix.flat)
        for _ in range(8):value = q1.positive_mul(value,value)
        bound = math.comb(8192,q)*len(bands)**q*rational((sum(value[:3],arb(0))*correction*arb(2)**(256*exponent)).upper())
        result.append((q,bound,dyadic_power(bound)))
    return result


def run(lower,upper,tag,verify=False):
    assert 2 <= lower <= upper <= 8192 and tag.isidentifier()
    directory = base.HERE/'generated'/f'larger_gap_{tag}'
    output = base.HERE/'generated'/f'larger_range_{tag}_outward.json'
    t,s,spectrum,kernel = maps.load('t64_s20')
    caps,cap_sources = cap_loader.latest_caps();bands = groups.BANDS
    assert sorted(w for band in bands for w in band) == list(base.WEIGHTS)
    assert caps[256] == 1 and tuple(bands[-1]) == (256,)
    files = [Path(__file__),Path(poly.__file__),Path(groups.__file__),Path(cap_loader.__file__),
             Path(scaled.__file__),Path(hull.__file__),Path(q1.__file__),
             base.HERE/'screen_constant_row_split.py',base.HERE/'certify_larger_state_range.py',
             base.HERE/'general_occupancy.py',base.HERE/'tightened_occupancy.py',
             base.HERE/'general_batch_certificate.py',base.HERE/'certify_q3_compact.py']+cap_sources
    hashes = {**maps.sources(),**{str(p.relative_to(base.HERE)):base.sha(p) for p in files}}
    ctx.prec = 512 if verify else 256
    @lru_cache(maxsize=16)
    def region(tilt):
        lam = (arb(tilt)/10).exp()
        return poly.regions(t,s,spectrum,kernel,(-lam).exp(),upper)
    best = {}; shards = []
    if verify:
        old = base.read(output)
        assert old['configuration'] == 't64_s20' and old['occupancy_range'] == [lower,upper]
        for source,digest in old['local_sha256'].items():assert base.sha(base.HERE/source) == digest
        for index,name in enumerate(old['shards']):
            saved = base.read(base.HERE/name)
            assert saved['local_sha256'] == hashes
            lo,hi = saved['occupancy_range']
            values = evaluate(region(saved['witness']['tilt']),saved['witness'],lo,hi,bands,caps)
            assert len(values) == len(saved['rows'])
            for (q,bound,_),row in zip(values,saved['rows']):
                assert row['occupation'] == q and bound <= as_fraction(row['upper_power_of_two'])
                power = row['upper_power_of_two']
                if power <= -70:best[q] = min(best.get(q,power),power)
            if index % 5 == 0:print('Replay shard',index,'covered',len(best),flush=True)
        assert sorted(best) == list(range(lower,upper+1))
        assert [dict(occupation=q,upper=base.encode(as_fraction(best[q]))) for q in sorted(best)] == old['rows']
        assert sum((as_fraction(v) for v in best.values()),F(0)) == base.decode(old['range_upper'])
        base.write_new(output.with_name(f'larger_range_{tag}_replay.json'),dict(
            status='LARGER_STATE_512_BIT_RANGE_REPLAY_PASSED',rows_checked=len(best),
            producer_sha256=base.sha(output),verifier_sha256=base.sha(Path(__file__))))
        print('Complete gap replay passed',lower,upper,flush=True);return
    assert not directory.exists() and not output.exists()
    directory.mkdir()
    trials = []
    while len(best) < upper-lower+1:
        anchor = next(q for q in range(lower,upper+1) if q not in best)
        prediction = round((-76+6*math.log2(anchor))/5)*5
        tilts = sorted(set(max(-100,min(15,prediction+delta)) for delta in (-5,0,5)))
        choices = [choose(region(tilt),anchor,tilt,bands,caps) for tilt in tilts]
        witness = max(choices,key=lambda row:row['margin_bits_diagnostic'])
        hi = min(upper,max(anchor,int(anchor*1.12)))
        values = evaluate(region(witness['tilt']),witness,anchor,hi,bands,caps)
        rows = [dict(occupation=q,upper_power_of_two=power,
                     margin_bits_diagnostic=math.log2(bound.denominator)-math.log2(bound.numerator)) for q,bound,power in values]
        shard = directory/f'shard{len(shards):03d}.json'
        base.write_new(shard,dict(status='OUTWARD_GAP_WITNESS',occupancy_range=[anchor,hi],witness=witness,
            discovery_trials=choices,rows=rows,local_sha256=hashes))
        shards.append(shard)
        for q,_,power in values:
            if power <= -70:best[q] = min(best.get(q,power),power)
        print('Gap',anchor,hi,'tilt',witness['tilt'],'anchor margin',round(rows[0]['margin_bits_diagnostic'],3),
              'covered',len(best),'of',upper-lower+1,flush=True)
        if anchor not in best:
            print('FAILED anchor retained; need a different witness search',anchor,flush=True)
            return
    total = sum((as_fraction(v) for v in best.values()),F(0))
    base.write_new(output,dict(status='LARGER_STATE_OUTWARD_RANGE',configuration='t64_s20',
        occupancy_range=[lower,upper],rows=[dict(occupation=q,upper=base.encode(as_fraction(best[q]))) for q in sorted(best)],
        range_upper=base.encode(total),range_below_2_to_minus_60=total<F(1,1<<60),all_occupations_certified=False,
        parameters=dict(message_bits=1<<20,output_bits=1<<21,step_bits=t,state_bits=s,distance_cutoff=base.CUTOFF),
        shards=[str(p.relative_to(base.HERE)) for p in shards],
        local_sha256={**hashes,**{str(p.relative_to(base.HERE)):base.sha(p) for p in shards}}))
    print('Complete gap produced',lower,upper,'shards',len(shards),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--lower',type=int,required=True);parser.add_argument('--upper',type=int,required=True)
    parser.add_argument('--tag',required=True);parser.add_argument('--verify',action='store_true')
    args = parser.parse_args();run(args.lower,args.upper,args.tag,args.verify)
