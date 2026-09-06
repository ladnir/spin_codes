"""Directed one-row tail bounds from cellwise brackets on active time.

Arb builds probability polynomials and Hoeffding enclosures. A positive,
unnormalized support-count DP uses outward binary64 operations. Exact final
aggregation gives true-tail brackets and a stronger deterministic M22 bound.
"""
from __future__ import annotations
import itertools
import json
import math
import sys
from fractions import Fraction
from pathlib import Path
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
import numpy as np
from flint import arb,ctx
from audit_bch_q1_full_arb import encode,decode,rational
from certify_bch_m25_closure import deterministic_caps,old_higher_occupations
from run_higher_endpoint_preflight import sha,write_new

ROOT=Path(__file__).resolve().parents[1]
GEN=ROOT/'generated'
OUTPUT=GEN/'bch256_q1_activity_cells_outward.json'
CELLS=16
LENGTH=8192//CELLS
CUTOFF=209716
MAX_CHARGE=(2*CUTOFF+25000)//LENGTH


def cell_polynomials(length,memory,mode,scalar):
    eps=scalar(1)/(1<<memory)
    stay=(1-eps)**length
    after=(1-eps)*(1-stay)/(length*eps)
    zero=[[[scalar(0) for _ in range(2)] for _ in range(2)] for _ in range(2)]
    one=[[[scalar(0) for _ in range(2)] for _ in range(2)] for _ in range(2)]
    zero[0][0][0]=scalar(1)
    if mode=='ceil':
        zero[1][0][1]=1-stay
        zero[1][1][1]=stay
        for start in range(2):
            one[start][0][1]=1-after
            one[start][1][1]=after
    else:
        assert mode=='floor'
        zero[1][0][0]=1-stay
        zero[1][1][1]=stay
        one[0][0][0]=one[1][0][0]=1-after
        one[0][1][0]=after
        one[1][1][0]=after-stay
        one[1][1][1]=stay
    return zero,one


def region_polynomials(length,cells,memory,mode,scalar):
    assert length%cells==0
    z,a=cell_polynomials(length//cells,memory,mode,scalar)
    zero=[[[scalar(0) for _ in range(cells+1)] for _ in range(2)] for _ in range(2)]
    one=[[[scalar(0) for _ in range(cells+1)] for _ in range(2)] for _ in range(2)]
    zero[0][0][0]=zero[1][1][0]=scalar(1)
    for n in range(cells):
        nz=[[[scalar(0) for _ in range(cells+1)] for _ in range(2)] for _ in range(2)]
        na=[[[scalar(0) for _ in range(cells+1)] for _ in range(2)] for _ in range(2)]
        for start,end,mid in itertools.product(range(2),repeat=3):
            for k in range(n+1):
                for charge in (0,1):
                    nz[start][end][k+charge]+=zero[start][mid][k]*z[mid][end][charge]
                    na[start][end][k+charge]+=(one[start][mid][k]*z[mid][end][charge]
                                            +zero[start][mid][k]*a[mid][end][charge])
        zero,one=nz,na
    return zero,[[[v/cells for v in row] for row in matrix] for matrix in one]


def scalar_endpoint(value,upper):
    endpoint=rational(value.upper() if upper else value.lower())
    endpoint=max(Fraction(0),endpoint)
    result=float(endpoint)
    if upper:
        if Fraction.from_float(result)<endpoint:
            result=np.nextafter(result,np.inf)
    elif Fraction.from_float(result)>endpoint:
        result=np.nextafter(result,0.)
    return float(result)


def arrays(poly,upper):
    return np.array([[[scalar_endpoint(v,upper) for v in row] for row in matrix] for matrix in poly])


def directed(value,upper):
    # All arguments are nonnegative. Downward rounding to zero is safe;
    # upward rounding repairs positive underflow and may safely inflate zero.
    return np.nextafter(value,np.inf if upper else 0.)


def counts(zero,one,positions,max_weight,max_charge,upper):
    cells=zero.shape[-1]-1
    current=np.zeros((max_weight+1,2,max_charge+1))
    current[0,0,0]=1
    for completed in range(positions):
        maximum=min(completed+1,max_weight)
        oldmax=min(completed,max_weight)
        width0=min(completed*cells,max_charge)+1
        updated=np.zeros_like(current)
        for start,end in itertools.product(range(2),repeat=2):
            for charge in range(cells+1):
                width=min(width0,max_charge+1-charge)
                if width<=0:
                    continue
                z=zero[start,end,charge]
                a=one[start,end,charge]
                if z:
                    dest=updated[:oldmax+1,end,charge:charge+width]
                    term=directed(current[:oldmax+1,start,:width]*z,upper)
                    dest[:]=directed(dest+term,upper)
                if a:
                    dest=updated[1:maximum+1,end,charge:charge+width]
                    term=directed(current[:maximum,start,:width]*a,upper)
                    dest[:]=directed(dest+term,upper)
        current=updated
    assert np.isfinite(current).all() and (current>=0).all()
    return current


def toy_check():
    # Independently enumerate all input routes and all state-reset coins at
    # bit granularity; unused reset coins are still sampled and marginalized.
    positions,length,cells,memory=2,4,2,2
    bins={mode:[[[Fraction(0) for _ in range(5)] for _ in range(2)] for _ in range(3)] for mode in ('ceil','floor')}
    exact_cdf=[Fraction(0) for _ in range(3)]
    cases=0
    for w in range(3):
        for support in itertools.combinations(range(positions),w):
            for locations in itertools.product(range(length),repeat=w):
                ones={r*length+j for r,j in zip(support,locations)}
                route=Fraction(1,math.comb(positions,w)*length**w)
                for mask in range(1<<(positions*length)):
                    active=0
                    actual=ceil=floor=0
                    probability=route
                    for cell in range(positions*cells):
                        entry=active
                        contains=False
                        no_reset=True
                        for j in range(length//cells):
                            i=cell*(length//cells)+j
                            bit=int(i in ones)
                            contains=contains or bool(bit)
                            survive=(mask>>i)&1
                            probability*=Fraction(3 if survive else 1,4)
                            no_reset=no_reset and bool(survive)
                            actual+=int(bool(active or bit))
                            active=int(bool(active or bit))*survive
                        ceil+=int(bool(entry or contains))
                        floor+=int(bool(entry and no_reset))
                    assert floor*(length//cells)<=actual<=ceil*(length//cells)
                    bins['ceil'][w][active][ceil]+=probability
                    bins['floor'][w][active][floor]+=probability
                    exact_cdf[w]+=probability*Fraction(sum(math.comb(actual,k) for k in range(min(2,actual)+1)),1<<actual)
                    cases+=1
    for mode in ('ceil','floor'):
        z,a=region_polynomials(length,cells,memory,mode,Fraction)
        rz,ra=([[[arb(v.numerator)/v.denominator for v in row] for row in matrix] for matrix in poly] for poly in (z,a))
        for upper in (False,True):
            result=counts(arrays(rz,upper),arrays(ra,upper),positions,2,4,upper)
            for w,end,k in itertools.product(range(3),range(2),range(5)):
                bound=Fraction.from_float(float(result[w,end,k]))/math.comb(positions,w)
                exact=bins[mode][w][end][k]
                assert bound>=exact if upper else bound<=exact
        for w in range(3):
            value=sum((bins[mode][w][end][k]*Fraction(sum(math.comb(2*k,j) for j in range(min(2,2*k)+1)),1<<(2*k))
                       for end in range(2) for k in range(5)),Fraction(0))
            assert value<=exact_cdf[w] if mode=='ceil' else value>=exact_cdf[w]
    tiny=np.nextafter(0.,np.inf)
    assert directed(np.array([tiny*.5]),True)[0]>0
    assert directed(np.array([tiny*.5]),False)[0]==0
    return dict(bit_level_paths_enumerated=cases,all_count_bins_bracketed=True,
                exact_binomial_tail_brackets_checked=True,underflow_directions_checked=True)


def hoeffding(count,lower):
    h=count*LENGTH
    if lower:
        if h<=CUTOFF:
            return arb(1)
        if h>=2*(CUTOFF+1):
            return arb(0)
        t=arb(CUTOFF+1)-arb(h)/2
        return 1-(-2*t*t/h).exp()
    if h<=2*CUTOFF:
        return arb(1)
    t=arb(h)/2-CUTOFF
    return (-2*t*t/h).exp()


def aggregate(cache):
    old=json.loads((GEN/'bch256_q1_full_arb_transfer.json').read_text())
    oa=json.loads((GEN/'oa21_closure_probe/audit.json').read_text())
    coefficients={int(w):decode(r['coefficient_upper']) for w,r in old['coefficient_rows'].items()}
    caps=deterministic_caps()
    for row in oa['shells']:
        caps[row['weight']]=caps[256-row['weight']]=row['cap']
    lower_factors=[max(Fraction(0),rational(hoeffding(k,True).lower())) for k in range(MAX_CHARGE+1)]
    upper_factors=[min(Fraction(1),rational(hoeffding(k,False).upper())) for k in range(MAX_CHARGE+1)]
    remainder=min(Fraction(1),rational(hoeffding(MAX_CHARGE+1,False).upper()))
    rows=[]
    for w in (38,40,42):
        def integrate(key,factors):
            return sum((Fraction.from_float(float(cache[key][w,end,k]))*factors[k]
                        for end in range(2) for k in range(MAX_CHARGE+1)),Fraction(0))/math.comb(256,w)
        lo=8192*integrate('ceil_lower',lower_factors)
        hi=8192*(integrate('floor_upper',upper_factors)+remainder)
        assert 0<lo<hi<coefficients[w]
        rows.append(dict(weight=w,coefficient_lower=encode(lo),coefficient_upper=encode(hi),
                         old_chernoff_coefficient_upper=encode(coefficients[w]),
                         recovered_bits_diagnostic=encode(coefficients[w])['log2_diagnostic']-encode(hi)['log2_diagnostic'],
                         maximum_further_tail_improvement_bits_diagnostic=encode(hi)['log2_diagnostic']-encode(lo)['log2_diagnostic']))
        coefficients[w]=hi
    q1=sum((caps[w]*v for w,v in coefficients.items()),Fraction(0))
    q2,q3,tail=old_higher_occupations()
    full=q1+Fraction(128,127)*(q2+q3+tail)
    h38=decode(oa['shells'][0]['physical_optimum'])
    obstruction=31*h38*decode(rows[0]['coefficient_lower'])
    assert obstruction>Fraction(1,1<<40)
    return dict(shells=rows,Q1_upper=encode(q1),full_first_moment_upper=encode(full),
                full_margin_bits_diagnostic=-encode(full)['log2_diagnostic'],
                exact_full_bound_below_2_to_minus_40=full<Fraction(1,1<<40),
                feasible_OA21_primal_weight38_true_tail_lower=encode(obstruction),
                OA21_relaxation_insufficient_even_with_exact_inner_tails=True,
                omitted_upper_count_mass_allowance=encode(remainder),original_M22_target_closed=False)


def sources():
    return [Path(__file__),ROOT/'code/audit_bch_q1_full_arb.py',ROOT/'code/certify_bch_m25_closure.py',
            GEN/'bch256_q1_full_arb_transfer.json',GEN/'oa21_closure_probe/audit.json',
            GEN/'bch256_closure_deterministic_envelope.json',GEN/'johnson_n256_w52_d38.json',
            GEN/'johnson_n256_w54_d38.json',GEN/'bch256_q2_positive_outward.json',
            GEN/'bch256_q3_capped_outward.json',GEN/'bch256_occupation_tail_outward.json']


def main():
    ctx.prec=192
    checking='--verify' in sys.argv
    if not checking:
        assert not OUTPUT.exists() and not OUTPUT.with_suffix('.npz').exists()
    toy=toy_check()
    result={}
    for mode,upper in (('ceil',False),('floor',True)):
        z,a=region_polynomials(8192,CELLS,22,mode,arb)
        rz,ra=arrays(z,upper),arrays(a,upper)
        label=mode+('_upper' if upper else '_lower')
        result[label]=counts(rz,ra,256,42,MAX_CHARGE,upper)
        result[label+'_zero']=rz
        result[label+'_one']=ra
        print('DIRECTED COUNTS COMPLETE',label,flush=True)
    receipt=dict(classification='Outward true one-row tail brackets using cell charges; deterministic M22 bound remains above target',
                 parameters=dict(memory_bits=22,outer_rows=8192,outer_regions=256,cutoff=CUTOFF,
                                 cells_per_region=CELLS,cell_length=LENGTH,maximum_charge=MAX_CHARGE),
                 toy_checks=toy,**aggregate(result),
                 source_sha256={str(p.relative_to(ROOT)):sha(p) for p in sources()})
    if checking:
        saved=np.load(OUTPUT.with_suffix('.npz'))
        assert set(saved.files)==set(result)
        for key in result:
            assert np.array_equal(saved[key],result[key])
        assert json.loads(OUTPUT.read_text())==receipt
    else:
        np.savez_compressed(OUTPUT.with_suffix('.npz'),**result)
        write_new(OUTPUT,receipt)
    print(json.dumps({k:receipt[k] for k in ('full_margin_bits_diagnostic','exact_full_bound_below_2_to_minus_40',
                      'OA21_relaxation_insufficient_even_with_exact_inner_tails','toy_checks')},indent=2))
    for row in receipt['shells']:
        print(row['weight'],row['coefficient_lower']['log2_diagnostic'],row['coefficient_upper']['log2_diagnostic'])


if __name__=='__main__':
    main()
