"""Directed run-mass envelope for a band of actual source sizes1..4096.

At common z,alpha, the caller supplies own-normalized regions R_(r,f),
V_(r,f) and positive scales s_r=Omega_r^(alpha/256). build_envelope
takes entrywise maxima of s_r times each region. Missing marginal indices
contribute zero. Sparse mappings support selected-cell pilots; band_join
requires complete original coverage before claiming a complete band.

For d=a+b, C_d=sum_r N_r(r-d), with invalid run counts zero. The result is
sum_d C_d max_(a+b=d) e0(Rbar_a Rbar_b)^127 Rbar_a Vbar_b.
There are256 regional factors, including one output-only terminal, so
the source scales have already contributed Omega_r^alpha. DO NOT apply
another source factor. The caller applies the common z^(-alpha*L) once.
Counts stay outside alpha. Ordered pairs retain both terminal orders.

Use ACTUAL full-pair indices, not holes. A hole-indexed calculation needs
a different d regrouping and is intentionally unsupported here. This
helper does not prove that arbitrary caller-supplied matrices/scales have
the required local meanings; it validates their shape and coverage.

APIs: build_envelope(...)->Envelope (dict subclass); band_counts(...)->
exact integer tuple indexed by d; band_join(...)->dict with bound,
degree_bounds,report. The caller owns time, memory and disk guards.
At most17 distinct weights, each<=4096. No work runs on import.
"""

from collections import Counter
from collections.abc import Mapping
from fractions import Fraction as F
import hashlib
import importlib.util
from math import comb
from pathlib import Path
import sys
import time

sys.dont_write_bytecode=True
JOIN_PIN="54f2ff06a1476ccf84fa4db1b757f1d5f8e4e814a1f0f733dc73ccf18b2ddabb"
POINTS,MAX_WEIGHTS=8192,17


class Envelope(dict):
    """Sparse region mapping with immutable original-coverage metadata."""
    __slots__=("weights","coverage","input_intersections")

    def __init__(self,weights,coverage,intersections=0):
        super().__init__()
        self.weights=tuple(weights)
        self.coverage=tuple((r,tuple(coverage[r])) for r in weights)
        self.input_intersections=intersections


def _weights(weights,points=POINTS):
    values=tuple(weights)
    if type(points) is not int or points<2 or points%2:
        raise ValueError("points must be a positive even integer")
    if not 1<=len(values)<=MAX_WEIGHTS or len(set(values))!=len(values):
        raise ValueError("require1..17 distinct actual source sizes")
    if any(type(r) is not int or not 1<=r<=min(4096,points//2) for r in values):
        raise ValueError("actual weights must lie in1..min(4096,points/2); no holes")
    return tuple(sorted(values))


def load_join():
    path=Path(__file__).with_name("cheap_lane_directed_run_mass.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=JOIN_PIN:
        raise RuntimeError("frozen run-mass helper mismatch")
    spec=importlib.util.spec_from_file_location("band_frozen_run_mass",path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def _region(value,arb):
    if len(value)!=2 or len(value[0])!=4 or len(value[1])!=2:
        raise ValueError("expected normalized flat(R4,V2), without log scales")
    rows=[]; intersections=0
    for row in value:
        clean=[]
        for x in row:
            if not isinstance(x,arb) or not x.is_finite() or x.upper()<0:
                raise ArithmeticError("regional entries must be finite nonnegative Arb bounds")
            if not x>=0: intersections+=1
            clean.append(x.nonnegative_part())
        rows.append(tuple(clean))
    return tuple(rows),intersections


def build_envelope(regions_by_weight,scales,guard):
    """Sparse {r:{f:(R4,V2)}} -> Envelope; scales contain s_r, not Omega_r.

    Every provided region is already divided by its OWN actual-r route
    normalizer. Endpoint fit proxies never change these integer indices.
    Missing indices are zero for the maximum, not a certificate that an
    omitted valid marginal vanishes. band_join rejects incomplete input.
    """
    from flint import arb
    if not isinstance(regions_by_weight,Mapping) or not isinstance(scales,Mapping):
        raise TypeError("regions and scales must be weight-indexed mappings")
    weights=_weights(regions_by_weight)
    if set(scales)!=set(weights):
        raise ValueError("source-scale weights must exactly match region weights")
    coverage={}; result=Envelope(weights,{r:() for r in weights})
    intersections=0
    zero=((arb(0),)*4,(arb(0),)*2)
    for r in weights:
        guard(); scale=scales[r]; regions=regions_by_weight[r]
        if not isinstance(scale,arb) or not scale.is_finite() or not scale>0:
            raise ArithmeticError("each source scale must be a finite positive Arb enclosure")
        if not isinstance(regions,Mapping):
            raise TypeError("each weight needs a sparse actual-full-index mapping")
        indices=tuple(sorted(regions))
        if any(type(f) is not int or not 0<=f<=r//2 for f in indices):
            raise ValueError("invalid actual full-pair index")
        coverage[r]=indices
        for f in indices:
            guard(); clean,number=_region(regions[f],arb); intersections+=number
            previous=result.get(f,zero)
            result[f]=tuple(tuple(old.max(scale*x).nonnegative_part()
                for old,x in zip(old_row,row)) for old_row,row in zip(previous,clean))
    result.coverage=tuple((r,coverage[r]) for r in weights)
    result.input_intersections=intersections
    guard()
    return result


def _band_counts(weights,points,join,guard):
    weights=_weights(weights,points); masses=[0]*max(weights)
    for r in weights:
        guard(); counts=join.run_counts(points=points,weight=r)
        if len(counts)!=r+1 or counts[0]!=0:
            raise ArithmeticError("unexpected positive-run count interface")
        for c in range(1,r+1):
            if c%64==0: guard()
            masses[r-c]+=counts[c]
    if sum(masses)!=sum(comb(points,r) for r in weights):
        raise ArithmeticError("band count checksum failed")
    return tuple(masses)


def band_counts(weights,points=POINTS):
    """Exact C_d vector; no source factor and no regional matrix work.

    Small even points are supported for exact toys. Production uses8192.
    The checksum includes all actual indexed supports, not necklaces.
    """
    return _band_counts(weights,points,load_join(),lambda:None)


def band_join(envelope,weights,join,guard,progress=None):
    """Complete band join. Returns bound, d-indexed degree_bounds, report.

    Input must come from build_envelope with every original valid index
    present. Partial envelopes are only for selected controls. No source,
    alpha, marker or z correction is introduced here.
    """
    from flint import arb
    started=time.monotonic(); weights=_weights(weights)
    if not isinstance(envelope,Envelope) or envelope.weights!=weights:
        raise ValueError("use build_envelope with exactly these actual weights")
    if any(indices!=tuple(range(r//2+1)) for r,indices in envelope.coverage):
        raise ValueError("incomplete per-weight marginals cannot enter a complete band join")
    full=max(weights)//2
    if set(envelope)!=set(range(full+1)):
        raise ValueError("missing full-index union coverage")
    clean=[]; intersections=0
    for f in range(full+1):
        guard(); region,number=_region(envelope[f],arb)
        clean.append(region); intersections+=number
    counts=_band_counts(weights,POINTS,join,guard)
    degree_bounds=[]; total=arb(0); ordered=unordered=0
    for d,number in enumerate(counts):
        guard(); maximum=arb(0)
        if number:
            lo,hi=max(0,d-full),min(full,d)
            if hi<lo: raise ArithmeticError("positive degree mass has no profile pair")
            for a in range(lo,min(hi,d//2)+1):
                b=d-a
                if unordered%128==0: guard()
                left,right=join.two_orientations(clean[a],clean[b],diagonal=a==b)
                maximum=maximum.max(left)
                if a!=b: maximum=maximum.max(right)
                ordered+=1 if a==b else 2; unordered+=1
        term=number*maximum
        if not term.is_finite() or term.upper()<0:
            raise ArithmeticError("invalid band degree contribution")
        term=term.nonnegative_part(); degree_bounds.append(term); total+=term
        if progress is not None and (d==0 or (d+1)%64==0 or d+1==len(counts)):
            progress("directed_band_run_mass",degree_sum=d,ordered_envelope_pairs=ordered,
                unordered_powers=unordered,seconds=time.monotonic()-started)
    expected=(full+1)**2-(1 if max(weights)%2==0 else 0)
    if ordered!=expected:
        raise ArithmeticError("complete ordered envelope-pair count mismatch")
    guard(); count_sum=sum(counts); logtwo=arb(2).log()
    report={"scope":"directed source-scaled regional envelope with exact indexed band run masses",
        "actual_weights":weights,"points":POINTS,"actual_full_indices_not_holes":True,
        "all_original_weight_marginals_present":True,"marginal_union_size":len(clean),
        "ordered_envelope_pairs":ordered,"unordered_powers":unordered,
        "degree_count_entries":len(counts),"all_band_counts_sum_to_binomials":True,
        "band_count_sum_sha256":hashlib.sha256(str(count_sum).encode()).hexdigest(),
        "source_scale_already_in_each_regional_factor":True,"new_source_factor_applied":False,
        "new_alpha_applied":False,"output_weight_correction_applied":False,
        "caller_final_correction":"multiply by common z^(-alpha*L) exactly once",
        "regional_order":"e0(Rbar_a Rbar_b)^127 Rbar_a Vbar_b; reverse retains Vbar_a",
        "regional_factors":256,"packet_updates":32767,
        "envelope_input_intersections":envelope.input_intersections,
        "join_input_intersections":intersections,
        "before_output_correction_log2_interval":str(total.log()/logtwo) if total>0 else "zero",
        "seconds":time.monotonic()-started,"whole_code_certificate":False}
    return {"bound":total,"degree_bounds":tuple(degree_bounds),"report":report}


def exact_tests(join,guard):
    """Tiny independent Fraction chains/counts; no million-profile test."""
    from flint import arb
    checks=masks=0
    def direct(ra,rb,bits):
        a,_=ra; b,terminal=rb; vector=terminal
        for step in range((1<<(bits+1))-1):
            matrix=a if step%2==0 else b
            vector=(matrix[0]*vector[0]+matrix[1]*vector[1],
                    matrix[2]*vector[0]+matrix[3]*vector[1])
        return vector[0]
    def toy(r,f):
        return ((F(r+f+1,31),F(0) if f==0 else F(f,29),F(2,23),F(r+2,37)),
                (F(r+f,41),F(0) if r==1 else F(f+1,43)))
    def fraction_envelope(regions,scales):
        last=max(regions)//2
        return {f:tuple(tuple(max((scales[r]*rows[f][i][j] for r,rows in regions.items()
                if f in rows),default=F(0)) for j in range(width))
                for i,width in enumerate((4,2))) for f in range(last+1)}
    def fraction_join(regions,counts,bits):
        result=F(0); full=max(regions)
        for d,number in enumerate(counts):
            if number:
                result+=number*max(direct(regions[a],regions[d-a],bits)
                    for a in range(max(0,d-full),min(full,d)+1))
        return result
    for points in (6,8,10):
        weights=tuple(range(1,min(4,points//2)+1))
        regions={r:{f:toy(r,f) for f in range(r//2+1)} for r in weights}
        scales={r:F(r+1,r+2) for r in weights}
        envelope=fraction_envelope(regions,scales)
        counts=_band_counts(weights,points,join,guard)
        assert sum(counts)==sum(comb(points,r) for r in weights); checks+=1
        empirical=Counter(); supports=[]
        for mask in range(1<<points):
            r=mask.bit_count()
            if r not in regions: continue
            bits0=tuple((mask>>i)&1 for i in range(points))
            c=sum(bits0[i] and not bits0[(i-1)%points] for i in range(points))
            a=sum(bits0[2*i]*bits0[2*i+1] for i in range(points//2))
            b=sum(bits0[2*i+1]*bits0[(2*i+2)%points] for i in range(points//2))
            assert a+b==r-c; empirical[a+b]+=1; supports.append((r,a,b)); masks+=1
        assert tuple(empirical[d] for d in range(len(counts)))==counts; checks+=1
        for bits in (1,2,3):
            guard(); factors=1<<(bits+1); individual=F(0)
            exact=sum((scales[r]**factors*direct(regions[r][a],regions[r][b],bits)
                       for r,a,b in supports),F(0))
            for r in weights:
                singleton=fraction_envelope({r:regions[r]},{r:scales[r]})
                one_counts=_band_counts((r,),points,join,guard)
                ordinary=F(0)
                for c,number in enumerate(join.run_counts(points=points,weight=r)):
                    if c:
                        d=r-c; maximum=F(0)
                        for a in range(max(0,d-r//2),min(r//2,d)+1):
                            b=d-a; forward,reverse=join.two_orientations(regions[r][a],regions[r][b],bits)
                            assert forward==direct(regions[r][a],regions[r][b],bits)
                            assert reverse==direct(regions[r][b],regions[r][a],bits)
                            maximum=max(maximum,forward); checks+=2
                        ordinary+=number*maximum
                singleton_bound=fraction_join(singleton,one_counts,bits)
                assert singleton_bound==scales[r]**factors*ordinary; checks+=1
                individual+=singleton_bound
            bound=fraction_join(envelope,counts,bits)
            assert exact<=individual<=bound; checks+=2
    # Exercise the public directed path on a small actual8192-cycle band.
    small_weights=(2,3); rational={r:{f:toy(r,f) for f in range(r//2+1)} for r in small_weights}
    scales={2:F(2,3),3:F(3,4)}
    toarb=lambda x:arb(x.numerator)/arb(x.denominator)
    original={r:{f:tuple(tuple(toarb(x) for x in row) for row in region)
                    for f,region in rows.items()} for r,rows in rational.items()}
    directed=build_envelope(original,{r:toarb(x) for r,x in scales.items()},guard)
    actual=band_join(directed,small_weights,join,guard)
    expected=fraction_join(fraction_envelope(rational,scales),
                           _band_counts(small_weights,POINTS,join,guard),7)
    assert actual["bound"].overlaps(toarb(expected)); checks+=1
    for r in small_weights:
        one=build_envelope({r:original[r]},{r:arb(1)},guard)
        band=band_join(one,(r,),join,guard)["bound"]
        ordinary=join.directed_join(tuple(original[r][f] for f in range(r//2+1)),guard,weight=r)["bound"]
        assert band.overlaps(ordinary); checks+=1
    # Odd-neighbor endpoint1024 is valid in2048/2049, invalid in2047.
    one=((arb(1),)*4,(arb(1),)*2)
    sparse=build_envelope({2047:{1023:one},2048:{1024:one},2049:{1024:one}},
                         {r:arb(1) for r in (2047,2048,2049)},guard)
    assert set(sparse)=={1023,1024}; checks+=1
    try: band_join(sparse,(2047,2048,2049),join,guard)
    except ValueError: checks+=1
    else: raise AssertionError("partial pilot accepted as complete band")
    edge_counts=_band_counts((2047,2048,2049),POINTS,join,guard)
    assert len(edge_counts)==2049 and edge_counts[2048]==POINTS; checks+=1
    try: build_envelope({2047:{1024:one}},{2047:arb(1)},guard)
    except ValueError: checks+=1
    else: raise AssertionError("invalid odd-weight endpoint accepted")
    return {"exact_fraction_and_interface_checks":checks,"exhaustive_small_cycle_supports":masks,
        "public_Arb_small_band_control":True,"singleton_band_reproduction":True,
        "odd_neighbor_endpoint_and_count_control":True,"actual_large_matrix_joins":0,
        "tests_require_parent_execution":True}
