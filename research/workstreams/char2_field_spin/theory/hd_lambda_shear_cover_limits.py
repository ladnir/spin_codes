"""Optimistic lower comparator for the disjoint-cover proof expression.

Fix a length-n code with b parity checks and a moment order p=2,3,4.
A witness consists of p disjoint independent parity-column sets, of
sizes s_i. The support-rank relaxation assigns the exact expression
 U(s)=q^(p*b-n+u) sum_k prod_i C(s_i,k_i)
                         (q-1)^(sum_i k_i-p*max_i k_i),
where u=n-sum_i s_i. This upper-bounds the actual likelihood moment
ONLY when the sizes come from an actual verified independence witness.

Here we instead minimize U(s) over abstract sorted size profiles. We
ignore whether the parity matrix admits those profiles. This minimum
is at most the expression produced by every actual cover witness;
it need not bound the actual likelihood moment in either direction.
It suffices to use sum_i s_i=min(n,p*b). To justify this restriction,
move an uncovered coordinate into any set below capacity b. For each
support S the maximum intersection size cannot decrease. Thus the
positive support sum cannot increase. This argument concerns abstract
sets, not a claim that a matrix allows the corresponding augmentation.

The size sum is evaluated by maximum level m. Put
 A_i(m)=sum_(k<=m) C(s_i,k)(q-1)^k and A_i(-1)=0. The total integer
weight of tuples whose maximum is m is prod_i A_i(m)-prod_i A_i(m-1).
Each difference is nonnegative. Dividing it by (q-1)^(p*m) gives an
exact rational sum without Cartesian enumeration or precision loss.

For each arity a=4..10, set n=32 and b=32-2a. Subtract the exact zero
atom q^(p*b-32) from the abstract minimum. The positive eight-pair tail
of degree at least five is monotone in the resulting nonzero atom.
Every actual cover-based tail for a fixed original A is at least this
optimistic tail. Consequently the raw group comparator is
 C(16,a)*256^(16a)/255^(8a)*tail^(1/p).
Original A is fixed through all eight pairs. The binomial factor occurs
once because this comparator is the same for every labelled A.

The frozen full-profile uniform reference turns this into a LOWER
comparator for the homogeneous r18 cover-method expression. It tests
only p2/3/4, mark691/700, alpha1/5, and the five-pair activity floor.
A lower endpoint above -40 rules out that prescribed cover relaxation,
not the construction or a proof with additional rank/input constraints.
A passing comparator does not certify a failure upper bound. This
module performs no parity search, moment inventory, or provider run.
"""
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from math import comb, prod
from pathlib import Path
from time import perf_counter

from flint import arb, ctx

import hd_lambda_shear_renyi_ideal as reference

DIRECTORY = Path(__file__).resolve().parent
ARITIES, ORDERS = tuple(range(4,11)), (2,3,4)
PINS = {'hd_lambda_shear_renyi_ideal.py':
    '7c2d5e9b720cd6cae09212de0fbf11ea0c19352f023706d08bd2919b5477c7a5'}


def dependencies():
    pins = reference.dependencies()
    for name,digest in PINS.items():
        if name in pins and pins[name] != digest:
            raise ValueError('Conflicting source pin: '+name)
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Frozen reference source changed: '+name)
        pins[name] = digest
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    return pins


def _parameters(width,checks,order,q):
    if (type(width) is not int or not 1 <= width <= 32
            or type(checks) is not int or not 0 <= checks <= width
            or type(order) is not int or order not in ORDERS
            or type(q) is not int or q not in (4,256)):
        raise ValueError('Width1..32, checks0..width, p2/3/4 and q4/256 required')


@lru_cache(maxsize=4096)
def size_sum(sizes,q=256):
    """Exact nonnegative max-level sum; no support or Cartesian loop."""
    if (type(sizes) is not tuple or len(sizes) not in ORDERS
            or any(type(s) is not int or not 0 <= s <= 32 for s in sizes)
            or type(q) is not int or q not in (4,256)):
        raise ValueError('A tuple of p=2/3/4 sizes0..32 and q4/256 required')
    p, base = len(sizes), q-1
    partial = [0]*p
    previous, power, result = 0, 1, F(0)
    for maximum in range(max(sizes)+1):
        for i,size in enumerate(sizes):
            if maximum <= size:
                partial[i] += comb(size,maximum)*power
        current = prod(partial)
        difference = current-previous
        if difference < 0:
            raise ArithmeticError('A maximum-level integer weight became negative')
        result += F(difference,power**p)
        previous = current
        power *= base
    return result


def abstract_cover_expression(width,checks,sizes,*,q=256):
    sizes = tuple(sizes)
    _parameters(width,checks,len(sizes),q)
    if any(type(s) is not int or not 0 <= s <= checks for s in sizes) or sum(sizes)>width:
        raise ValueError('Abstract sizes must fit the parity and coordinate capacities')
    zero = F(q)**(len(sizes)*checks-width)
    value = zero*F(q)**(width-sum(sizes))*size_sum(tuple(sorted(sizes)),q)
    if value < zero:
        raise ArithmeticError('The abstract expression excluded its zero term')
    return value


def size_profiles(width,checks,order):
    """Every sorted maximal-capacity profile, including empty sets."""
    _parameters(width,checks,order,256)
    target = min(width,order*checks)
    def visit(prefix,minimum,remaining,slots):
        if not slots:
            if remaining == 0:
                yield prefix
            return
        low = max(minimum,remaining-checks*(slots-1))
        high = min(checks,remaining//slots)
        for value in range(low,high+1):
            yield from visit(prefix+(value,),value,remaining-value,slots-1)
    yield from visit((),0,target,order)


def optimistic_cover(width,checks,order,*,q=256,check_time=None):
    """Minimum abstract expression, not a bound on an actual moment."""
    _parameters(width,checks,order,q)
    best, best_sizes, count, ties = None,None,0,0
    for sizes in size_profiles(width,checks,order):
        if check_time is not None:
            check_time()
        value = abstract_cover_expression(width,checks,sizes,q=q)
        count += 1
        if best is None or value < best:
            best,best_sizes,ties = value,sizes,1
        elif value == best:
            ties += 1
    if best is None:
        raise ArithmeticError('No maximal-capacity abstract profile was enumerated')
    zero = F(q)**(order*checks-width)
    return dict(width=width,checks=checks,moment_order=order,field_order=q,
        abstract_profiles=count,covered_coordinates=min(width,order*checks),
        minimizing_sizes=list(best_sizes),tied_minima=ties,
        zero_pair_moment=str(zero),minimum_abstract_expression=str(best),
        optimistic_nonzero_atom=str(best-zero),actual_witness_claimed=False)


def raw_group_comparator(arity,order,*,check_time=None):
    """Arb context required; complete labelled homogeneous-arity sum."""
    if type(arity) is not int or arity not in ARITIES:
        raise ValueError('The prescribed comparison covers arities4..10 only')
    diagnostic = optimistic_cover(32,32-2*arity,order,check_time=check_time)
    zero,nonzero = F(diagnostic['zero_pair_moment']),F(diagnostic['optimistic_nonzero_atom'])
    tail = reference.positive_tail(zero,nonzero)
    beta = F(256)**(16*arity)/F(255)**(8*arity)
    raw = (arb(comb(16,arity))*reference.A.uniform.rational(beta)
           *reference.A.uniform.rational(tail)**reference.A.uniform.rational(F(1,order)))
    if not raw > 0 or not raw.is_finite():
        raise ArithmeticError('Invalid positive group comparator')
    return raw,dict(diagnostic,arity=arity,labelled_rowsets=comb(16,arity),beta=str(beta),
        optimistic_activity_tail=str(tail),raw_group_log2_enclosure=str(reference.A.log2(raw)))


def controls_from_reference(uniform_log,*,check_time=None):
    if not uniform_log.is_finite():
        raise ValueError('A finite directed uniform-reference enclosure is required')
    controls,selected = [],[]
    for arity in ARITIES:
        candidates = []
        for order in ORDERS:
            if check_time is not None:
                check_time()
            raw,diagnostic = raw_group_comparator(arity,order,check_time=check_time)
            value = reference.method_expression(raw,order,uniform_log)
            row = dict(diagnostic,**reference._interval(value),
                outer_theta=str(F(1,5)*F(order,order-1)),
                interpretation='lower-on-prescribed-cover-method-expression-only')
            controls.append(row)
            candidates.append(row)
        low = min(candidates,key=lambda row:F(row['method_expression_log2_lower']))
        high = min(candidates,key=lambda row:F(row['method_expression_log2_upper']))
        lower = F(low['method_expression_log2_lower'])
        selected.append(dict(arity=arity,selected_lower_order=low['moment_order'],
            selected_upper_order=high['moment_order'],method_expression_log2_lower=str(lower),
            method_expression_log2_upper=high['method_expression_log2_upper'],
            this_prescribed_cover_method_cannot_reach_40_bits=bool(lower > -40)))
    return dict(controls=controls,selected_homogeneous_arities=selected)


def receipt(*,precision=256,max_seconds=120,progress=None):
    """Root-owned comparator; no actual matrix, inventory, or provider."""
    if type(precision) is not int or precision not in (256,384):
        raise ValueError('Use the matching cached256/384-bit reference')
    if type(max_seconds) not in (int,float) or not 0 < max_seconds <= 600:
        raise ValueError('A finite bounded budget in(0,600] seconds is required')
    pins,started = dependencies(),perf_counter()
    def check():
        if perf_counter()-started > max_seconds:
            raise TimeoutError('Optimistic cover comparator exceeded its budget')
    with ctx.workprec(precision):
        uniform_log,bank = reference.checked_reference(precision)
        if progress:
            progress(dict(stage='abstract-cover-profiles',arities=list(ARITIES),orders=list(ORDERS)))
        data = controls_from_reference(uniform_log,check_time=check)
        reference_text = str(uniform_log)
    check()
    if dependencies() != pins:
        raise ValueError('Frozen sources changed during the comparator')
    if sha256((DIRECTORY/bank['file']).read_bytes()).hexdigest() != bank['sha256']:
        raise ValueError('Cached uniform reference changed during the comparator')
    return dict(data,complete=True,construction='scheduled-lambda-bc-staggered-column-shear-final-D2',
        active_groups=18,homogeneous_arities=list(ARITIES),orders=list(ORDERS),
        source_profiles=99,pair_width=32,pair_count=8,minimum_active_pairs=5,
        all_labelled_original_rowsets=True,original_A_fixed_across_pairs=True,
        exact_original_row_floors_retained=False,restricted_density_not_renormalized=True,
        source_count_outside_theta=True,rowset_count_inside_theta=True,
        projective_normalizations=1,additional_beta_factor=False,
        mark='691/700',uniform_alpha='1/5',forbidden_weight=734003,
        comparator_arithmetic_enclosed=True,precision_bits=precision,
        cached_uniform_log2_enclosure=reference_text,cached_reference=bank,
        abstract_sizes_ignore_actual_parity=True,actual_independence_witnesses=False,
        comparator_is_lower_on_prescribed_method=True,
        comparator_is_failure_probability_lower_bound=False,
        failure_probability_upper_bound=False,actual_pair_moment_bound=False,
        construction_rejected=False,mixed_arities_covered=False,other_group_counts_covered=False,
        whole_code_certificate=False,fresh_regional_banks=0,parity_searches=0,
        exact_pair_inventory_runs=0,native_operations=0,reviewed_source_sha256=pins,
        elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--precision',type=int,choices=(256,384),default=256)
    parser.add_argument('--max-seconds',type=float,default=120)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite evidence')
    result = receipt(precision=args.precision,max_seconds=args.max_seconds,
                     progress=lambda row:print(json.dumps(row),flush=True))
    with args.output.open('x',encoding='utf-8') as handle:
        json.dump(result,handle,indent=2,allow_nan=False)
        handle.write('\n')
    print(json.dumps(dict(selected_homogeneous_arities=result['selected_homogeneous_arities'],
        interpretation='cover-method-lower-comparator-only',elapsed_seconds=result['elapsed_seconds'])))
