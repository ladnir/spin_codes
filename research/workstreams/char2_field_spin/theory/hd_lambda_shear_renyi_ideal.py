"""Ideal paired-block Renyi comparator, NOT a failure upper bound.

For one original row set A of size a, an invertible two-column map has
2a input dimensions and32 output coordinates. Include the final
independent nonzero output diagonal. Relative to uniform U32, the
likelihood's zero p-moment is exactly z=q**(p*(32-2a)-32).
Its nonzero probability mass is1-q**(-2a). Jensen gives the universal
LOWER bound n=(1-q**(-2a))**p/(1-q**(-32))**(p-1) for its nonzero moment.

Every active RS[16,8] row has at least9 nonzero original columns.
Therefore at least5 of the8 original column pairs are nonzero. Keep
this pair-activity restriction, not the stronger individual-row floor.
For independent pairs with fixed original A, replacing each nonzero
moment by n lowers the positive tail sum_{k>=5} C(8,k) z**(8-k)n**k.
The resulting raw group comparator is
  B_a=C(16,a)*q**(16a)/255**(8a)*tail**(1/p).
Pairing does NOT double the existing dimensional/RS normalization.

This file evaluates the ideal expression under the SAME fixed uniform
regional majorant used by the frozen full99-source-profile receipt.
Uniform U32 equals two independent U16 columns, so the uniform inner
reference and its physical256-region chronology are unchanged. For
alpha=1/5, theta=alpha*p/(p-1), the homogeneous expression is
  theta*(18*log2(B_a)-log2(255))-alpha*734003*log2(691/700)+uniform_log.
The cached uniform_log already sums indexed source counts, outside theta.

The output encloses a lower comparator for THIS moment/activity/reference
method. It is neither a failure-probability bound nor a lower bound on
the construction's distance. A failing comparator can reject this fixed
proof relaxation; a passing one only leaves room for an actual upper
moment bound. No exact pair inventory, provider, or native job is run.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb
from pathlib import Path
from time import perf_counter

from flint import arb, ctx

import hd_lambda_renyi_limits as limits

DIRECTORY = Path(__file__).resolve().parent
A = limits.A
ARITIES, ORDERS, R = (9,10,11,12), (2,3,4), 18
ALPHA, MARK = F(1,5), F(691,700)
PINS = {'hd_lambda_renyi_limits.py':
        'ae285af26a43c32e8bbc6119f4631cefea8f75d1fb7c2fd394408abd57106a2f'}
BANKS = {
    256: ('hd_lambda_renyi_limits_256.json', 'f3d433f3e55d369741f93e53370d57e1423ab5b12ff033328c168e817ddf71be'),
    384: ('hd_lambda_renyi_limits_384.json', '0a425a1acb6466880226c08abcd131ab910b4f76f3f9ed88a939a2b99644f31e'),
}


def dependencies():
    for name, digest in PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Frozen reference helper changed: '+name)
    return dict(limits.dependencies(), **PINS,
        **{Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()})


def checked_reference(precision):
    """Read only the pinned full99-profile UNIFORM reference enclosure."""
    if type(precision) is not int or precision not in BANKS:
        raise ValueError('Use the matching cached256/384-bit reference')
    name, digest = BANKS[precision]
    raw = (DIRECTORY/name).read_bytes()
    if sha256(raw).hexdigest() != digest:
        raise ValueError('Frozen full99-profile reference receipt changed')
    data = json.loads(raw)
    expected = dict(active_groups=R, source_profiles=99, field_order=256,
        source_dimension=1048576, output_length=2097152, forbidden_weight=734003,
        mark=str(MARK), uniform_alpha=str(ALPHA), precision_bits=precision,
        source_count_outside_theta=True, projective_normalizations=1,
        arithmetic_enclosures_directed=True, reviewed_source_sha256=limits.dependencies())
    for key, value in expected.items():
        if type(data.get(key)) is not type(value) or data[key] != value:
            raise ValueError('Wrong cached uniform reference scope: '+key)
    text = data.get('uniform_source_sum_log2_enclosure')
    if type(text) is not str:
        raise ValueError('An explicit Arb uniform-log enclosure is required')
    reference = arb(text)
    if not reference.is_finite():
        raise ValueError('The cached uniform reference is not finite')
    return reference, dict(file=name, sha256=digest, precision_bits=precision,
        imported_field='uniform_source_sum_log2_enclosure', profile_count=99)


def ideal_atoms(dimension, width, order, *, q=256):
    """Exact zero atom and Jensen nonzero LOWER atom, for tiny tests too."""
    if (type(width) is not int or not 1 <= width <= 32
            or type(dimension) is not int or not 1 <= dimension <= width
            or type(order) is not int or order not in ORDERS
            or type(q) is not int or q not in (4,256)):
        raise ValueError('Dimension1..width<=32, q4/256, and order2/3/4 required')
    zero = F(q)**(order*(width-dimension)-width)
    nonzero = (1-F(q)**(-dimension))**order/(1-F(q)**(-width))**(order-1)
    return zero, nonzero


def positive_tail(zero, nonzero, *, blocks=8, minimum=5):
    """Exact positive binomial tail; no total-minus-low cancellation."""
    if (type(blocks) is not int or not 1 <= blocks <= 8
            or type(minimum) is not int or not 0 <= minimum <= blocks
            or not isinstance(zero,(int,F)) or not isinstance(nonzero,(int,F))
            or zero < 0 or nonzero < 0):
        raise ValueError('Nonnegative rational atoms and a bounded feasible floor required')
    return sum((F(comb(blocks,k))*zero**(blocks-k)*nonzero**k
                for k in range(minimum,blocks+1)), F(0))


def raw_group_comparator(arity, order):
    """Call inside an Arb context. Includes rowsets and ONE old beta_a."""
    if type(arity) is not int or arity not in ARITIES:
        raise ValueError('The prescribed comparison covers arities9..12 only')
    zero, nonzero = ideal_atoms(2*arity,32,order)
    tail = positive_tail(zero,nonzero)
    beta = F(256)**(16*arity)/F(255)**(8*arity)
    raw = (arb(comb(16,arity))*A.uniform.rational(beta)
           * A.uniform.rational(tail)**A.uniform.rational(F(1,order)))
    if not raw > 0 or not raw.is_finite():
        raise ArithmeticError('Invalid ideal group comparator')
    return raw, dict(arity=arity,moment_order=order,zero_pair_moment=str(zero),
        ideal_nonzero_pair_moment=str(nonzero),ideal_activity_tail=str(tail),
        labelled_rowsets=comb(16,arity),beta=str(beta),raw_group_log2_enclosure=str(A.log2(raw)))


def method_expression(raw_weight, order, uniform_log):
    """The reference already includes the source union; no added beta."""
    if type(order) is not int or order not in ORDERS or not raw_weight > 0:
        raise ValueError('A positive raw weight and order2/3/4 are required')
    theta = ALPHA*F(order,order-1)
    return (A.uniform.rational(theta)*(R*A.log2(raw_weight)-A.log2(arb(255)))
        -A.uniform.rational(ALPHA)*734003*A.log2(A.uniform.rational(MARK))+uniform_log)


def _interval(value):
    lower = -A.uniform.rational_ceiling(-value)
    upper = A.uniform.rational_ceiling(value)
    return dict(method_expression_log2_lower=str(lower),
                method_expression_log2_upper=str(upper))


def controls_from_reference(uniform_log):
    """Twelve prescribed controls; select only among completed p2/3/4."""
    rows, selected = [], []
    for arity in ARITIES:
        candidates = []
        for order in ORDERS:
            raw, diagnostic = raw_group_comparator(arity,order)
            expression = method_expression(raw,order,uniform_log)
            row = dict(diagnostic, **_interval(expression),
                outer_theta=str(ALPHA*F(order,order-1)),
                interpretation='ideal-lower-on-this-fixed-method-not-failure-bound')
            candidates.append(row)
            rows.append(row)
        low = min(candidates,key=lambda x:F(x['method_expression_log2_lower']))
        high = min(candidates,key=lambda x:F(x['method_expression_log2_upper']))
        lower = F(low['method_expression_log2_lower'])
        selected.append(dict(arity=arity, selected_lower_order=low['moment_order'],
            selected_upper_order=high['moment_order'],
            method_expression_log2_lower=str(lower),
            method_expression_log2_upper=high['method_expression_log2_upper'],
            this_fixed_method_cannot_reach_40_bits=bool(lower > -40)))
    return dict(controls=rows,selected_homogeneous_arities=selected)


def receipt(precision=256):
    """Root owns real runs. No provider invocation or moment inventory."""
    if type(precision) is not int or precision not in BANKS:
        raise ValueError('Use cached256/384-bit reference precision')
    pins, started = dependencies(), perf_counter()
    with ctx.workprec(precision):
        reference, bank = checked_reference(precision)
        data = controls_from_reference(reference)
        reference_text = str(reference)
    # Recheck both the complete dependency graph and the cache itself.
    if dependencies() != pins:
        raise ValueError('Frozen sources changed during the comparator')
    if sha256((DIRECTORY/bank['file']).read_bytes()).hexdigest() != bank['sha256']:
        raise ValueError('Cached uniform reference changed during the comparator')
    return dict(data, construction='scheduled-lambda-bc-staggered-column-shear-final-D2',
        reference='uniform-GF256^32-per-pair-equals-two-independent-U16-columns',
        active_groups=R, homogeneous_arities=list(ARITIES), orders=list(ORDERS),
        mixed_arities_covered=False, other_group_counts_covered=False,
        all_labelled_original_rowsets=True, original_A_fixed_across_pairs=True,
        minimum_original_nonzero_columns=9, minimum_active_pairs=5, pair_count=8,
        exact_original_row_floors_retained=False, restricted_density_not_renormalized=True,
        full_independent_D1_comparison=True, permitted_gauges_per_pair=1,
        final_independent_D2_after_shear=True, source_profiles=99,
        source_count_outside_theta=True, rowset_sum_inside_theta=True,
        projective_normalizations=1, additional_beta_factor=False,
        mark=str(MARK),uniform_alpha=str(ALPHA),forbidden_weight=734003,
        cached_uniform_log2_enclosure=reference_text,cached_reference=bank,
        ideal_is_failure_bound=False, actual_pair_moment_upper_bound=False,
        construction_rejected=False, whole_code_certificate=False,
        comparator_arithmetic_enclosed=True,precision_bits=precision,
        fresh_regional_banks=0,exact_pair_inventory_runs=0,
        native_operations=0,reviewed_source_sha256=pins,
        elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--precision',type=int,choices=(256,384),default=256)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite evidence')
    result = receipt(args.precision)
    with args.output.open('x',encoding='utf-8') as handle:
        json.dump(result,handle,indent=2,allow_nan=False)
        handle.write('\n')
    print(json.dumps(dict(selected_homogeneous_arities=result['selected_homogeneous_arities'],
        interpretation='method-lower-comparator-only',elapsed_seconds=result['elapsed_seconds'])))
