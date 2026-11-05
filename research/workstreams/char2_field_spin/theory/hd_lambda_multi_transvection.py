"""Cheap rejection screens for commuting quiet transvections of Hlambda.

Fix distinct binary weight-four vectors v with XOR of their labels zero.
Require even pairwise intersections. Then Hv=v and all v^T w vanish, so
Hnew=H+sum(v*v^T) is a symmetric involution. Its entries remain lambda or
lambda+1. Independence modulo RM(1,4) is checked separately: disjoint
planes need not increase the five-dimensional syndrome rank separately.

For a fixed complement O, write the full parity matrix Hnew[O,:]=B+lambda J.
Each column is nonzero. Its binary pattern determines its projective class:
equal patterns coincide; the all-zero and all-one patterns also coincide.
No other classes merge when lambda is outside {0,1}. This statement uses
the full parity matrix, INCLUDING its diagonal identity entries.

Let the class sizes be k_i and b=|O|. Every nonempty support inside one class
has rank one. Its exact forced-nonzero codeword count is
Z_w=((q-1)^w+(q-1)(-1)^w)/q. Thus the normalized D2-symmetrized pth moment
has the following nonnegative, exactly counted part:

  q^(p*b-16) * (1 + sum_i sum_(w=1..k_i)
      binom(k_i,w)*Z_w^p/(q-1)^((p-1)*w)).

The 1 is the zero codeword. This is a LOWER bound, not a replacement upper
bound. Large classes can reject a candidate; small classes do not certify
its full moments, scheduled word bound, or code distance. The optional
bounded screen never computes full support tables or regional providers.

Circuit counts reuse the original p[0..7], q[0..3] sums. They describe an
explicit scalar circuit, not a native implementation or timing result.
No prior construction receipt is transferred to any configuration here.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, islice
from math import comb, isfinite
from pathlib import Path
from time import perf_counter

import hd_lambda_outer as original


DIRECTORY = Path(__file__).resolve().parent
DEFAULT_LAMBDA = 0xbc
FULL = 65535
PINS = {
    'hd_lambda_outer.py': 'f0c063843e295a160a1186c6c084a07dd1c2a3745005ed68057fde59913c99e7',
    'hd_outer_pair_layers.py': '017c5ba20747ae8e25acf02393e9702a03b086038e1814ac2424cdfe49d66f60',
}
BASE_SYNDROMES = (FULL,) + tuple(sum(((i >> bit) & 1) << i for i in range(16)) for bit in range(4))


def _support(support):
    support = tuple(support)
    if (len(support) != 4 or len(set(support)) != 4
            or any(type(i) is not int or not 0 <= i < 16 for i in support)):
        raise ValueError('A support has four distinct labels in0..15')
    a, b, c, d = support
    if a ^ b ^ c ^ d:
        raise ValueError('A quiet binary support must have XOR of labels zero')
    return tuple(sorted(support))


def configuration(supports):
    supports = tuple(sorted(_support(s) for s in supports))
    if len(supports) > 3 or len(set(supports)) != len(supports):
        raise ValueError('Use at most three distinct quiet transvections')
    masks = tuple(sum(1 << i for i in support) for support in supports)
    if any((a & b).bit_count() & 1 for a, b in combinations(masks, 2)):
        raise ValueError('Transvection supports must have even pairwise intersections')
    return supports


def construction(supports):
    supports = configuration(supports)
    suffix = '-'.join(format(sum(1 << i for i in s), '04x') for s in supports) or 'none'
    return 'g16-hlambda-bc-commuting-quiet-transvections-'+suffix+'-research'


def binary_rank(vectors):
    basis = {}
    for value in vectors:
        if type(value) is not int or value < 0:
            raise ValueError('Nonnegative packed binary rows required')
        while value:
            pivot = value.bit_length()-1
            if pivot not in basis:
                basis[pivot] = value
                break
            value ^= basis[pivot]
    return len(basis)


def syndrome_extension(supports):
    supports = configuration(supports)
    rank = binary_rank(BASE_SYNDROMES + tuple(sum(1 << i for i in s) for s in supports))
    independent = rank == 5+len(supports)
    return dict(extended_binary_span_rank=rank, independent_modulo_RM1=independent,
                correction_rank_if_independent=rank if independent else None,
                quiet_dimension_if_independent=16-rank if independent else None)


def quiet_supports():
    return tuple(s for s in combinations(range(16), 4) if not (s[0] ^ s[1] ^ s[2] ^ s[3]))


def concrete_candidates():
    """A fixed cost-oriented list, not a completed search or endorsement."""
    return (
        ((0, 1, 2, 3), (4, 5, 6, 7)),  # dependent modulo RM(1,4): control
        ((0, 1, 2, 3), (0, 1, 4, 5)),
        ((0, 1, 2, 3), (0, 1, 4, 5), (0, 1, 8, 9)),
        ((0, 1, 2, 3), (0, 1, 4, 5), (0, 2, 4, 6)),
    )


def binary_rows(supports):
    """Packed full B in Hnew=B+lambda J, before the physical output schedule."""
    supports = configuration(supports)
    masks = tuple(sum(1 << i for i in s) for s in supports)
    rows = []
    for i in range(16):
        row = sum((int(i == j) ^ ((i & j).bit_count() & 1)) << j for j in range(16))
        for mask in masks:
            if mask & (1 << i):
                row ^= mask
        rows.append(row)
    return tuple(rows)


def matrix(supports, lambda_value=DEFAULT_LAMBDA, q=256):
    original._parameters(lambda_value, q)
    return tuple(tuple(lambda_value ^ ((row >> j) & 1) for j in range(16))
                 for row in binary_rows(supports))


def _circuit(x, supports, add, multiply, lambda_value):
    p = tuple(add(x[2*i], x[2*i+1]) for i in range(8))
    a0 = add(add(add(x[1], x[3]), add(x[5], x[7])), add(add(x[9], x[11]), add(x[13], x[15])))
    a1 = add(add(p[1], p[3]), add(p[5], p[7]))
    q = tuple(add(p[2*i], p[2*i+1]) for i in range(4))
    a2, a3 = add(q[1], q[3]), add(q[2], q[3])
    total = multiply(lambda_value, add(add(q[0], q[1]), a3))
    correction = [total]*16
    syndromes = (a0, a1, a2, a3)
    for i in range(1, 16):
        low = i & -i
        correction[i] = add(correction[i ^ low], syndromes[low.bit_length()-1])
    output = [add(x[i], correction[i]) for i in range(16)]
    # Hv=v, so these parities may be taken from the original input.
    for support in supports:
        if support == tuple(range(support[0], support[0]+4)) and support[0] % 4 == 0:
            parity = q[support[0]//4]
        elif support[0] % 2 == 0 and support[1] == support[0]+1 and support[2] % 2 == 0 and support[3] == support[2]+1:
            parity = add(p[support[0]//2], p[support[2]//2])
        else:
            parity = add(add(x[support[0]], x[support[1]]), add(x[support[2]], x[support[3]]))
        for i in support:
            output[i] = add(output[i], parity)
    return tuple(output)


def circuit(values, supports, lambda_value=DEFAULT_LAMBDA, q=256):
    original._parameters(lambda_value, q)
    values, supports = tuple(values), configuration(supports)
    if len(values) != 16 or any(type(x) is not int or not 0 <= x < q for x in values):
        raise ValueError('Sixteen field values required')
    return _circuit(values, supports, lambda a, b: a ^ b,
                    lambda a, b: original._multiply(a, b, q), lambda_value)


def primitive_counts(supports):
    supports = configuration(supports)
    counts = Counter()
    def add(a, b):
        counts['binary_xors'] += 1
        return a ^ b
    def multiply(a, b):
        counts['fixed_products'] += 1
        return original._multiply(a, b, 256)
    _circuit((0,)*16, supports, add, multiply, DEFAULT_LAMBDA)
    return dict(construction=construction(supports), supports=supports,
        binary_xors=counts['binary_xors'], fixed_products=counts['fixed_products'], variable_products=0,
        extra_xors_over_Hlambda=counts['binary_xors']-57, **syndrome_extension(supports),
        optimized_minimum_cost_claim=False, native_implementation_exists=False, native_performance_claim=False)


def projective_key(pattern, height):
    """Constant0/1 patterns merge; each nonconstant binary pattern is its own ray."""
    if (type(height) is not int or not 1 <= height <= 16 or type(pattern) is not int
            or not 0 <= pattern < 1 << height):
        raise ValueError('A packed column of a nonempty binary parity matrix is required')
    return 0 if pattern == (1 << height)-1 else pattern


def projective_classes(parity_rows, *, width=16):
    rows = tuple(parity_rows)
    if (type(width) is not int or not 1 <= width <= 16 or not 1 <= len(rows) <= 16
            or any(type(row) is not int or not 0 <= row < 1 << width for row in rows)):
        raise ValueError('A nonempty binary parity matrix of width1..16 is required')
    classes = {}
    for i in range(width):
        pattern = sum(((row >> i) & 1) << j for j, row in enumerate(rows))
        classes.setdefault(projective_key(pattern, len(rows)), []).append(i)
    return tuple((key, tuple(members)) for key, members in sorted(classes.items()))


def zero_sum_nonzero_count(weight, q=256):
    if type(weight) is not int or weight < 0 or q not in (4, 256):
        raise ValueError('Nonnegative weight and the explicit field order required')
    return ((q-1)**weight+(q-1)*(-1)**weight)//q


@lru_cache(maxsize=2048)
def rank_one_moment_part(class_sizes, codimension, order, *, q=256, width=16):
    """Exact zero+rank1 part; the caller must establish full parity-row rank."""
    if (type(width) is not int or not 1 <= width <= 16
            or type(codimension) is not int or not 1 <= codimension <= width
            or type(order) is not int or order not in (2, 3, 4) or q not in (4, 256)
            or not isinstance(class_sizes, tuple) or sum(class_sizes) != width
            or any(type(k) is not int or k <= 0 for k in class_sizes)):
        raise ValueError('Valid positive class sizes, full codimension, and order2..4 required')
    rank_one = F(0)
    for size in class_sizes:
        for w in range(1, size+1):
            count = zero_sum_nonzero_count(w, q)
            rank_one += F(comb(size, w)*count**order, (q-1)**((order-1)*w))
    scale = F(q)**(order*codimension-width)
    return scale*(1+rank_one)


def pair_profile(supports, *, lambda_value=DEFAULT_LAMBDA, q=256):
    supports = configuration(supports)
    h, ratios = matrix(supports, lambda_value, q), original.exceptional_ratios(lambda_value, q)
    rows, histogram, minimum, first_worst = [], Counter(), 16, None
    for a, b in combinations(range(16), 2):
        masks = tuple(sum(bool(original._multiply(row[a], ratio, q) ^ row[b]) << i
                          for i, row in enumerate(h)) for ratio in ratios)
        branches = ((0, 0, 1), (1, FULL, q-1), (1, FULL, q-1),
                    *((2, mask, q-1) for mask in masks), (2, FULL, (q-1)*(q-4)))
        for index, mask in enumerate(masks):
            weight = mask.bit_count()
            histogram[weight] += 1
            if weight < minimum:
                minimum, first_worst = weight, dict(source=(a, b), exceptional_ratio_index=index, output_mask=mask)
        rows.append(dict(source=(a, b), branches=branches))
    return dict(construction=construction(supports), pair_count=120, singleton_output_weight=16,
        field_order=q, lambda_value=lambda_value, toy_field_only=q != 256,
        two_active_minimum_output_weight=minimum, first_worst=first_worst,
        exceptional_ratio_histogram=dict(sorted(histogram.items())), pairs=rows,
        old_lowarity_receipts_transferred=False)


def checked_dependencies():
    for name, digest in PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Reviewed dependency changed: '+name)
    return dict(PINS, **{Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()})


def complement_screen(supports, *, complements=None, max_seconds=10, orders=(3, 4)):
    """Bounded LOWER-bound screen; default scope is all4368 size-five complements.

    Passing an explicit subset is a diagnostic only. Time exhaustion preserves
    partial witnesses and never claims that every complement was checked.
    """
    supports, orders = configuration(supports), tuple(orders)
    if (not orders or len(set(orders)) != len(orders)
            or any(type(p) is not int or p not in (2, 3, 4) for p in orders)):
        raise ValueError('Distinct orders from2,3,4 required')
    if (type(max_seconds) not in (int, float) or not isfinite(max_seconds)
            or not 0 < max_seconds <= 30):
        raise ValueError('Use a positive screen budget at most30 seconds')
    full_scope = complements is None
    plan = tuple(combinations(range(16), 5)) if full_scope else tuple(tuple(o) for o in islice(complements, 4369))
    if (len(plan) > 4368 or len(set(plan)) != len(plan) or any(len(o) != 5 or tuple(sorted(set(o))) != o
            or any(type(i) is not int or not 0 <= i < 16 for i in o) for o in plan)):
        raise ValueError('Distinct sorted size-five complements required')
    pins, started = checked_dependencies(), perf_counter()
    rows = binary_rows(supports)
    pairs = pair_profile(supports)
    maximum, max_witness, histogram, examined = 0, None, Counter(), 0
    moment_max = {p: F(0) for p in orders}
    moment_witness = {p: None for p in orders}
    stop = None
    for complement in plan:
        if perf_counter()-started >= max_seconds:
            stop = 'Time budget; partial screen only'
            break
        classes = projective_classes(tuple(rows[o] for o in complement))
        sizes = tuple(sorted(len(indices) for _, indices in classes))
        largest = max(sizes)
        histogram[largest] += 1
        if largest > maximum:
            maximum = largest
            max_witness = dict(complement=complement, classes=classes)
        for p in orders:
            value = rank_one_moment_part(sizes, 5, p)
            if value > moment_max[p]:
                moment_max[p], moment_witness[p] = value, complement
        examined += 1
    if checked_dependencies() != pins:
        raise ValueError('Sources changed during the screen')
    union = set().union(*(set(s) for s in supports))
    return dict(construction=construction(supports), supports=supports, **syndrome_extension(supports),
        checked_complements=examined, requested_complements=len(plan), complete_requested_scope=examined == len(plan),
        complete_all4368=full_scope and examined == 4368, stop_reason=stop,
        maximum_projective_class=maximum, maximum_class_witness=max_witness,
        largest_class_histogram=dict(sorted(histogram.items())),
        maximum_zero_plus_rank1_moment_parts={str(p): str(v) for p, v in moment_max.items()},
        moment_part_witnesses={str(p): v for p, v in moment_witness.items()},
        bounds_are_lower_not_upper=True, class_minimization_is_not_certificate=True,
        unchanged_disjoint_complements=comb(16-len(union), 5) if len(union) <= 11 else 0,
        circuit_counts=primitive_counts(supports), reviewed_source_sha256=pins,
        exact_two_source_profile={k: v for k, v in pairs.items() if k != 'pairs'},
        regional_provider_runs=0, full_support_tables=0, old_receipts_transferred=False,
        whole_code_certificate=False, elapsed_seconds=perf_counter()-started)


def screen_candidates(candidates, *, max_candidates=16, max_seconds=10):
    """At most16 supplied configurations, one global budget; no implicit search."""
    if type(max_candidates) is not int or not 1 <= max_candidates <= 16:
        raise ValueError('Candidate cap must be1..16')
    if (type(max_seconds) not in (int, float) or not isfinite(max_seconds)
            or not 0 < max_seconds <= 30):
        raise ValueError('Use a positive global budget at most30 seconds')
    candidates = tuple(configuration(s) for s in islice(candidates, max_candidates+1))
    if len(candidates) > max_candidates or len(set(candidates)) != len(candidates):
        raise ValueError('Distinct supplied candidates must fit the explicit cap')
    pins, started, results = checked_dependencies(), perf_counter(), []
    for supports in candidates:
        remaining = max_seconds-(perf_counter()-started)
        if remaining <= 0:
            break
        row = complement_screen(supports, max_seconds=remaining)
        results.append(row)
        if not row['complete_all4368']:
            break
    if checked_dependencies() != pins:
        raise ValueError('Sources changed during the candidate screen')
    return dict(candidates=results, requested_candidates=len(candidates),
        complete=len(results) == len(candidates) and all(r['complete_all4368'] for r in results),
        lower_bound_rejection_screen_only=True, no_success_certificate=True,
        reviewed_source_sha256=pins, elapsed_seconds=perf_counter()-started)
