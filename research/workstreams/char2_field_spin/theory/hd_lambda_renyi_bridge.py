"""Directed, homogeneous arity12 comparison through a third density moment.

Fix an indexed set G of r=18 active outer superblocks. Every block has
exactly twelve active RS rows, but their labelled row sets may differ.
For a fixed such row set A and column j, normalize the raw allowed-input
column count to a probability P_(A,j), including zero inputs. Independent
D2 scalars make its likelihood relative to uniform U on GF256^16 equal
L(y)=256^4*n_C(supp(y))/255^|supp(y)|. Let I3=E_U[L^3].

If every arity12 column has I3<=M, the product likelihood for a fixed
row-set tuple has L3 norm at most M^(16*r/3). The original RS counting
majorant is 256^(192*r)/255^(96*r) times that product probability.
Sum all 1820^r row-set tuples before taking any fractional power.

For fixed route/inner setup omega, Holder with exponents3 and3/2 bounds
the bad-event probability under the product probability by its L3 norm
times the uniform bad-event probability to power2/3. Markov's inequality
and min(1,x)<=x^(3/10) leave uniform exponent1/5. Thus the existing
uniform regional provider at Z=691/700 and alpha=1/5 applies unchanged.
The single projective factor255 is raised to3/10, not1/5. Indexed G
counts are unioned outside that power. No saved numerical bank is used.

This closes only the homogeneous arity12/r18 class if the bound passes.
It does not cover mixed arities, other r, or the whole code. A complete
exact moment inventory is required; a representative does not suffice.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from math import comb, isfinite
from pathlib import Path
from time import perf_counter

from flint import arb, ctx

import hd_lambda_high_arity_certificate as binding


DIRECTORY = Path(__file__).resolve().parent
A = binding.arithmetic
R, ARITY, MOMENT_ORDER = 18, 12, 3
MARK, ALPHA, THETA = F(691, 700), F(1, 5), F(3, 10)
ADAPTER_SHA256 = '156ca16d6f18370eb2484e6a1bf74ef83cbede670b63b2ba03db27f465d2e036'
MOMENT_PINS = {
    'hd_lambda_renyi_support.py': 'c4662bb016f248faa7ffcc1972e01349bc66306c40762e12df9245fc04a3c03e',
    'hd_binary_outer_relations.py': 'ec1c13cd4a685760983b33ebd8e722efc3735f6335ed55bfe444aea01b1a20f2',
}
MOMENT_SCOPE = dict(field_order=256, field_polynomial='0x11b', lambda_aes_byte=0xbc,
    arity=12, moment_order=3, exact_rational_arithmetic=True,
    code='ker(Hlambda[complement,:])', full_parity_includes_identity=True,
    normalized_likelihood='256^4*n_C(S)/255^|S|',
    output_diagonal='independent uniform nonzero field multipliers',
    coordinate_permutations_leave_moment_unchanged=True)


def dependencies():
    adapter = DIRECTORY/'hd_lambda_high_arity_certificate.py'
    if sha256(adapter.read_bytes()).hexdigest() != ADAPTER_SHA256:
        raise ValueError('Reviewed construction adapter changed')
    for name, digest in MOMENT_PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Reviewed moment helper changed: '+name)
    return dict(binding.checked_dependencies(), **{
        adapter.name: ADAPTER_SHA256,
        Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()}, **MOMENT_PINS)


def checked_moment_inventory(data):
    """Require all labelled four-coordinate complements, exactly once."""
    expected = set(combinations(range(16), 4))
    seen, maximum = set(), F(1)
    for row in data.get('rows', []):
        complement = row.get('complement')
        if (not isinstance(complement, list) or len(complement) != 4
                or any(type(i) is not int for i in complement)):
            raise ValueError('Four ordered integer complement labels required')
        key = tuple(complement)
        if key not in expected or key in seen:
            raise ValueError('Unexpected or duplicate complement')
        if type(row.get('moment')) is not str:
            raise ValueError('An exact rational moment is required')
        moment = F(row['moment'])
        if moment < 1:
            raise ValueError('Normalized density moment cannot be below one')
        maximum = max(maximum, moment)
        seen.add(key)
    if seen != expected:
        raise ValueError('Incomplete1820-complement inventory')
    for key, value in dict(MOMENT_SCOPE, complete=True).items():
        if type(data.get(key)) is not type(value) or data[key] != value:
            raise ValueError('Wrong exact moment scope: '+key)
    if (data.get('reviewed_source_sha256') != MOMENT_PINS
            or data.get('source_sha256') != MOMENT_PINS['hd_lambda_renyi_support.py']
            or data.get('source_pins_unchanged') is not True):
        raise ValueError('Exact inventory is not bound to the reviewed helper')
    return maximum


def moment_inventory(*, progress=None, max_seconds=360):
    """Serial exact enumeration of all mapped complements; no regional runs."""
    import hd_lambda_renyi_support as local
    if (type(max_seconds) not in (int, float) or not isfinite(max_seconds)
            or not 0 < max_seconds <= 600):
        raise ValueError('Use a bounded exact inventory run, at most600 seconds')
    helper = Path(local.__file__)
    dependencies()
    pins = dict(MOMENT_PINS)
    started, rows, maximum, stop = perf_counter(), [], F(1), None
    for complement in combinations(range(16), 4):
        if perf_counter()-started >= max_seconds:
            stop = 'Time budget: partial inventory is not a certificate'
            break
        data = local.renyi_moment(complement, q=256, order=3)
        moment = data['moment']
        rows.append(dict(complement=list(complement), moment=str(moment)))
        maximum = max(maximum, moment)
        if progress and (len(rows) == 1 or len(rows) % 100 == 0):
            progress(dict(completed=len(rows), total=1820,
                maximum_moment_float=float(maximum), elapsed_seconds=perf_counter()-started))
    after = {name: sha256((DIRECTORY/name).read_bytes()).hexdigest() for name in pins}
    if pins != after:
        raise ValueError('Exact moment sources changed during enumeration')
    return dict(MOMENT_SCOPE, complete=len(rows) == 1820,
        rows=rows, maximum_moment=str(maximum),
        stop_reason=stop, reviewed_source_sha256=pins, source_sha256=pins[helper.name],
        source_pins_unchanged=True, elapsed_seconds=perf_counter()-started)


def prefix(maximum_moment):
    """Call inside an Arb precision context; all powers use exact rationals."""
    maximum_moment = F(maximum_moment)
    if maximum_moment < 1:
        raise ValueError('A normalized third moment upper bound at least one is required')
    log_block = (16*ARITY*8-8*ARITY*A.log2(arb(255))
                 + A.log2(arb(comb(16, ARITY)))
                 + A.uniform.rational(F(16, 3))*A.log2(A.uniform.rational(maximum_moment)))
    return (A.uniform.rational(THETA)*(R*log_block-A.log2(arb(255)))
            - A.uniform.rational(ALPHA)*binding.SCOPE['forbidden_weight']
              * A.log2(A.uniform.rational(MARK)))


def receipt(inventory, precision=256):
    if type(precision) is not int or precision not in (256, 384):
        raise ValueError('Use a fresh256/384-bit replay')
    pins = dependencies()
    maximum = checked_moment_inventory(inventory)
    started = perf_counter()
    with ctx.workprec(precision):
        regions = A.uniform.regional_directed(MARK, ALPHA, R, 60)
        terms = A.source_profile_terms(R, regions)
        common = prefix(maximum)
        rows = [dict(even_full=e, odd_full=o, indexed_source_count_upper=n,
            failure_log2_upper=str(A.uniform.rational_ceiling(common+A.log2(term))))
            for e, o, n, term in terms]
        value = common+A.log2(sum((term for _, _, _, term in terms), arb(0)))
        ceiling = A.uniform.rational_ceiling(value)
    if dependencies() != pins:
        raise ValueError('Reviewed dependencies changed during replay')
    return dict(construction=binding.SCOPE['construction'], field_order=256,
        active_groups=R, every_group_arity=ARITY, all1820_labelled_rowsets=True,
        mixed_arities_covered=False, other_active_group_weights_covered=False,
        all99_indexed_source_profiles=True, profile_count=len(rows), profiles=rows,
        original_A_fixed_across_columns=True, moment_order=3,
        maximum_column_moment=str(maximum), moment_comparison='normalized-L3-before-fractionalization',
        mark=str(MARK), uniform_alpha=str(ALPHA), outer_theta=str(THETA),
        rowset_choices_per_group=comb(16, ARITY), rowset_union_inside_theta=True,
        indexed_source_union_outside_theta=True, projective_normalizations=1,
        source_dimension=1048576, output_length=2097152, forbidden_weight=734003,
        precision_bits=precision, fresh_regional_banks=1, directed=True,
        whole_code_certificate=False, failure_log2_upper=str(ceiling),
        failure_log2_enclosure=str(value), selected_class_reaches_60_bits=bool(ceiling <= -60),
        reviewed_source_sha256=pins, exact_moment_helper_sha256=inventory.get('source_sha256'),
        elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--inventory', type=Path)
    modes.add_argument('--compute-inventory', action='store_true')
    parser.add_argument('--precision', type=int, choices=(256, 384), default=256)
    parser.add_argument('--max-seconds', type=float, default=360)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite evidence')
    if args.compute_inventory:
        result = moment_inventory(max_seconds=args.max_seconds,
            progress=lambda row: print(json.dumps(row), flush=True))
    else:
        raw = args.inventory.read_bytes()
        result = receipt(json.loads(raw), args.precision)
        result['inventory_sha256'] = sha256(raw).hexdigest()
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('rows', 'profiles', 'reviewed_source_sha256')}))
