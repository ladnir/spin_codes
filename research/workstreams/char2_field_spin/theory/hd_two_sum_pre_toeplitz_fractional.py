"""Fixed-A pre-H pair comparison, composed before the fractional power.

Exact rational mask masses; FLOAT powers and scalar summaries. No regional
provider, setup-failure claim, row-set resampling, or native work. The only
production scope is a=8, u=0..4, alpha=1/5. Group8 supports tiny controls.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import comb, factorial, isfinite, lcm, log2

from hd_two_sum_input_weight import exact_support_count, marked_count
from hd_two_sum_pre_toeplitz_pairs import pre_pair_support_masses
from hd_two_sum_toeplitz_pairs import checked_matching


ALPHA = F(1, 5)
MODES = ('iid-majorant', 'exact-support', 'diagonal')


def checked_scope(group, arity, internal):
    if (group not in (8, 16) or any(type(x) is not int for x in (group, arity, internal))
            or not 0 <= arity <= group or not max(0, arity-group//2) <= internal <= arity//2):
        raise ValueError('Feasible g8/g16 fixed row set and matching are required')
    return arity-2*internal, group//2-arity+internal


def matching_for_internal(arity, internal, group=16):
    """A is always the FIRST a coordinates; only the fixed matching changes."""
    crossing, empty = checked_scope(group, arity, internal)
    return (tuple((2*i, 2*i+1) for i in range(internal))
            +tuple((2*internal+i, arity+i) for i in range(crossing))
            +tuple((arity+crossing+2*i, arity+crossing+2*i+1) for i in range(empty)))


def compositions(total, parts):
    if parts == 1:
        yield (total,)
    else:
        for first in range(total+1):
            for suffix in compositions(total-first, parts-1):
                yield (first,)+suffix


def multinomial(counts):
    result = factorial(sum(counts))
    for count in counts:
        result //= factorial(count)
    return result


@lru_cache(None)
def target_types(arity, internal, group=16):
    """One representative per joint A/T pair type, exact labelled T count.

    AA/empty singleton targets have two orientations. Crossing pairs keep
    their A leg and non-A leg distinct. No C(group,a) row-set factor.
    """
    crossing, empty = checked_scope(group, arity, internal)
    matching = matching_for_internal(arity, internal, group)
    rows = []
    for inside in compositions(internal, 3):
        for cross in compositions(crossing, 4):
            for outside in compositions(empty, 3):
                masks = (tuple(m for m, count in zip((0, 1, 3), inside) for _ in range(count))
                         +tuple(m for m, count in enumerate(cross) for _ in range(count))
                         +tuple(m for m, count in zip((0, 1, 3), outside) for _ in range(count)))
                target = sum(((mask & 1) << left)+((mask >> 1) << right)
                             for mask, (left, right) in zip(masks, matching))
                multiplicity = (multinomial(inside)*multinomial(cross)*multinomial(outside)
                                *2**(inside[1]+outside[1]))
                rows.append((target, multiplicity, inside, cross, outside))
    expected = comb(internal+2, 2)*comb(empty+2, 2)*comb(crossing+3, 3)
    if len(rows) != expected or sum(row[1] for row in rows) != 2**group:
        raise ArithmeticError('Target types must partition all labelled masks exactly')
    return tuple(rows)


@lru_cache(None)
def scaled_pair_terms(active, target, marker, mode):
    """Integer local polynomial plus common denominator, cached by pair type."""
    masses = pre_pair_support_masses(active, input_marker=marker, mode=mode)
    values = tuple(mass/255**b.bit_count() for b, mass in enumerate(masses))
    denominator = lcm(*(value.denominator for value in values))
    terms = {}
    for b, value in enumerate(values):
        if value:
            key = (b.bit_count(), (b & target).bit_count())
            terms[key] = terms.get(key, 0)+value.numerator*(denominator//value.denominator)
    return denominator, tuple((k, overlap, value) for (k, overlap), value in terms.items())


def fixed_mask_mass(active, matching, target, *, group=16, input_marker=F(1), mode='iid-majorant'):
    """C_A,M(T;v), after the ENTIRE positive composition, as a Fraction.

    Sum_B L_A,M(B;v)/255^|B| * N_NZ(|B|,|B intersect T|,|T|).
    No RS normalization, alpha, source count, or inner operator is applied.
    """
    matching = checked_matching(matching, group)
    marker = F(input_marker)
    if (group not in (8, 16) or any(type(x) is not int or not 0 <= x < 1 << group
                                   for x in (active, target))
            or marker < 1 or mode not in MODES):
        raise ValueError('Fixed labelled masks, v>=1, and declared comparison required')
    table, denominator = {(0, 0): 1}, 1
    for left, right in matching:
        a = ((active >> left) & 1) | (((active >> right) & 1) << 1)
        t = ((target >> left) & 1) | (((target >> right) & 1) << 1)
        scale, branches = scaled_pair_terms(a, t, marker, mode)
        updated = {}
        for (k, overlap), count in table.items():
            for dk, dl, weight in branches:
                key = (k+dk, overlap+dl)
                updated[key] = updated.get(key, 0)+count*weight
        table, denominator = updated, denominator*scale
    weight = target.bit_count()
    numerator = sum(count*exact_support_count(k, overlap, weight, group)
                    for (k, overlap), count in table.items())
    return F(numerator, denominator)


def exact_type_masses(arity, internal, *, group=16, input_marker=F(1), mode='iid-majorant', check_time=None):
    """Compact exact mask table; old control uses the SAME A and input mark."""
    checked_scope(group, arity, internal)
    matching, active = matching_for_internal(arity, internal, group), (1 << arity)-1
    rows = []
    for target, multiplicity, *_ in target_types(arity, internal, group):
        if check_time:
            check_time()
        new = fixed_mask_mass(active, matching, target, group=group, input_marker=input_marker, mode=mode)
        old = marked_count(arity, (active & target).bit_count(), target.bit_count(), input_marker, group)
        rows.append((target, multiplicity, new, old))
    return tuple(rows)


def logadd(left, right):
    if left == -float('inf'):
        return right
    if right == -float('inf'):
        return left
    top, bottom = max(left, right), min(left, right)
    return top+log2(1+2**(bottom-top))


def fractional_coefficients(rows, *, group=16, alpha=ALPHA, columns=16):
    """Log2 sum_{|T|=w} C(T)^(columns*alpha), multiplicity OUTSIDE power."""
    alpha = F(alpha)
    if not 0 < alpha <= 1 or type(columns) is not int or columns < 1 or group not in (8, 16):
        raise ValueError('Positive fractional exponent, columns, and supported group required')
    result = [[-float('inf')]*(group+1) for _ in range(2)]
    for target, multiplicity, new, old in rows:
        if type(multiplicity) is not int or multiplicity <= 0 or min(new, old) < 0:
            raise ValueError('Positive multiplicities and nonnegative exact masses required')
        for destination, value in zip(result, (F(new), F(old))):
            if value:
                term = log2(multiplicity)+float(columns*alpha)*(log2(value.numerator)-log2(value.denominator))
                weight = target.bit_count()
                destination[weight] = logadd(destination[weight], term)
    return tuple(result[0]), tuple(result[1])


def scalar_summary(log_coefficients, output_marker, log2_activity=0.):
    """FLOAT scalar diagnostic sum_w F(w) z^w y^[w>0]; not a regional bound."""
    z = F(output_marker)
    if not 0 < z <= 1 or not isfinite(log2_activity):
        raise ValueError('A scalar z in(0,1] and finite log activity mark required')
    total = -float('inf')
    for weight, value in enumerate(log_coefficients):
        total = logadd(total, value+weight*log2(float(z))+(log2_activity if weight else 0.))
    return total


def saved_logs(values):
    return [None if x == -float('inf') else x for x in values]


def scope_metadata():
    return dict(group=16, arity=8, fixed_active_set=255, internal_pair_counts=list(range(5)),
                alpha=str(ALPHA), power='16*alpha', target_multiplicity_outside_power=True,
                positive_mask_composition_before_power=True, row_set_resampled=False,
                matching_fixed_per_column=True, matching_cases_are_separate_controls=True,
                complete_sixteen_column_schedule_evaluated=False,
                comparison_mode='iid-majorant', independent_output_diagonal_retained=True,
                artificial_input_diagonal_measure_bridge=True,
                exact_rational_mask_masses=True, fractional_powers_are_FLOAT=True,
                regional_banks_built=0, witness_fits=0, source_count_applied=False,
                RS_normalization_applied=False, RS_floor_applied=False,
                directed=False, whole_code_certificate=False)


FROZEN_BANK = 'hd_two_sum_gf16_source_strata_anchors_float.json'
FROZEN_SHA = '5050033694bb6a1ef7852aa309ee3abca252c9740fc713a83aa8945a380fb1cf'


def frozen_scalar_marks(directory=None):
    """Read only literal metadata, never construct a regional operator."""
    from hashlib import sha256
    import json
    from pathlib import Path
    directory = Path(__file__).resolve().parent if directory is None else Path(directory)
    raw = (directory/FROZEN_BANK).read_bytes()
    if sha256(raw).hexdigest() != FROZEN_SHA:
        raise ValueError('Frozen bank bytes changed')
    source = json.loads(raw)
    expected = ((18, '691/700'), (24, '172/175'), (32, '171/175'))
    if len(source.get('banks', [])) != 3 or source.get('directed') is not False:
        raise ValueError('The frozen three-anchor FLOAT bank is required')
    rows = []
    for bank, (r, z) in zip(source['banks'], expected):
        if (bank['active_groups'], bank['mark'], bank['alpha']) != (r, z, str(ALPHA)):
            raise ValueError('Frozen anchor changed')
        templates = [row for row in bank['type_templates'] if row['name'].startswith('old-')]
        if len(templates) != 3:
            raise ValueError('Exactly the three original activity templates are required')
        for row in templates:
            if row['family'] != 'nonempty-column-floor' or len(row['markers']) != 16:
                raise ValueError('Complete fixed arity activity markers are required')
            log_activity = F(row['markers'][7])
            if log_activity < 0:
                raise ValueError('An activity lower-floor marker must be at least1')
            rows.append(dict(active_groups=r, output_marker=z, template=row['name'],
                             log2_activity_marker=str(log_activity)))
    rows.append(dict(active_groups=None, output_marker='1/2', template='standalone-half', log2_activity_marker='0'))
    return rows


def subset_summary(coefficients, first, last):
    value = -float('inf')
    for coefficient in coefficients[first:last+1]:
        value = logadd(value, coefficient)
    return None if value == -float('inf') else value


def run(max_seconds=30, progress=None):
    """Only five fixed matching cases; v1, no fit or regional provider."""
    from time import perf_counter
    if type(max_seconds) not in (int, float) or not isfinite(max_seconds) or not 0 < max_seconds <= 30:
        raise ValueError('An explicit budget of at most30 seconds is required')
    started = perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('The bounded per-mask control exhausted its budget')
    result = scope_metadata()
    marks = frozen_scalar_marks()
    result.update(input_marker='1', frozen_bank=FROZEN_BANK, frozen_bank_sha256=FROZEN_SHA,
                  scalar_marks=marks, scalar_mark_interpretation=
                  'sum_w F(w)*z^w*2^(16*ell*1[w>0]); diagnostic only, not old W_a',
                  compute_budget_seconds=max_seconds, complete_requested_cases=False, results=[])
    try:
        for internal in range(5):
            rows = exact_type_masses(8, internal, input_marker=F(1), check_time=check_time)
            zero = [row for row in rows if row[0] == 0]
            if len(zero) != 1 or zero[0][1:] != (1, F(1), F(1)):
                raise ArithmeticError('Both exact all-zero mask masses must equal1')
            new, old = fractional_coefficients(rows)
            summaries = []
            for marks_row in marks:
                # ell is already a post-alpha activity log in the frozen bank.
                # This diagnostic uses its16th-power value; no floor is enforced.
                ell = 16*float(F(marks_row['log2_activity_marker']))
                new_sum = scalar_summary(new, F(marks_row['output_marker']), ell)
                old_sum = scalar_summary(old, F(marks_row['output_marker']), ell)
                summaries.append(dict(marks_row, log2_applied_activity=ell,
                                      new_log2_sum=new_sum, old_log2_sum=old_sum,
                                      old_minus_new_bits=old_sum-new_sum))
            subsets = []
            for first, last in ((0, 0), (1, 4), (5, 8), (9, 16), (1, 16)):
                new_sum, old_sum = subset_summary(new, first, last), subset_summary(old, first, last)
                subsets.append(dict(weights=[first, last], new_log2_sum=new_sum, old_log2_sum=old_sum,
                    old_minus_new_bits=None if new_sum is None or old_sum is None else old_sum-new_sum))
            row = dict(internal_pairs=internal, crossing_pairs=8-2*internal, empty_pairs=internal,
                       fixed_A=255, matching=[list(pair) for pair in matching_for_internal(8, internal)],
                       target_types=len(rows), labelled_target_masks=sum(x[1] for x in rows),
                       zero_mass_new='1', zero_mass_old='1',
                       new_log2_fractional_coefficients=saved_logs(new),
                       old_log2_fractional_coefficients=saved_logs(old),
                       weight_subsets=subsets, scalar_summaries=summaries)
            result['results'].append(row)
            if progress:
                progress(dict(internal_pairs=internal, target_types=len(rows), weight_subsets=subsets,
                              elapsed_seconds=perf_counter()-started))
        check_time()
        result.update(status='completed', complete_requested_cases=True)
    except TimeoutError as error:
        result.update(status='compute-budget-exceeded', reason=str(error))
    result['elapsed_seconds'] = perf_counter()-started
    return result


if __name__ == '__main__':
    import argparse
    from hashlib import sha256
    import json
    from pathlib import Path
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-seconds', type=float, default=30)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output must be a new artifact')
    result = run(args.max_seconds, progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    result['dependency_pins'] = {name: sha256((Path(__file__).resolve().parent/name).read_bytes()).hexdigest()
        for name in ('hd_two_sum_pre_toeplitz_pairs.py', 'hd_two_sum_toeplitz_pairs.py',
                     'hd_two_sum_input_weight.py', 'hd_two_sum_closed_kernel.py')}
    with args.output.open('x', encoding='utf-8', newline='\n') as destination:
        destination.write(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(output=str(args.output), sha256=sha256(args.output.read_bytes()).hexdigest(),
                         status=result['status'], elapsed_seconds=result['elapsed_seconds'])))
