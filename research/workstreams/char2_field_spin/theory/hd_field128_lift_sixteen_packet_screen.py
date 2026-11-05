"""Nested positive sixteen-packet FLOAT refinement of the literal ensemble.

Within the frozen T8 comparator, select clean W^8 arrangements whose total
full-pair count is at most small_full. For the default small_full=2 there
are53 ordered raw records. Write T8=Csmall8+R8 as a positive partition;
T16=Csmall16+Csmall8*R8+R8*T8 <= T8^2. No arithmetic here is directed.
"""

import ast
from fractions import Fraction as F
from hashlib import sha256
from math import isfinite, log2
from pathlib import Path

import numpy as np

import hd_field128_lift_eight_packet_screen as eight


four, old = eight.four, eight.old
CONSTRUCTION = eight.CONSTRUCTION


def small_octets(z, q, maximum, budget, max_v_full=2, small_full=2):
    if type(small_full) is not int or not 0 <= small_full <= 3:
        raise ValueError('This bounded provider supports small-eight full counts0..3')
    cutoff = min(maximum, 2*small_full)
    quartets = eight.clean_quartets(z, q, cutoff, budget, max_v_full)
    return eight.combine_raw(quartets, cutoff, budget)


def eight_partition(z, alpha, maximum, q, budget, max_v_full=2, small_full=2):
    """Return Csmall8 and R8, retaining every frozen T8 case positively."""
    if type(small_full) is not int or not 0 <= small_full <= 3:
        raise ValueError('Expected small-eight full count0..3')
    cvv, r2, cv4 = four.positive_parts(z, alpha, maximum, q, budget, max_v_full)
    two = np.logaddexp2(cvv, r2)
    cross4 = np.logaddexp2(old.polynomial_product(cvv, r2, maximum, budget),
                          old.polynomial_product(r2, two, maximum, budget))
    clean4, other4 = cv4.copy(), cv4.copy()
    clean4[:, 1:, :], other4[:, :1, :] = old.NEG, old.NEG
    r4 = np.logaddexp2(other4, cross4)
    four_step = np.logaddexp2(clean4, r4)
    quartets = eight.clean_quartets(z, q, maximum, budget, max_v_full)
    clean8 = four.four_from_pairs(quartets, alpha, maximum, budget)
    selected, others = clean8.copy(), clean8.copy()
    selected[:, :, small_full+1:], others[:, :, :small_full+1] = old.NEG, old.NEG
    cross8 = np.logaddexp2(old.polynomial_product(clean4, r4, maximum, budget),
                          old.polynomial_product(r4, four_step, maximum, budget))
    return selected, np.logaddexp2(others, cross8)


def sixteen_operator(z, alpha, maximum, q, budget, max_v_full=2, small_full=2):
    selected, remainder = eight_partition(z, alpha, maximum, q, budget, max_v_full, small_full)
    baseline = np.logaddexp2(selected, remainder)
    raw_octets = small_octets(z, q, maximum, budget, max_v_full, small_full)
    composed = four.four_from_pairs(raw_octets, alpha, maximum, budget)
    refined = np.logaddexp2(composed, old.polynomial_product(selected, remainder, maximum, budget))
    refined = np.logaddexp2(refined, old.polynomial_product(remainder, baseline, maximum, budget))
    return refined, baseline


def regional_screen(z, alpha, maximum, *, q=old.FIELD_ORDER, max_v_full=2, small_full=2,
                    max_seconds=60., _budget=None):
    old.checked_q(q, construction=True)
    if (type(maximum) is not int or not 1 <= maximum <= old.MAXIMUM
            or not isfinite(z) or not isfinite(alpha) or not 0 < z < 1 or not 0 < alpha <= 1
            or type(max_v_full) is not int or not 0 <= max_v_full <= 2
            or type(small_full) is not int or not 0 <= small_full <= 3):
        raise ValueError('Expected cutoff1..64, V cutoff0..2, small count0..3, and finite valid marks')
    budget = _budget or old.Budget(max_seconds)
    if alpha == 1:
        return old.regional_screen(z, alpha, maximum, q=q, paired=False, _budget=budget)
    operator, _ = sixteen_operator(z, alpha, maximum, q, budget, max_v_full, small_full)
    run = old.polynomial_power(operator, 8, maximum, budget)
    regions = {}
    for r in range(maximum+1):
        for f in range(r//2+1):
            s = r-2*f
            values = run[:, s, f]-old.log_multinomial((old.PAIRS-s-f, s, f))
            regions[r, f] = tuple(map(tuple, values.reshape((2, 2))))
    budget.check()
    return regions


def weight_screen(r, z, alpha, *, q=old.FIELD_ORDER, distance=F(7, 20), max_v_full=2,
                  small_full=2, outer_measure='invertible-toeplitz', profile=None, max_seconds=60.):
    bits = old.checked_q(q, construction=True)
    distance = F(distance)
    if not 0 < distance < F(1, 2):
        raise ValueError('Expected distance strictly between0and1/2')
    if profile is not None and (len(profile) != 2 or any(type(x) is not int or not 0 <= x <= r//2 for x in profile)):
        raise ValueError('Invalid matching profile')
    budget = old.Budget(max_seconds)
    regions = regional_screen(z, alpha, r, q=q, max_v_full=max_v_full, small_full=small_full, _budget=budget)
    forbidden = int(distance*old.OUTPUT_LENGTH)
    prefix = alpha*(r*old.outer_log_density(q, outer_measure)-(bits+old.log1mexp2(-bits))-forbidden*log2(z))
    candidates = ((a, b) for a in range(r//2+1) for b in range(r//2+1)) if profile is None else (tuple(profile),)
    profiles = []
    for a, b in candidates:
        budget.check()
        count = old.profile_count_upper(old.PAIRS, r, a, b)
        if count:
            run = old.power(old.multiply(regions[r, a], regions[r, b]), old.REGIONS//2)
            profiles.append((float(log2(count)+prefix+old.total(run[0])), a, b))
    if not profiles:
        raise ValueError('The requested profile has no support')
    value = float(old.total(row[0] for row in profiles))
    if not isfinite(value):
        raise ArithmeticError('Nonfinite selected-scope FLOAT result')
    budget.check()
    return dict(construction=CONSTRUCTION, active_groups=r, z=z, alpha=alpha, field_order=q, field_bits=bits,
                target_distance=str(distance), forbidden_weight=forbidden, outer_measure=outer_measure,
                composition_steps=16, pair_swaps=True, max_v_full_pairs=max_v_full,
                small_eight_total_full_pairs=small_full,
                positive_partition='Csmall16 + Csmall8*R8 + R8*T8',
                log2_upper_screen=value, worst_profile=max(profiles), profile_count=len(profiles),
                selected_profile=None if profile is None else list(profile),
                covers_all_supports_at_selected_weight=profile is None,
                covers_all_even_odd_splits=profile is None, encoder_changed=False,
                coefficient_method='positive-truncated-bivariate-log2-FLOAT',
                arithmetic='FLOAT', directed=False, whole_code_certificate=False,
                scope='selected weight only' if profile is None else 'one matching-profile term only')


def source_pins():
    directory = Path(__file__).resolve().parent
    pending, seen = [Path(__file__).name], {}
    while pending:
        name = pending.pop()
        if name in seen:
            continue
        raw = (directory/name).read_bytes()
        seen[name] = sha256(raw).hexdigest()
        for node in ast.walk(ast.parse(raw)):
            modules = ([alias.name for alias in node.names] if isinstance(node, ast.Import)
                       else [node.module] if isinstance(node, ast.ImportFrom) and node.module else [])
            for module in modules:
                target = module.split('.')[0]+'.py'
                if (directory/target).is_file() and target not in seen:
                    pending.append(target)
    return dict(sorted(seen.items()))


def main():
    import argparse, json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--r', type=int, required=True)
    parser.add_argument('--z', type=float, required=True)
    parser.add_argument('--alpha', type=float, required=True)
    parser.add_argument('--distance', default='7/20')
    parser.add_argument('--field-bits', type=int, default=128)
    parser.add_argument('--max-v-full', type=int, default=2)
    parser.add_argument('--small-full', type=int, default=2)
    parser.add_argument('--profile', type=int, nargs=2)
    parser.add_argument('--max-seconds', type=float, default=60.)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.output and args.output.exists():
        parser.error('Output must be a new file')
    if not 8 <= args.field_bits <= 128 or args.field_bits%8:
        parser.error('Field bits must be a multiple of8in8..128')
    result = weight_screen(args.r, args.z, args.alpha, q=1 << args.field_bits, distance=F(args.distance),
                           max_v_full=args.max_v_full, small_full=args.small_full,
                           profile=args.profile, max_seconds=args.max_seconds)
    result['source_sha256'] = source_pins()
    if args.output:
        with args.output.open('x', encoding='utf-8') as stream:
            json.dump(result, stream, indent=2, allow_nan=False)
            stream.write('\n')
    print(json.dumps(result, allow_nan=False))


if __name__ == '__main__':
    main()
