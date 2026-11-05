"""Selected-weight FLOAT regional Jensen cap for the literal construction.

Fix a support histogram h and a region route pi. Let P_h(pi) be the ordered
product of the region's 128 raw nonnegative packet matrices. The common
regional target is Q_h=E_pi[P_h(pi)^circ(alpha)]. Concavity gives
Q_h <= J_h=(E_pi[P_h(pi)])^circ(alpha). The exact histogram normalizer
belongs inside this power. Route probabilities are not powered separately.

The frozen eight-packet cap also bounds Q_h, so their entrywise minimum is
valid. Neither cap need dominate the other. Regional routes are independent;
the resulting two-state matrices compose without resetting the actual state.
This changes only the proof calculation. All output remains nondirected FLOAT.
"""

from fractions import Fraction as F
from math import isfinite, log2
from pathlib import Path

import numpy as np

import hd_field128_lift_eight_packet_screen as eight


old = eight.old
CONSTRUCTION = eight.CONSTRUCTION
MODES = ('min8', 'jensen')


def regional_tables(z, alpha, maximum, *, q=old.FIELD_ORDER, mode='min8', max_v_full=2,
                    max_seconds=60., _budget=None):
    """Return (selected cap, pure Jensen cap), both as normalized log tables."""
    old.checked_q(q, construction=True)
    if (type(maximum) is not int or not 1 <= maximum <= old.MAXIMUM
            or not isfinite(z) or not isfinite(alpha) or not 0 < z < 1 or not 0 < alpha <= 1
            or mode not in MODES or type(max_v_full) is not int or not 0 <= max_v_full <= 2):
        raise ValueError('Expected cutoff1..64, finite valid marks, V cutoff0..2, and mode min8 or jensen')
    budget = _budget or old.Budget(max_seconds)
    # This table already divides by the regional multinomial and averages
    # independent fair singleton swaps. Only then apply entrywise alpha.
    mean = old.regional_screen(z, 1., maximum, q=q, paired=False, _budget=budget)
    jensen = {key: tuple(map(tuple, alpha * np.asarray(value))) for key, value in mean.items()}
    if mode == 'jensen' or alpha == 1:
        return jensen, jensen
    baseline = eight.regional_screen(z, alpha, maximum, q=q, max_v_full=max_v_full, _budget=budget)
    selected = {key: tuple(map(tuple, np.minimum(jensen[key], baseline[key]))) for key in jensen}
    budget.check()
    return selected, jensen


def regional_screen(z, alpha, maximum, *, q=old.FIELD_ORDER, mode='min8', max_v_full=2,
                    max_seconds=60., _budget=None):
    return regional_tables(z, alpha, maximum, q=q, mode=mode, max_v_full=max_v_full,
                           max_seconds=max_seconds, _budget=_budget)[0]


def _profile_values(regions, candidates, r, prefix, budget):
    values = []
    for a, b in candidates:
        budget.check()
        count = old.profile_count_upper(old.PAIRS, r, a, b)
        if count:
            run = old.power(old.multiply(regions[r, a], regions[r, b]), old.REGIONS // 2)
            values.append((float(log2(count) + prefix + old.total(run[0])), a, b))
    if not values:
        raise ValueError('The requested profile has no support')
    total = float(old.total(row[0] for row in values))
    if not isfinite(total):
        raise ArithmeticError('Nonfinite selected-scope FLOAT result')
    return total, max(values), len(values)


def weight_screen(r, z, alpha, *, q=old.FIELD_ORDER, distance=F(7, 20), mode='min8', max_v_full=2,
                  outer_measure='invertible-toeplitz', profile=None, max_seconds=60.):
    """All profiles at r, or one explicitly scoped term; retain pure Jensen."""
    bits = old.checked_q(q, construction=True)
    distance = F(distance)
    if not 0 < distance < F(1, 2):
        raise ValueError('Expected distance strictly between0and1/2')
    if profile is not None and (len(profile) != 2 or any(type(x) is not int or not 0 <= x <= r//2 for x in profile)):
        raise ValueError('Invalid matching profile')
    budget = old.Budget(max_seconds)
    selected, jensen = regional_tables(z, alpha, r, q=q, mode=mode, max_v_full=max_v_full, _budget=budget)
    forbidden = int(distance * old.OUTPUT_LENGTH)
    prefix = alpha * (r * old.outer_log_density(q, outer_measure)
                      - (bits + old.log1mexp2(-bits)) - forbidden * log2(z))
    candidates = (tuple((a, b) for a in range(r//2 + 1) for b in range(r//2 + 1))
                  if profile is None else (tuple(profile),))
    value, worst, count = _profile_values(selected, candidates, r, prefix, budget)
    pure, pure_worst, _ = ((value, worst, count) if selected is jensen else
                           _profile_values(jensen, candidates, r, prefix, budget))
    improved = sum(int(np.any(np.asarray(selected[key]) < np.asarray(jensen[key]))) for key in selected)
    budget.check()
    return dict(construction=CONSTRUCTION, active_groups=r, z=z, alpha=alpha,
                field_order=q, field_bits=bits, source_dimension=old.DIMENSION, output_length=old.OUTPUT_LENGTH,
                outer_groups=old.GROUPS, target_distance=str(distance), forbidden_weight=forbidden,
                outer_measure=outer_measure, projective_divisor=str(q - 1), pair_swaps=True,
                state_persists_across_regions=True, full_field_random_coefficients=True,
                combination_mode=mode, conditional_region_average=True,
                region_packets=old.STEPS, normalizer_inside_fractional_power=True,
                min8_max_v_full_pairs=max_v_full, tables_improved_over_pure_jensen=improved,
                log2_upper_screen=value, pure_jensen_log2_upper_screen=pure,
                worst_profile=worst, pure_jensen_worst_profile=pure_worst, profile_count=count,
                selected_profile=None if profile is None else list(profile),
                covers_all_supports_at_selected_weight=profile is None,
                covers_all_even_odd_splits=profile is None, encoder_changed=False,
                coefficient_method='alpha1-positive-truncated-region-mean-then-entrywise-alpha-FLOAT',
                arithmetic='FLOAT', directed=False, whole_code_certificate=False,
                scope='selected weight only' if profile is None else 'one matching-profile term only')


def source_pins():
    # Static closure includes the optional min8 branch even in pure mode.
    import ast
    from hashlib import sha256
    directory = Path(__file__).resolve().parent
    pending, seen = [Path(__file__).name], {}
    while pending:
        name = pending.pop()
        if name in seen:
            continue
        raw = (directory / name).read_bytes()
        seen[name] = sha256(raw).hexdigest()
        for node in ast.walk(ast.parse(raw)):
            modules = ([item.name for item in node.names] if isinstance(node, ast.Import)
                       else [node.module] if isinstance(node, ast.ImportFrom) and node.module else [])
            for module in modules:
                target = module.split('.')[0] + '.py'
                if (directory / target).is_file() and target not in seen:
                    pending.append(target)
    return dict(sorted(seen.items()))


def main():
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--r', type=int, required=True)
    parser.add_argument('--z', type=float, required=True)
    parser.add_argument('--alpha', type=float, required=True)
    parser.add_argument('--distance', default='7/20')
    parser.add_argument('--field-bits', type=int, default=128)
    parser.add_argument('--mode', choices=MODES, default='min8')
    parser.add_argument('--max-v-full', type=int, default=2)
    parser.add_argument('--profile', type=int, nargs=2)
    parser.add_argument('--max-seconds', type=float, default=60.)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.output and args.output.exists():
        parser.error('Output must be a new file')
    if not 8 <= args.field_bits <= 128 or args.field_bits % 8:
        parser.error('Field bits must be a multiple of8in8..128')
    result = weight_screen(args.r, args.z, args.alpha, q=1 << args.field_bits, distance=F(args.distance),
                           mode=args.mode, max_v_full=args.max_v_full, profile=args.profile,
                           max_seconds=args.max_seconds)
    result['source_sha256'] = source_pins()
    if args.output:
        with args.output.open('x', encoding='utf-8') as stream:
            json.dump(result, stream, indent=2, allow_nan=False)
            stream.write('\n')
    print(json.dumps(result, allow_nan=False))


if __name__ == '__main__':
    main()
