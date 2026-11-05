"""Nested positive eight-packet FLOAT refinement for the literal construction.

Keep the complete frozen four-packet comparator. Refine only W^8, where W
is the clean subset (s=0) of V (s<=1,f<=2). There are four W types,256 raw
W^4 records, and at most65536 ordered W^8 products. The positive partition
T4=C_W4+R4 gives T8=C_W8+C_W4*R4+R4*T4 <= T4^2 coefficientwise.
No numerical result from this provider is directed or a whole-code claim.
"""

import ast
from fractions import Fraction as F
from hashlib import sha256
from math import isfinite, log2
from pathlib import Path

import numpy as np

import hd_field128_lift_four_packet_screen as four


old = four.old
CONSTRUCTION = four.CONSTRUCTION


def combine_raw(records, maximum, budget):
    """Ordered concatenation of raw records without alpha or aggregation."""
    singles, fulls, weights, logs = records
    next_s, next_f, next_w, next_logs = [], [], [], []
    for i in range(len(singles)):
        budget.check()
        ss, ff = singles[i]+singles, fulls[i]+fulls
        keep = ss+2*ff <= maximum
        next_s.extend(ss[keep])
        next_f.extend(ff[keep])
        next_w.extend(weights[i]+weights[keep])
        next_logs.extend(four.compose(logs[i], logs[keep]))
    return (np.asarray(next_s, dtype=int), np.asarray(next_f, dtype=int),
            np.asarray(next_w), np.asarray(next_logs).reshape((-1, 2, 2)))


def clean_quartets(z, q, maximum, budget, max_v_full=2):
    if type(max_v_full) is not int or not 0 <= max_v_full <= 2:
        raise ValueError('Expected full-pair cutoff0..2')
    s, f, counts, indices = old.good_geometry()
    keep = (s == 0) & (f <= max_v_full) & (s+2*f <= maximum)
    records = (s[keep], f[keep], np.log2(counts[keep]), old.local_logs(z, q, old.OCC[indices[keep]]))
    return combine_raw(combine_raw(records, maximum, budget), maximum, budget)


def eight_operator(z, alpha, maximum, q, budget, max_v_full=2):
    """Return the refined eight-step and complete original four-step caps."""
    cvv, remainder, cv4 = four.positive_parts(z, alpha, maximum, q, budget, max_v_full)
    two = np.logaddexp2(cvv, remainder)
    cross = np.logaddexp2(old.polynomial_product(cvv, remainder, maximum, budget),
                         old.polynomial_product(remainder, two, maximum, budget))
    clean, nonclean = cv4.copy(), cv4.copy()
    clean[:, 1:, :], nonclean[:, :1, :] = old.NEG, old.NEG
    # Nonnegative singleton degrees mean s=0 selects exactly W^4, not a superset.
    r4 = np.logaddexp2(nonclean, cross)
    baseline = np.logaddexp2(clean, r4)
    raw_quartets = clean_quartets(z, q, maximum, budget, max_v_full)
    clean8 = four.four_from_pairs(raw_quartets, alpha, maximum, budget)
    refined = np.logaddexp2(clean8, old.polynomial_product(clean, r4, maximum, budget))
    refined = np.logaddexp2(refined, old.polynomial_product(r4, baseline, maximum, budget))
    return refined, baseline


def regional_screen(z, alpha, maximum, *, q=old.FIELD_ORDER, max_v_full=2,
                    max_seconds=60., _budget=None):
    old.checked_q(q, construction=True)
    if (type(maximum) is not int or not 1 <= maximum <= old.MAXIMUM
            or not isfinite(z) or not isfinite(alpha) or not 0 < z < 1 or not 0 < alpha <= 1
            or type(max_v_full) is not int or not 0 <= max_v_full <= 2):
        raise ValueError('Expected cutoff1..64, V cutoff0..2, and finite valid marks')
    budget = _budget or old.Budget(max_seconds)
    if alpha == 1:
        return old.regional_screen(z, alpha, maximum, q=q, paired=False, _budget=budget)
    operator, _ = eight_operator(z, alpha, maximum, q, budget, max_v_full)
    run = old.polynomial_power(operator, 16, maximum, budget)
    regions = {}
    for r in range(maximum+1):
        for f in range(r//2+1):
            s = r-2*f
            values = run[:, s, f]-old.log_multinomial((old.PAIRS-s-f, s, f))
            regions[r, f] = tuple(map(tuple, values.reshape((2, 2))))
    budget.check()
    return regions


def weight_screen(r, z, alpha, *, q=old.FIELD_ORDER, distance=F(7, 20), max_v_full=2,
                  outer_measure='invertible-toeplitz', profile=None, max_seconds=60.):
    bits = old.checked_q(q, construction=True)
    distance = F(distance)
    if not 0 < distance < F(1, 2):
        raise ValueError('Expected distance strictly between0and1/2')
    if profile is not None and (len(profile) != 2 or any(type(x) is not int or not 0 <= x <= r//2 for x in profile)):
        raise ValueError('Invalid matching profile')
    budget = old.Budget(max_seconds)
    regions = regional_screen(z, alpha, r, q=q, max_v_full=max_v_full, _budget=budget)
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
                composition_steps=8, pair_swaps=True, max_v_full_pairs=max_v_full,
                eight_refinement_subset='clean arrangements s=0 inside V',
                positive_partition='C_W8 + C_W4*R4 + R4*T4',
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
    parser.add_argument('--profile', type=int, nargs=2)
    parser.add_argument('--max-seconds', type=float, default=60.)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.output and args.output.exists():
        parser.error('Output must be a new file')
    if not 8 <= args.field_bits <= 128 or args.field_bits%8:
        parser.error('Field bits must be a multiple of8in8..128')
    result = weight_screen(args.r, args.z, args.alpha, q=1 << args.field_bits, distance=F(args.distance),
                           max_v_full=args.max_v_full, profile=args.profile, max_seconds=args.max_seconds)
    result['source_sha256'] = source_pins()
    if args.output:
        with args.output.open('x', encoding='utf-8') as stream:
            json.dump(result, stream, indent=2, allow_nan=False)
            stream.write('\n')
    print(json.dumps(result, allow_nan=False))


if __name__ == '__main__':
    main()
