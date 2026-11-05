"""Bounded FLOAT four-packet refinement for the frozen literal construction.

V consists of good arrangements with at most max_v_full full pairs. The
positive two-packet partition is T2=C_VV+R. Compose four raw V matrices
before alpha, then use T4=C_VVVV+C_VV*R+R*(C_VV+R). This equals the
four-term positive partition and never subtracts rounded upper bounds.
No encoder changes, numerical certificates, or whole-code claims occur here.
"""

import ast
from fractions import Fraction as F
from hashlib import sha256
from math import isfinite, log2
from pathlib import Path

import numpy as np

import hd_field128_lift_screen as old


CONSTRUCTION = 'toeplitz-fullfield-h4-alternating-swapped-pairs-k20-rate-half-v1'
NEG = old.NEG


def compose(left, right):
    """Ordered log2 product; either operand may have leading batch axes."""
    return np.logaddexp2(left[..., :, 0, None]+right[..., None, 0, :],
                         left[..., :, 1, None]+right[..., None, 1, :])


def accumulate(result, singles, fulls, values):
    # at handles repeated histogram indices without losing positive terms.
    flat = values.reshape((-1, 4))
    for entry in range(4):
        np.logaddexp2.at(result[entry], (singles, fulls), flat[:, entry])


def pair_partition(singles, fulls, counts, logs, alpha, maximum, budget, max_v_full=2):
    """Partition good-good geometries directly; retain raw ordered V pairs."""
    if type(max_v_full) is not int or not 0 <= max_v_full <= 2:
        raise ValueError('This bounded provider supports V full-pair cutoffs0..2')
    v = fulls <= max_v_full
    weights = np.log2(counts)
    selected, remainder = old._empty(maximum), old._empty(maximum)
    pair_s, pair_f, pair_w, pair_logs = [], [], [], []
    for i in range(len(singles)):
        budget.check()
        ss, ff = singles[i]+singles, fulls[i]+fulls
        keep = ss+2*ff <= maximum
        both_v = keep & v[i] & v
        others = keep & ~both_v
        raw = compose(logs[i], logs)
        weighted = weights[i]+weights
        values = weighted[:, None, None]+alpha*raw
        accumulate(selected, ss[both_v], ff[both_v], values[both_v])
        accumulate(remainder, ss[others], ff[others], values[others])
        pair_s.extend(ss[both_v])
        pair_f.extend(ff[both_v])
        pair_w.extend(weighted[both_v])
        pair_logs.extend(raw[both_v])
    raw_pairs = (np.asarray(pair_s, dtype=int), np.asarray(pair_f, dtype=int),
                 np.asarray(pair_w), np.asarray(pair_logs).reshape((-1, 2, 2)))
    return selected, remainder, raw_pairs


def four_from_pairs(raw_pairs, alpha, maximum, budget):
    """C_VVVV from ordered raw pair products, before the fractional power."""
    singles, fulls, weights, logs = raw_pairs
    result = old._empty(maximum)
    for i in range(len(singles)):
        budget.check()
        ss, ff = singles[i]+singles, fulls[i]+fulls
        keep = ss+2*ff <= maximum
        raw = compose(logs[i], logs[keep])
        values = weights[i]+weights[keep, None, None]+alpha*raw
        accumulate(result, ss[keep], ff[keep], values)
    return result


def positive_parts(z, alpha, maximum, q, budget, max_v_full=2):
    """Return C_VV, R, C_VVVV without any coefficient subtraction."""
    if type(max_v_full) is not int or not 0 <= max_v_full <= 2:
        raise ValueError('This bounded provider supports V full-pair cutoffs0..2')
    single = old.single_coefficients(z, alpha, maximum, q, budget)
    good, dirty = single.copy(), single.copy()
    good[:, 2:, :], dirty[:, :2, :] = NEG, NEG
    s, f, counts, indices = old.good_geometry()
    keep = s+2*f <= maximum
    s, f, counts, indices = s[keep], f[keep], counts[keep], indices[keep]
    logs = old.local_logs(z, q, old.OCC[indices])
    selected, remainder, raw_pairs = pair_partition(s, f, counts, logs, alpha, maximum, budget, max_v_full)
    for left, right in ((good, dirty), (dirty, good), (dirty, dirty)):
        remainder = np.logaddexp2(remainder, old.polynomial_product(left, right, maximum, budget))
    return selected, remainder, four_from_pairs(raw_pairs, alpha, maximum, budget)


def four_operator(selected, remainder, composed_four, maximum, budget):
    """Positive evaluation of C4+C*R+R*C+R*R using two convolutions."""
    two = np.logaddexp2(selected, remainder)
    four = np.logaddexp2(composed_four, old.polynomial_product(selected, remainder, maximum, budget))
    four = np.logaddexp2(four, old.polynomial_product(remainder, two, maximum, budget))
    return four, two


def regional_screen(z, alpha, maximum, *, q=old.FIELD_ORDER, max_v_full=2,
                    max_seconds=60., _budget=None):
    old.checked_q(q, construction=True)
    if (type(maximum) is not int or not 1 <= maximum <= old.MAXIMUM
            or not isfinite(z) or not isfinite(alpha) or not 0 < z < 1 or not 0 < alpha <= 1
            or type(max_v_full) is not int or not 0 <= max_v_full <= 2):
        raise ValueError('Expected bounded cutoff1..64, V cutoff0..2, and finite valid marks')
    budget = _budget or old.Budget(max_seconds)
    if alpha == 1:
        # Alpha1 has no fractional-path loss and the frozen single power is equal.
        return old.regional_screen(z, alpha, maximum, q=q, paired=False, _budget=budget)
    selected, remainder, four = positive_parts(z, alpha, maximum, q, budget, max_v_full)
    operator, _ = four_operator(selected, remainder, four, maximum, budget)
    run = old.polynomial_power(operator, 32, maximum, budget)
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
    """All profiles at r, or one explicitly scoped term when profile is supplied."""
    bits = old.checked_q(q, construction=True)
    distance = F(distance)
    if not 0 < distance < F(1, 2):
        raise ValueError('Expected target distance strictly between0and1/2')
    if profile is not None and (len(profile) != 2 or any(type(x) is not int or not 0 <= x <= r//2 for x in profile)):
        raise ValueError('Invalid selected matching profile')
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
                composition_steps=4, pair_swaps=True, max_v_full_pairs=max_v_full,
                positive_partition='C_VVVV + C_VV*R + R*(C_VV+R)',
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
