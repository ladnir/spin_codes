"""Isolated FLOAT full-RS[width,state] sparse screens, default RS[128,64].

The fixed maps use GF256 evaluation points0..width-1 and interpolation
points128..127+state, with C=A^T. For power-of-two width and2state<=width,
the additive evaluation set gives CA=0. Independent nonzero full-field
diagonals surround an MDS Lagrange H_state. No principal-nullity table is
used. Frozen polynomial, route-count, and outer-density helpers are reused
through explicit calls; their modules and constants are never mutated.

Outer four-RS[64,32] groups, full4x4 Toeplitz randomizers,4096 pair slots
per region,256 regions, and persistent state are unchanged. A128-symbol
packet has64 pair slots, so each8192-symbol region has64 packets and32
composed packet pairs. This is not an encoder or a distance certificate.
"""

import ast
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from math import comb, isfinite, log2
from pathlib import Path

import numpy as np

import hd_field128_full_rs_state_screen as frozen


old = frozen.old
FIELD_ORDER, NEG = old.FIELD_ORDER, old.NEG
WIDTHS = (64, 128)
DEFAULT_WIDTH, DEFAULT_STATE = 128, 64


def parameters(q, n, state):
    bits = old.checked_q(q)
    if (type(n) is not int or type(state) is not int or not 1 <= state <= 64
            or not 2*state <= n <= 128 or n&(n-1) or q < n):
        raise ValueError('Use power-of-two width<=128, state<=width/2, and q>=width')
    return bits


def construction(width=DEFAULT_WIDTH, state=DEFAULT_STATE):
    parameters(256, width, state)
    if width not in WIDTHS:
        raise ValueError('Public packet widths are64and128')
    return f'toeplitz-fullfield-full-rs{width}-{state}-alternating-swapped-pairs-k20-rate-half-v1'


def routing(width=DEFAULT_WIDTH):
    if type(width) is not int or width not in WIDTHS:
        raise ValueError('Public packet widths are64and128')
    symbols = 2*old.PAIRS
    return dict(pair_slots_per_packet=width//2, pair_slots_per_region=old.PAIRS,
                symbols_per_region=symbols, packets_per_region=symbols//width,
                composed_pairs_per_region=symbols//(2*width), regions=old.REGIONS)


@lru_cache(128)
def fixed_maps(width=DEFAULT_WIDTH, state=DEFAULT_STATE):
    construction(width, state)
    a = frozen.lagrange(range(width), range(128, 128+state))
    refresh = frozen.lagrange(range(state, 2*state), range(state))
    return a, tuple(zip(*a)), refresh


@lru_cache(256)
def exact_moments(j, z, q=FIELD_ORDER, n=DEFAULT_WIDTH, state=DEFAULT_STATE, joint=True):
    """Positive rational quiet, birth, return cap, and uniform-live moments."""
    parameters(q, n, state)
    if isinstance(z, float):
        raise ValueError('Exact local oracle requires a rational mark')
    z = F(z)
    if type(j) is not int or not 0 <= j <= n or not 0 < z < 1:
        raise ValueError('Invalid occupancy or exact mark')
    w, t, rho = 1-z, z+(1-z)/q, min(state, j)
    quiet = sum((F(comb(j, c), q**(c+min(state, j-c)))*z**(j-c)*w**c
                 for c in range(j+1)), F(0))
    birth = sum((F(comb(j, c), q**c)*z**(j-c)*w**c*(1-F(1, q**min(state, j-c)))
                 for c in range(j)), F(0))
    inactive = n-j
    excess = sum((F(comb(inactive, e), q**e)*z**(inactive-e)*w**e*(1-F(1, q**(state-e)))
                  for e in range(min(state, inactive+1))), F(0))
    live = t**j*excess/(1-F(1, q**state))
    ret = min(birth, (1-F(1, q**rho))*z**max(0, n-state+1-j), (q**state-1)*live)
    if joint:
        bound = F(0)
        for c in range(j+1):
            for e in range(min(state, inactive+1)):
                probability = min(1-F(1, q**rho),
                                  F(q**min(rho, state-e)-1, q**(j-max(0, j-c-state))))
                bound += comb(j, c)*comb(inactive, e)*z**(n-c-e)*w**(c+e)*probability
        ret = min(ret, bound)
    return quiet, birth, ret, live


def exact_packet(j, z, q=FIELD_ORDER, n=DEFAULT_WIDTH, state=DEFAULT_STATE, joint=True):
    quiet, birth, ret, live = exact_moments(j, z, q, n, state, joint)
    d = q**state-1
    density = F(d, (q-1)**state)
    return quiet, density*birth, ret/d, min(density*live, F(z)**(n-state+1)) if not j else density*live


@lru_cache(64)
def local_logs(z, q=FIELD_ORDER, n=DEFAULT_WIDTH, state=DEFAULT_STATE, joint=True):
    """Positive expansions in log2 form; no q**state binary64 conversion."""
    bits = parameters(q, n, state)
    if not isfinite(z) or not 0 < z < 1:
        raise ValueError('Expected finite0<z<1')
    lz, lw = log2(z), log2(1-z)
    tilt = float(np.logaddexp2(lz, lw-bits))
    nonzero = old.log1mexp2(-state*bits)
    log_d, density = state*bits+nonzero, nonzero-state*old.log1mexp2(-bits)
    quiet, birth, excess = [], [], []
    for j in range(n+1):
        quiet.append(float(np.logaddexp2.reduce([
            log2(comb(j, c))+(j-c)*lz+c*lw-(c+min(state, j-c))*bits for c in range(j+1)])))
        birth.append(float(np.logaddexp2.reduce([
            log2(comb(j, c))+(j-c)*lz+c*lw-c*bits+old.log1mexp2(-min(state, j-c)*bits)
            for c in range(j)])) if j else NEG)
        excess.append(float(np.logaddexp2.reduce([
            log2(comb(j, e))+(j-e)*lz+e*lw-e*bits+old.log1mexp2(-(state-e)*bits)
            for e in range(min(state, j+1))])))
    result = []
    for j in range(n+1):
        rho, inactive = min(state, j), n-j
        live = j*tilt+excess[inactive]-nonzero
        ret = min(birth[j], old.log1mexp2(-rho*bits)+max(0, n-state+1-j)*lz, log_d+live)
        if joint and j:
            terms = []
            for c in range(j+1):
                for e in range(min(state, inactive+1)):
                    dimension = min(rho, state-e)
                    count = (dimension+max(0, j-c-state)-j)*bits+old.log1mexp2(-dimension*bits)
                    probability = min(old.log1mexp2(-rho*bits), count)
                    terms.append(log2(comb(j, c)*comb(inactive, e))+(n-c-e)*lz+(c+e)*lw+probability)
            ret = min(ret, float(np.logaddexp2.reduce(terms)))
        live += density
        if not j:
            live = min(live, (n-state+1)*lz)
        result.append((quiet[j], density+birth[j], ret-log_d, live))
    entries = np.array(result).reshape((n+1, 2, 2))
    if np.any(np.isnan(entries)) or np.any(np.isposinf(entries)):
        raise ArithmeticError('Invalid width/state local logarithm')
    entries.setflags(write=False)
    return entries


@lru_cache(256)
def geometry(maximum, *, width=DEFAULT_WIDTH, good=False):
    """Fair orientation swaps cancel2**single; occupancy is single+2full."""
    pairs = routing(width)['pair_slots_per_packet']
    if type(maximum) is not int or not 0 <= maximum <= 2*width:
        raise ValueError('Expected a nonnegative bounded geometry cutoff')
    return tuple((s, f, comb(pairs, f)*comb(pairs-f, s))
                 for f in range(min(pairs, maximum//2)+1)
                 for s in range(min(1 if good else pairs, pairs-f, maximum-2*f)+1))


def single_coefficients(z, alpha, maximum, q, width, state, budget, joint=True):
    budget.check()
    logs = local_logs(z, q, width, state, joint).reshape((width+1, 4))
    result = old._empty(maximum)
    for s, f, count in geometry(maximum, width=width):
        budget.check()
        result[:, s, f] = log2(count)+alpha*logs[s+2*f]
    return result


def good_pair_coefficients(z, alpha, maximum, q, width, state, budget, joint=True):
    logs = local_logs(z, q, width, state, joint)
    good = geometry(maximum, width=width, good=True)
    result = old._empty(maximum)
    for s, f, count in good:
        budget.check()
        left = logs[s+2*f]
        for t, g, other in good:
            if s+t+2*(f+g) > maximum:
                continue
            right = logs[t+2*g]
            composed = np.logaddexp2(left[:, 0, None]+right[None, 0, :], left[:, 1, None]+right[None, 1, :])
            result[:, s+t, f+g] = np.logaddexp2(result[:, s+t, f+g], log2(count*other)+alpha*composed.ravel())
    return result


def regional_screen(z, alpha, maximum, *, q=FIELD_ORDER, width=DEFAULT_WIDTH, state=DEFAULT_STATE,
                    paired=True, joint=True, max_seconds=60., _budget=None):
    construction(width, state)
    old.checked_q(q, construction=True)
    if (type(maximum) is not int or not 1 <= maximum <= 64 or not isfinite(z) or not isfinite(alpha)
            or not 0 < z < 1 or not 0 < alpha <= 1):
        raise ValueError('Expected coefficient cutoff1..64 and finite valid marks')
    budget = _budget or old.Budget(max_seconds)
    route = routing(width)
    single = single_coefficients(z, alpha, maximum, q, width, state, budget, joint)
    if paired and alpha != 1:
        good, dirty = single.copy(), single.copy()
        good[:, 2:, :], dirty[:, :2, :] = NEG, NEG
        operator = good_pair_coefficients(z, alpha, maximum, q, width, state, budget, joint)
        for left, right in ((good, dirty), (dirty, good), (dirty, dirty)):
            operator = np.logaddexp2(operator, old.polynomial_product(left, right, maximum, budget))
        run = old.polynomial_power(operator, route['composed_pairs_per_region'], maximum, budget)
    else:
        run = old.polynomial_power(single, route['packets_per_region'], maximum, budget)
    regions = {}
    for r in range(maximum+1):
        for f in range(r//2+1):
            s = r-2*f
            values = run[:, s, f]-old.log_multinomial((old.PAIRS-s-f, s, f))
            regions[r, f] = tuple(map(tuple, values.reshape((2, 2))))
    budget.check()
    return regions


def metadata(width=DEFAULT_WIDTH, state=DEFAULT_STATE):
    name = construction(width, state)
    route = routing(width)
    return dict(construction=name, state_dimension=state, inner_width=width,
                emission_minimum=width-state+1, inner_emission=f'single-RS[{width},{state}]',
                feedback='A transpose', CA_zero=True, fixed_coefficient_subfield=256,
                fixed_nodes=f'GF256 additive subspace0..{width-1}', state_nodes=f'128..{127+state}',
                refresh=f'D2 H_{state} D1; H fixed MDS Lagrange from0..{state-1} to{state}..{2*state-1}',
                refresh_diagonal_fullfield_products_per_packet=2*state,
                refresh_diagonal_fullfield_products_per_codeword=2*state*(old.OUTPUT_LENGTH//width),
                refresh_diagonal_key_field_elements_per_codeword=2*state*(old.OUTPUT_LENGTH//width),
                outer_groups=old.GROUPS, source_dimension=old.DIMENSION, output_length=old.OUTPUT_LENGTH,
                outer_emission='four-RS[64,32];128-source to256-output symbols',
                outer_randomizer='independent full-field invertible4x4 Toeplitz per outer group and output column',
                outer_randomizers_per_group=64,
                pair_swaps=True, state_persists_across_regions=True, full_field_random_coefficients=True,
                **route)


def weight_screen(r, z, alpha, *, q=FIELD_ORDER, width=DEFAULT_WIDTH, state=DEFAULT_STATE,
                  distance=F(12, 25), paired=True, joint=True, outer_measure='invertible-toeplitz',
                  max_seconds=60.):
    details, bits = metadata(width, state), old.checked_q(q, construction=True)
    distance = F(distance)
    if not 0 < distance < F(1, 2):
        raise ValueError('Expected distance strictly between0and1/2')
    budget = old.Budget(max_seconds)
    regions = regional_screen(z, alpha, r, q=q, width=width, state=state,
                              paired=paired, joint=joint, _budget=budget)
    forbidden = int(distance*old.OUTPUT_LENGTH)
    prefix = alpha*(r*old.outer_log_density(q, outer_measure)-(bits+old.log1mexp2(-bits))-forbidden*log2(z))
    profiles = []
    for a in range(r//2+1):
        budget.check()
        for b in range(r//2+1):
            count = old.profile_count_upper(old.PAIRS, r, a, b)
            if count:
                run = old.power(old.multiply(regions[r, a], regions[r, b]), old.REGIONS//2)
                profiles.append((float(log2(count)+prefix+old.total(run[0])), a, b))
    value = float(old.total(row[0] for row in profiles))
    if not isfinite(value):
        raise ArithmeticError('Nonfinite selected-weight FLOAT result')
    budget.check()
    return dict(**details, active_groups=r, z=z, alpha=alpha, field_order=q, field_bits=bits,
                target_distance=str(distance), forbidden_weight=forbidden, outer_measure=outer_measure,
                composition_steps=2 if paired else 1, joint_rank_return_cap=joint,
                log2_upper_screen=value, worst_profile=max(profiles), profile_count=len(profiles),
                covers_all_supports_at_selected_weight=True, covers_all_even_odd_splits=True,
                directed=False, whole_code_certificate=False, arithmetic='FLOAT',
                scope='new width/state construction; selected weight only, no encoder or timing claim')


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
    parser.add_argument('--width', type=int, choices=WIDTHS, default=DEFAULT_WIDTH)
    parser.add_argument('--state', type=int, default=DEFAULT_STATE)
    parser.add_argument('--distance', default='12/25')
    parser.add_argument('--field-bits', type=int, default=128)
    parser.add_argument('--simple-return', action='store_true')
    parser.add_argument('--unpaired', action='store_true')
    parser.add_argument('--max-seconds', type=float, default=60.)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.output and args.output.exists():
        parser.error('Output must be a new file')
    if not 8 <= args.field_bits <= 128 or args.field_bits%8:
        parser.error('Field bits must be a multiple of8in8..128')
    result = weight_screen(args.r, args.z, args.alpha, q=1 << args.field_bits,
                           width=args.width, state=args.state, distance=F(args.distance),
                           paired=not args.unpaired, joint=not args.simple_return, max_seconds=args.max_seconds)
    result['source_sha256'] = source_pins()
    if args.output:
        with args.output.open('x', encoding='utf-8') as stream:
            json.dump(result, stream, indent=2, allow_nan=False)
            stream.write('\n')
    print(json.dumps(result, allow_nan=False))


if __name__ == '__main__':
    main()
