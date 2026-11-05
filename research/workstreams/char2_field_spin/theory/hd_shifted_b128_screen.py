"""Isolated FLOAT gate for 64-input/128-output grouped Toeplitz outers.

Four RS[32,16] constituents have independent invertible Toeplitz4 column
mixers. The 16384 groups feed 128 alternating, independently swapped pair
regions. Each region contains 256 unchanged 64-input/16-state/h4 packets.
Its regional polynomial is T**128; the ordered regional pair repeats64
times. The outer density is 256**80/255**16. This is a new construction,
not a changed default or an inherited certificate for the B256 encoder.
"""
from fractions import Fraction as F
from copy import deepcopy
from math import comb, isfinite, lgamma, log, log2, sqrt
from pathlib import Path

import numpy as np

from hd_grouped_packet_screen import LIMIT, log_multinomial
from hd_randomized_outer_state_profile import multiply, power, total
from hd_shifted_swap_sparse import (swap_single_screen, good_pair_screen,
                                    product_screen, power_screen)
from hd_shifted_good_pair_screen import (good_step_coefficients, good_pair_coefficients,
                                        paired_frequency)
from hd_shifted_swap_screen import swap_coefficients
from hd_shifted_unoriented_screen import low_degree_frequency
from hd_shifted_potential_band_neighbor import source_witness
from hd_shifted_run_count_screen import run_log_screen_float


CONSTRUCTION = 'toeplitz-gf256-h4-alternating-swapped-pairs-b128-k20-rate-half'
GROUPS, REGIONS, PAIRS, PACKETS_PER_REGION = 16384, 128, 8192, 256
LOG_DENSITY = 640-16*log2(255)


def geometry():
    return dict(construction=CONSTRUCTION, source_dimension=1 << 20, output_length=1 << 21,
                input_symbols_per_group=64, output_symbols_per_group=128, outer_groups=GROUPS,
                regions=REGIONS, pairs_per_region=PAIRS, packets_per_region=PACKETS_PER_REGION,
                paired_blocks_per_region=PACKETS_PER_REGION//2, regional_pair_repetitions=REGIONS//2,
                outer_density_formula='256^80/255^16', field_order=256, forbidden_weight=LIMIT,
                pair_swaps=True, directed=False, whole_code_certificate=False)


def sparse_profile_counts(groups, weight):
    """Exact joint matching counts for all indexed weight-one/two supports."""
    if type(groups) is not int or groups < 4 or groups % 2 or weight not in (1, 2):
        raise ValueError('An even cycle of at least four groups and weight one or two are required')
    return {(0, 0): groups} if weight == 1 else {
        (0, 0): comb(groups, 2)-groups, (1, 0): groups//2, (0, 1): groups//2}


def sparse_regions(z, alpha, maximum, *, pairs=PAIRS, packet_steps=PACKETS_PER_REGION, paired=True):
    """Coefficient-exact truncation to weight<=2, with FLOAT arithmetic only.

    Explicit smaller dimensions support toy/regression checks. This helper
    never changes the existing sparse or regional proof constants.
    """
    if (maximum not in (1, 2) or not 0 < z < 1 or not 0 < alpha <= 1
            or type(packet_steps) is not int or packet_steps < 2 or packet_steps % 2
            or pairs != 32*packet_steps):
        raise ValueError('Valid weight-one/two witnesses and 32 pairs per packet are required')
    marker = maximum/pairs
    degree = np.arange(maximum+1)[:, None]+np.arange(maximum//2+1)
    single = swap_single_screen(z, alpha, maximum).copy()*marker**degree
    if paired:
        good, dirty = single.copy(), single.copy()
        good[:, 2:, :] = 0.
        dirty[:, :2, :] = 0.
        composed = good_pair_screen(z, alpha, maximum)*marker**degree
        base = (composed+product_screen(good, dirty, maximum)
                +product_screen(dirty, good, maximum)+product_screen(dirty, dirty, maximum))
        run = power_screen(base, packet_steps//2, maximum)
    else:
        run = power_screen(single, packet_steps, maximum)
    result = {}
    for weight in range(maximum+1):
        for full in range(weight//2+1):
            singleton = weight-2*full
            values = run[:, singleton, full]
            if np.any(values < 0) or not np.all(np.isfinite(values)):
                raise ArithmeticError('Invalid sparse coefficient')
            with np.errstate(divide='ignore'):
                logs = np.log2(values).reshape(2, 2)
            normalizer = (log_multinomial((pairs-singleton-full, singleton, full))
                          +(singleton+full)*log2(marker))
            result[weight, full] = tuple(map(tuple, logs-normalizer))
    return result


def sparse_weight_screen(weight, z, alpha=1.):
    """Sum every support at weight one/two, retaining both ordered matchings."""
    if weight not in (1, 2):
        raise ValueError('This bounded structural gate covers only weights one and two')
    regions = sparse_regions(z, alpha, weight)
    prefix = alpha*(weight*LOG_DENSITY-log2(255)-LIMIT*log2(z))
    profiles = []
    for (a, b), count in sparse_profile_counts(GROUPS, weight).items():
        ordered = power(multiply(regions[weight, a], regions[weight, b]), REGIONS//2)
        value = float(log2(count)+prefix+total(ordered[0]))
        profiles.append(dict(a=a, b=b, support_count=count, log2_upper_screen=value))
    value = float(total(row['log2_upper_screen'] for row in profiles))
    return dict(geometry(), active_groups=weight, z=str(F(z)), alpha=str(F(alpha)),
                log2_upper_screen=value, profiles=profiles, all_indexed_supports=True,
                coefficient_method='positive-weight-two-truncation', witness_fits=0,
                passes_40_bits=value <= -40.)


def sparse_gate():
    """Two predetermined density-scaled witnesses; no search or fitting."""
    old_single_z = .9991588324466113
    rows = [sparse_weight_screen(1, sqrt(old_single_z)),
            sparse_weight_screen(2, old_single_z)]
    return dict(geometry(), method='bounded-sparse-gate', results=rows,
                all_pass_40_bits=all(row['passes_40_bits'] for row in rows),
                witness_rule='old B256 singleton witness; halve log(z) for B128 weight one',
                new_fits=0)


def regional_alias_log_table(z, alpha, odds, *, grid=128, paired_blocks=128):
    """FLOAT alias proposal with an explicit regional two-step exponent.

    The 1e-12 normalized cushion is a numerical screening safeguard, not a
    directed error estimate. No bound in this module is a certificate.
    """
    if (grid not in (64, 128) or paired_blocks not in (64, 128)
            or not 0 < z < 1 or not 0 < alpha <= 1 or len(odds) != 2
            or not all(map(isfinite, odds))):
        raise ValueError('A bounded fixed-witness alias screen is required')
    dirty = swap_coefficients(z, alpha).copy()
    dirty[:, :2, :] = 0.
    dirty *= np.exp2(np.arange(33)[None, :, None]*odds[0]
                     +np.arange(33)[None, None, :]*odds[1])
    dirty = np.fft.fftn(dirty, s=(grid, grid), axes=(1, 2))
    good = low_degree_frequency(good_step_coefficients(z, alpha), odds, grid, (64, 64))
    composed = low_degree_frequency(good_pair_coefficients(z, alpha), odds, grid, (64, 64))
    values = paired_frequency(dirty, good, composed)
    peak = float(np.max(values[:, 0, 0].real))
    if not 0 < peak < float('inf') or not np.all(np.isfinite(values)):
        raise ArithmeticError('Invalid positive matrix evaluation')
    values /= peak
    a, b, c, d = values[:, 0, 0].real
    rho = (a+d+np.sqrt((a-d)**2+4*b*c))/2
    vector = np.array([b, rho-a] if rho-a >= rho-d else [rho-d, c])
    if not rho > 0 or not np.all(vector > 0) or not np.all(np.isfinite(vector)):
        raise ArithmeticError('Invalid positive similarity scaling')
    vector /= max(vector)
    values *= (vector[None, :]/vector[:, None]/rho).ravel()[:, None, None]
    for _ in range(paired_blocks.bit_length()-1):
        a, b, c, d = values
        trace, bc = a+d, b*c
        values = np.stack((a*a+bc, b*trace, c*trace, d*d+bc))
    aliases = np.fft.ifftn(values, axes=(1, 2))
    imaginary = float(np.max(abs(aliases.imag)))
    minimum = float(np.min(aliases.real))
    if imaginary > 1.e-10 or minimum < -1.e-10:
        raise ArithmeticError('FLOAT alias diagnostic failed')
    logs = np.log2(np.maximum(aliases.real, 0.)+1.e-12)
    logs += paired_blocks*(log2(peak)+log2(rho))
    logs += np.log2(vector[:, None]/vector[None, :]).ravel()[:, None, None]
    logs.setflags(write=False)
    return logs, dict(maximum_imaginary_alias=imaginary, minimum_real_alias=minimum,
                      normalized_cushion=1.e-12, paired_blocks=paired_blocks)


def medium_weight_screen(artifact, *, progress=None):
    """One all-profile weight2048 screen, reusing the fixed old1024 witness."""
    source, grid = source_witness(artifact)
    if source.get('active_groups') != 1024 or grid != 128:
        raise ValueError('This gate permits only the frozen old1024 bank at grid128')
    z, alpha, potential = (F(source[key]) for key in ('z', 'alpha', 'potential'))
    markers = tuple(dict.fromkeys(tuple(F(value) for value in row.split())
                                 for row in source['markers']))
    if (not 0 < z < 1 or not 0 < alpha <= 1 or potential <= 0
            or not markers or len(markers) > 100
            or any(len(row) != 2 or min(row) <= 0 for row in markers)):
        raise ValueError('Admissible exact frozen witness literals are required')
    z, alpha, ratio = float(z), float(alpha), log2(potential)
    weight = 2048
    full = np.arange(weight//2+1)
    singleton, empty = weight-2*full, PAIRS-weight+full
    factorial = np.array([lgamma(n+1)/log(2) for n in range(PAIRS+1)])
    multinomial = factorial[PAIRS]-factorial[empty]-factorial[singleton]-factorial[full]
    best = np.full((4, len(full)), np.inf)
    diagnostics = []
    for index, row in enumerate(markers):
        odds = tuple(log2(value) for value in row)
        table, diagnostic = regional_alias_log_table(z, alpha, odds, grid=grid,
                                                     paired_blocks=PACKETS_PER_REGION//2)
        normalized = table[:, singleton % grid, full % grid]-multinomial-singleton*odds[0]-full*odds[1]
        np.minimum(best, normalized, out=best)
        diagnostics.append(diagnostic)
        if progress:
            progress(dict(stage='b128-medium-atlas', completed=index+1, total=len(markers)))
    if not np.all(np.isfinite(best)):
        raise ArithmeticError('Invalid normalized regional screen')
    regions = {int(f): best[:, index].reshape(2, 2) for index, f in enumerate(full)}
    support_sum = run_log_screen_float(PAIRS, weight, regions, ratio, repetitions=REGIONS)
    prefix = alpha*(weight*LOG_DENSITY-log2(255)-LIMIT*log2(z))
    value = float(support_sum+prefix)
    return dict(geometry(), active_groups=weight, log2_upper_screen=value,
                log2_support_sum_screen=support_sum, log2_outer_prefix=prefix,
                regional_histograms=len(regions), covers_all_matching_profiles=True,
                covers_all_even_odd_splits=True, count_method='FLOAT-log-factorial-cyclic-run-counts',
                coefficient_method='FLOAT-fixed-alias-bank', grid=grid, atlases_built=len(markers),
                source_weight=1024, exact_frozen_witness=deepcopy(source),
                marker_fits=0, witness_fits=0, passes_40_bits=value <= -40.,
                diagnostics=diagnostics)


if __name__ == '__main__':
    import argparse, json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--medium-source', type=Path,
                        help='Run only the bounded2048 gate with this frozen old1024 proposal')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output exists; structural artifacts are never overwritten')
    if args.medium_source:
        from hashlib import sha256
        raw = args.medium_source.read_bytes()
        result = medium_weight_screen(json.loads(raw),
            progress=lambda row: print(json.dumps(row), flush=True))
        result['source_artifact_sha256'] = sha256(raw).hexdigest()
    else:
        result = sparse_gate()
    encoded = json.dumps(result, separators=(',', ':'), allow_nan=False)
    with args.output.open('x', encoding='utf-8') as target:
        target.write(encoded+'\n')
    summary = (dict(active_groups=result['active_groups'], log2_upper_screen=result['log2_upper_screen'],
                    passes_40_bits=result['passes_40_bits'], atlases_built=result['atlases_built'])
               if args.medium_source else dict(all_pass_40_bits=result['all_pass_40_bits'],
                    results=[{key: row[key] for key in ('active_groups', 'log2_upper_screen')}
                             for row in result['results']]))
    print(json.dumps(dict(path=str(args.output.resolve()), **summary, directed=False)), flush=True)
