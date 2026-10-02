"""Diagnostic deletion of paths from one fixed physical transfer envelope.

The retained recurrence emits x+A*a, then sets a_next=M*a+C*x. Fix the
prepared maps, tilt, and complete physical envelopes before any deletion.
Coordinate zero represents zero state; every other coordinate represents
live mass or a live-mass envelope. Birth-row selection is performed once.

Each variant changes only the stated physical entries. Macro composition,
regional placement, and chronological composition across regions are then
rebuilt. The initial coordinate remains zero and the terminal vector is all
ones. Every result retains C(L,q)*beta**q*exp(tilt*cutoff).

These are path sums of modified ENVELOPES, not probabilities for modified
codes and not upper or lower bounds on the original code. In particular,
deleting live returns does not measure their true probability. The scalar
zero_only variant projects the same envelope onto its zero coordinate; it
does not silently substitute an exact physical law. Even baseline is emitted
under a diagnostic-only schema, never as a certificate component.

Arb arithmetic avoids floating-point underflow. A zero result is reported
only when the modified transfer evaluates to exact zero. Ratios compare
saved diagnostic endpoints, not probabilities.
"""
from __future__ import annotations

from dataclasses import asdict
from fractions import Fraction as Q
from math import comb
from pathlib import Path
from time import monotonic

from flint import arb, arb_mat, ctx
import packet_rs_state_ladder as ladder
import packet_regional_power as polynomial
import packet_uniform_tail as uniform

q1, kernel = ladder.q1, ladder.q1.kernel_t64
SCHEMA = 'packet-physical-return-ablation-diagnostic-1'
DEFAULT_VARIANTS = ('baseline', 'no_occupied_zero_stay', 'no_pair_zero_stay',
                    'no_live_return', 'zero_only')
DESCRIPTIONS = {
    'baseline': 'unchanged fixed physical envelope',
    'no_occupied_zero_stay': 'set physical L_j[0,0]=0 for every j>0',
    'no_pair_zero_stay': 'set only physical L_2[0,0]=0',
    'no_live_return': 'set physical L_j[i,0]=0 for every i>0 and every j',
    'zero_only': 'project every physical L_j onto its scalar [0,0] entry',
}


def _variants(variants):
    if isinstance(variants, (str, bytes)):
        raise ValueError('explicit distinct variant sequence required')
    variants = tuple(variants)
    if (not variants or variants[0] != 'baseline' or len(set(variants)) != len(variants)
            or any(name not in DESCRIPTIONS for name in variants)):
        raise ValueError('baseline first, followed by distinct recognized variants required')
    return variants


def alter_physical(operators, variant):
    """Copy and alter physical entries; never mutate or reselect source rows."""
    if variant not in DESCRIPTIONS:
        raise ValueError('recognized diagnostic variant required')
    operators = tuple(operators)
    if not 2 <= len(operators) <= 17:
        raise ValueError('complete physical occupancy family with 1..16 slots required')
    n = operators[0].nrows()
    if n < 1 or any(m.nrows() != n or m.ncols() != n for m in operators):
        raise ValueError('matching nonempty square physical matrices required')
    if any(not m[i, k].is_finite() or not m[i, k] >= 0
           for m in operators for i in range(n) for k in range(n)):
        raise ValueError('finite nonnegative physical envelopes required')
    result = []
    for j, original in enumerate(operators):
        if variant == 'zero_only':
            result.append(arb_mat([[original[0, 0]]]))
            continue
        value = arb_mat(original)
        if ((variant == 'no_occupied_zero_stay' and j > 0)
                or (variant == 'no_pair_zero_stay' and j == 2)):
            value[0, 0] = 0
        if variant == 'no_live_return':
            for i in range(1, n):
                value[i, 0] = 0
        result.append(value)
    return result


def _endpoint(value):
    value = kernel.up(value)
    if not value.is_finite() or not value >= 0:
        raise ArithmeticError('finite nonnegative diagnostic value required')
    return [int(v) for v in value.man_exp()]


def _arb_endpoint(pair):
    return arb(pair[0]) * arb(2) ** pair[1]


def evaluate_physical(operators, *, group_count, regions, occupancy, tilt,
                      beta, threshold, variants=DEFAULT_VARIANTS):
    """Compare toy or real physical envelopes, without authenticating maps.

    This lower-level interface is for tests and diagnostic composition only.
    It returns JSON-ready noncertificate values, not proof components.
    """
    variants = _variants(variants)
    operators = tuple(operators)
    # Validate the frozen source even when only a subset of variants is used.
    source = alter_physical(operators, 'baseline')
    windows = len(source)-1
    if (type(group_count) is not int or group_count < 2*windows
            or group_count % (2*windows) or type(regions) is not int or regions < 1
            or type(occupancy) is not int or not 0 <= occupancy <= group_count
            or type(threshold) is not int or not 0 <= threshold <= 4*group_count*regions):
        raise ValueError('integral complete macro geometry, occupancy, and cutoff required')
    tilt, beta = Q(tilt), Q(beta)
    if tilt <= 0 or beta <= 0:
        raise ValueError('positive exact tilt and beta required')
    prefactor = (arb(comb(group_count, occupancy)) * kernel.aq(beta)**occupancy *
                 (kernel.aq(tilt)*threshold).exp())
    results = {}
    for variant in variants:
        start = monotonic()
        physical = alter_physical(source, variant)
        macro = kernel.convolve(physical)
        regional = polynomial.placement_power(macro, epochs=group_count//(2*windows),
            windows=2*windows, rounding=q1.rounded, maximum_groups=occupancy)
        mixed = uniform.regional_uniform(regional, occupancy)
        complete = mixed**regions
        moment = sum((complete[0, j] for j in range(complete.ncols())), arb(0))
        endpoint = _endpoint(prefactor*moment)
        value = _arb_endpoint(endpoint)
        results[variant] = dict(description=DESCRIPTIONS[variant],
            diagnostic_mass_dyadic=endpoint, transfer_moment_dyadic=_endpoint(moment),
            negative_log2_diagnostic_mass=None if value == 0 else str(-value.log()/arb(2).log()),
            exact_zero_in_modified_transfer=bool(value == 0),
            transfer_dimension=physical[0].nrows(), elapsed_seconds=monotonic()-start,
            is_probability_bound=False, whole_code_certificate=False)
    baseline = _arb_endpoint(results['baseline']['diagnostic_mass_dyadic'])
    for result in results.values():
        value = _arb_endpoint(result['diagnostic_mass_dyadic'])
        result['ratio_to_baseline_endpoint'] = None if baseline == 0 else str(value/baseline)
        result['log2_endpoint_reduction'] = (None if value == 0 or baseline == 0
                                            else str((baseline/value).log()/arb(2).log()))
    return results


def run(data, map_record, *, K, occupancy, tilt, precision=192, activity=Q(1, 2),
        output=None, variants=DEFAULT_VARIANTS):
    """Use one already prepared state-17..22 map at a single fixed witness.

    No census is performed here. Obtain both inputs from the same adapter's
    prepare() call. The caller owns any expensive preparation or numerical run.
    """
    if data is None or not isinstance(map_record, dict):
        raise ValueError('paired fresh prepared data and map record required')
    variants = _variants(variants)
    geometry = ladder.lengths.geometry(K)
    bits = map_record.get('s')
    adapter = ladder.adapter_for(bits)
    tilts, output = ladder.lengths.checked_options((tilt,), precision, output)
    if type(occupancy) is not int or not 3 <= occupancy <= geometry.group_count:
        raise ValueError('one tail occupancy in 3..L required')
    activity, tilt = Q(activity), Q(tilts[0])
    if not 0 <= activity <= 1:
        raise ValueError('birth-row selection activity in [0,1] required')
    if output is not None and Path(output).exists():
        raise ValueError('fresh diagnostic output path required')
    start = monotonic()
    data, map_record, before = ladder._prepare(adapter, bits, precision, data, map_record)
    physical_data = kernel.authenticate(data)
    z = (-kernel.aq(tilt)).exp()
    physical = kernel.sparse_kernel.outward_at_z(physical_data, z)
    physical = kernel.kernel_birth_density.refine_local(physical_data, physical, z, activity)
    sources = q1.source_snapshot()
    if any(sources.get(path) != digest for path, digest in before.items()):
        raise RuntimeError('source changed while constructing the fixed physical family')
    beta, counts = ladder.lengths.exact_outer()
    results = evaluate_physical(physical, group_count=geometry.group_count,
        regions=geometry.regions, occupancy=occupancy, tilt=tilt, beta=beta,
        threshold=geometry.N//10, variants=variants)
    adapter.authenticate(data, map_record)
    if ctx.prec != precision or q1.source_snapshot() != sources:
        raise RuntimeError('precision or sources changed during diagnostic computation')
    record = dict(schema=SCHEMA, diagnostic_only=True, is_probability_bound=False,
        whole_code_certificate=False, all_occupancies_covered=False,
        K=K, N=geometry.N, geometry=asdict(geometry), state_bits=bits,
        physical_t=64, macro_t=128, physical_steps=geometry.N//64,
        regions=geometry.regions, macros_per_region=geometry.macros_per_region,
        occupancy=occupancy, tilt=str(tilt), threshold=geometry.N//10, distance='1/10',
        beta=str(beta), count_sha256=ladder.lengths.count_hash(counts),
        multiplicity=f'C({geometry.group_count},{occupancy})', precision=precision,
        placement_backend='binary', birth_density='capped', activity=str(activity),
        row_selection='once before all deletions; no variant-specific reselection',
        zero_initial_state=True, final_flush=False,
        terminal='sum of all retained mass coordinates; zero_only is an exact coordinate projection',
        state_continuity='retained between every physical step and region',
        map_record=map_record, source_sha256=sources, variants=list(variants), results=results,
        elapsed_seconds=monotonic()-start,
        scope='Fixed-witness path sums of altered physical envelopes. Not probabilities, '
              'not bounds for the original or a modified code, and not certificate inputs.')
    if output is not None:
        ladder.arithmetic._write_new(output, record)
    return record
