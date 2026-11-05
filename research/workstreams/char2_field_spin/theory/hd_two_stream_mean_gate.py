"""Matched one/two-stream whole-a8/F1 mean gate; no parameter search.

The joint provider must extract both streams' TOTAL occupancy coefficient
before applying the single regional normalizer. A common product potential
then bounds each chronological matrix word. ONE prescribed scalar secant
gives a positive moment measure; its total mass is charged explicitly.
Full fixed-A Holder mask weights and the complete positive word tail are
unchanged. This FLOAT diagnostic is not a distance certificate.
"""
from fractions import Fraction
from hashlib import sha256
import json
from math import log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_two_sum_gf16_mean_column import column_polynomial, raw_holder_weights
from hd_two_sum_late_alpha_live import load_inputs
from hd_two_sum_late_alpha_quiet import log_power
from hd_two_sum_matrix_moment import load_saved as load_one_stream_reference
from hd_xor_multigroup import LIMIT


DIRECTORY = Path(__file__).resolve().parent
GROUPS, ARITY, COLUMNS, FLOOR = 18, 8, 16, 1296
MARK, LOG_P = Fraction(691, 700), -64.
SOURCE_PINS = {
    'hd_two_stream_regions.py': 'd4a562cd4497d1d64747eb53b6a084d22b5b8d75c93e8ba944c13435e6d59d0d',
    'hd_two_sum_gf16_mean_column.py': 'fa6377002fa52fbfcb980ddd4b72d30df16b6c12c7c6db24eed54554e61a1593',
    'hd_two_sum_closed_kernel.py': '2ee0fe57de750a19466155317d595589695252de450518b436eab6c5c0fbf979',
}


def _logs(values):
    values = np.asarray(values, dtype=float)
    if np.any(np.isnan(values)) or np.any(np.isposinf(values)):
        raise ValueError('Finite-or-zero log entries required')
    return values


def product_potential(streams, log_p=LOG_P, *, lumped=False):
    if (streams not in (1, 2) or not np.isfinite(log_p) or log_p > 0
            or (lumped and streams != 2)):
        raise ValueError('One/two streams and one common potential 0<p<=1 required')
    if streams == 1:
        return np.array([0., log_p])
    return np.array([0., log_p, 2*log_p] if lumped else [0., log_p, log_p, 2*log_p])


def potential_curve(bank, maximum, log_v, *, full_max=1):
    """Maximize the weighted row sum over actual f, not each entry separately."""
    log_v = _logs(log_v)
    if (log_v.ndim != 1 or np.any(~np.isfinite(log_v)) or log_v[0] != 0.
            or type(maximum) is not int or maximum < 1 or type(full_max) is not int
            or full_max < 0):
        raise ValueError('Finite positive state potential and complete occupancy range required')
    expected = {(m, f) for m in range(maximum+1)
                for f in range(min(full_max, m//2)+1)}
    if set(bank) != expected:
        raise ValueError('Every requested (m,f), and no other histogram, is required')
    curve, contacts = [], []
    for m in range(maximum+1):
        rows = []
        for f in range(min(full_max, m//2)+1):
            matrix = _logs(bank[m, f])
            if matrix.shape != (len(log_v), len(log_v)):
                raise ValueError('Provider state dimension does not match its potential')
            contractions = np.logaddexp2.reduce(matrix+log_v[None, :], axis=1)-log_v
            rows.extend((float(value), f, state) for state, value in enumerate(contractions))
        value = max(row[0] for row in rows)
        curve.append(value)
        contacts.append([dict(full_pairs=f, initial_state=state)
                         for candidate, f, state in rows if abs(candidate-value) <= 1e-10])
    if not np.isfinite(curve).all():
        raise ValueError('Every occupancy needs a positive finite contraction')
    # The empty input is diagonal, fixes zero, and contracts live mass.
    # Do not silently round a computed positive empty contraction to one.
    if curve[0] > 0.:
        raise ArithmeticError('The provider did not preserve the exact empty contraction <=1')
    curve[0] = 0.
    return np.asarray(curve), contacts


def secant_measure(curve, *, endpoints=(2, 18)):
    """One finite positive measure, with total mass explicitly retained.

    If eta<=1, mu=(1-eta)delta0+eta deltaB. Otherwise mu=eta deltaB.
    The latter deliberately relaxes the zero moment to eta>1. No alternate
    slope is fitted if the prescribed secant produces this weaker case.
    """
    curve = _logs(curve)
    left, right = endpoints
    if (curve.ndim != 1 or np.any(~np.isfinite(curve)) or curve[0] != 0.
            or not 1 <= left < right < len(curve)):
        raise ValueError('Finite scalar curve c0=1 and one valid secant required')
    slope = float((curve[right]-curve[left])/(right-left))
    degrees = np.arange(len(curve))
    intercept = float(np.max(curve[1:]-slope*degrees[1:]))
    mass = max(0., intercept)
    log_probability = intercept-mass
    envelope = intercept+slope*degrees
    envelope[0] = mass
    if np.any(envelope < curve-1e-10):
        raise ArithmeticError('Prescribed positive measure misses an occupancy')
    return dict(log2_atom=slope, log2_nonzero_weight=intercept,
        log2_total_mass=mass, log2_probability=log_probability,
        secant_endpoints=list(endpoints), contacts=np.flatnonzero(envelope-curve < 1e-9).tolist(),
        minimum_slack=float(np.min(envelope-curve)),
        zero_moment_relaxed=bool(mass > 0), exact_float_parameters={
            'atom_log': str(Fraction(slope)), 'nonzero_weight_log': str(Fraction(intercept))})


def contraction_screen(curve, streams, log_source_count, *, weights=None,
                       log_p=LOG_P, check_time=lambda: None):
    """All outer masks/states in the declared homogeneous r18/a8/F1 class."""
    log_v = product_potential(streams, log_p)
    if not np.isfinite(log_source_count):
        raise ValueError('One finite exact-source-count logarithm required')
    weights = raw_holder_weights() if weights is None else _logs(weights)
    if weights.shape != (17,):
        raise ValueError('The complete sixteen-region raw arity-eight table is required')
    measure = secant_measure(curve)
    # Normalized scalar law is used by this existing exact positive mixture.
    # Its removed total mass occurs once for each of the256 regional factors.
    column = column_polynomial(weights, GROUPS, measure['log2_probability'],
        measure['log2_atom'], arity=ARITY, check_time=check_time)
    word = log_power(column, COLUMNS, check_time)
    tail = float(np.logaddexp2.reduce(word[FLOOR:]))
    endpoint = float(log_v[0]-np.min(log_v))
    mass = COLUMNS*COLUMNS*measure['log2_total_mass']
    prefix = log_source_count-LIMIT*log2(float(MARK))-log2(255)-8*ARITY*GROUPS*log2(255)
    return dict(streams=streams, measure=measure, raw_word_tail_log2=tail,
        endpoint_factor_log2=endpoint, measure_mass_log2=mass,
        prefix_log2=prefix, log2_screen=prefix+tail+endpoint+mass,
        log2_potential=log_v.tolist(), peak_tail_degree=FLOOR+int(np.argmax(word[FLOOR:])))


def run(*, max_seconds=30, progress=None):
    """The provider import is lazy: tiny tests never build a regional bank."""
    if type(max_seconds) not in (int, float) or not 0 < max_seconds <= 30:
        raise ValueError('One explicit budget at most30 seconds required')
    # Provider module is supplied by the independent geometry implementation.
    from hd_two_stream_regions import raw_region_controls
    for name, expected in SOURCE_PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != expected:
            raise ValueError('Reviewed dependency changed: '+name)
    _, _, _, source_count = load_inputs()
    _, reference_curves = load_one_stream_reference()
    started = perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Matched two-stream control exhausted its budget')
    result = dict(active_groups=GROUPS, arity=ARITY, source_full_pair_max=1,
        source_count=str(source_count), mark=str(MARK), alpha=1,
        scope='complete homogeneous r18/arity8/Fmax1 class; proposed two-stream inner',
        existing_one_stream_control=True, paired_route_normalizers_per_region=1,
        regional_alpha=1, regional_potential_log2=LOG_P, source_count_once=True,
        rowsets_fixed_across_columns=True, raw_holder_exponent=16,
        extra_rowset_factor=False, rs_normalizer_exponent=1152,
        aggregate_surrogate_floor=FLOOR, surrogate='min(output_weight,8)',
        all_outer_moments=True, all_boundary_state_paths=True,
        potential_applied_before_mask_sum=True, terminal_factor_once=True,
        regional_masks_assumed_independent=False, parameter_fits=0,
        scalar_secant=[2, 18], raw_packet_banks_built=0,
        directed=False, whole_code_certificate=False,
        reviewed_source_sha256=SOURCE_PINS,
        compute_budget_seconds=max_seconds, complete_requested_controls=False, controls=[])
    try:
        check_time()
        banks = raw_region_controls(float(MARK), GROUPS, full_max=1)
        result['raw_packet_banks_built'] = 1
        result['provider_metadata'] = banks['metadata']
        reproduced = np.stack([np.maximum.reduce([banks['one_stream'][m, f]
            for f in range(min(1, m//2)+1)]).ravel() for m in range(GROUPS+1)]).T
        if (not np.array_equal(np.isfinite(reproduced), np.isfinite(reference_curves))
                or np.max(np.abs(reproduced[np.isfinite(reproduced)]
                                 -reference_curves[np.isfinite(reference_curves)])) > 1e-8):
            raise ArithmeticError('Matched one-stream provider does not reproduce the frozen full F1 bank')
        result['one_stream_reference_maximum_log2_difference'] = float(np.max(np.abs(
            reproduced[np.isfinite(reproduced)]-reference_curves[np.isfinite(reference_curves)])))
        check_time()
        weights = raw_holder_weights()
        for streams, key in ((1, 'one_stream'), (2, 'two_stream')):
            potential = product_potential(streams)
            curve, contacts = potential_curve(banks[key], GROUPS, potential)
            if streams == 2:
                lumped, _ = potential_curve(banks['two_stream_lumped'], GROUPS,
                                            product_potential(2, lumped=True))
                if not np.allclose(curve, lumped, atol=1e-10, rtol=0):
                    raise ArithmeticError('Three-state quotient did not preserve the common potential')
            row = contraction_screen(curve, streams, log2(source_count),
                                     weights=weights, check_time=check_time)
            row.update(regional_contraction_log2=curve.tolist(), regional_contacts=contacts,
                state_order=['0', '1'] if streams == 1 else ['00', '01', '10', '11'])
            result['controls'].append(row)
            if progress:
                progress(dict(streams=streams, log2_screen=row['log2_screen'], measure=row['measure']))
        result['two_minus_one_log2'] = result['controls'][1]['log2_screen']-result['controls'][0]['log2_screen']
        check_time()
        result['complete_requested_controls'] = True
    except TimeoutError as error:
        result['incomplete_reason'] = str(error)
    result['elapsed_seconds'] = perf_counter()-started
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--max-seconds', type=float, default=30)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a frozen diagnostic')
    answer = run(max_seconds=args.max_seconds, progress=lambda row: print(json.dumps(row), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(answer, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(dict(complete=answer['complete_requested_controls'], elapsed_seconds=answer['elapsed_seconds'])))
