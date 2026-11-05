"""One FLOAT regional diagnostic for a common-source good-route event.

No regional provider is rebuilt. Exact route probabilities and an exact
two-state good-route cap are compared with one externally pinned saved R1
curve. This is not a column bound, a directed receipt, or a distance proof.
The good event is common to every actual support U contained in fixed G.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isfinite, log2, log1p
from pathlib import Path
from time import perf_counter

from hd_shifted_common_source_partition import good_probability, matrix_product, IDENTITY


DIRECTORY = Path(__file__).resolve().parent
ARTIFACT = 'hd_two_sum_late_alpha_quiet_gate_float.json'
ARTIFACT_SHA = '1f9bcb13d4a5120a130e3b23dfad264ffe25a3fc42f9b5c2dc4e8488732d4521'
SOURCE_PINS = {
    'hd_two_sum_late_alpha_quiet.py': '8a753036a9be44759206b5beb3baa044f28d80b265a9304da81bdde6e8f11ae3',
    'hd_shifted_common_source_partition.py': 'fd4bfc68dbb14f53851e00be6d162705f3c9f10cbbc1edffd6417b96918f48f3',
    'hd_xor_pair_route.py': 'dec1f2e2ba064f60529c31675a3779ef53dc6021035bcd024366e5614b3861bc',
    'hd_two_sum_closed_kernel.py': '2ee0fe57de750a19466155317d595589695252de450518b436eab6c5c0fbf979',
    'hd_two_sum_gf16_factorized.py': 'd34f054a6a1a07cec5ae97c75acdf475763a5c5e5ccab85d8be05cdd281e5da2',
}
CONSTRUCTION = 'two-sum-gf16-coset-g16-rs16x8-h4-swapped-pairs-k20-rate-half'
MARK, ALPHA, GROUPS, ACTUAL_OCCUPANCY = F(691, 700), F(1, 5), 18, 9
SOURCE_HISTOGRAMS = ((18, 0), (16, 1))


def checked_artifact(value):
    """Validate the frozen input's meaning separately from its byte pin."""
    expected = dict(construction=CONSTRUCTION, field_order=256,
        active_groups=GROUPS, arity=8, columns=16, constituent_dimension=8,
        forbidden_weight=734003, mark=str(MARK), alpha=str(ALPHA),
        fresh_regional_alpha=1, packet_refresh_keys_already_averaged=True,
        regional_envelope_entry='00',
        regional_geometry_envelope='max over every feasible full-pair count',
        complete_requested_control=True, directed=False,
        whole_code_certificate=False)
    if any(value.get(key) != entry for key, entry in expected.items()):
        raise ValueError('Saved curve has a different proof scope or witness')
    if value.get('source_sha256') != SOURCE_PINS['hd_two_sum_late_alpha_quiet.py']:
        raise ValueError('Saved curve provider identity mismatch')
    for name in ('hd_two_sum_closed_kernel.py', 'hd_two_sum_gf16_factorized.py'):
        if value.get('direct_dependency_sha256', {}).get(name) != SOURCE_PINS[name]:
            raise ValueError('Saved curve dependency identity mismatch')
    envelope = value.get('regional_envelope', {})
    curves = envelope.get('log2_curves', [])
    if (len(curves) != GROUPS+1 or curves[0] != 0
            or envelope.get('null_log_curve_means_exact_zero') is not True
            or any(x is not None and (type(x) not in (int, float) or not isfinite(x)) for x in curves)
            or curves[ACTUAL_OCCUPANCY] is None
            or not isfinite(envelope.get('log2_slope', float('nan')))
            or envelope.get('log2_intercept') != 0):
        raise ValueError('Incomplete or invalid saved unpowered curve')
    return float(curves[ACTUAL_OCCUPANCY]), float(envelope['log2_slope'])*ACTUAL_OCCUPANCY


def load_saved_curve(path=DIRECTORY/ARTIFACT):
    data = Path(path).read_bytes()
    if sha256(data).hexdigest() != ARTIFACT_SHA:
        raise ValueError('Frozen saved curve byte pin mismatch')
    for name, expected in SOURCE_PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != expected:
            raise ValueError('Reviewed source pin mismatch: '+name)
    return checked_artifact(json.loads(data))


def good_region_cap(z, actual_occupancy, *, q=256, state=16, emission_distance=13):
    """Exact generic cone cap, conditional on every actual block having <=4.

    A=(0,beta;1/(q^state-1),beta). Its second row is a density-reference
    bound, not a probability row. For each fixed current input, at most one
    nonzero reference state returns to zero. Fresh MDS D2 H D1 refresh has
    every nonzero atom <=(q-1)^(-state), hence density beta.
    Empty packet diag(1,beta*z^emission_distance) must be <=I before deletion.
    Toy q/state values test this algebra; production uses only the defaults.
    """
    z = F(z)
    if (not 0 < z <= 1 or type(actual_occupancy) is not int or actual_occupancy < 1
            or any(type(x) is not int or x < 1 for x in (q, state, emission_distance)) or q < 2):
        raise ValueError('Positive occupancy and valid field/distance/mark required')
    dimension = q**state-1
    beta = F(dimension, (q-1)**state)
    empty_live = beta*z**emission_distance
    if empty_live > 1:
        raise ValueError('Empty packets are not certified contractions at this mark')
    packet = (F(0), beta, F(1, dimension), beta)
    matrix, powers = IDENTITY, []
    for k in range(1, actual_occupancy+1):
        matrix = matrix_product(matrix, packet)
        powers.append(matrix[0])
    return dict(beta=beta, flat_reference_dimension=dimension,
                return_entry=F(1, dimension), empty_live_cap=empty_live,
                packet=packet, power00=powers, good00=max(powers))


def _log_fraction(value):
    value = F(value)
    if value < 0:
        raise ValueError('Nonnegative rational required')
    if not value:
        return float('-inf')
    approximation = float(value)
    return log2(approximation) if approximation else log2(value.numerator)-log2(value.denominator)


def _logadd(left, right):
    if left == float('-inf'):
        return right
    if right == float('-inf'):
        return left
    high, low = max(left, right), min(left, right)
    return high+log1p(2.**(low-high))/log1p(1.)


def conditional_comparison(probability, good_cap, log2_unconditional, alpha):
    """p*Cgood^alpha +(1-p)^(1-alpha)*Rall^alpha, with exact p.

    Bad conditional mean <=Rall/(1-p) uses positivity only. This may be
    worse than unconditional Jensen; the minimum is also a valid bound.
    No mean is powered before division by its conditioning probability.
    """
    probability, good_cap, alpha = F(probability), F(good_cap), F(alpha)
    if (not 0 <= probability <= 1 or good_cap < 0 or not 0 < alpha <= 1
            or not isfinite(log2_unconditional)):
        raise ValueError('Valid probability, cap, finite mean log and alpha required')
    good = (_log_fraction(probability)+float(alpha)*_log_fraction(good_cap)
            if probability and good_cap else float('-inf'))
    bad_probability = 1-probability
    bad = ((1-float(alpha))*_log_fraction(bad_probability)+float(alpha)*log2_unconditional
           if bad_probability else float('-inf'))
    total = _logadd(good, bad)
    baseline = float(alpha)*log2_unconditional
    serial_log = lambda value: value if isfinite(value) else None
    return dict(good_component_log2=serial_log(good), bad_component_log2=serial_log(bad),
                good_plus_bad_log2=serial_log(total), unconditional_jensen_log2=baseline,
                change_from_unconditional_bits=total-baseline,
                minimum_with_unconditional_log2=min(total, baseline))


def run():
    start = perf_counter()
    curve, line = load_saved_curve()
    caps = good_region_cap(MARK, ACTUAL_OCCUPANCY)
    rows = []
    for singletons, full_pairs in SOURCE_HISTOGRAMS:
        probability = good_probability(singletons, full_pairs)
        rows.append(dict(potential_singletons=singletons, potential_full_pairs=full_pairs,
            potential_active_groups=singletons+2*full_pairs,
            actual_occupancy=ACTUAL_OCCUPANCY, actual_full_pairs_not_fixed=True,
            good_probability=str(probability), bad_probability=str(1-probability),
            bad_probability_log2=_log_fraction(1-probability),
            saved_curve_comparison=conditional_comparison(probability, caps['good00'], curve, ALPHA),
            saved_line_comparison=conditional_comparison(probability, caps['good00'], line, ALPHA)))
    return dict(construction=CONSTRUCTION, field_order=256, mark=str(MARK), alpha=str(ALPHA),
        active_groups=GROUPS, actual_occupancy=ACTUAL_OCCUPANCY,
        scope='regional00 at one fixed actual mask U contained in fixed source G',
        source_good_event_common_to_all_masks=True, source_geometry_resampled=False,
        regional00_zero_claim=False, actual_histogram_specialization=False,
        saved_curve_is_unconditional_max_over_f=True,
        whole_column_bound=False, quiet_event_bound=False, whole_code_certificate=False,
        directed=False, numerical_status='FLOAT diagnostic with exact rational event/cap inputs',
        provider_rebuilds=0, parameter_fits=0, source_count_charged=False,
        probability_outside_fractional_power=True, normalizer_is_flat_nonzero_state=True,
        saved_curve_log2=curve, saved_line_log2=line,
        beta=str(caps['beta']), flat_reference_dimension=str(caps['flat_reference_dimension']),
        nonempty_packet_cap=[str(x) for x in caps['packet']],
        empty_live_cap=str(caps['empty_live_cap']), empty_cap_at_most_identity=True,
        nonempty_packet_power00=[dict(k=k, cap=str(value), log2_cap=_log_fraction(value) if value else None)
                               for k, value in enumerate(caps['power00'], 1)],
        good_region00_cap=str(caps['good00']), good_region00_log2=_log_fraction(caps['good00']),
        cases=rows, saved_artifact=ARTIFACT, saved_artifact_sha256=ARTIFACT_SHA,
        reviewed_source_sha256=SOURCE_PINS,
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        elapsed_seconds=perf_counter()-start)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a frozen diagnostic')
    answer = run()
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(answer, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(dict(elapsed_seconds=answer['elapsed_seconds'], cases=answer['cases']), allow_nan=False))
