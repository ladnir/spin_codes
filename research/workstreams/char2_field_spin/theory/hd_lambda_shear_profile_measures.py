"""Full r18 NZ regional bank and common measures for fixed source profiles.

Fix r source groups with F full pairs in one matching. An active m-subset
has f full pairs only if max(0,m-r+F)<=f<=min(F,floor(m/2)). Maximizing
R[m,f] over this range gives a regional envelope for that fixed source.
For alternating profile(a,b), one measure dominates both envelopes. It
is fixed across all row sets and regional positions. Unordered(a,b)
keys reuse identical envelopes; they do not identify source supports.

The provider command saves all100 original-one-state matrices at9927/10000.
The profiles command reads that pinned bank, fits prescribed secants,
and checks every feasible matrix entry against the actual five-atom law.
All scalar atoms must lie in[0,1]. A failed fit is recorded, not repaired
by a witness search. Small intercept roundoff is increased and rechecked.

The99 indexed profile counts come from profile_count_upper: they are exact
integer UPPER counts, not an exact partition. They remain outside moments.
Neither command contracts an outer kernel or produces a distance bound.
Everything numerical here is FLOAT, not a directed certificate.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isfinite, log2
from pathlib import Path
from time import perf_counter

import numpy as np

import hd_lambda_column_shear_control as frozen
import hd_mds_exact_empty_regions as provider
from hd_coupled_state_moment_gate import entry_measure, MomentFitError
from hd_pair_state_kernel import common_moment_atoms
from hd_two_sum_gf16_high_arity_certificate import profile_counts


DIRECTORY = Path(__file__).resolve().parent
MARK, MAXIMUM, LIMIT = F(9927, 10000), 18, 734003
SCHEMA = 'lambda-shear-NZ-profile-measures-v1'
PINS = {
    'hd_lambda_column_shear_control.py': 'a1ede595d659ddd4c20d1b5ea33a98158e7c9484e34a958d2d682f985214de96',
    'hd_mds_exact_empty_regions.py': '120fa611191903574b4a090e1ad0e0861bf02680df389890c70ea65e3adc6153',
    'hd_coupled_state_moment_gate.py': 'd2258e48f19bf3c8b7abea6a3d95024c10ce4d7c7286f49d009f52d8a4d7e89a',
    'hd_pair_state_kernel.py': 'b28142e19e51a906ae742128325d1e7557dfb425624c6da45778a2a75753e195',
    'hd_two_sum_gf16_high_arity_certificate.py': '9c9c2e4d1bfc95f2952d6bc743b9047489128d5cc26a43f5186ea721b1da4027',
    'hd_shifted_pairs_count.py': 'bb1700278361343325007cf558751af2db6df19dbf6c2090502941f1a95397be',
}
BANK_SCOPE = dict(schema=SCHEMA, kind='full-raw-bank', mark=str(MARK), maximum=MAXIMUM,
    expected_matrices=100, regional_alpha=1, regional_amplitudes='independent-NZ-conditional-on-full-mask',
    inner='original-one-state-h4-width64', final_outer_D2_required=True,
    shifted_pair_permutations=True, independent_pair_swaps=True,
    matrix_entry_order=['00', '01', '10', '11'], empty11='exact-inverse-refresh-weight-cap',
    source_counts_in_local_matrices=False, arithmetic='positive-FLOAT-log2-not-directed',
    directed=False, whole_code_certificate=False)


def dependencies():
    pins = frozen.dependencies()
    for name, expected in PINS.items():
        if (name in pins and pins[name] != expected
                or sha256((DIRECTORY/name).read_bytes()).hexdigest() != expected):
            raise ValueError('Reviewed numerical dependency changed: '+name)
        pins[name] = expected
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    return pins


def feasible_full_counts(r, m, full):
    if (any(type(x) is not int for x in (r, m, full)) or not 1 <= r <= 8192
            or not 0 <= m <= r or not 0 <= full <= r//2):
        raise ValueError('Valid total source size, active subset size and full-pair count required')
    return tuple(range(max(0, m-r+full), min(full, m//2)+1))


def full_keys(maximum):
    return tuple((m, f) for m in range(maximum+1) for f in range(m//2+1))


def checked_regions(regions, maximum):
    if set(regions) != set(full_keys(maximum)):
        raise ValueError('Every raw(m,f) entry is required exactly once')
    for key, value in regions.items():
        if any(type(x) is not int for x in key):
            raise ValueError('Integer raw histogram indices required')
        values = np.asarray(value, dtype=float)
        if values.shape != (2, 2) or np.any(np.isnan(values)) or np.any(np.isposinf(values)):
            raise ValueError('Raw matrices must contain finite logarithms or exact zeros')
    return regions


def _serial_matrix(matrix):
    return [[None if x == -np.inf else float(x) for x in row] for row in matrix]


def encode_regions(regions, maximum):
    checked_regions(regions, maximum)
    return [dict(m=m, f=f, log2_matrix=_serial_matrix(regions[m, f])) for m, f in full_keys(maximum)]


def decode_regions(rows, maximum):
    regions = {}
    for row in rows:
        key = row.get('m'), row.get('f')
        if any(type(x) is not int for x in key) or key in regions:
            raise ValueError('Duplicate or noninteger raw histogram')
        values = row.get('log2_matrix')
        if (not isinstance(values, list) or len(values) != 2
                or any(not isinstance(line, list) or len(line) != 2 for line in values)
                or any(x is not None and (type(x) not in (int, float) or not isfinite(x))
                       for line in values for x in line)):
            raise ValueError('Malformed raw logarithmic matrix')
        regions[key] = np.array([[-np.inf if x is None else x for x in line] for line in values])
    return checked_regions(regions, maximum)


def make_bank(regions, metadata, pins, *, elapsed_seconds=0.):
    if (metadata.get('mark') != str(MARK) or metadata.get('maximum') != MAXIMUM
            or metadata.get('reference_cone_unchanged') is not True
            or metadata.get('nonempty_packet_caps_unchanged') is not True):
        raise ValueError('Provider metadata changed the original regional interface')
    cap = metadata['exact_cap']
    if not 0 < F(cap) <= 1 or sha256(cap.encode()).hexdigest() != metadata['exact_cap_sha256']:
        raise ValueError('Invalid exact empty-cap provenance')
    return dict(BANK_SCOPE, complete=True, rows=encode_regions(regions, MAXIMUM),
        provider_metadata=metadata, reviewed_source_sha256=pins, source_pins_unchanged=True,
        raw_regional_banks_built=1, elapsed_seconds=elapsed_seconds, timeout=None)


def validate_bank(data, pins):
    if any(data.get(key) != value for key, value in BANK_SCOPE.items()):
        raise ValueError('Wrong raw regional bank scope')
    if (data.get('complete') is not True or data.get('reviewed_source_sha256') != pins
            or data.get('source_pins_unchanged') is not True or data.get('raw_regional_banks_built') != 1):
        raise ValueError('Incomplete or unpinned raw bank')
    regions = decode_regions(data['rows'], MAXIMUM)
    make_bank(regions, data['provider_metadata'], pins)
    return regions


def load_bank(path, expected_sha256):
    pins = dependencies()
    raw = Path(path).read_bytes()
    if sha256(raw).hexdigest() != expected_sha256:
        raise ValueError('Full raw bank differs from its explicit receipt pin')
    data = json.loads(raw)
    return data, validate_bank(data, pins)


def build_bank(*, max_seconds=120):
    if type(max_seconds) not in (int, float) or not isfinite(max_seconds) or not 0 < max_seconds <= 180:
        raise ValueError('Positive provider budget at most180 seconds required')
    pins, started = dependencies(), perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Full raw regional bank exhausted its budget')
    try:
        regions, metadata = provider.forced_nonzero_regions(MARK, MAXIMUM, check_time=check_time)
        check_time()
        result = make_bank(regions, metadata, pins, elapsed_seconds=perf_counter()-started)
    except TimeoutError as exception:
        result = dict(BANK_SCOPE, complete=False, rows=[], timeout=str(exception),
            reviewed_source_sha256=pins, source_pins_unchanged=True,
            raw_regional_banks_built=0, elapsed_seconds=perf_counter()-started)
    if dependencies() != pins: raise ValueError('Numerical sources changed during raw bank construction')
    result['max_seconds'] = max_seconds
    return result


def profile_curves(regions, r, full_counts):
    if type(r) is not int or not 2 <= r <= MAXIMUM:
        raise ValueError('Profile fitting supports2<=r<=18')
    full_counts = tuple(sorted(set(full_counts)))
    if not full_counts or any(type(f) is not int or not 0 <= f <= r//2 for f in full_counts):
        raise ValueError('Nonempty in-scope fixed-source full-pair counts required')
    rows, feasible = [], []
    for m in range(r+1):
        counts = sorted({f for full in full_counts for f in feasible_full_counts(r, m, full)})
        matrices = np.stack([np.asarray(regions[m, f], dtype=float) for f in counts])
        if matrices.shape != (len(counts), 2, 2) or np.any(np.isnan(matrices)) or np.any(np.isposinf(matrices)):
            raise ValueError('Every feasible positive regional matrix is required')
        rows.append(np.max(matrices, axis=0).ravel())
        feasible.append(counts)
    return np.array(rows).T, feasible


def moment_curves(measure, r):
    """Evaluate the emitted five-atom law, including0^0=1 at degree zero."""
    atoms, matrices = common_moment_atoms(measure)
    if np.any(atoms > 0): raise MomentFitError('Shared-ratio support comparison requires every atom<=1')
    result = []
    for m in range(r+1):
        powers = np.zeros(len(atoms)) if m == 0 else m*atoms
        result.append(np.logaddexp2.reduce(matrices+powers[:, None, None], axis=0).ravel())
    return np.array(result).T


def fit_profile_measure(regions, r, profile, *, zero_endpoints=None, cross_endpoints=None):
    if len(profile) != 2 or any(type(f) is not int for f in profile):
        raise ValueError('An ordered pair of fixed matching full-pair counts is required')
    key = tuple(sorted(profile))
    curves, feasible = profile_curves(regions, r, key)
    zero_endpoints = (min(5, r-1), r) if zero_endpoints is None else tuple(zero_endpoints)
    cross_endpoints = (min(2, r-1), r) if cross_endpoints is None else tuple(cross_endpoints)
    measure = entry_measure(curves, zero_endpoints=zero_endpoints, cross_endpoints=cross_endpoints)
    intercepts = ('log2_probability', 'log2_scale01', 'log2_scale10', 'log2_scale11')
    # Compensate only floating evaluation/intercept rounding; slopes never change.
    adjustments = 0
    for _ in range(8):
        bound = moment_curves(measure, r)
        deficits = np.max(curves-bound, axis=1)
        if np.all(deficits <= 0): break
        for index, deficit in enumerate(deficits):
            if deficit > 0:
                name = intercepts[index]
                measure[name] = float(np.nextafter(measure[name]+deficit, np.inf))
                adjustments += 1
        if measure['log2_probability'] > 0:
            raise MomentFitError('Roundoff repair would make the00 probability exceed one')
    bound = moment_curves(measure, r)
    if np.any(bound < curves): raise ArithmeticError('Five-atom law misses a feasible regional entry')
    slacks, checks = [], 0
    for m, values in enumerate(feasible):
        for f in values:
            slack = bound[:, m]-np.asarray(regions[m, f]).ravel()
            if np.any(slack < 0): raise ArithmeticError('A raw feasible state entry escaped curve domination')
            slacks.extend(slack[np.isfinite(slack)])
            checks += 4
    mapping = dict(probability_log='log2_probability', atom00_log='log2_atom00',
        atom01_log='log2_atom01', atom10_log='log2_atom10', scale01_log='log2_scale01',
        scale10_log='log2_scale10', scale11_log='log2_scale11')
    measure['exact_float_parameters'] = {name: str(F(measure[field])) for name, field in mapping.items()}
    measure['minimum_slack'] = float(np.min(bound-curves))
    measure['contacts'] = [np.flatnonzero(row < 1e-9).tolist() for row in bound-curves]
    return dict(profile_key=list(key), active_groups=r, state_measure=measure,
        regional_log2_curves=[[None if x == -np.inf else float(x) for x in row] for row in curves],
        feasible_full_counts_by_occupancy=feasible, every_feasible_state_checked=True,
        feasible_matrix_entries_checked=checks, minimum_raw_slack=min(slacks),
        all_scalar_atoms_at_most_one=True, floating_intercept_adjustments=adjustments,
        common_measure_all_rowsets=True, same_measure_all_regions=True,
        alternating_profile_order_not_quotiented=True, source_count_inside_measure=False)


def profile_inventory(path, expected_sha256, *, zero_endpoints=None, cross_endpoints=None, max_seconds=30):
    if type(max_seconds) not in (int, float) or not isfinite(max_seconds) or not 0 < max_seconds <= 60:
        raise ValueError('Positive fitting budget at most60 seconds required')
    pins, started = dependencies(), perf_counter()
    _, regions = load_bank(path, expected_sha256)
    counts = profile_counts(MAXIMUM)
    if len(counts) != 99: raise ValueError('The r18 profile inventory must have99 indexed entries')
    keys = sorted({tuple(sorted(key)) for key in counts})
    measures, failures = [], []
    for key in keys:
        if perf_counter()-started >= max_seconds: break
        try:
            measures.append(fit_profile_measure(regions, MAXIMUM, key,
                zero_endpoints=zero_endpoints, cross_endpoints=cross_endpoints))
        except (ValueError, ArithmeticError) as exception:
            failures.append(dict(profile_key=list(key), reason=str(exception)))
    if dependencies() != pins or sha256(Path(path).read_bytes()).hexdigest() != expected_sha256:
        raise ValueError('Bank or source changed during profile fitting')
    examined = len(measures)+len(failures)
    return dict(schema=SCHEMA, kind='profile-measures', complete=examined == len(keys),
        all_profiles_have_valid_measures=examined == len(keys) and not failures,
        mark=str(MARK), active_groups=MAXIMUM, profile_entries=len(counts),
        unordered_measure_keys=len(keys), examined_measure_keys=examined,
        profiles=[dict(even_full=a, odd_full=b, indexed_source_count_upper=str(count),
            measure_key=sorted((a, b)), source_and_projective_prefix_log2=
                log2(count)-LIMIT*log2(MARK)-log2(255)) for (a, b), count in counts.items()],
        profile_count_method='profile_count_upper-exact-integer-upper-not-exact-partition',
        source_count_outside_moments_once=True, extra_source_factor=False,
        measures=measures, failed_measure_keys=failures,
        input_bank=str(Path(path).resolve()), input_bank_sha256=expected_sha256,
        reviewed_source_sha256=pins, source_pins_unchanged=True, raw_regional_banks_built=0,
        arithmetic='positive-FLOAT-log2-not-directed', directed=False, whole_code_certificate=False,
        elapsed_seconds=perf_counter()-started, max_seconds=max_seconds)


def load_profile_measures(path, expected_sha256, *, bank_path=None):
    """Return(metadata, by_ordered_profile); failed fits have no usable entry.

    Each entry contains state_measure, integer source_count (an upper count),
    common_prefix_log2, and its unordered measure key. No provider or refit.
    The raw bank is reloaded so domination is checked independently of saved
    contact/slack claims. All99 indexed count records are required.
    """
    pins = dependencies()
    raw = Path(path).read_bytes()
    if sha256(raw).hexdigest() != expected_sha256:
        raise ValueError('Profile-measure receipt differs from its explicit pin')
    data = json.loads(raw)
    required = dict(schema=SCHEMA, kind='profile-measures', complete=True, mark=str(MARK),
        active_groups=MAXIMUM, profile_entries=99,
        profile_count_method='profile_count_upper-exact-integer-upper-not-exact-partition',
        source_count_outside_moments_once=True, extra_source_factor=False,
        reviewed_source_sha256=pins, source_pins_unchanged=True, raw_regional_banks_built=0,
        arithmetic='positive-FLOAT-log2-not-directed', directed=False, whole_code_certificate=False)
    if any(data.get(key) != value for key, value in required.items()):
        raise ValueError('Profile receipt is incomplete or has wrong scope/provenance')
    actual_bank_path = Path(data['input_bank']) if bank_path is None else Path(bank_path)
    _, regions = load_bank(actual_bank_path, data['input_bank_sha256'])
    counts = profile_counts(MAXIMUM)
    indexed, measures, failed = {}, {}, set()
    for row in data['profiles']:
        key = row.get('even_full'), row.get('odd_full')
        if (any(type(x) is not int for x in key) or key not in counts or key in indexed
                or row.get('indexed_source_count_upper') != str(counts[key])
                or row.get('measure_key') != sorted(key)):
            raise ValueError('Wrong, duplicated, or relabelled indexed profile count')
        prefix = log2(counts[key])-LIMIT*log2(MARK)-log2(255)
        saved_prefix = row.get('source_and_projective_prefix_log2')
        if type(saved_prefix) not in (int, float) or not isfinite(saved_prefix) or abs(saved_prefix-prefix) > 1e-9:
            raise ValueError('Profile prefix must charge its count and projective factor once')
        indexed[key] = prefix
    if set(indexed) != set(counts): raise ValueError('Missing indexed source profile')
    expected_keys = {tuple(sorted(key)) for key in counts}
    for row in data['measures']:
        key = tuple(row['profile_key'])
        if any(type(x) is not int for x in key) or key not in expected_keys or key in measures:
            raise ValueError('Wrong or duplicate fitted profile key')
        curves, feasible = profile_curves(regions, MAXIMUM, key)
        encoded = [[None if x == -np.inf else float(x) for x in line] for line in curves]
        if (row.get('active_groups') != MAXIMUM or row.get('feasible_full_counts_by_occupancy') != feasible
                or row.get('regional_log2_curves') != encoded):
            raise ValueError('Saved curves do not come from the pinned feasible raw entries')
        bound = moment_curves(row['state_measure'], MAXIMUM)
        if np.any(bound < curves): raise ValueError('Saved measure misses a feasible state entry')
        measures[key] = row['state_measure']
    for row in data['failed_measure_keys']:
        key = tuple(row['profile_key'])
        if (any(type(x) is not int for x in key) or key not in expected_keys
                or key in failed or key in measures):
            raise ValueError('Wrong or duplicate failed measure key')
        failed.add(key)
    if (set(measures) | failed) != expected_keys or data.get('examined_measure_keys') != len(expected_keys):
        raise ValueError('Not every unordered profile was examined')
    if (data.get('unordered_measure_keys') != len(expected_keys)
            or data.get('all_profiles_have_valid_measures') is not (not failed)):
        raise ValueError('False profile-coverage claim')
    by_profile = {key: dict(state_measure=measures[tuple(sorted(key))], source_count=counts[key],
        source_count_is_upper=True, common_prefix_log2=indexed[key], profile_key=list(sorted(key)),
        active_groups=MAXIMUM, even_full=key[0], odd_full=key[1], mark=str(MARK))
        for key in counts if tuple(sorted(key)) in measures}
    if (dependencies() != pins or sha256(Path(path).read_bytes()).hexdigest() != expected_sha256
            or sha256(actual_bank_path.read_bytes()).hexdigest() != data['input_bank_sha256']):
        raise ValueError('Profile receipt, raw bank, or sources changed during validation')
    return data, by_profile


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    bank = commands.add_parser('bank', help='Build one full raw bank; no outer gate or measure fit')
    bank.add_argument('--max-seconds', type=float, default=120)
    bank.add_argument('--output', type=Path, required=True)
    profiles = commands.add_parser('profiles', help='Fit profiles from one explicitly pinned saved bank')
    profiles.add_argument('--bank', type=Path, required=True)
    profiles.add_argument('--bank-sha256', required=True)
    profiles.add_argument('--zero-endpoints', type=int, nargs=2)
    profiles.add_argument('--cross-endpoints', type=int, nargs=2)
    profiles.add_argument('--max-seconds', type=float, default=30)
    profiles.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists(): parser.error('Refuse to overwrite evidence')
    result = build_bank(max_seconds=args.max_seconds) if args.command == 'bank' else profile_inventory(
        args.bank, args.bank_sha256, zero_endpoints=args.zero_endpoints,
        cross_endpoints=args.cross_endpoints, max_seconds=args.max_seconds)
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({key: result[key] for key in ('kind', 'complete', 'elapsed_seconds')}, allow_nan=False), flush=True)
