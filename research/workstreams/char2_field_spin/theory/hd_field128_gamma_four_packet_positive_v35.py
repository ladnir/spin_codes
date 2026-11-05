"""One-profile FLOAT test of delayed alpha; no FFT or distance certificate.

Keep every routed histogram: selected occupancy quartets are composed before
alpha, and all omitted pairs remain in a directly accumulated positive remainder.
The state32 caps are conditional on the same global setup event as frozen v34.
"""
import argparse
import ast
from contextlib import contextmanager
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, isfinite, log2
from pathlib import Path
import shutil
import sys
from time import monotonic

import hd_field128_gamma_state_ladder_v34 as previous

HERE = Path(__file__).resolve().parent
SOURCE_SHA = '938de22bbf7e1b25fc343813c404bb50fb73564017b010c4c8efc223250d3aa9'
PREVIOUS_SHA = '3eaba3ddd877f0f38e6efe16f9f55c1e05a5a3c9ae790478551a8c1f405fb265'
OUTPUT_LIMIT, RESERVE = 64 << 10, 20 << 30
P, R, DISTANCE = 4096, 4096, F(9, 20)


def source_pins():
    pending, pins = [Path(__file__).name], {}
    while pending:
        name = pending.pop()
        if name in pins:
            continue
        raw = (HERE/name).read_bytes()
        pins[name] = sha256(raw).hexdigest()
        for node in ast.walk(ast.parse(raw)):
            modules = ([item.name for item in node.names] if isinstance(node, ast.Import)
                       else [node.module] if isinstance(node, ast.ImportFrom) and node.module else [])
            for module in modules:
                child = module.split('.')[0]+'.py'
                if (HERE/child).is_file() and child not in pins:
                    pending.append(child)
    for name, digest in previous.source_pins().items():
        if name in pins and pins[name] != digest:
            raise RuntimeError('Conflicting dependency bytes')
        pins[name] = digest
    if pins.get(Path(previous.__file__).name) != PREVIOUS_SHA:
        raise RuntimeError('Frozen v34 identity changed')
    note = '../HD_FIELD128_ALL_PAIR_COMPOSITION.md'
    pins[note] = sha256((HERE/note).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def validate_source(artifact):
    if (artifact.get('schema') != 'gamma-state-ladder-conditional-FLOAT-v34'
            or artifact.get('active_groups') != R or F(artifact['target_distance']) != DISTANCE):
        raise ValueError('Only the pinned original-r4096/d45 comparison is supported')
    if artifact.get('source_sha256') != previous.source_pins():
        raise ValueError('V34 source closure differs')
    row = artifact['best_by_state']['32']
    for key, value in previous.state_metadata(32).items():
        if row.get(key) != value:
            raise ValueError('Conditional state32 metadata mismatch: '+key)
    if (row.get('active_groups') != R or F(row['target_distance']) != DISTANCE
            or row.get('profile_pairs_evaluated') != 2049**2
            or row.get('complete_selected_weight_screen') is not True
            or row.get('all_profiles_retained') is not True):
        raise ValueError('Incomplete or wrong original-r input')
    witness = row['witness']
    z, alpha = F(witness['z']), F(witness['alpha'])
    if not 0 < z < 1 or not 0 < alpha <= 1 or F(float(z)) != z or F(float(alpha)) != alpha:
        raise ValueError('Exact positive binary64 z/alpha required')
    markers = witness['markers']
    if markers != artifact.get('fixed_positive_marker_literals') or not 1 <= len(markers) <= 128:
        raise ValueError('Fixed marker lineage mismatch')
    for pair in markers:
        if not isinstance(pair, list) or len(pair) != 2:
            raise ValueError('Two positive literals per marker required')
        if any(F(v) <= 0 or not isfinite(float(F(v))) or float(F(v)) <= 0 for v in pair):
            raise ValueError('Invalid positive marker')
    profile = row['largest_profile']
    for key in ('even_full', 'odd_full'):
        if type(profile[key]) is not int or not 0 <= profile[key] <= R//2:
            raise ValueError('Invalid original-r profile')
    for value in (profile['log2_count_upper'], row['log2_outer_prefix'], row['log2_upper_screen']):
        if not isfinite(value):
            raise ValueError('Nonfinite saved diagnostic')
    return row


def compose(left, right):
    """Chronological 2x2 log product, including leading batch axes."""
    import numpy as np
    return np.logaddexp2(left[..., :, 0, None]+right[..., None, 0, :],
                        left[..., :, 1, None]+right[..., None, 1, :])


def power(matrix, exponent):
    import numpy as np
    if type(exponent) is not int or exponent < 0:
        raise ValueError('Nonnegative integer exponent required')
    result = np.array([[0., -np.inf], [-np.inf, 0.]])
    while exponent:
        if exponent & 1:
            result = compose(result, matrix)
        matrix = compose(matrix, matrix)
        exponent >>= 1
    return result


def matrix_sum(values):
    import numpy as np
    return np.logaddexp2.reduce(np.asarray(values).reshape((-1, 2, 2)), axis=0,
                               initial=-np.inf)


def occupancy_weights(odds, slots=32):
    """g_j(x,y): exact histogram multiplicities outside alpha; swaps average to1."""
    import numpy as np
    if type(slots) is not int or not 1 <= slots <= 32:
        raise ValueError('Invalid slot count')
    weights = np.full(2*slots+1, -np.inf)
    for singles in range(slots+1):
        for full in range(slots-singles+1):
            value = (log2(comb(slots, singles)*comb(slots-singles, full))
                     + singles*odds[0]+full*odds[1])
            j = singles+2*full
            weights[j] = np.logaddexp2(weights[j], value)
    return weights


def pair_terms(raw_pairs, weights, alpha):
    return weights[:, None, None, None]+weights[None, :, None, None]+alpha*raw_pairs


def choose_occupancies(terms, limit):
    """Top tilted pair marginals; selection changes cost only, never support."""
    import numpy as np
    total = matrix_sum(terms)
    shifted = np.full_like(terms, -np.inf)
    np.subtract(terms, total, out=shifted, where=np.isfinite(total)[None, None])
    with np.errstate(under='ignore'):
        mass = np.exp2(shifted)
    marginal = mass.sum(axis=(1, 2, 3))+mass.sum(axis=(0, 2, 3))
    return tuple(sorted(map(int, np.argsort(-marginal, kind='stable')[:limit])))


def positive_partition(raw_pairs, weights, alpha, selected, check=lambda: None):
    """C4+C2R2+R2C2+R2^2. R2 is accumulated from omitted pairs, not subtracted."""
    import numpy as np
    size = len(weights)
    if (raw_pairs.shape != (size, size, 2, 2) or not 0 < alpha <= 1
            or not 1 <= len(selected) <= 16 or len(set(selected)) != len(selected)
            or any(type(j) is not int or not 0 <= j < size for j in selected)):
        raise ValueError('Invalid raw pair/selection scope')
    mask = np.zeros(size, dtype=bool)
    mask[list(selected)] = True
    both = mask[:, None] & mask[None, :]
    terms = pair_terms(raw_pairs, weights, alpha)
    c2, rest, t2 = matrix_sum(terms[both]), matrix_sum(terms[~both]), matrix_sum(terms)
    small = raw_pairs[both]
    pair_weights = (weights[:, None]+weights[None, :])[both]
    c4 = np.full((2, 2), -np.inf)
    for start in range(0, len(small), 32):
        check()
        raw = compose(small[start:start+32, None], small[None])
        values = (pair_weights[start:start+32, None, None, None]
                  +pair_weights[None, :, None, None]+alpha*raw)
        c4 = np.logaddexp2(c4, matrix_sum(values))
    old4 = compose(t2, t2)
    new4 = matrix_sum([c4, compose(c2, rest), compose(rest, c2), compose(rest, rest)])
    finite = np.isfinite(old4)
    tolerance = 2e-10*max(1., float(np.max(np.abs(old4[finite]), initial=0.)))
    if (np.any(new4 > old4+tolerance)
            or np.any(new4 < old4-(1-alpha)-tolerance)):
        raise ArithmeticError('Fractional boundary inequalities failed')
    retained = np.zeros((2, 2))
    with np.errstate(under='ignore'):
        retained[finite] = np.exp2(c2[finite]-t2[finite])
    return dict(two=t2, four=new4, selected_two=c2, remainder_two=rest,
                selected_four=c4, old_four=old4, retained_pair_entry_fraction=retained,
                quartets_evaluated=len(small)**2, selected_occupancies=list(selected))


def regional(matrix, blocks, weight, full, odds, pairs=P):
    singles, empty = weight-2*full, pairs-weight+full
    if min(singles, empty, full) < 0:
        raise ValueError('Impossible histogram')
    normalizer = log2(comb(pairs, singles)*comb(pairs-singles, full))
    return power(matrix, blocks)-normalizer-singles*odds[0]-full*odds[1]


def profile_transfer(even, odd, repetitions=128):
    import numpy as np
    run = power(compose(even, odd), repetitions)
    return float(np.logaddexp2(run[0, 0], run[0, 1]))


def write_result(path, result):
    data = (json.dumps(previous.json_ready(result), sort_keys=True, indent=2, allow_nan=False)+'\n').encode()
    if len(data) > OUTPUT_LIMIT or shutil.disk_usage(HERE).free-len(data) < RESERVE:
        raise RuntimeError('64KiB output or20GiB planned reserve exceeded')
    with path.open('xb') as stream:
        stream.write(data)
    if path.read_bytes() != data:
        raise RuntimeError('Exclusive output changed')
    return dict(path=path.relative_to(HERE).as_posix(), sha256=sha256(data).hexdigest(), bytes=len(data))


@contextmanager
def remote_limits():
    if sys.platform != 'linux':
        raise RuntimeError('Production requires guarded Linux')
    import resource
    import signal
    resource.setrlimit(resource.RLIMIT_AS, (512 << 20, 512 << 20))
    resource.setrlimit(resource.RLIMIT_CPU, (50, 55))
    resource.setrlimit(resource.RLIMIT_FSIZE, (OUTPUT_LIMIT, OUTPUT_LIMIT))
    signal.alarm(60)
    try:
        yield
    finally:
        signal.alarm(0)


def _run(source_name, expected_sha256, output_name, *, max_markers=4, max_occupancies=16, max_seconds=50.):
    if (type(max_markers) is not int or not 1 <= max_markers <= 4
            or type(max_occupancies) is not int or not 1 <= max_occupancies <= 16
            or not isfinite(max_seconds) or not 0 < max_seconds <= 50
            or expected_sha256 != SOURCE_SHA):
        raise ValueError('Use approved source and at most4markers/16occupancies/50seconds')
    source, output = previous.inside(source_name), previous.inside(output_name)
    if output.exists() or source.stat().st_size > 1 << 20:
        raise ValueError('Exclusive output and input<=1MiB required')
    if shutil.disk_usage(HERE).free-OUTPUT_LIMIT < RESERVE:
        raise RuntimeError('20GiB reserve after output budget required')
    raw, pins = source.read_bytes(), source_pins()
    if sha256(raw).hexdigest() != expected_sha256:
        raise ValueError('Pinned source mismatch')
    row = validate_source(json.loads(raw))
    def check_inputs():
        if source_pins() != pins or source.read_bytes() != raw:
            raise RuntimeError('Input or full source closure changed')
    started = monotonic()
    def check():
        if monotonic()-started > max_seconds:
            raise TimeoutError('Bounded one-profile diagnostic timed out; no complete output')
    try:
        import numpy as np
        import hd_field128_full_rs_gamma_screen as gamma
        z, alpha = float(F(row['witness']['z'])), float(F(row['witness']['alpha']))
        logs = gamma.local_logs(z, previous.Q, state=32)
        raw_pairs = compose(logs[:, None], logs[None])
        profile = row['largest_profile']
        even, odd = profile['even_full'], profile['odd_full']
        markers = row['witness']['markers']
        ranked = []
        for index, pair in enumerate(markers):
            check()
            odds = tuple(log2(float(F(v))) for v in pair)
            weights = occupancy_weights(odds)
            two = matrix_sum(pair_terms(raw_pairs, weights, alpha))
            a, b = regional(two, 64, R, even, odds), regional(two, 64, R, odd, odds)
            ranked.append((profile_transfer(a, b), index, odds, weights, a, b))
        selected = sorted(ranked, key=lambda item: (item[0], item[1]))[:max_markers]
        old_even = np.full((2, 2), np.inf)
        old_odd, new_even, new_odd = old_even.copy(), old_even.copy(), old_even.copy()
        details = []
        for _, index, odds, weights, a, b in selected:
            check()
            terms = pair_terms(raw_pairs, weights, alpha)
            chosen = choose_occupancies(terms, max_occupancies)
            parts = positive_partition(raw_pairs, weights, alpha, chosen, check)
            aa = regional(parts['four'], 32, R, even, odds)
            bb = regional(parts['four'], 32, R, odd, odds)
            old_even, old_odd = np.minimum(old_even, a), np.minimum(old_odd, b)
            new_even, new_odd = np.minimum(new_even, aa), np.minimum(new_odd, bb)
            details.append(dict(marker_index=index, marker_literals=markers[index],
                selected_occupancies=parts['selected_occupancies'],
                retained_pair_entry_fraction=parts['retained_pair_entry_fraction'].tolist(),
                quartets_evaluated=parts['quartets_evaluated'],
                log2_old_four=parts['old_four'].tolist(), log2_refined_four=parts['four'].tolist(),
                log2_old_even=a.tolist(), log2_new_even=aa.tolist(),
                log2_old_odd=b.tolist(), log2_new_odd=bb.tolist()))
        check()
        prefix = previous.prefactor(z, alpha, DISTANCE, R, gamma.old.outer_log_density(previous.Q, 'invertible-toeplitz'))
        if abs(prefix-row['log2_outer_prefix']) > 1e-8*max(1., abs(prefix)):
            raise ArithmeticError('Original-r prefactor differs')
        count_log = profile['log2_count_upper']
        old_score = prefix+count_log+profile_transfer(old_even, old_odd)
        new_score = prefix+count_log+profile_transfer(new_even, new_odd)
        if not isfinite(new_score) or new_score > old_score+1e-6:
            raise ArithmeticError('Refined positive profile bound did not dominate baseline')
        result = dict(schema='gamma-conditional32-four-packet-positive-diagnostic-v35',
            source=dict(path=source_name, sha256=expected_sha256), source_sha256=pins,
            conditional_construction=previous.state_metadata(32), active_groups=R, target_distance=str(DISTANCE),
            selected_profile=[even, odd], witness_z=row['witness']['z'], witness_alpha=row['witness']['alpha'],
            marker_candidates=len(markers), selected_markers=len(selected), max_occupancies=max_occupancies,
            marker_selection='lowest same-marker positive bound for the pinned profile; no optimization',
            occupancy_selection='largest pair-entry normalized marginals; all omitted pairs retained',
            positive_partition='C4(V^4)+C2*R2+R2*C2+R2*R2; R2 summed directly',
            pair_slots_per_packet=32, packets_per_region=128, four_blocks_per_region=32, regions=256,
            saved_log2_count_upper=count_log, count_source='copied pinned FLOAT profile cap; not recomputed',
            log2_outer_prefix=prefix, histogram_normalizers_recomputed=True,
            log2_old_positive_profile=old_score, log2_new_positive_profile=new_score,
            gain_bits=old_score-new_score, full_four_packet_gain_ceiling_bits=8192*(1-alpha),
            saved_alias_profile_reference=prefix+profile['log2_term'],
            saved_all_profile_reference=row['log2_upper_screen'],
            comparison_scope='old/new positive markers only; not improvement relative to saved aliases',
            details=details, all_other_route_cases_retained=True, counts_outside_alpha=True,
            chronology_preserved=True, selected_profile_terms_evaluated=1, all_profiles_replayed=False,
            coefficient_atlas_built=False, optimizer_run=False, directed=False,
            selected_weight_certificate=False, whole_code_certificate=False,
            setup_risk_added_to_score=False, arithmetic='FLOAT diagnostics only',
            elapsed_seconds=monotonic()-started,
            limits=dict(wall_seconds=60, address_space_bytes=512 << 20, output_bytes=OUTPUT_LIMIT,
                        minimum_free_after_write=RESERVE))
        check_inputs()
        receipt = write_result(output, result)
        return dict(receipt=receipt, gain_bits=result['gain_bits'],
                    log2_new_positive_profile=new_score, directed=False, all_profiles_replayed=False)
    finally:
        check_inputs()


def run(*args, **kwargs):
    with remote_limits():
        return _run(*args, **kwargs)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--sha256', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--max-markers', type=int, default=4)
    parser.add_argument('--max-occupancies', type=int, default=16)
    parser.add_argument('--max-seconds', type=float, default=50.)
    args = parser.parse_args()
    print(json.dumps(run(args.source, args.sha256, args.output, max_markers=args.max_markers,
        max_occupancies=args.max_occupancies, max_seconds=args.max_seconds), allow_nan=False), flush=True)


if __name__ == '__main__':
    main()
