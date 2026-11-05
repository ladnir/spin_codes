"""Optimistic return ablation at pinned marks; NOT an upper bound/certificate.

The original local caps are retained except that every live-to-zero entry
is set to exact zero BEFORE packet composition and the fractional power.
This diagnoses the headroom of return-only changes in the positive-marker
method at these fixed marks. It is not an alias replay, optimized witness,
construction, or claim about the actual minimum distance.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import log2
from pathlib import Path
import shutil
import sys
from time import monotonic

import hd_field128_gamma_four_packet_positive_v35 as base

LIMIT = 16 << 10


def ablate_returns(logs):
    import numpy as np
    result = np.asarray(logs, dtype=float).copy()
    if result.ndim != 3 or result.shape[1:] != (2, 2):
        raise ValueError('Expected occupancy-indexed 2x2 logarithms')
    if np.any(np.isnan(result)) or np.any(np.isposinf(result)):
        raise ValueError('Invalid log entries')
    result[:, 1, 0] = -np.inf
    return result


def evaluate(logs, alpha, markers, weight, even, odd, *, pairs=4096,
             blocks=64, repetitions=128):
    import numpy as np
    raw_pairs = base.compose(logs[:, None], logs[None])
    best_even = np.full((2, 2), np.inf)
    best_odd = best_even.copy()
    scores = []
    for index, pair in enumerate(markers):
        odds = tuple(log2(float(F(v))) for v in pair)
        weights = base.occupancy_weights(odds, (len(logs)-1)//2)
        two = base.matrix_sum(base.pair_terms(raw_pairs, weights, alpha))
        a = base.regional(two, blocks, weight, even, odds, pairs=pairs)
        b = base.regional(two, blocks, weight, odd, odds, pairs=pairs)
        np.minimum(best_even, a, out=best_even)
        np.minimum(best_odd, b, out=best_odd)
        scores.append((base.profile_transfer(a, b, repetitions), index))
    return dict(transfer=base.profile_transfer(best_even, best_odd, repetitions),
                best_single_marker=min(scores), even=best_even, odd=best_odd)


def finite_json(value):
    """Exact structural log-zero is written as a string, never a JSON NaN."""
    import numpy as np
    if isinstance(value, np.ndarray):
        return finite_json(value.tolist())
    if isinstance(value, dict):
        return {k: finite_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [finite_json(v) for v in value]
    if isinstance(value, (float, np.floating)):
        if np.isneginf(value):
            return '-inf'
        if not np.isfinite(value):
            raise ValueError('Unexpected nonfinite output')
        return float(value)
    return value


def run(source_name, output_name):
    if sys.platform != 'linux':
        raise RuntimeError('Guarded remote execution only')
    import resource
    import signal
    resource.setrlimit(resource.RLIMIT_AS, (512 << 20, 512 << 20))
    resource.setrlimit(resource.RLIMIT_CPU, (25, 30))
    resource.setrlimit(resource.RLIMIT_FSIZE, (LIMIT, LIMIT))
    signal.alarm(35)
    here = Path(__file__).resolve().parent
    source, output = base.previous.inside(source_name), base.previous.inside(output_name)
    if output.exists() or source.stat().st_size > 1 << 20:
        raise ValueError('Exclusive output; input at most1MiB')
    if shutil.disk_usage(here).free-LIMIT < 20 << 30:
        raise RuntimeError('Insufficient space after bounded output')
    raw = source.read_bytes()
    if sha256(raw).hexdigest() != base.SOURCE_SHA:
        raise ValueError('Pinned v34 result required')
    pins = base.source_pins()
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    row = base.validate_source(json.loads(raw))
    started = monotonic()
    import hd_field128_full_rs_gamma_screen as gamma
    z, alpha = float(F(row['witness']['z'])), float(F(row['witness']['alpha']))
    profile = row['largest_profile']
    original = gamma.local_logs(z, base.previous.Q, state=32)
    rows = {}
    prefix = base.previous.prefactor(z, alpha, base.DISTANCE, base.R,
        gamma.old.outer_log_density(base.previous.Q, 'invertible-toeplitz'))
    if abs(prefix-row['log2_outer_prefix']) > 1e-8:
        raise ArithmeticError('Prefix changed')
    for name, logs in [('baseline', original), ('return_zero', ablate_returns(original))]:
        item = evaluate(logs, alpha, row['witness']['markers'], base.R,
                        profile['even_full'], profile['odd_full'])
        item['log2_positive_profile'] = prefix+profile['log2_count_upper']+item['transfer']
        rows[name] = item
    if rows['return_zero']['transfer'] > rows['baseline']['transfer']+1e-8:
        raise ArithmeticError('Ablation increased transfer')
    result = dict(schema='optimistic-return-ablation-v36',
        source=dict(path=source_name, sha256=base.SOURCE_SHA), source_sha256=pins,
        conditional_construction=base.previous.state_metadata(32),
        active_groups=base.R, target_distance=str(base.DISTANCE),
        selected_profile=[profile['even_full'], profile['odd_full']],
        z=row['witness']['z'], alpha=row['witness']['alpha'],
        marker_count=len(row['witness']['markers']),
        count_source='pinned FLOAT profile cap; not recomputed',
        count_log2=profile['log2_count_upper'], prefix_log2=prefix,
        rows=rows, gain_bits=rows['baseline']['transfer']-rows['return_zero']['transfer'],
        saved_alias_profile_reference=prefix+profile['log2_term'],
        actual_probability_upper_bound=False, directed=False,
        all_profiles_replayed=False, coefficient_atlas_built=False,
        optimizer_run=False, whole_code_certificate=False,
        interpretation='Headroom of return-only changes at fixed marks in positive-marker method only; not alias replay or impossibility',
        elapsed_seconds=monotonic()-started)
    if base.source_pins() != {k: v for k, v in pins.items() if k != Path(__file__).name}:
        raise RuntimeError('Frozen sources changed')
    if source.read_bytes() != raw or sha256(Path(__file__).read_bytes()).hexdigest() != pins[Path(__file__).name]:
        raise RuntimeError('Input/source changed')
    data = (json.dumps(finite_json(result), indent=2, sort_keys=True, allow_nan=False)+'\n').encode()
    if len(data) > LIMIT or shutil.disk_usage(here).free-len(data) < 20 << 30:
        raise RuntimeError('Output storage limit')
    with output.open('xb') as stream:
        stream.write(data)
    signal.alarm(0)
    return dict(path=output_name, sha256=sha256(data).hexdigest(), bytes=len(data),
        baseline=rows['baseline']['log2_positive_profile'],
        return_zero=rows['return_zero']['log2_positive_profile'],
        gain_bits=result['gain_bits'], elapsed_seconds=result['elapsed_seconds'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.source, args.output), allow_nan=False), flush=True)
