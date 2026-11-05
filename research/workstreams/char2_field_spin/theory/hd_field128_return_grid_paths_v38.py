"""Fixed-profile optimistic J=0 grid/path diagnosis, never a certificate.
Rebuild v37 at128; only its winning markers continue at256/512. No search.
"""
import argparse
from contextlib import contextmanager
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isfinite, log2
from pathlib import Path
import shutil
import sys
from time import monotonic

import hd_field128_return_alias_ablation_v37 as old

HERE = Path(__file__).resolve().parent
OLD_SHA = '86f15a8be2e71a8015443a5fa6f30cfda47d227793ec535c29603a898e48b98c'
EXPECTED = (4069.075989561272, 3561.912788278656)
LIMIT, RESERVE = 32 << 10, 20 << 30


def source_pins():
    pins = old.source_pins()
    if pins.get(Path(old.__file__).name) != OLD_SHA:
        raise RuntimeError('Frozen v37 changed')
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


class Best:
    """First strict winning marker; no mixing of different z/alpha."""
    def __init__(self):
        import numpy as np
        self.logs = np.full((2, 2), np.inf)
        self.winners = np.full((2, 2), -1, dtype=int)
        self.reliable = np.zeros((2, 2), dtype=bool)

    def add(self, index, logs, reliable):
        import numpy as np
        logs, reliable = np.asarray(logs), np.asarray(reliable).reshape(2, 2)
        if logs.shape != (2, 2) or np.any(np.isnan(logs)) or np.any(np.isposinf(logs)):
            raise ValueError('Invalid selected regional matrix')
        better = logs < self.logs
        self.winners[better], self.reliable[better] = index, reliable[better]
        np.minimum(self.logs, logs, out=self.logs)

    def winning_subset(self):
        import numpy as np
        if not np.isneginf(self.logs[1, 0]):
            raise ValueError('Expected algebraic return zero')
        indices = sorted({int(self.winners[i, j]) for i, j in ((0, 0), (0, 1), (1, 1))})
        if not 1 <= len(indices) <= 8 or indices[0] < 0:
            raise ValueError('Invalid winning-marker subset')
        return indices


def paths(logs, repetitions=256):
    import numpy as np
    logs = np.asarray(logs)
    if (logs.shape != (2, 2) or not np.isneginf(logs[1, 0])
            or not np.isfinite(logs[[0, 1], [0, 1]]).all()
            or np.isnan(logs[0, 1]) or np.isposinf(logs[0, 1])
            or type(repetitions) is not int or repetitions < 1):
        raise ValueError('Finite diagonal/triangular matrix and positive exponent required')
    powered = old.prior.power(logs, repetitions)
    return dict(log2_Rn00=float(powered[0, 0]), log2_Rn01=float(powered[0, 1]),
                log2_Rn_row0=float(np.logaddexp2(powered[0, 0], powered[0, 1])), repetitions=repetitions)


def snapshot(best, prefix_count):
    import numpy as np
    result = paths(best.logs)
    for name, key in (('all_quiet_score', 'log2_Rn00'), ('one_birth_score', 'log2_Rn01'),
                      ('total_score', 'log2_Rn_row0')):
        result[name] = prefix_count+result[key]
    result.update(regional_log2=[[None if np.isneginf(x) else float(x) for x in row] for row in best.logs],
        regional_log2_binary64_literals=[[None if np.isneginf(x) else str(F(float(x))) for x in row] for row in best.logs],
        winning_marker_indices=[[None if (i, j) == (1, 0) else int(best.winners[i, j]) for j in range(2)] for i in range(2)],
        selected_reliable=best.reliable.tolist(), structural_zero_entry=[1, 0])
    return result


@contextmanager
def remote_limits():
    if sys.platform != 'linux':
        raise RuntimeError('Production requires Linux')
    import resource
    import signal
    resource.setrlimit(resource.RLIMIT_AS, (1 << 30, 1 << 30))
    resource.setrlimit(resource.RLIMIT_CPU, (90, 90))
    resource.setrlimit(resource.RLIMIT_FSIZE, (LIMIT, LIMIT))
    signal.alarm(96)
    try:
        yield
    finally:
        signal.alarm(0)


def run(source_name, output_name):
    source, output = old.prior.previous.inside(source_name), old.prior.previous.inside(output_name)
    if output.exists() or source.stat().st_size > 1 << 20:
        raise ValueError('Exclusive output / bounded source required')
    if shutil.disk_usage(HERE).free-LIMIT < RESERVE:
        raise RuntimeError('20GiB planned reserve required')
    raw, pins = source.read_bytes(), source_pins()
    if sha256(raw).hexdigest() != old.prior.SOURCE_SHA:
        raise ValueError('Pinned v34 source required')
    row = old.prior.validate_source(json.loads(raw))
    witness, profile = row['witness'], row['largest_profile']
    if (len(witness['markers']) != 79 or witness['grid'] != 128
            or (profile['even_full'], profile['odd_full']) != (2026, 2026)):
        raise ValueError('Expected79 markers/grid128/profile2026,2026')
    import numpy as np
    import hd_field128_full_rs_gamma_screen as gamma
    import hd_field128_lift_alias_screen as alias
    started, budget = monotonic(), gamma.old.Budget(88.)
    z, alpha = float(F(witness['z'])), float(F(witness['alpha']))
    local = gamma.local_logs(z, old.prior.previous.Q, state=32)
    base = old.Coefficients(local, z, alpha, budget)
    zero = old.Coefficients(old.ablate(local), z, alpha, budget)
    prefix = old.prior.previous.prefactor(z, alpha, F(9, 20), 4096,
        gamma.old.outer_log_density(old.prior.previous.Q, 'invertible-toeplitz'))
    if abs(prefix-row['log2_outer_prefix']) > 1e-7:
        raise ArithmeticError('Original-r prefix mismatch')
    offset = prefix+profile['log2_count_upper']
    baseline, initial = Best(), Best()
    cached, underflows = {}, []
    def triangle(index, grid):
        odds = tuple(log2(float(F(x))) for x in witness['markers'][index])
        table = old.triangular_atlas(zero.composed, odds, base.pair_logs(odds, _budget=budget),
                                    grid=grid, check=budget.check)
        value, reliable = old.selected(table, odds)
        return value, reliable, table.diagnostic['input_underflows']
    for index, literal in enumerate(witness['markers']):
        budget.check()
        odds = tuple(log2(float(F(x))) for x in literal)
        table = alias.alias_table(z, alpha, odds, _coefficients=base, _budget=budget)
        value, reliable = old.selected(table, odds)
        baseline.add(index, value, reliable)
        lost = table.diagnostic['normalized_input_underflow_terms']
        del table
        value, reliable, triangular_lost = triangle(index, 128)
        initial.add(index, value, reliable)
        cached[index] = (value, reliable)
        underflows.append([index, lost, triangular_lost])
    baseline_score = offset+old.prior.profile_transfer(baseline.logs, baseline.logs)
    first = snapshot(initial, offset)
    if (not isfinite(baseline_score) or abs(baseline_score-EXPECTED[0]) > 1e-4
            or not isfinite(first['total_score']) or abs(first['total_score']-EXPECTED[1]) > 1e-4):
        raise ArithmeticError('V37 baseline/J0 reproduction failed')
    indices = initial.winning_subset()
    control = Best()
    for index in indices:
        control.add(index, *cached[index])
    if not np.array_equal(control.logs, initial.logs):
        raise ArithmeticError('Winning-subset128 control differs')
    results = [dict(grid=128, marker_indices=indices, reused128_selected_entries=True, **snapshot(control, offset))]
    for grid in (256, 512):
        selected, diagnostics = Best(), []
        for index in indices:
            budget.check()
            value, reliable, lost = triangle(index, grid)
            selected.add(index, value, reliable)
            diagnostics.append([index, lost])
        result = snapshot(selected, offset)
        if not isfinite(result['total_score']):
            raise ArithmeticError('Nonfinite grid score')
        results.append(dict(grid=grid, marker_indices=indices, input_underflows=diagnostics, **result))
    budget.check()
    result = dict(schema='optimistic-return-grid-paths-v38', source=dict(path=source_name, sha256=old.prior.SOURCE_SHA),
        source_sha256=pins, conditional_construction=old.prior.previous.state_metadata(32), witness=witness,
        active_groups=4096, selected_profile=[2026, 2026], target_distance='9/20',
        log2_outer_prefix=prefix, log2_profile_count=profile['log2_count_upper'],
        matched_baseline_score=baseline_score, baseline_winning_indices=baseline.winners.tolist(),
        baseline_selected_reliable=baseline.reliable.tolist(), all79_grid128=first,
        grid128_input_underflows=underflows, subset_grids=results,
        matrix_values='FLOAT log2; rational literals encode binary64 only; None denotes exact log2(0)',
        paths='R^256[00] allquiet; [01] exactly one birth; original prefix/count once per contribution',
        no_monotonic_grid_gain_assumed=True, all_profiles_replayed=False, optimizer_run=False,
        optimistic_ablation=True, valid_distance_upper_bound=False, directed=False,
        selected_weight_certificate=False, whole_code_certificate=False, setup_risk_added=False,
        limits=dict(address_space=1 << 30, cpu_seconds=90, wall_seconds=96, output_bytes=LIMIT),
        scope='Fixed markers/profile only; neither code failure nor proof impossibility', elapsed_seconds=monotonic()-started)
    if source.read_bytes() != raw or source_pins() != pins:
        raise RuntimeError('Input/source pins changed')
    data = (json.dumps(result, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n').encode()
    if len(data) > LIMIT or shutil.disk_usage(HERE).free-len(data) < RESERVE:
        raise RuntimeError('32KiB output or20GiB planned reserve exceeded')
    with output.open('xb') as stream:
        stream.write(data)
    if output.read_bytes() != data or source.read_bytes() != raw or source_pins() != pins:
        raise RuntimeError('Closing pin check failed')
    return dict(path=output_name, sha256=sha256(data).hexdigest(), scores=[r['total_score'] for r in results])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    with remote_limits():
        print(json.dumps(run(args.source, args.output), allow_nan=False), flush=True)
