"""Optimistic J=0 alias diagnostic; NOT a valid distance upper bound."""
import argparse
from contextlib import contextmanager
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, isfinite, log2
from pathlib import Path
import shutil
import sys
from time import monotonic

import hd_field128_gamma_four_packet_positive_v35 as prior

HERE = Path(__file__).resolve().parent
PRIOR_SHA = '291dd2aaf6bd2502eb3c7aad444c2daf731980eeaf5827919d3dd1ab130358d1'
LIMIT, RESERVE = 32 << 10, 20 << 30
EXPECTED_BASELINE = 4069.07598956


def source_pins():
    import hd_field128_gamma_all_pair_screen as pairs
    pins = prior.source_pins()
    for name, digest in pairs.source_pins().items():
        if name in pins and pins[name] != digest:
            raise RuntimeError('Source conflict')
        pins[name] = digest
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    if pins.get(Path(prior.__file__).name) != PRIOR_SHA:
        raise RuntimeError('v35 changed')
    return dict(sorted(pins.items()))


def ablate(logs):
    import numpy as np
    result = np.array(logs, dtype=float, copy=True)
    if result.ndim != 3 or result.shape[1:] != (2, 2):
        raise ValueError('Expected occupancy-indexed2x2')
    result[:, 1, 0] = -np.inf
    return result


class Coefficients:
    def __init__(self, logs, z, alpha, budget):
        import numpy as np
        import hd_field128_gamma_all_pair_screen as pairs
        self.z, self.alpha, self.q = z, alpha, prior.previous.Q
        self.composed = pairs.all_pair_coefficients(logs, alpha, _budget=budget)
        self.single = np.full((4, 33, 33), -np.inf)
        self.packet = self

    def pair_logs(self, odds, **kwargs):
        import hd_field128_gamma_all_pair_screen as pairs
        return pairs.AliasCoefficients.pair_logs(self, odds, **kwargs)


def triangular_atlas(polynomial, odds, baseline_pair, *, grid=128, blocks=64, check=lambda: None):
    """Triangular FFT with baseline scaling."""
    import numpy as np
    import hd_field128_lift_alias_screen as alias
    polynomial = np.asarray(polynomial)
    if (polynomial.ndim != 3 or polynomial.shape[0] != 4
            or not np.all(np.isneginf(polynomial[2]))
            or type(blocks) is not int or blocks < 1 or blocks & (blocks-1)):
        raise ValueError('Expected triangular/power-of-two')
    rho, ratio = alias.perron_logs(baseline_pair)
    vector = np.array((0., ratio))
    similarity = (vector[None, :]-vector[:, None]).ravel()
    values, lost = alias._frequency(polynomial, odds, similarity-rho, grid)
    row_sums = values[:, 0, 0].real.reshape(2, 2).sum(axis=1)
    if not np.all(np.isfinite(values)) or np.any(row_sums < -1e-12) or np.any(row_sums > 1+1e-9):
        raise ArithmeticError('Baseline row cap failed')
    count = blocks
    while count > 1:
        check()
        values = alias._matrix_product(values, values)
        count //= 2
    residues = np.fft.ifftn(values, axes=(1, 2))
    if (not np.all(np.isfinite(residues)) or np.max(abs(residues.imag)) > 1e-10
            or np.min(residues.real) < -1e-10 or np.any(residues[2] != 0)):
        raise ArithmeticError('FFT diagnostics failed')
    reliable = residues.real >= 1e-10
    logs = np.log2(np.maximum(residues.real, 0)+1e-12)
    logs += blocks*rho-similarity[:, None, None]
    degree = np.arange(polynomial.shape[1])[:, None]*odds[0]
    degree = degree+np.arange(polynomial.shape[2])[None, :]*odds[1]
    point = np.logaddexp2.reduce((polynomial+degree).reshape(4, -1), axis=1).reshape(2, 2)
    positive = prior.power(point, blocks).ravel()
    # Algebraic zero.
    logs[2], reliable[2], positive[2] = -np.inf, True, -np.inf
    check()
    return alias.AliasTable(logs, reliable, positive, dict(input_underflows=lost))


def selected(table, odds, *, full=2026):
    import numpy as np
    single, pairs = 4096-2*full, 4096
    norm = log2(comb(pairs, single)*comb(pairs-single, full))
    grid = table.logs.shape[1]
    reliable = table.reliable[:, single % grid, full % grid]
    numerator = np.where(reliable, table.logs[:, single % grid, full % grid], table.positive_numerator)
    return (numerator-norm-single*odds[0]-full*odds[1]).reshape(2, 2), reliable


@contextmanager
def remote_limits():
    if sys.platform != 'linux':
        raise RuntimeError('Production requires Linux')
    import resource
    import signal
    resource.setrlimit(resource.RLIMIT_AS, (1 << 30, 1 << 30))
    resource.setrlimit(resource.RLIMIT_CPU, (55, 60))
    resource.setrlimit(resource.RLIMIT_FSIZE, (LIMIT, LIMIT))
    signal.alarm(60)
    try:
        yield
    finally:
        signal.alarm(0)


def run(source_name, output_name):
    source, output = prior.previous.inside(source_name), prior.previous.inside(output_name)
    if output.exists() or source.stat().st_size > 1 << 20:
        raise ValueError('Output exists / oversized source')
    if shutil.disk_usage(HERE).free-LIMIT < RESERVE:
        raise RuntimeError('20GiB planned reserve required')
    raw, pins = source.read_bytes(), source_pins()
    if sha256(raw).hexdigest() != prior.SOURCE_SHA:
        raise ValueError('Pinned v34 result required')
    row = prior.validate_source(json.loads(raw))
    witness, profile = row['witness'], row['largest_profile']
    if (len(witness['markers']) != 79 or witness['grid'] != 128
            or (profile['even_full'], profile['odd_full']) != (2026, 2026)):
        raise ValueError('Expected79 markers/grid128/profile2026,2026')
    import numpy as np
    import hd_field128_full_rs_gamma_screen as gamma
    import hd_field128_lift_alias_screen as alias
    started, budget = monotonic(), gamma.old.Budget(50.)
    z, alpha = float(F(witness['z'])), float(F(witness['alpha']))
    logs = gamma.local_logs(z, prior.previous.Q, state=32)
    baseline = Coefficients(logs, z, alpha, budget)
    zero = Coefficients(ablate(logs), z, alpha, budget)
    best = [np.full((2, 2), np.inf), np.full((2, 2), np.inf)]
    reliability = [np.zeros((2, 2), bool), np.zeros((2, 2), bool)]
    details = []
    for index, literal in enumerate(witness['markers']):
        budget.check()
        odds = tuple(log2(float(F(x))) for x in literal)
        base = alias.alias_table(z, alpha, odds, _coefficients=baseline, _budget=budget)
        tri = triangular_atlas(zero.composed, odds, baseline.pair_logs(odds, _budget=budget), check=budget.check)
        for position, table in enumerate((base, tri)):
            value, reliable = selected(table, odds)
            improve = value < best[position]
            reliability[position][improve] = reliable.reshape(2, 2)[improve]
            best[position] = np.minimum(best[position], value)
        details.append([index, base.diagnostic['normalized_input_underflow_terms'], tri.diagnostic['input_underflows']])
    prefix = prior.previous.prefactor(z, alpha, prior.DISTANCE, 4096,
        gamma.old.outer_log_density(prior.previous.Q, 'invertible-toeplitz'))
    if abs(prefix-row['log2_outer_prefix']) > 1e-7:
        raise ArithmeticError('Original-r prefix mismatch')
    scores = [prefix+profile['log2_count_upper']+prior.profile_transfer(value, value) for value in best]
    if not all(map(isfinite, scores)) or abs(scores[0]-EXPECTED_BASELINE) > 1e-4:
        raise ArithmeticError('Matched baseline mismatch')
    result = dict(schema='optimistic-return-alias-ablation-v37', source=dict(path=source_name, sha256=prior.SOURCE_SHA),
        source_sha256=pins, conditional_construction=prior.previous.state_metadata(32),
        active_groups=4096, selected_profile=[2026, 2026], target_distance='9/20',
        witness=witness, log2_profile_count=profile['log2_count_upper'], log2_outer_prefix=prefix,
        log2_baseline_profile=scores[0], log2_optimistic_J0_profile=scores[1], gain_bits=scores[0]-scores[1],
        selected_reliable=[x.tolist() for x in reliability],
        marker_underflow_counts=details, structural_zero_entry=[1, 0], all79_markers_replayed=True,
        all_profiles_replayed=False, arithmetic='FLOAT; positive fallback; original normalization',
        optimistic_ablation=True, valid_distance_upper_bound=False, directed=False,
        whole_code_certificate=False,
        scope='Fixed witness only; no code-failure or optimized-impossibility claim',
        elapsed_seconds=monotonic()-started)
    if source.read_bytes() != raw or source_pins() != pins:
        raise RuntimeError('Source closure/input changed')
    data = (json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n').encode()
    if len(data) > LIMIT or shutil.disk_usage(HERE).free-len(data) < RESERVE:
        raise RuntimeError('32KiB output or20GiB reserve exceeded')
    with output.open('xb') as stream:
        stream.write(data)
    if output.read_bytes() != data or source.read_bytes() != raw or source_pins() != pins:
        raise RuntimeError('Closing pin check failed')
    return dict(path=output_name, sha256=sha256(data).hexdigest(), scores=scores, directed=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    with remote_limits():
        print(json.dumps(run(args.source, args.output), allow_nan=False), flush=True)
