"""V34 serialization repair of the frozen v33 FLOAT state ladder.

State32 uses one independent global gamma uniform outside GF(2^64). Every
local cap is conditional on all active principal minors being nonsingular.
Its setup-failure bound is recorded separately, never charged per profile.
Frozen coefficient, alias, midpoint-count and transfer arithmetic is reused.
"""
import argparse
import ast
from contextlib import contextmanager
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isfinite, log2
from pathlib import Path, PurePosixPath
import shutil
import sys
from time import monotonic

import hd_field128_gamma_state_ladder_v33 as frozen_v33

FROZEN_V33_SHA256 = '17d9b935d2cdad4219c59fb7279053899268b3ee94f44dc9cac546b43765efe6'
HERE = Path(__file__).resolve().parent
Q, N, P = 1 << 128, 1 << 21, 4096
RESERVE, OUTPUT_LIMIT = 20 << 30, 1 << 20
STATES = (30, 31, 32)
NOTES = ('HD_FIELD128_GAMMA_RETURN.md', 'HD_FIELD128_GAMMA31.md',
         'HD_FIELD128_RANDOM_GAMMA32.md', 'HD_FIELD128_MIDPOINT_COUNTS.md')
FIXED30 = 'toeplitz-fullfield-full-rs64-30-gamma16-feedback-alternating-swapped-pairs-k20-rate-half-v1'


def state_metadata(state):
    if type(state) is not int or state not in STATES:
        raise ValueError('State must be30,31,32')
    conditional = state == 32
    bad = F(16*((1 << 62)-1), Q-(1 << 64)) if conditional else F(0)
    return dict(state_dimension=state, field_order=Q, inner_width=64,
        construction=(f'toeplitz-fullfield-full-rs64-{state}-'
            + ('random-global-degree16-gamma-conditional' if conditional else 'gamma16-feedback')
            + '-alternating-swapped-pairs-k20-rate-half-v1'),
        conditional_on_global_principal_nonsingularity=conditional,
        gamma_rule=('one uniform global element of F\\GF(2^64), independent of outer/route/refresh'
                    if conditional else 'fixed degree16 element over GF256; towerY permitted'),
        setup_failure_probability_upper_rational=str(bad),
        setup_failure_charged_once_at_whole_setup=conditional,
        setup_failure_included_in_selected_FLOAT_score=False,
        principal_nonsingularity=('hypothesis G(gamma); random-global setup bound'
                                 if conditional else 'deterministic degree-halving lemma'),
        emission_minimum=65-state, state_nodes=[128, 127+state], output_nodes=[0, 63],
        refresh_input_nodes=[0, state-1], refresh_output_nodes=[state, 2*state-1],
        feedback='gamma A transpose', refresh='D2 H_state D1; independent nonzero full-field diagonals',
        CA_zero=True, initial_state='zero', final_state='unrestricted',
        gamma_sampled_or_verified_in_this_pilot=False,
        pair_swaps=True, state_persists_across_regions=True,
        full_field_random_coefficients=True, outer_measure='invertible-toeplitz',
        counts_outside_alpha=True, all_pair_before_alpha=True,
        directed=False, selected_weight_certificate=False, whole_code_certificate=False)


def source_pins():
    pending, pins = [Path(__file__).name], {}
    while pending:
        name = pending.pop()
        if name in pins:
            continue
        raw = (HERE/name).read_bytes()
        pins[name] = sha256(raw).hexdigest()
        for node in ast.walk(ast.parse(raw)):
            modules = ([x.name for x in node.names] if isinstance(node, ast.Import)
                       else [node.module] if isinstance(node, ast.ImportFrom) and node.module else [])
            for module in modules:
                child = module.split('.')[0]+'.py'
                if (HERE/child).is_file() and child not in pins:
                    pending.append(child)
    for name in NOTES:
        pins['../'+name] = sha256((HERE.parent/name).read_bytes()).hexdigest()
    if pins.get(Path(frozen_v33.__file__).name) != FROZEN_V33_SHA256:
        raise RuntimeError('Frozen v33 source changed')
    return dict(sorted(pins.items()))


def check_pins(pins):
    if source_pins() != pins:
        raise RuntimeError('Pilot source/proof-note closure changed')


def inside(name):
    if not isinstance(name, str) or not name or '\\' in name or ':' in name:
        raise ValueError('Use a theory-relative POSIX path')
    relative = PurePosixPath(name)
    if relative.is_absolute() or '..' in relative.parts or relative == PurePosixPath('.'):
        raise ValueError('Path outside theory')
    result = (HERE/name).resolve()
    if not result.is_relative_to(HERE) or result == HERE:
        raise ValueError('Resolved path outside theory')
    return result


def validate_options(weight, states, distance, mode, search_state, evaluations, seconds):
    if type(weight) is not int or not 2 <= weight < 8192:
        raise ValueError('Original weight2..8191 required')
    if (not 1 <= len(states) <= 3 or len(set(states)) != len(states)
            or any(type(s) is not int or s not in STATES for s in states)):
        raise ValueError('Distinct states from30,31,32 required')
    if not F(2, 5) < distance < F(1, 2):
        raise ValueError('Target strictly between40% and50% required')
    if mode not in ('compare', 'search') or type(search_state) is not int or search_state not in states:
        raise ValueError('Search state must be among compared states')
    if type(evaluations) is not int or not 1 <= evaluations <= 24:
        raise ValueError('At most24 complete search objectives')
    if not isfinite(seconds) or not 0 < seconds <= 150:
        raise ValueError('At most150 cooperative seconds')


def project_seed(saved, weight, distance):
    row = saved.get('best', saved)
    if not isinstance(row, dict) or row.get('construction') != FIXED30 or row.get('state_dimension') != 30:
        raise ValueError('Expected a pinned gamma30 literal source')
    if row.get('field_order', Q) != Q:
        raise ValueError('Wrong source field')
    witness = row['witness']
    z, alpha = F(witness['z']), F(witness['alpha'])
    if (not 0 < z < 1 or not 0 < alpha <= 1 or F(float(z)) != z or F(float(alpha)) != alpha
            or not 2.**-1022 <= float(z) <= 2.**-1e-8 or not 2.**-32 <= float(alpha) <= 1):
        raise ValueError('Marks must be exact binary64 values within search chart')
    key = 'exact_positive_marker_proposals' if 'exact_positive_marker_proposals' in witness else 'markers'
    markers = witness.get(key)
    if (not isinstance(markers, list) or not 1 <= len(markers) <= 128
            or any(not isinstance(pair, (list, tuple)) or len(pair) != 2 for pair in markers)):
        raise ValueError('Expected1..128 explicit positive marker pairs')
    markers = [tuple(F(v) for v in pair) for pair in markers]
    if any(v <= 0 or not isfinite(float(v)) or float(v) <= 0 for pair in markers for v in pair):
        raise ValueError('Positive finite binary64 marker conversions required')
    grid = witness.get('grid', 128)
    if type(grid) is not int or grid not in (64, 128, 256):
        raise ValueError('Unsupported alias grid')
    return dict(active_groups=weight, target_distance=str(distance), field_order=Q,
        witness=dict(z=str(z), alpha=str(alpha), markers=[list(map(str, p)) for p in markers], grid=grid),
        source_marker_key=key, literal_only=True, old_numerical_bounds_reused=False,
        source_weight=row.get('active_groups'), source_distance=row.get('target_distance'))


def prefactor(z, alpha, distance, weight, density):
    # Exactly the generic bank expression: projective divisor once, count outside alpha.
    from hd_field128_lift_alias_screen import log1mexp2
    return alpha*(density*weight-128-log1mexp2(-128)-int(distance*N)*log2(z))


@contextmanager
def remote_limits():
    if sys.platform != 'linux':
        raise RuntimeError('Production run requires guarded Linux execution')
    import resource
    import signal
    resource.setrlimit(resource.RLIMIT_AS, (4 << 30, 4 << 30))
    resource.setrlimit(resource.RLIMIT_CPU, (160, 165))
    resource.setrlimit(resource.RLIMIT_FSIZE, (OUTPUT_LIMIT, OUTPUT_LIMIT))
    signal.alarm(180)
    try:
        yield
    finally:
        signal.alarm(0)


def json_ready(value):
    """Convert scalar diagnostics only; arrays remain rejected by json.dumps."""
    import numpy as np
    def convert(item):
        if isinstance(item, np.generic):
            return item.item()
        if isinstance(item, dict):
            return {convert(key): convert(child) for key, child in item.items()}
        if isinstance(item, (list, tuple)):
            return [convert(child) for child in item]
        return item
    return convert(value)


class Writer:
    def __init__(self, directory):
        if shutil.disk_usage(HERE).free-OUTPUT_LIMIT < RESERVE:
            raise RuntimeError('20GiB reserve after entire output budget required')
        directory.mkdir(exist_ok=False)
        self.directory, self.written, self.pins = directory, 0, {}

    def write(self, name, value):
        if Path(name).name != name:
            raise ValueError('Output basename required')
        data = (json.dumps(json_ready(value), sort_keys=True, indent=2, allow_nan=False)+'\n').encode()
        if self.written+len(data) > OUTPUT_LIMIT or shutil.disk_usage(self.directory).free-len(data) < RESERVE:
            raise RuntimeError('Output budget or20GiB planned-write reserve')
        path = self.directory/name
        self.written += len(data)
        with path.open('xb') as stream:
            stream.write(data)
        self.pins[name] = sha256(data).hexdigest()
        return dict(path=path.relative_to(HERE).as_posix(), sha256=self.pins[name], bytes=len(data))

    def check(self):
        if any(sha256((self.directory/name).read_bytes()).hexdigest() != digest for name, digest in self.pins.items()):
            raise RuntimeError('Exclusive output bytes changed')


def _run(source_name, expected_sha256, weight, output_name, *, states=(30, 31, 32),
         distance=F(9, 20), mode='compare', search_state=32, max_evaluations=24, max_seconds=150.):
    distance, states = F(distance), tuple(states)
    validate_options(weight, states, distance, mode, search_state, max_evaluations, max_seconds)
    path, output = inside(source_name), inside(output_name)
    if path.stat().st_size > OUTPUT_LIMIT or output.exists():
        raise ValueError('Source cap1MiB / exclusive new output required')
    raw, pins = path.read_bytes(), source_pins()
    if sha256(raw).hexdigest() != expected_sha256:
        raise ValueError('Source SHA256 mismatch')
    seed = project_seed(json.loads(raw), weight, distance)
    writer = Writer(output)
    def check():
        check_pins(pins)
        if path.read_bytes() != raw:
            raise RuntimeError('Literal source changed')
        writer.check()
    try:
        writer.write('seed.json', dict(seed, source=dict(path=source_name, sha256=expected_sha256)))
        import numpy as np
        import scipy
        from scipy.optimize import minimize
        import hd_field128_full_rs_gamma_screen as gamma
        import hd_field128_gamma_all_pair_screen as pairs
        import hd_field128_gamma_midpoint_counts as counts_module
        started, budget = monotonic(), pairs.old.Budget(max_seconds)
        full = list(pairs.alias.full_counts(P, weight))
        counts = counts_module.MidpointCounts(P, weight, max_singletons=64, budget=budget)
        witness = seed['witness']
        z0, alpha0 = float(F(witness['z'])), float(F(witness['alpha']))
        odds = [tuple(log2(float(F(v))) for v in pair) for pair in witness['markers']]
        history, cache, best = [], {}, {}
        class Coefficients:
            pair_logs = pairs.AliasCoefficients.pair_logs
            def __init__(self, z, alpha, state):
                self.z, self.alpha, self.q, self.state = z, alpha, Q, state
                # Generic formula's hypothesis is G(gamma); no public whitelist is modified.
                logs = gamma.local_logs(z, Q, n=64, state=state)
                self.composed = pairs.all_pair_coefficients(logs, alpha, _budget=budget)
                self.single = np.full((4, 33, 33), pairs.NEG)
                self.packet = self
        def objective(state, z, alpha):
            budget.check()
            key = state, z, alpha
            if key in cache:
                return cache[key]
            coefficient = Coefficients(z, alpha, state)
            bank = pairs.alias.AliasBankScreen(z, alpha, odds, 0., q=Q, distance=distance,
                grid=witness['grid'], floor_policy='positive-fallback', max_table_bytes=384 << 20,
                _coefficients=coefficient, _budget=budget)
            regions = bank.regional(weight)
            if regions.shape != (len(full), 2, 2):
                raise ArithmeticError('Original-r regional domain mismatch')
            row = counts_module.support_sum(P, weight, regions, counts=counts, budget=budget)
            if row['profile_pairs_evaluated'] != len(full)**2:
                raise ArithmeticError('Incomplete profile domain')
            prefix = bank.constant+bank.alpha*bank.density*weight
            expected = prefactor(z, alpha, distance, weight, bank.density)
            if abs(prefix-expected) > 1e-8*max(1., abs(prefix)):
                raise ArithmeticError('Prefactor normalization mismatch')
            score = prefix+row['log2_support_sum_screen']
            if not isfinite(score):
                raise ArithmeticError('Nonfinite full objective')
            row.update(state_metadata(state), active_groups=weight, target_distance=str(distance),
                forbidden_weight=int(distance*N), log2_outer_prefix=prefix, log2_upper_screen=score,
                feasible_full_range=[full[0], full[-1]], count_complement_only=True,
                original_r_used_for_coefficients_and_prefix=True, complete_selected_weight_screen=True,
                witness=dict(z=str(F(z)), alpha=str(F(alpha)), markers=witness['markers'], grid=witness['grid']),
                regional_alias_diagnostics=dict(bank.regional_diagnostics),
                regional_logs_sha256=sha256(regions.tobytes()).hexdigest(),
                coefficient_table_sha256=sha256(coefficient.composed.tobytes()).hexdigest(),
                common_potential_used=False, old_numerical_bounds_reused=False,
                elapsed_seconds=monotonic()-started)
            cache[key] = score
            history.append(dict(state=state, z=str(F(z)), alpha=str(F(alpha)), score=score,
                                largest_profile=row['largest_profile']))
            if state not in best or score < best[state]['log2_upper_screen']:
                best[state] = row
            del regions, bank, coefficient
            return score
        status, error, optimizer = 'complete', None, None
        comparisons = []
        try:
            for state in states:
                objective(state, z0, alpha0)
                comparisons.append(best[state])
            if mode == 'search':
                initial = np.array((-log2(z0), log2(alpha0)))
                simplex = np.array((initial, (initial[0]+.5 if initial[0] <= 1021.5 else initial[0]-.5, initial[1]),
                                    (initial[0], initial[1]+.25 if initial[1] <= -.25 else initial[1]-.25)))
                def target(point):
                    ell, exponent = map(float, point)
                    if not 1e-8 <= ell <= 1022 or not -32 <= exponent <= 0:
                        budget.check()
                        return 1e300
                    return objective(search_state, z0 if ell == initial[0] else 2.**-ell,
                                     alpha0 if exponent == initial[1] else 2.**exponent)
                fitted = minimize(target, initial, method='Nelder-Mead', options=dict(maxfev=max_evaluations,
                    initial_simplex=simplex, xatol=1e-7, fatol=1e-5))
                optimizer = dict(success=bool(fitted.success), message=str(fitted.message), calls=int(fitted.nfev))
        except (TimeoutError, MemoryError, ArithmeticError) as exception:
            status, error = 'partial', dict(type=type(exception).__name__, message=str(exception))
        result = dict(schema='gamma-state-ladder-conditional-FLOAT-v34', source_sha256=pins,
            source=dict(path=source_name, sha256=expected_sha256), active_groups=weight,
            target_distance=str(distance), requested_states=list(states), comparisons=comparisons,
            completed_comparison_states=[row['state_dimension'] for row in comparisons],
            best_by_state={str(s): row for s, row in best.items()}, history=history, status=status,
            failure=error, optimizer=optimizer, search_state=search_state, mode=mode,
            max_evaluations=max_evaluations, max_seconds=max_seconds, counts_built=1,
            fixed_positive_marker_literals=witness['markers'],
            runtime=dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__,
                         platform=sys.platform, observation_only=True),
            limits=dict(wall_seconds=180, address_space_bytes=4 << 30, total_json_bytes=OUTPUT_LIMIT,
                        minimum_free_after_write=RESERVE),
            directed=False, whole_code_certificate=False, selected_weight_certificate=False,
            scope='FLOAT selected-weight diagnostics; state32 conditional; setup risk separate; no timing claim')
        check()
        receipt = writer.write('result.json', result)
        return dict(receipt=receipt, status=status, best_by_state={str(s): row['log2_upper_screen'] for s, row in best.items()},
                    completed_trials=len(history), directed=False, whole_code_certificate=False)
    finally:
        check()


def run(*args, **kwargs):
    with remote_limits():
        return _run(*args, **kwargs)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--sha256', required=True)
    parser.add_argument('--r', required=True, type=int)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--states', nargs='+', type=int, default=list(STATES))
    parser.add_argument('--distance', default='9/20')
    parser.add_argument('--mode', choices=('compare', 'search'), default='compare')
    parser.add_argument('--search-state', type=int, default=32)
    parser.add_argument('--max-evaluations', type=int, default=24)
    parser.add_argument('--max-seconds', type=float, default=150.)
    args = parser.parse_args()
    print(json.dumps(run(args.source, args.sha256, args.r, args.output_dir, states=args.states,
        distance=F(args.distance), mode=args.mode, search_state=args.search_state,
        max_evaluations=args.max_evaluations, max_seconds=args.max_seconds), allow_nan=False), flush=True)


if __name__ == '__main__':
    main()
