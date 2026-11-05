"""One-profile FLOAT input-weight tilt; unchanged conditional gamma32 code.

Each active four-RS[64,32] outer group has >=33 nonzero columns. Its
invertible column mixers preserve zero/nonzero columns. Consequently its
routed symbol weight W is >=33r. For eta>=1, retain this necessary condition
by 1_allowed <= eta**(W-33r) in the outer counting comparison BEFORE alpha.
The uniform product input law is not changed or renormalized. The prefix
loses 33*r*alpha*log2(eta); projective divisor and support counts stay once.

For an occupied set S of size j, X is uniform on F**S independently of the
incoming weighted-state measure. Quiet/birth are K(eta*z), T(eta*z)-K(eta*z).
For uniform nonzero Q, put U=AQ. Averaging each active X coordinate gives
u0=(1+(q-1)*eta*z)/q when U_i=0, or
u1=(z+eta+(q-2)*eta*z)/q otherwise. The exact nonzero MDS support enumerator
therefore gives V_eta below. Its support counts are exact integers before
taking logarithms; the live sum uses only positive terms.

The incoming cone represents zero mass Z and nonzero atoms <=W/(q**s-1).
Return fixes Q=gamma*C*X, so use J_eta <= min(eta**j*J_old,(q**s-1)*V_eta).
Do NOT substitute eta*z in J_old: the return map need not preserve input
weight. The old refresh density multiplies only birth/live. The j=0 live
distance cap stays unchanged. Condition on the same global good-gamma
event, whose failure probability is added once only in a future full proof.

No directed arithmetic, complete-profile sum, or distance claim is produced.
Optional J=0 is an optimistic invalid upper-bound ablation, labeled separately.
"""
import argparse
from contextlib import contextmanager
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from math import comb, isfinite, log2
from pathlib import Path
import shutil
import sys
from time import monotonic

import hd_field128_return_alias_ablation_v37 as old

HERE = Path(__file__).resolve().parent
OLD_SHA = '86f15a8be2e71a8015443a5fa6f30cfda47d227793ec535c29603a898e48b98c'
Q, R, LIMIT, RESERVE = 1 << 128, 4096, 32 << 10, 20 << 30
DEFAULT_ETAS = (0, 1, 4, 8, 16, 32, 48, 64)


def source_pins():
    pins = old.source_pins()
    if pins.get(Path(old.__file__).name) != OLD_SHA:
        raise RuntimeError('Frozen v37 changed')
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    note = '../HD_FIELD128_LIFT_PROOF.md'
    digest = sha256((HERE/note).read_bytes()).hexdigest()
    if note in pins and pins[note] != digest:
        raise RuntimeError('Conflicting proof-note pin')
    pins[note] = digest
    return dict(sorted(pins.items()))


def validate_local(z, eta, q, n, state):
    import hd_field128_full_rs_gamma_screen as gamma
    gamma.frozen.parameters(q, n, state)
    if not 0 < z < 1 or not 1 <= eta or not eta*z < 1:
        raise ValueError('Require 0<z<1, eta>=1, eta*z<1')


@lru_cache(8)
def support_counts(q, n, state):
    """Nonzero MDS words on ONE specified exact support, including a_0=0."""
    import hd_field128_full_rs_gamma_screen as gamma
    gamma.frozen.parameters(q, n, state)
    distance, values = n-state+1, [0]*(n+1)
    for h in range(distance, n+1):
        values[h] = (q-1)*sum((-1)**i*comb(h-1, i)*q**(h-distance-i)
                              for i in range(h-distance+1))
    if any(v < 0 for v in values) or sum(comb(n, h)*v for h, v in enumerate(values)) != q**state-1:
        raise ArithmeticError('Exact MDS support mass mismatch')
    return tuple(values)


def exact_live(j, z, eta, *, q=Q, n=64, state=32):
    """Tiny exact oracle; production uses log_live, not huge rational powers."""
    if isinstance(z, float) or isinstance(eta, float):
        raise ValueError('Exact oracle requires rational marks')
    z, eta = F(z), F(eta)
    validate_local(z, eta, q, n, state)
    if type(j) is not int or not 0 <= j <= n:
        raise ValueError('Invalid occupancy')
    counts = support_counts(q, n, state)
    u0, u1 = (1+(q-1)*eta*z)/q, (z+eta+(q-2)*eta*z)/q
    return sum((comb(j, b)*comb(n-j, c)*counts[b+c]
                *u1**b*u0**(j-b)*z**c
                for b in range(j+1) for c in range(n-j+1)), F(0))/(q**state-1)


def exact_packet(j, z, eta, *, q=Q, n=64, state=32):
    """Tiny conditional-gamma oracle for the four new nonnegative caps."""
    import hd_field128_full_rs_gamma_screen as gamma
    if isinstance(z, float) or isinstance(eta, float):
        raise ValueError('Exact oracle requires rational marks')
    z, eta = F(z), F(eta)
    validate_local(z, eta, q, n, state)
    quiet, birth, _, _ = gamma.frozen.exact_moments(j, eta*z, q, n, state, False)
    old_return = gamma.exact_moments(j, z, q, n, state)[2]
    live, d = exact_live(j, z, eta, q=q, n=n, state=state), q**state-1
    density = F(d, (q-1)**state)
    ret = min(eta**j*old_return, d*live)
    return quiet, density*birth, ret/d, min(density*live, z**(n-state+1)) if not j else density*live


@lru_cache(8)
def live_geometry(q, n, state):
    import numpy as np
    counts = support_counts(q, n, state)
    rows = []
    for j in range(n+1):
        terms = [(b, c, log2(comb(j, b)*comb(n-j, c)*counts[b+c]))
                 for b in range(j+1) for c in range(n-j+1) if counts[b+c]]
        values = np.array(terms, dtype=float).T
        values.setflags(write=False)
        rows.append(values)
    return tuple(rows)


def log_live(z, eta, *, q=Q, n=64, state=32):
    import numpy as np
    validate_local(z, eta, q, n, state)
    lz, le, bits = log2(z), log2(eta), log2(q)
    u0 = float(np.logaddexp2(0., log2(q-1)+le+lz))-bits
    terms = [lz, le] + ([log2(q-2)+le+lz] if q > 2 else [])
    u1 = float(np.logaddexp2.reduce(terms))-bits
    logd = log2(q**state-1)
    result = np.empty(n+1)
    for j, (b, c, masses) in enumerate(live_geometry(q, n, state)):
        result[j] = np.logaddexp2.reduce(masses+b*u1+(j-b)*u0+c*lz)-logd
    return result


def tilted_logs(z, eta, *, q=Q, n=64, state=32):
    import numpy as np
    import hd_field128_full_rs_gamma_screen as gamma
    if not isfinite(z) or not isfinite(eta):
        raise ValueError('Finite marks required')
    validate_local(z, eta, q, n, state)
    original = gamma.local_logs(z, q, n, state)
    if eta == 1:
        # Exact control identity; independent live identity is tested separately.
        return original.copy()
    tilted = gamma.frozen.local_logs(eta*z, q, n, state, False).copy()
    live = log_live(z, eta, q=q, n=n, state=state)
    density = log2(q**state-1)-state*log2(q-1)
    tilted[:, 1, 0] = np.minimum(original[:, 1, 0]+np.arange(n+1)*log2(eta), live)
    tilted[:, 1, 1] = density+live
    tilted[0] = original[0]  # Empty input has eta^wt(X)=1, including distance cap.
    if np.any(np.isnan(tilted)) or np.any(np.isposinf(tilted)):
        raise ArithmeticError('Nonfinite tilted local cap')
    return tilted


def prefix_correction(alpha, eta, weight=R):
    if (not isfinite(alpha) or not 0 < alpha <= 1 or not isfinite(eta) or eta < 1
            or type(weight) is not int or not 1 <= weight <= 8192):
        raise ValueError('Invalid tilted prefix')
    return -33*weight*alpha*log2(eta)


def validate_options(exponents, seconds):
    if (not 1 <= len(exponents) <= 8 or exponents[0] != 0
            or any(type(e) is not int or not 0 <= e <= 64 for e in exponents)
            or tuple(sorted(set(exponents))) != tuple(exponents)
            or not isfinite(seconds) or not 0 < seconds <= 110):
        raise ValueError('Use up to8 increasing integer eta exponents0..64 starting0, <=110seconds')


def matrix_json(values):
    return [[None if x == -float('inf') else float(x) for x in row] for row in values]


def selected_eta(local, z, alpha, markers, grid, budget, include_j0):
    import numpy as np
    import hd_field128_lift_alias_screen as alias
    coefficients = old.Coefficients(local, z, alpha, budget)
    zero = old.Coefficients(old.ablate(local), z, alpha, budget) if include_j0 else None
    best = [np.full((2, 2), np.inf) for _ in range(1+int(include_j0))]
    winners = [np.full((2, 2), -1, dtype=int) for _ in best]
    reliable = [np.zeros((2, 2), bool) for _ in best]
    underflows = [0]*len(best)
    for index, literal in enumerate(markers):
        budget.check()
        odds = tuple(log2(float(F(x))) for x in literal)
        table = alias.alias_table(z, alpha, odds, grid=grid, _coefficients=coefficients, _budget=budget)
        tables = [table]
        if include_j0:
            tables.append(old.triangular_atlas(zero.composed, odds, coefficients.pair_logs(odds, _budget=budget),
                                               grid=grid, check=budget.check))
        for k, table in enumerate(tables):
            values, good = old.selected(table, odds, full=2026)
            better = values < best[k]
            winners[k][better], reliable[k][better] = index, good.reshape(2, 2)[better]
            np.minimum(best[k], values, out=best[k])
            key = 'input_underflows' if k else 'normalized_input_underflow_terms'
            underflows[k] += int(table.diagnostic[key])
    budget.check()
    return best, [v.tolist() for v in winners], [v.tolist() for v in reliable], underflows


@contextmanager
def remote_limits():
    if sys.platform != 'linux':
        raise RuntimeError('Production requires Linux')
    import resource
    import signal
    resource.setrlimit(resource.RLIMIT_AS, (1 << 30, 1 << 30))
    resource.setrlimit(resource.RLIMIT_CPU, (115, 120))
    resource.setrlimit(resource.RLIMIT_FSIZE, (LIMIT, LIMIT))
    signal.alarm(120)
    try:
        yield
    finally:
        signal.alarm(0)


def run(source_name, output_name, *, eta_log2=DEFAULT_ETAS, include_j0=False, max_seconds=110.):
    with remote_limits():
        return _run(source_name, output_name, eta_log2=eta_log2, include_j0=include_j0, max_seconds=max_seconds)


def _run(source_name, output_name, *, eta_log2=DEFAULT_ETAS, include_j0=False, max_seconds=110.):
    validate_options(eta_log2, max_seconds)
    if type(include_j0) is not bool:
        raise ValueError('Explicit optional J0 flag required')
    source, output = old.prior.previous.inside(source_name), old.prior.previous.inside(output_name)
    if output.exists() or source.stat().st_size > 1 << 20:
        raise ValueError('Exclusive output / <=1MiB input required')
    if shutil.disk_usage(HERE).free-LIMIT < RESERVE:
        raise RuntimeError('20GiB planned reserve required')
    raw, pins = source.read_bytes(), source_pins()
    if sha256(raw).hexdigest() != old.prior.SOURCE_SHA:
        raise ValueError('Pinned v34 result required')
    row = old.prior.validate_source(json.loads(raw))
    witness, profile = row['witness'], row['largest_profile']
    if (len(witness['markers']) != 79 or witness['grid'] != 128
            or (profile['even_full'], profile['odd_full']) != (2026, 2026)):
        raise ValueError('Expected79 markers/grid128/profile2026,2026')
    z, alpha = float(F(witness['z'])), float(F(witness['alpha']))
    if any(not (2.**e)*z < 1 for e in eta_log2):
        raise ValueError('Every requested eta*z must be strictly below1')
    def check_inputs():
        if source.read_bytes() != raw or source_pins() != pins:
            raise RuntimeError('Pinned input or source changed')
    started, results, failure = monotonic(), [], None
    try:
        import hd_field128_full_rs_gamma_screen as gamma
        import numpy as np
        budget = gamma.old.Budget(max_seconds)
        prefix = old.prior.previous.prefactor(z, alpha, F(9, 20), R,
            gamma.old.outer_log_density(Q, 'invertible-toeplitz'))
        if abs(prefix-row['log2_outer_prefix']) > 1e-7:
            raise ArithmeticError('Original-r prefix mismatch')
        # Check formula independently of the exact-control branch.
        original = gamma.local_logs(z, Q, state=32)
        density = log2(Q**32-1)-32*log2(Q-1)
        if np.max(abs(log_live(z, 1.)[1:]+density-original[1:, 1, 1])) > 2e-9:
            raise ArithmeticError('Eta1 exact live formula does not match frozen caps')
        try:
            for exponent in eta_log2:
                budget.check()
                eta = 2.**exponent
                matrices, winners, reliability, underflows = selected_eta(
                    tilted_logs(z, eta), z, alpha, witness['markers'], 128, budget, include_j0)
                adjusted = prefix+prefix_correction(alpha, eta)
                scores = [adjusted+profile['log2_count_upper']+old.prior.profile_transfer(v, v) for v in matrices]
                if not all(map(isfinite, scores)):
                    raise ArithmeticError('Nonfinite completed profile score')
                if exponent == 0 and abs(scores[0]-old.EXPECTED_BASELINE) > 1e-4:
                    raise ArithmeticError('Frozen eta1 profile baseline mismatch')
                results.append(dict(eta_log2=exponent, eta_literal=str(F(eta)),
                    log2_outer_prefix=adjusted, log2_prefix_correction=prefix_correction(alpha, eta),
                    full_return_profile_FLOAT=scores[0],
                    optimistic_J0_profile_FLOAT=scores[1] if include_j0 else None,
                    regional_log2=[matrix_json(v) for v in matrices],
                    winning_marker_indices=winners, selected_reliable=reliability,
                    normalized_input_underflow_terms=underflows, markers_evaluated=79,
                    profile_complete=True, elapsed_seconds=monotonic()-started))
        except (TimeoutError, MemoryError) as error:
            failure = dict(type=type(error).__name__, message=str(error))
        result = dict(schema='gamma32-input-weight-tilt-selected-profile-FLOAT-v39',
            source=dict(path=source_name, sha256=old.prior.SOURCE_SHA), source_sha256=pins,
            conditional_construction=old.prior.previous.state_metadata(32), witness=witness,
            active_groups=R, selected_profile=[2026, 2026], target_distance='9/20',
            eta_log2_requested=list(eta_log2), optional_J0=include_j0, results=results,
            status='complete' if len(results) == len(eta_log2) else 'partial', failure=failure,
            completed_eta_count=len(results), baseline_reproduced=bool(results),
            original_log2_outer_prefix=prefix, log2_profile_count=profile['log2_count_upper'],
            count_source='pinned FLOAT midpoint count cap; not recomputed',
            input_weight_constraint='W_in>=33r from nonzero RS column support and invertible mixers',
            tilt_measure='uniform product input unchanged; eta^W inside expectation, correction inside alpha',
            local_return='min(eta^j J_old(z), (q^32-1) V_eta); no tilted-birth shortcut',
            full_return_formula_valid_conditional=True, J0_is_valid_upper_bound=False,
            all_profiles_replayed=False, setup_risk_added_to_score=False, optimizer_run=False,
            arithmetic='FLOAT with reliability fallback; no directed rounding',
            directed=False, selected_weight_certificate=False, whole_code_certificate=False,
            scope='one support-profile contribution only; no code-failure/impossibility claim',
            elapsed_seconds=monotonic()-started,
            limits=dict(wall_seconds=120, address_space_bytes=1 << 30, output_bytes=LIMIT,
                        minimum_free_after_write=RESERVE))
        check_inputs()
        data = (json.dumps(result, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n').encode()
        if len(data) > LIMIT or shutil.disk_usage(HERE).free-len(data) < RESERVE:
            raise RuntimeError('32KiB output /20GiB planned reserve exceeded')
        with output.open('xb') as stream:
            stream.write(data)
        if output.read_bytes() != data:
            raise RuntimeError('Exclusive output changed')
        return dict(path=output_name, sha256=sha256(data).hexdigest(), status=result['status'],
                    completed_etas=len(results), scores=[r['full_return_profile_FLOAT'] for r in results])
    finally:
        check_inputs()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--eta-log2', type=int, nargs='+', default=list(DEFAULT_ETAS))
    parser.add_argument('--include-j0', action='store_true')
    parser.add_argument('--max-seconds', type=float, default=110.)
    args = parser.parse_args()
    print(json.dumps(run(args.source, args.output, eta_log2=args.eta_log2,
        include_j0=args.include_j0, max_seconds=args.max_seconds), allow_nan=False), flush=True)
