"""Bounded FLOAT profile-wise minimum of complete, separately marked bounds.

For each fixed matching profile, finish each witness's chronological transfer
and its own outer prefix before taking the minimum. The common integer count
cap is outside every fractional power. No regional entry is mixed across
witnesses. This screen is not directed evidence or a certificate.
"""
import argparse
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isfinite, log2
from pathlib import Path, PurePosixPath
import re
from time import monotonic

import hd_field128_gamma_portable_neighbors_v27 as portable

HERE = Path(__file__).resolve().parent
NAME = Path(__file__).name
PORTABLE_SHA256 = '169ee5c7719873ce4528ec9e9a5eea1e87e8c866ae32e713abde0d3678339a94'
CONSTRUCTION = 'toeplitz-fullfield-full-rs64-30-gamma16-feedback-alternating-swapped-pairs-k20-rate-half-v1'
FIELD_ORDER, PAIRS, REPETITIONS = 1 << 128, 4096, 256
MAX_INPUT_BYTES, SCREEN_SECONDS = 1 << 20, 90.


def validate_target(weight):
    if type(weight) is not int or not 2 <= weight <= 4095:
        raise ValueError('Expected one original weight in 2..4095')
    return weight


def proposal_name(value):
    if (not isinstance(value, str) or '\\' in value
            or not value.endswith('.json') or PurePosixPath(value).is_absolute()
            or any(re.fullmatch(r'[A-Za-z0-9_][A-Za-z0-9_.-]*', p) is None
                   or p.endswith('.') for p in value.split('/'))):
        raise ValueError('Expected a portable relative JSON path without traversal')
    return value


def exact_literal(value):
    if not isinstance(value, str):
        raise ValueError('Expected canonical rational strings')
    try:
        number = F(value)
        if str(number) != value or not isfinite(float(number)) or float(number) <= 0:
            raise ValueError('Literal must be canonical, positive, and finite in binary64')
    except (ZeroDivisionError, OverflowError) as error:
        raise ValueError('Invalid positive rational literal') from error
    return number


def project_literals(saved, weight):
    """Read literals only, deliberately ignoring ancestral bounds and links."""
    validate_target(weight)
    if (not isinstance(saved, dict) or type(saved.get('active_groups')) is not int
            or saved['active_groups'] != weight or saved.get('construction') != CONSTRUCTION
            or type(saved.get('state_dimension')) is not int or saved['state_dimension'] != 30
            or saved.get('target_distance') != '2/5'):
        raise ValueError('Proposal must have the same original weight and gamma30 construction')
    witness = saved.get('witness')
    if (not isinstance(witness, dict) or type(witness.get('grid')) is not int or witness['grid'] != 128
            or type(witness.get('state_dimension')) is not int or witness['state_dimension'] != 30):
        raise ValueError('Expected a grid-128 state-30 witness')
    for obj in (saved, witness):
        if ('field_order' in obj and (type(obj['field_order']) is not int or obj['field_order'] != FIELD_ORDER)
                or 'inner_width' in obj and (type(obj['inner_width']) is not int or obj['inner_width'] != 64)):
            raise ValueError('Wrong field or packet width')
    z, alpha = (exact_literal(witness.get(key)) for key in ('z', 'alpha'))
    if not 0 < z < 1 or not 0 < float(z) < 1 or not 0 < alpha <= 1:
        raise ValueError('Expected 0 < z < 1 and 0 < alpha <= 1')
    key = 'exact_positive_marker_proposals' if 'exact_positive_marker_proposals' in witness else 'markers'
    rows = witness.get(key)
    if (not isinstance(rows, list) or not 1 <= len(rows) <= 128
            or any(not isinstance(row, list) or len(row) != 2 for row in rows)):
        raise ValueError('Expected 1..128 exact positive marker pairs')
    markers = [tuple(map(exact_literal, row)) for row in rows]
    if len(set(markers)) != len(markers):
        raise ValueError('Duplicate marker literals')
    return dict(z=str(z), alpha=str(alpha), markers=[list(map(str, p)) for p in markers],
                grid=128, state_dimension=30), key


def load_proposals(specifications, weight):
    validate_target(weight)
    if not 2 <= len(specifications) <= 4:
        raise ValueError('Expected 2..4 explicitly pinned proposals')
    result, seen = [], set()
    for name, digest in specifications:
        name = proposal_name(name)
        if not isinstance(digest, str) or re.fullmatch(r'[0-9a-f]{64}', digest) is None:
            raise ValueError('Expected a lowercase SHA-256 digest')
        path = HERE/name
        if not path.resolve().is_relative_to(HERE.resolve()) or name in seen:
            raise ValueError('Duplicate or escaping proposal path')
        seen.add(name)
        if path.stat().st_size > MAX_INPUT_BYTES:
            raise ValueError('Proposal exceeds the 1 MiB input limit')
        raw = path.read_bytes()
        if sha256(raw).hexdigest() != digest:
            raise ValueError('Proposal byte pin mismatch: '+name)
        witness, key = project_literals(json.loads(raw), weight)
        result.append(dict(source=dict(path=name, sha256=digest), witness=witness,
                           source_marker_key=key))
    return result


def helpers():
    if sha256((HERE/portable.NAME).read_bytes()).hexdigest() != PORTABLE_SHA256:
        raise RuntimeError('Frozen portable wrapper changed')
    return portable.helpers()


def source_pins(loaded):
    pins = portable.source_pins(loaded)
    pins[portable.NAME] = PORTABLE_SHA256
    pins[NAME] = sha256((HERE/NAME).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def combine_profiles(regions, prefixes, full, counts, transfer, budget, *, repetitions=256, chunk=32):
    """Pure array layer; a one-witness control is allowed for tiny tests.

    Every ordered (a,b) class occurs once. Counts are shared; transfer and
    prefix belong to the same witness. Only finished class terms are minimized.
    """
    import numpy as np
    full = np.asarray(full)
    regions = tuple(np.asarray(row, dtype=float) for row in regions)
    if (not 1 <= len(regions) <= 4 or len(prefixes) != len(regions)
            or full.ndim != 1 or len(full) == 0 or not np.issubdtype(full.dtype, np.integer)
            or np.any(np.diff(full) != 1) or not all(isfinite(p) for p in prefixes)
            or type(repetitions) is not int or repetitions < 2 or repetitions % 2
            or type(chunk) is not int or not 1 <= chunk <= 128
            or any(row.shape != (len(full), 2, 2) or np.any(np.isnan(row))
                   or np.any(np.isposinf(row)) for row in regions)):
        raise ValueError('Invalid full profile domain, regional matrices, prefixes, or batching')
    totals = np.full(len(regions), -np.inf)
    total, largest = -np.inf, None
    wins = np.zeros(len(regions), dtype=np.int64)
    evaluated = positive = 0
    for start in range(0, len(full), chunk):
        budget.check()
        logs = counts.log_counts(full[start:start+chunk, None], full[None, :], budget=budget)
        shape = (min(chunk, len(full)-start), len(full))
        if logs.shape != shape or np.any(np.isnan(logs)) or np.any(np.isposinf(logs)):
            raise ArithmeticError('Invalid common count cap array')
        best, selected = np.full(shape, np.inf), np.zeros(shape, dtype=np.int64)
        for index, (row, prefix) in enumerate(zip(regions, prefixes)):
            budget.check()
            product = transfer.multiply(row[start:start+chunk, None], row[None])
            repeated = transfer.power(product, repetitions//2)
            output = np.logaddexp2(repeated[..., 0, 0], repeated[..., 0, 1])
            terms = logs+output+prefix
            if np.any(np.isnan(terms)) or np.any(np.isposinf(terms)):
                raise ArithmeticError('Nonfinite completed class term')
            totals[index] = np.logaddexp2(totals[index], np.logaddexp2.reduce(terms.ravel()))
            improved = terms < best
            best[improved], selected[improved] = terms[improved], index
        valid = np.isfinite(logs)
        wins += np.bincount(selected[valid], minlength=len(regions))
        positive += int(valid.sum())
        evaluated += best.size
        total = float(np.logaddexp2(total, np.logaddexp2.reduce(best.ravel())))
        position = np.unravel_index(np.argmax(best), best.shape)
        if largest is None or best[position] > largest['log2_term']:
            largest = dict(even_full=int(full[start+position[0]]), odd_full=int(full[position[1]]),
                winning_witness_index=int(selected[position]), log2_term=float(best[position]),
                log2_count_upper=float(logs[position]))
    budget.check()
    if not isfinite(total) or not np.all(np.isfinite(totals)):
        raise ArithmeticError('Zero or nonfinite FLOAT union is not accepted as evidence')
    if total > float(min(totals))+1e-9:
        raise ArithmeticError('Profile minimum exceeded a single-witness union')
    return dict(single_log2_upper_screens=list(map(float, totals)), combined_log2_upper_screen=total,
        largest_profile=largest, witness_winning_profile_counts=list(map(int, wins)),
        zero_count_profiles=evaluated-positive, positive_count_profiles=positive,
        winner_tie_rule='first supplied witness; zero-count profiles excluded from wins',
        profile_pairs_evaluated=evaluated, all_profiles_retained=True,
        minimum_stage='each complete chronological transfer plus its own outer prefix and common count',
        entrywise_cross_witness_minimum=False, common_count_cap_outside_alpha=True,
        directed=False, selected_weight_certificate=False, whole_code_certificate=False)


def build_regions(proposal, weight, provider, budget):
    """One fresh bank at a time; only the O(r) regional array escapes."""
    witness = proposal['witness']
    z, alpha = float(F(witness['z'])), float(F(witness['alpha']))
    odds = [tuple(log2(float(F(v))) for v in pair) for pair in witness['markers']]
    coeff = provider.AliasCoefficients(z, alpha, state=30, q=FIELD_ORDER, _budget=budget)
    bank = provider.AliasBankScreen(z, alpha, odds, 0., state=30, q=FIELD_ORDER, distance=F(2, 5),
        grid=128, floor_policy='positive-fallback', _coefficients=coeff, _budget=budget)
    regions = bank.regional(weight).copy()
    prefix = bank.constant+bank.alpha*bank.density*weight
    metadata = dict(source=deepcopy(proposal['source']), witness=deepcopy(witness),
        source_marker_key=proposal['source_marker_key'], log2_outer_prefix=prefix,
        regional_logs_sha256=sha256(regions.tobytes()).hexdigest(),
        coefficient_table_sha256=sha256(coeff.composed.tobytes()).hexdigest(),
        regional_alias_diagnostics=dict(bank.regional_diagnostics), bank_metadata=bank.metadata(),
        marker_conversion='Exact positive rationals retained; FLOAT log2(float(rational)) only')
    return regions, prefix, metadata


def run(specifications, weight, output_dir):
    validate_target(weight)
    with portable.remote_limits():
        loaded = helpers()
        validator, normalizer, neighbors = loaded
        proposals = load_proposals(specifications, weight)
        pins = source_pins(loaded)
        def check():
            for proposal in proposals:
                portable.check_pins(pins, proposal['source'])
        check()
        destination = HERE/validator.relative_name(output_dir)
        portable.reserve_guard(HERE, portable.MAX_OUTPUT_BYTES)
        destination.mkdir(exist_ok=False)
        outputs = portable.Outputs(destination, normalizer.normalize_new_metadata)
        rows, regions, prefixes = [], [], []
        stage, result, status = 'bank-build', None, 'incomplete'
        started = monotonic()
        try:
            outputs.write('driver.json', dict(schema='gamma-portable-profile-min-driver-v28',
                source_sha256=pins, sources=[p['source'] for p in proposals], active_groups=weight,
                hard_wall_seconds=portable.WALL_SECONDS, cooperative_seconds=SCREEN_SECONDS,
                memory_bytes=portable.MEMORY_BYTES, minimum_free_bytes=portable.RESERVE,
                maximum_output_bytes=portable.MAX_OUTPUT_BYTES, maximum_simultaneous_banks=1,
                output_budget_scope='All new JSON; caller must separately bound console logs',
                directed=False, whole_code_certificate=False))
            provider, mode = neighbors.provider, neighbors.mode
            budget = provider.old.Budget(SCREEN_SECONDS)
            try:
                for proposal in proposals:
                    budget.check()
                    matrix, prefix, metadata = build_regions(proposal, weight, provider, budget)
                    regions.append(matrix); prefixes.append(prefix); rows.append(metadata)
                stage = 'common-counts'
                counts = mode.MidpointCounts(PAIRS, weight, max_singletons=64, budget=budget)
                stage = 'complete-profile-union'
                transfer = provider.profiles.transfer
                result = combine_profiles(regions, prefixes, transfer.tail.full_counts(PAIRS, weight),
                    counts, transfer, budget, repetitions=REPETITIONS)
                if result['profile_pairs_evaluated'] != (weight//2+1)**2:
                    raise ArithmeticError('Incomplete original-weight profile domain')
                result['count_provenance'] = counts.diagnostics()
                status, stage = 'complete', None
            except (TimeoutError, MemoryError) as error:
                status = 'timeout' if isinstance(error, TimeoutError) else 'memory-error'
            check(); outputs.check()
            summary = dict(schema='gamma-portable-profile-min-FLOAT-v28', status=status,
                active_groups=weight, construction=CONSTRUCTION, state_dimension=30,
                field_order=FIELD_ORDER, target_distance='2/5', source_sha256=pins,
                proposals=proposals, completed_banks=rows, pending_stage=stage,
                result=result, complete_selected_weight_screen=status == 'complete',
                elapsed_seconds=monotonic()-started, original_weight_used_for_all_histograms_and_prefixes=True,
                old_numerical_bounds_reused=False, historical_origin_verified=False,
                no_weight_interpolation=True, common_potential_used=False,
                directed=False, selected_weight_certificate=False, whole_code_certificate=False,
                scope='One original-weight FLOAT screen; class-wise witness choice only; no directed bound',
                arithmetic='Frozen FLOAT coefficient/alias/transfer arithmetic; no rounding certificate')
            pin = outputs.write('summary.json', summary)
            check(); outputs.check()
            return dict(status=status, summary=pin, active_groups=weight, result=result,
                        whole_code_certificate=False, directed=False)
        finally:
            check(); outputs.check()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--proposal', nargs=2, action='append', metavar=('RELATIVE_JSON', 'SHA256'), required=True)
    parser.add_argument('--target-weight', type=int, required=True)
    parser.add_argument('--output-dir', required=True)
    args = parser.parse_args()
    result = run(args.proposal, args.target_weight, args.output_dir)
    print(json.dumps(result, allow_nan=False), flush=True)
    return 0 if result['status'] == 'complete' else 2


if __name__ == '__main__':
    raise SystemExit(main())
