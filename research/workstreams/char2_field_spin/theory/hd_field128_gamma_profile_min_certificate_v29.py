"""Directed selected-weight union of complete, separately marked profile caps.

Each witness has its own fractional exponent, normalized regional matrices,
and outer prefactor. Finish its chronological transfer before taking the
minimum for a fixed ordered matching profile. The common exact integer count
cap stays outside every fractional power. No matrix entry is mixed between
witnesses. This producer has a new receipt identity: rows ALREADY include the
outer prefactors and must not receive another common prefix in an audit.

Only exact proposal literals are imported. Historical calibration bounds and
ancestry are deliberately not verified or used. Numerical helpers are frozen.
Actual execution is remote-only; tiny arithmetic controls need no bank.
"""
import argparse
import ast
from contextlib import contextmanager
from copy import deepcopy
from fractions import Fraction as F
import gc
from hashlib import sha256
import importlib
import json
import os
from pathlib import Path
import signal
import sys
from time import monotonic

import hd_field128_gamma_portable_profile_min_v28 as literals

HERE, NAME = Path(__file__).resolve().parent, Path(__file__).name
LITERALS_SHA = '029cecfc78f697b6b660b3d241165324ee5f09eae6de9e2512e0ee92e84679bf'
CERT_NAME = 'hd_field128_gamma_midpoint_certificate'
CERT_SHA = 'ff8b6bb70597dc06ec678d7b9b1dda935e0e2c6e64d15f25d4f1767fb9a0f8cb'
HARD_SECONDS, COMPUTE_SECONDS, TABLE_BYTES = 600, 540., 64 << 20
PRECISION, PACKET_BITS, AGGREGATE_PRECISION = 256, 80, 384
SCHEMA = 'gamma-profile-min-directed-selected-weight-v29'
portable = literals.portable


def helpers():
    for name, expected in ((literals.NAME, LITERALS_SHA), (CERT_NAME+'.py', CERT_SHA),
                           (portable.NAME, literals.PORTABLE_SHA256)):
        if sha256((HERE/name).read_bytes()).hexdigest() != expected:
            raise RuntimeError('Frozen helper changed: '+name)
    return importlib.import_module(CERT_NAME)


def source_pins(cert):
    """Static local import closure plus the frozen numerical/proof-note closure."""
    pins, pending = dict(cert.source_pins()), [NAME]
    visited = set()
    while pending:
        name = pending.pop()
        if name in visited:
            continue
        visited.add(name)
        raw = (HERE/name).read_bytes()
        digest = sha256(raw).hexdigest()
        if name in pins and pins[name] != digest:
            raise RuntimeError('Conflicting source digest: '+name)
        pins[name] = digest
        for node in ast.walk(ast.parse(raw)):
            modules = ([item.name for item in node.names] if isinstance(node, ast.Import)
                       else [node.module] if isinstance(node, ast.ImportFrom) and node.module else [])
            for module in modules:
                target = module.split('.')[0]+'.py'
                if (HERE/target).is_file() and target not in visited:
                    pending.append(target)
    return dict(sorted(pins.items()))


def check_pins(pins, proposals):
    for proposal in proposals:
        portable.check_pins(pins, proposal['source'])


@contextmanager
def remote_limits():
    """Install hard limits before numerical imports; no local bank execution."""
    if not sys.platform.startswith('linux'):
        raise RuntimeError('Actual directed production requires Linux')
    import resource
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    ceiling = min([portable.MEMORY_BYTES]+[n for n in (soft, hard) if n != resource.RLIM_INFINITY])
    resource.setrlimit(resource.RLIMIT_AS, (ceiling, ceiling))
    old = signal.getsignal(signal.SIGALRM)
    signal.signal(signal.SIGALRM, lambda unused_signum, unused_frame: os._exit(124))
    signal.setitimer(signal.ITIMER_REAL, HARD_SECONDS)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, old)


def remaining(budget):
    budget.check()
    return max(.001, budget.max_seconds-(monotonic()-budget.started))


def build_regions(proposal, weight, cert, budget):
    """Only O(r) exact upper endpoints and metadata escape a fresh actual bank."""
    from flint import arb, ctx
    witness = proposal['witness']
    bank = None
    try:
        bank = cert.provider.AliasBank(F(witness['z']), F(witness['alpha']),
            [tuple(map(F, row)) for row in witness['markers']], state=30, q=1 << 128,
            distance=F(2, 5), grid=128, precision=PRECISION, bits=PACKET_BITS,
            outer_measure='invertible-toeplitz', max_seconds=remaining(budget),
            max_table_bytes=TABLE_BYTES)
        if type(bank) is not cert.provider.AliasBank:
            raise TypeError('Expected the actual frozen all-pair AliasBank')
        bank.check_pins()
        with ctx.workprec(PRECISION):
            regions = cert.frozen.regional_matrices(bank, weight, budget)
            regions = {key: tuple(cert.base.upper(x) for x in row) for key, row in regions.items()}
            exponent = bank.constant+bank.slope*weight
            prefactor = cert.base.upper(arb(2)**exponent)
            if not prefactor > 0:
                raise ArithmeticError('Expected a finite positive outer prefactor')
            metadata = dict(bank.metadata(), proposal_source=deepcopy(proposal['source']),
                source_marker_key=proposal['source_marker_key'], original_weight=weight,
                exact_input_witness=deepcopy(witness), fresh_bank_rebuilt=True,
                normalized_regional_full_range=[min(regions), max(regions)],
                normalized_regional_matrix_count=len(regions),
                outer_prefactor=cert.frozen.frozen.support_log(prefactor),
                outer_prefactor_formula='exp2(bank.constant + bank.slope * original_weight)',
                historical_origin_verified=False, old_numerical_bounds_reused=False)
        budget.check()
        return regions, prefactor, metadata
    finally:
        if bank is not None:
            bank.check_pins()
        # Frozen RegionalCache has an lru_cache/bound-method self-cycle. Drop
        # the last local bank reference and reclaim that unreachable cycle
        # before the next bank is built; only detached Arb endpoints escape.
        bank = None
        gc.collect()


def directed_rows(region_banks, prefactors, pairs, weight, counts, cert, *,
                  repetitions=256, budget=None, checkpoint=None):
    """Tiny generic arithmetic interface; production supplies checked midpoint counts.

    Upper endpoints are exact dyadics, so their minimum is an ordered comparison,
    not an uncertain interval comparison. Only completed rows survive a stop.
    """
    from flint import arb
    full = cert.tail.full_counts(pairs, weight)
    if (not 1 <= len(region_banks) <= 4 or len(prefactors) != len(region_banks)
            or type(repetitions) is not int or repetitions < 2 or repetitions % 2
            or any(set(rows) != set(full) or any(len(row) != 4 for row in rows.values())
                   for rows in region_banks)):
        raise ValueError('Complete matching regions, prefactors and even repetitions required')
    upper = cert.base.upper
    region_banks = [{f: tuple(upper(x) for x in rows[f]) for f in full} for rows in region_banks]
    prefactors = tuple(upper(x) for x in prefactors)
    if any(not x > 0 for x in prefactors):
        raise ValueError('Prefactors must be strictly positive')
    budget = budget or cert.base.Budget()
    total, rows, positive = arb(0), [], 0
    winners = [0]*len(region_banks)
    try:
        for first in full:
            budget.check()
            row_sum, nonzero = arb(0), 0
            row_wins = [0]*len(region_banks)
            left = [bank[first] for bank in region_banks]
            for position, second in enumerate(full):
                if not position % 32:
                    budget.check()
                count = counts.count(first, second)
                if type(count) is not int or count < 0:
                    raise ArithmeticError('Expected an exact nonnegative integer count cap')
                if not count:
                    continue
                candidates = []
                for index, bank in enumerate(region_banks):
                    if repetitions == 256:
                        transfer = cert.frozen.frozen.pair_row_sum128(left[index], bank[second])
                    else:
                        product = cert.base.matrix_product(left[index], bank[second])
                        power = cert.base.matrix_power(product, repetitions//2)
                        transfer = upper(power[0]+power[1])
                    candidates.append(upper(count*prefactors[index]*transfer))
                winner = min(range(len(candidates)), key=candidates.__getitem__)
                row_sum += candidates[winner]
                row_wins[winner] += 1
                nonzero += 1
            row_sum = upper(row_sum)
            budget.check()
            row = dict(even_full=first, odd_full_range=[full.start, full.stop-1],
                profile_pairs_evaluated=len(full), profiles_with_positive_count=nonzero,
                witness_winning_profile_counts=row_wins,
                weighted_row_bound=cert.frozen.frozen.support_log(row_sum),
                outer_prefactors_already_included=True, complete_even_profile_row=True,
                covers_all_odd_profiles_for_this_even_profile=True,
                complete_selected_weight=False, whole_code_certificate=False)
            if checkpoint:
                checkpoint(row)
            rows.append(row)
            total = upper(total+row_sum)
            positive += nonzero
            winners = [a+b for a, b in zip(winners, row_wins)]
    except (TimeoutError, ArithmeticError) as error:
        raise cert.IncompleteTransfer(error, dict(rows=rows, completed_even_rows=len(rows),
            pending_even_full_range=[full.start+len(rows), full.stop-1],
            profile_pairs_evaluated=len(rows)*len(full), profiles_with_positive_count=positive,
            witness_winning_profile_counts=winners)) from error
    return total, dict(rows=rows, completed_even_rows=len(rows), pending_even_full_range=None,
        profile_pairs_evaluated=len(full)**2, profiles_with_positive_count=positive,
        witness_winning_profile_counts=winners)


def aggregate_rows(rows, full, cert):
    """Fresh 384-bit union of saved row ceilings; no outer prefix is added."""
    from flint import arb, ctx
    if len(rows) != len(full):
        raise ValueError('A selected-weight union needs every complete even-profile row')
    with ctx.workprec(AGGREGATE_PRECISION):
        total = arb(0)
        for expected, row in zip(full, rows):
            if (row.get('even_full') != expected or row.get('odd_full_range') != [full.start, full.stop-1]
                    or row.get('profile_pairs_evaluated') != len(full)
                    or row.get('complete_even_profile_row') is not True
                    or row.get('covers_all_odd_profiles_for_this_even_profile') is not True
                    or row.get('outer_prefactors_already_included') is not True):
                raise ValueError('Incomplete, duplicate or wrongly scoped weighted row')
            bound = row['weighted_row_bound']
            if bound.get('zero') is True:
                if bound.get('log2_upper_rational') is not None:
                    raise ValueError('Zero row has a spurious logarithmic ceiling')
                continue
            if bound.get('zero') is not False:
                raise ValueError('Explicit zero semantics required')
            text = bound['log2_upper_rational']
            if not isinstance(text, str) or str(F(text)) != text:
                raise ValueError('Expected a canonical rational row ceiling')
            total = cert.base.upper(total+arb(2)**cert.base.rational(F(text)))
        return cert.frozen.frozen.support_log(total)


def run(specifications, weight, output_dir, *, failure_bits=54):
    literals.validate_target(weight)
    if type(failure_bits) is not int or not 0 <= failure_bits <= 4096:
        raise ValueError('Invalid selected-weight failure target')
    with remote_limits():
        cert = helpers()
        from flint import arb, ctx
        proposals = literals.load_proposals(specifications, weight)
        pins = source_pins(cert)
        check = lambda: check_pins(pins, proposals)
        check()
        if ('/' in output_dir or not output_dir or output_dir.endswith('.')
                or any(c not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.-' for c in output_dir)
                or output_dir in ('.', '..')):
            raise ValueError('Expected a new portable output-directory basename')
        destination = HERE/output_dir
        portable.reserve_guard(HERE, portable.MAX_OUTPUT_BYTES)
        destination.mkdir(exist_ok=False)
        outputs = portable.Outputs(destination, lambda value, unused_root: value)
        budget = cert.base.Budget(COMPUTE_SECONDS)
        regions, prefactors, banks, checkpoints = [], [], [], []
        full = cert.tail.full_counts(literals.PAIRS, weight)
        partial = dict(rows=[], completed_even_rows=0, pending_even_full_range=[full.start, full.stop-1],
                       profile_pairs_evaluated=0, profiles_with_positive_count=0,
                       witness_winning_profile_counts=[0]*len(proposals))
        stage, failure, result, context_hash = 'fresh-directed-banks', None, None, None
        try:
            outputs.write('driver.json', dict(schema='gamma-profile-min-directed-driver-v29',
                source_sha256=pins, proposals=proposals, active_groups=weight,
                hard_wall_seconds=HARD_SECONDS, cooperative_seconds=COMPUTE_SECONDS,
                memory_bytes=portable.MEMORY_BYTES, table_bytes=TABLE_BYTES,
                minimum_free_bytes=portable.RESERVE, maximum_output_bytes=portable.MAX_OUTPUT_BYTES,
                maximum_simultaneous_banks=1, fresh_banks_required=True,
                directed=False, whole_code_certificate=False))
            try:
                for index, proposal in enumerate(proposals):
                    matrix, prefix, metadata = build_regions(proposal, weight, cert, budget)
                    check()
                    regions.append(matrix); prefactors.append(prefix)
                    pin = outputs.write(f'bank_{index:02d}.json', metadata)
                    banks.append(dict(index=index, metadata=pin, proposal=deepcopy(proposal)))
                stage = 'common-counts'
                counts = cert.get_counts(literals.PAIRS, weight, None, 64, budget)
                context = dict(schema=SCHEMA, construction=literals.CONSTRUCTION,
                    field_order=literals.FIELD_ORDER, state_dimension=30, target_distance='2/5',
                    active_groups=weight, source_sha256=pins, banks=banks,
                    precision=PRECISION, packet_rounding_bits=PACKET_BITS, grid=128,
                    expected_even_rows=len(full), expected_profile_pairs=len(full)**2,
                    repetitions=256, paired_region_power=128,
                    count_provenance=counts.diagnostics(), common_count_cap_outside_alpha=True,
                    own_outer_prefactor_inside_each_completed_profile_term=True,
                    entrywise_cross_witness_minimum=False, historical_origin_verified=False,
                    old_numerical_bounds_reused=False, original_weight_used_for_all_histograms_and_prefixes=True,
                    directed=True, whole_code_certificate=False,
                    scope='All matching profiles at one weight only if complete; no common outer prefix')
                context_hash = sha256(json.dumps(context, sort_keys=True, separators=(',', ':'),
                                                 allow_nan=False).encode()).hexdigest()
                outputs.write('context.json', dict(context=context, context_sha256=context_hash,
                    complete_selected_weight=False, whole_code_certificate=False))

                def save(row):
                    checkpoints.append(outputs.write(f"row_{len(checkpoints):06d}.json",
                        dict(row, context_sha256=context_hash, active_groups=weight,
                             directed=True, scope='One fully weighted row only; not a selected-weight bound')))

                stage = 'all-profile-directed-union'
                with ctx.workprec(PRECISION):
                    total, partial = directed_rows(regions, prefactors, literals.PAIRS, weight,
                        counts, cert, budget=budget, checkpoint=save)
                    raw_bound = cert.frozen.frozen.support_log(total)
                stage = 'fresh-384-saved-row-union'
                aggregate = aggregate_rows(partial['rows'], full, cert)
                if raw_bound['zero'] != aggregate['zero']:
                    raise ArithmeticError('Raw and saved-row zero semantics disagree')
                ceiling = (None if aggregate['zero'] else str(max(F(raw_bound['log2_upper_rational']),
                                                                 F(aggregate['log2_upper_rational']))))
                result = dict(zero=aggregate['zero'], log2_upper_rational=ceiling,
                    direct_union_256=raw_bound, saved_row_union_384=aggregate,
                    aggregation_precision=AGGREGATE_PRECISION,
                    meets_requested_selected_weight_target=aggregate['zero'] or F(ceiling) <= -failure_bits)
                budget.check()
                stage = None
            except cert.IncompleteTransfer as error:
                partial = error.partial
                failure = dict(error_type=type(error.cause).__name__, message=str(error.cause))
            except (TimeoutError, MemoryError, ArithmeticError) as error:
                failure = dict(error_type=type(error).__name__, message=str(error))
            complete = failure is None and result is not None and partial['completed_even_rows'] == len(full)
            check(); outputs.check()
            summary = dict(schema=SCHEMA, construction=literals.CONSTRUCTION, active_groups=weight,
                field_order=literals.FIELD_ORDER, state_dimension=30, target_distance='2/5',
                status='complete' if complete else 'incomplete', failure=failure, pending_stage=stage,
                source_sha256=pins, proposals=proposals, banks=banks, context_sha256=context_hash,
                **partial, row_checkpoints=checkpoints, result=result if complete else None,
                log2_upper_rational=result['log2_upper_rational'] if complete else None,
                complete_selected_weight=complete, covers_all_supports_at_weight=complete,
                covers_all_even_odd_splits=complete, requested_failure_bits=failure_bits,
                directed=True, whole_code_certificate=False, historical_origin_verified=False,
                old_numerical_bounds_reused=False, outer_prefactors_already_included=True,
                elapsed_seconds=monotonic()-budget.started,
                scope='Fresh directed selected-weight producer; no whole-code or ancestry-authentication claim')
            pin = outputs.write('summary.json', summary)
            check(); outputs.check()
            return dict(status=summary['status'], summary=pin, active_groups=weight,
                        result=summary['result'], whole_code_certificate=False)
        finally:
            check(); outputs.check()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--proposal', nargs=2, action='append', required=True, metavar=('RELATIVE_JSON', 'SHA256'))
    parser.add_argument('--target-weight', type=int, required=True)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--failure-bits', type=int, default=54)
    args = parser.parse_args()
    result = run(args.proposal, args.target_weight, args.output_dir, failure_bits=args.failure_bits)
    print(json.dumps(result, allow_nan=False), flush=True)
    return 0 if result['status'] == 'complete' else 2


if __name__ == '__main__':
    raise SystemExit(main())
