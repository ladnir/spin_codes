"""Resume the frozen all-A5 FLOAT control in XOR-translation-orbit order.

For P_j(i)=inverse(i xor j), P_j(A xor d)=P_(j xor d)(A).
The sixteen translates of an odd-cardinality A are distinct. They therefore
reuse the same sixteen pre-Q fibre tables in the existing 512-entry cache.
This changes execution order only. Every original labelled A still receives
its own sixteen ordered columns, output permutations, and whole-word tail.
No orbit multiplier or symmetry quotient is used.

The caller supplies an external SHA256 pin for an observed seed receipt.
We validate its scope, current source/input pins, and recorded root algebra.
Hashes bind cached evidence; they do not independently verify FLOAT arithmetic.
Missing A are evaluated through the unchanged frozen controller. The exact
fallback budget is 128 attempts across the seed and this continuation.

Only complete, unique coverage of all 4368 labelled A can expose the selected
class bound. The result remains FLOAT, not a directed or whole-code certificate.
No provider, native work, or proof control runs when this module is imported.
"""
from copy import deepcopy
from hashlib import sha256
import json
from math import isfinite
from pathlib import Path
from time import perf_counter

import hd_lambda_arity5_cover_control as base


DIRECTORY = Path(__file__).resolve().parent
BASE_SHA256 = '3ad62578884e1849babfaba71c049e55a662d3d0b924e6b978c9e6285d3316b9'
SOURCE_FILES = base.SOURCE_FILES + ('hd_lambda_arity5_resume_control.py',)
MAX_EXACT_ATTEMPTS = 128
SCOPE = dict(
    construction='scheduled-Hlambda-bc-no-cross', lambda_aes_byte=0xbc,
    schedule='input-inverse-output-xor', groups=18, arity=5,
    source_class='Fmax<=1', original_row_arity_per_active_group=5,
    heterogeneous_labelled_five_row_choices_permitted=True,
    arbitrary_arity_mixtures_with_total_arity90_covered=False,
    original_A_fixed_across_columns=True, full_independent_D1=True,
    raw_field_order=256, rs_root_normalization='255**(-40)',
    aggregate_input_activity_floor=810, retained_deficit=630,
    original_per_row_activity_floors_relaxed=True,
    minimum_only_after_complete_fixed_A_word_tail=True,
    coefficientwise_minimum_between_exact_and_surrogate=False,
    extra_rowset_multiplier=False, source_projective_prefix_included_once=True,
    exact_fallback_trigger_log2=base.FALLBACK_THRESHOLD,
    target_selected_arity_log2=base.TARGET_LOG2,
    refinement_cap_does_not_stop_marked_coverage=True,
    arithmetic='positive-FLOAT-log2-not-directed', null_log2_means_exact_zero=True,
    directed=False, whole_code_certificate=False, regional_provider_runs=0)


def translation_orbits():
    """Return 273 disjoint sixteen-set orbits, retaining every original A."""
    seen, orbits = set(), []
    for source in base.original_sources():
        mask = base.source_mask(source)
        if mask in seen:
            continue
        orbit = tuple(sorted(tuple(sorted(i ^ d for i in source)) for d in range(16)))
        masks = {base.source_mask(row) for row in orbit}
        if len(masks) != 16 or seen.intersection(masks):
            raise AssertionError('Odd-arity translation orbits must be disjoint and free')
        seen.update(masks)
        orbits.append(orbit)
    if len(orbits) != 273 or len(seen) != base.SOURCE_COUNT:
        raise AssertionError('Incomplete labelled orbit inventory')
    return tuple(orbits)


def remaining_sources(records):
    """Filter covered labels without changing any source or column labels."""
    covered = {base.source_mask(row['original_source']) for row in records}
    return tuple(source for orbit in translation_orbits() for source in orbit
                 if base.source_mask(source) not in covered)


def _pins():
    return {name: sha256((DIRECTORY/name).read_bytes()).hexdigest() for name in SOURCE_FILES}


def _restore_logs(value):
    """Revive exact-zero log values without changing non-log null fields."""
    if isinstance(value, dict):
        return {key: (float('-inf') if item is None and key.endswith('_log2')
                      else _restore_logs(item)) for key, item in value.items()}
    if isinstance(value, list):
        return [_restore_logs(item) for item in value]
    return value


def _same(actual, expected, label):
    if actual != expected:
        raise ValueError(f'Seed {label} does not match the frozen control')


def validate_seed(seed, pins, saved, measure, prefix):
    """Validate observed cached rows; do not treat hashes as arithmetic proof.

    `pins` contains the currently reviewed source files. A resumed seed may
    additionally pin this wrapper; a first-run seed pins exactly the base set.
    All present source pins must match, and every base pin is mandatory.
    """
    seed = _restore_logs(deepcopy(seed))
    for key, expected in SCOPE.items():
        _same(seed.get(key), expected, key)
    for key, expected in (('input_artifact', base.INPUT), ('input_sha256', base.INPUT_SHA),
                          ('raw_artifact_sha256_after', base.INPUT_SHA),
                          ('source_pins_unchanged', True), ('mark', saved['mark']),
                          ('source_count', saved['source_count']), ('common_measure', measure),
                          ('common_prefix_log2', prefix)):
        _same(seed.get(key), expected, key)
    before, after = seed.get('source_pins_before'), seed.get('source_pins_after')
    if not isinstance(before, dict) or not isinstance(after, dict):
        raise ValueError('Seed source pins are missing')
    if set(before) not in (set(base.SOURCE_FILES), set(SOURCE_FILES)) or before != after:
        raise ValueError('Seed source inventory or before/after pins changed')
    if any(pins.get(name) != digest for name, digest in before.items()):
        raise ValueError('Seed source bytes differ from the current reviewed files')
    _same(before.get('hd_lambda_arity5_cover_control.py'), BASE_SHA256, 'base helper pin')
    if (type(seed.get('exact_fallback_limit')) is not int
            or not 0 <= seed['exact_fallback_limit'] <= MAX_EXACT_ATTEMPTS):
        raise ValueError('Invalid seed fallback budget')
    records, attempts = seed.get('records'), seed.get('exact_fallbacks')
    if not isinstance(records, list) or not isinstance(attempts, list):
        raise ValueError('Seed must contain rows and its exact-attempt ledger')
    summary = base.coverage_summary(records, prefix)
    if summary['duplicate_original_masks'] or summary['invalid_records']:
        raise ValueError('Seed has duplicate, malformed, or non-minimal recorded roots')
    if len(attempts) > seed['exact_fallback_limit']:
        raise ValueError('Seed exact attempts exceed its declared cumulative budget')
    by_mask = {row['original_mask']: row for row in records}
    attempted = set()
    for row in records:
        marked = base._log_value(row['marked_tail_log2'])/base.GROUPS-base.RS_ROOT_CHARGE
        _same(row['marked_normalized_root_log2'], marked, 'marked tail normalization')
        chosen = row['normalized_fixed_A_root_log2']
        _same(row['fixed_A_repeated_r_bound_log2'], prefix+base.GROUPS*chosen,
              'fixed-A source/RS normalization')
        exact_present = 'exact_normalized_root_log2' in row
        if exact_present:
            exact = base._log_value(row['exact_tail_log2'])/base.GROUPS-base.RS_ROOT_CHARGE
            _same(row['exact_normalized_root_log2'], exact, 'exact tail normalization')
            _same(row.get('exact_fallback_status'), 'complete', 'completed exact status')
        method = ('exact-GF4' if exact_present and row['exact_normalized_root_log2'] < marked
                  else 'marked-fibre')
        _same(row.get('selected_method'), method, 'selected complete method')
        if row.get('exact_fallback_status') not in (
                'not_needed', 'attempt_cap_skipped', 'complete', 'timeout_marked_retained'):
            raise ValueError('Unknown seed fallback status')
    for attempt in attempts:
        mask = base.source_mask(attempt['original_source'])
        if mask != attempt['original_mask'] or mask not in by_mask or mask in attempted:
            raise ValueError('Seed fallback ledger has a duplicate or foreign source')
        attempted.add(mask)
        row = by_mask[mask]
        status = attempt.get('status')
        if status not in ('complete', 'timeout_marked_retained'):
            raise ValueError('Seed has an unfinished or unrecognized exact attempt')
        _same(row['exact_fallback_status'], status, 'fallback row/ledger status')
        _same(attempt['marked_fixed_A_bound_log2'],
              prefix+base.GROUPS*row['marked_normalized_root_log2'], 'fallback marked bound')
        if status == 'complete':
            for key in ('exact_tail_log2', 'exact_normalized_root_log2', 'selected_method'):
                _same(attempt.get(key), row.get(key), f'fallback {key}')
            _same(attempt['selected_fixed_A_bound_log2'], row['fixed_A_repeated_r_bound_log2'],
                  'fallback selected bound')
        elif 'exact_normalized_root_log2' in row:
            raise ValueError('Timed-out fallback cannot supply an exact root')
    required = {row['original_mask'] for row in records
                if row['exact_fallback_status'] in ('complete', 'timeout_marked_retained')}
    if attempted != required:
        raise ValueError('Seed fallback ledger omits a charged attempt')
    skipped = sum(row['exact_fallback_status'] == 'attempt_cap_skipped' for row in records)
    _same(seed.get('exact_fallback_cap_skipped'), skipped, 'fallback skipped count')
    # The seed's cached coverage and target claims are ignored: recompute them.
    return seed, summary


def run(seed_path, *, seed_sha256, batch_size=16, max_seconds=360, progress=None):
    """One bounded continuation, with no new witnesses or numerical provider."""
    if (type(batch_size) is not int or batch_size not in (16, 32)
            or not isfinite(max_seconds) or not 0 < max_seconds <= 360):
        raise ValueError('Use batch16/32 and a wall budget in (0,360] seconds')
    if (not isinstance(seed_sha256, str) or len(seed_sha256) != 64
            or any(c not in '0123456789abcdef' for c in seed_sha256)):
        raise ValueError('An explicit lowercase SHA256 pin for the observed seed is required')
    start = perf_counter()
    def check_time():
        if perf_counter()-start > max_seconds:
            raise TimeoutError('All-A5 resume exceeded its wall-time budget')
    seed_path = Path(seed_path).resolve()
    seed_bytes = seed_path.read_bytes()
    if sha256(seed_bytes).hexdigest() != seed_sha256:
        raise ValueError('Seed bytes do not match the supplied external pin')
    pins_before, raw_before = _pins(), base._input_hash()
    _same(pins_before['hd_lambda_arity5_cover_control.py'], BASE_SHA256, 'current base helper')
    _same(raw_before, base.INPUT_SHA, 'current regional artifact')
    saved, measure, prefix = base.load_measure()
    seed, seed_summary = validate_seed(json.loads(seed_bytes), pins_before, saved, measure, prefix)
    sources = remaining_sources(seed['records'])
    prior_attempts = len(seed['exact_fallbacks'])
    available_attempts = MAX_EXACT_ATTEMPTS-prior_attempts
    work = dict(records=[], exact_fallbacks=[], exact_fallback_cap_skipped=0,
                pending_preparations_without_tail=0, stop_reason=None)
    def report(row):
        if progress:
            row = dict(row)
            row['new_covered'] = row.get('covered', 0)
            row['covered'] = len(seed['records'])+row['new_covered']
            row['exact_attempts'] = prior_attempts+row.get('exact_attempts', 0)
            progress(row)
    try:
        check_time()
        if sources:
            table = base.moment_table(measure, base.GROUPS)
            table2 = base.moment_table(measure, 2*base.GROUPS)
            check_time()
            work = base._process_sources(sources, table, table2, prefix,
                batch_size=batch_size, fallback_limit=available_attempts,
                check_time=check_time, progress=report)
    except TimeoutError as error:
        work['stop_reason'] = dict(type=type(error).__name__, message=str(error))
    records = seed['records']+work['records']
    attempts = seed['exact_fallbacks']+work['exact_fallbacks']
    if len(attempts) > MAX_EXACT_ATTEMPTS:
        raise AssertionError('Cumulative exact-fallback cap exceeded')
    pins_after, raw_after = _pins(), base._input_hash()
    seed_after = sha256(seed_path.read_bytes()).hexdigest()
    unchanged = (pins_before == pins_after and raw_before == raw_after == base.INPUT_SHA
                 and seed_after == seed_sha256)
    coverage = base.coverage_summary(records, prefix)
    complete = (coverage['all_original_A_covered'] and work['stop_reason'] is None and unchanged)
    bound = coverage['full_selected_arity_bound_log2'] if complete else None
    coverage['full_selected_arity_bound_log2'] = bound
    skipped = seed['exact_fallback_cap_skipped']+work['exact_fallback_cap_skipped']
    result = dict(seed)
    result.update(records=records, exact_fallbacks=attempts,
        exact_fallback_cap_skipped=skipped, exact_fallback_limit=MAX_EXACT_ATTEMPTS,
        fallback_budget_exhausted=bool(skipped), stop_reason=work['stop_reason'],
        pending_preparations_without_tail=work['pending_preparations_without_tail'],
        coverage=coverage, complete=bool(complete),
        target_met=bool(complete and bound <= base.TARGET_LOG2),
        batch_size=batch_size, max_seconds=max_seconds,
        source_pins_before=pins_before, source_pins_after=pins_after,
        source_pins_unchanged=unchanged, raw_artifact_sha256_after=raw_after,
        resume_version=1, seed_artifact=str(seed_path), seed_sha256=seed_sha256,
        seed_sha256_after=seed_after,
        cached_seed_trust='externally pinned observed FLOAT receipt; hashes do not verify arithmetic',
        seed_recorded_original_A=seed_summary['unique_original_A'],
        new_recorded_original_A=len(work['records']), remaining_original_A_at_start=len(sources),
        seed_exact_attempts=prior_attempts, new_exact_attempts=len(work['exact_fallbacks']),
        available_exact_attempts_at_start=available_attempts,
        execution_order='XOR-translation orbits, all labelled A retained',
        orbit_quotient_used=False, orbit_multiplier_used=False,
        fibre_cache_modified=False, frozen_process_sources_reused=True,
        elapsed_seconds=perf_counter()-start)
    if progress:
        progress(dict(stage='finished', covered=coverage['unique_original_A'],
            total=base.SOURCE_COUNT, complete=bool(complete), target_met=result['target_met'],
            full_selected_arity_bound_log2=bound, exact_attempts=len(attempts),
            elapsed_seconds=result['elapsed_seconds']))
    return base._jsonable(result)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed', type=Path, required=True)
    parser.add_argument('--seed-sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--batch-size', type=int, choices=(16, 32), default=16)
    parser.add_argument('--max-seconds', type=float, default=360)
    args = parser.parse_args()
    if args.output.exists() or args.output.resolve() == args.seed.resolve():
        parser.error('Refuse to overwrite seed or existing evidence')
    result = run(args.seed, seed_sha256=args.seed_sha256, batch_size=args.batch_size,
        max_seconds=args.max_seconds,
        progress=lambda row: print(json.dumps(base._jsonable(row), allow_nan=False), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
