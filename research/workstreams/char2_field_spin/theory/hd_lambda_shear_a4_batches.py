"""Bounded batches and strict assembly of the 1820 labelled a4 SHEAR sets.

Ordinals index itertools.combinations(range(16), 4) in lexicographic order.
Each batch evaluates at most384 original sets with the frozen shared-ratio
kernel and one prepared common measure. Original A stays fixed through all
eight chronological pairs. The completed-A root is tail/18-32log2(255),
where the tail retains aggregate original input degree at least648.

The assembler checks every completed row against its ordinal, normalization,
common measure, and source pins. Completed rows may not overlap. Requested
ranges may overlap after a timeout, but a resume must omit completed rows.
A full bound is prefix+18log2(sum_A 2^root_A) only after all1820 labels occur.
Missing labels produce a partial receipt with no full-class bound. Neither
stage applies another row-set factor. This is selected r18/Fmax<=1 FLOAT
evidence, not a directed or whole-code certificate. No new provider is run.
"""
from hashlib import sha256
from itertools import combinations
import json
from math import isfinite, log2
from pathlib import Path
from time import perf_counter

import numpy as np

import hd_lambda_shear_shared_control as core


DIRECTORY = Path(__file__).resolve().parent
CORE_SHA = '933a8a423af5b65a00c6a070d68b66eeb14171b5c00328a492eb0311fbf9b38d'
SOURCES = tuple(combinations(range(16), 4))
GROUPS, ARITY, FLOOR, DEFICIT = 18, 4, 648, 504
MAX_BATCH, MAX_SECONDS = 384, 180
SCHEMA = 'lambda-shear-shared-a4-batch-v1'


def dependencies():
    if sha256((DIRECTORY/'hd_lambda_shear_shared_control.py').read_bytes()).hexdigest() != CORE_SHA:
        raise ValueError('Frozen selected shared-ratio controller changed')
    return {**core.dependencies(), Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()}


def select_ordinals(*, begin=None, end=None, ordinals=None):
    if ordinals is None:
        if type(begin) is not int or type(end) is not int or not 0 <= begin < end <= len(SOURCES):
            raise ValueError('Use a half-open ordinal interval inside [0,1820)')
        ordinals = tuple(range(begin, end))
    else:
        if begin is not None or end is not None:
            raise ValueError('Choose an interval or an explicit ordinal list, not both')
        ordinals = tuple(ordinals)
    if (not 1 <= len(ordinals) <= MAX_BATCH
            or any(type(i) is not int or not 0 <= i < len(SOURCES) for i in ordinals)
            or tuple(sorted(set(ordinals))) != ordinals):
        raise ValueError('Require at most384 distinct increasing original-set ordinals')
    return ordinals


def common_context(bank, measure, prefix):
    if (bank['groups'] != GROUPS or bank['source_class'] != 'Fmax<=1'
            or bank['mark'] != '9927/10000' or not isfinite(prefix)):
        raise ValueError('Wrong common regional scope or prefix')
    return dict(construction='scheduled-lambda-bc-staggered-column-shear-final-D2',
        arithmetic='positive-FLOAT-log2-not-directed', groups=GROUPS, arity=ARITY,
        source_class=bank['source_class'], source_count=bank['source_count'], mark=bank['mark'],
        common_measure=measure, common_prefix_log2=prefix, schedule=core.frozen.SCHEDULE, normal=1,
        common_measure_receipt=core.frozen.INPUT, common_measure_sha256=core.frozen.INPUT_SHA,
        aggregate_original_input_floor=FLOOR, retained_deficit=DEFICIT,
        RS_factor_per_group='255^(-32)', RS_factor_whole_word='255^(-576)',
        original_A_fixed_across8pairs=True, shared_ratio_weights=['254/255', '1/255'],
        shared_ratio_weights_outside_norm_roots=True, raw_count_products_outside_norm_roots=True,
        helper_result_rooted_again=False, entrywise_18_power_before_chronological_matrix_products=True,
        independent_column_support_assumption=False, source_and_projective_prefix_once=True,
        extra_rowset_factor_applied=False, final_D2_after_shear=True,
        full_independent_D1_comparison=True, permitted_gauges_per_pair=1,
        restored_D1_products_per_superblock=8, regional_provider_runs=0,
        directed=False, whole_code_certificate=False)


def validate_row(row, prefix):
    ordinal = row.get('ordinal')
    if (type(ordinal) is not int or not 0 <= ordinal < len(SOURCES)
            or row.get('original_source') != list(SOURCES[ordinal])
            or any(type(i) is not int for i in row['original_source'])
            or row.get('arity') != ARITY or row.get('status') != 'complete'
            or row.get('aggregate_original_input_floor') != FLOOR):
        raise ValueError('Foreign, reordered, or incomplete original-A record')
    tail, root, bound = (row.get(name) for name in
        ('word_tail_log2', 'normalized_group_root_log2', 'fixed_A_repeated18_bound_log2'))
    if any(type(x) not in (int, float) or not isfinite(x) for x in (tail, root, bound)):
        raise ValueError('Completed records require finite FLOAT logarithms')
    if abs(root-(tail/GROUPS-8*ARITY*log2(255))) > 1e-9:
        raise ValueError('Arity-four RS normalization changed')
    if abs(bound-(prefix+GROUPS*root)) > 1e-7:
        raise ValueError('Fixed-A bound does not use the common prefix exactly once')
    return ordinal, root


def validate_batch(saved, context, pins):
    if saved.get('schema') != SCHEMA or saved.get('kind') != 'batch':
        raise ValueError('Not a shared-ratio a4 batch receipt')
    if any(saved.get(key) != value for key, value in context.items()):
        raise ValueError('Batch construction, measure, prefix, or proof scope differs')
    if saved.get('reviewed_source_sha256') != pins or saved.get('source_pins_unchanged') is not True:
        raise ValueError('Batch did not retain the current reviewed source pins')
    requested = select_ordinals(ordinals=saved['requested_ordinals'])
    if saved.get('requested_original_sources') != [list(SOURCES[i]) for i in requested]:
        raise ValueError('Requested labels do not match the explicit ordinal list')
    rows = saved['rows']
    if not isinstance(rows, list) or len(rows) > len(requested):
        raise ValueError('Malformed batch rows')
    seen = tuple(validate_row(row, context['common_prefix_log2'])[0] for row in rows)
    if seen != requested[:len(rows)]:
        raise ValueError('Batch must contain the completed prefix of its requested ordinal list')
    if (saved.get('completed_sources') != len(rows)
            or saved.get('batch_complete') is not (len(rows) == len(requested))
            or saved.get('full_a4_class_covered') is not False
            or saved.get('all_a4_bound_log2') is not None):
        raise ValueError('Inconsistent batch completeness or illicit full-class claim')
    return rows


def run_batch(*, begin=None, end=None, ordinals=None, max_seconds=120, progress=None):
    requested = select_ordinals(begin=begin, end=end, ordinals=ordinals)
    if type(max_seconds) not in (int, float) or not isfinite(max_seconds) or not 0 < max_seconds <= MAX_SECONDS:
        raise ValueError('Use a positive wall budget no greater than180 seconds')
    pins, started = dependencies(), perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Shared-ratio a4 batch exhausted its wall budget')
    bank, measure, prefix = core.frozen.load_measure()
    context = common_context(bank, measure, prefix)
    rows, error, preparations = [], None, 0
    try:
        check_time()
        prepared = core.ratio.prepare_measure(measure, groups=GROUPS)
        preparations = 1
        for ordinal in requested:
            check_time()
            source = SOURCES[ordinal]
            columns = core.prepare_source(source, prepared, check_time=check_time)
            tail = float(core.tail_backend.positive_deficit_tail(columns, GROUPS, FLOOR, check_time=check_time))
            root = tail/GROUPS-8*ARITY*log2(255)
            row = dict(ordinal=ordinal, original_source=list(source), arity=ARITY, status='complete',
                aggregate_original_input_floor=FLOOR, word_tail_log2=tail,
                normalized_group_root_log2=root, fixed_A_repeated18_bound_log2=prefix+GROUPS*root)
            validate_row(row, prefix)
            rows.append(row)
            if progress and (len(rows) % 32 == 0 or len(rows) == len(requested)):
                progress(dict(stage='completed-original-rowsets', completed=len(rows), requested=len(requested),
                    last_ordinal=ordinal, last_original_source=list(source),
                    largest_fixed_A_bound_log2=max(x['fixed_A_repeated18_bound_log2'] for x in rows),
                    elapsed_seconds=perf_counter()-started))
    except TimeoutError as exception:
        error = str(exception)
    if dependencies() != pins:
        raise ValueError('Reviewed sources or common bank changed during the batch')
    result = dict(schema=SCHEMA, kind='batch', **context,
        requested_ordinals=list(requested), requested_original_sources=[list(SOURCES[i]) for i in requested],
        ordinal_begin=requested[0], ordinal_end_exclusive=requested[-1]+1,
        contiguous_ordinals=requested == tuple(range(requested[0], requested[-1]+1)),
        rows=rows, completed_sources=len(rows), batch_complete=len(rows) == len(requested), timeout=error,
        full_a4_class_covered=False, all_a4_bound_log2=None,
        reviewed_source_sha256=pins, source_pins_unchanged=True, measure_preparations=preparations,
        elapsed_seconds=perf_counter()-started, max_seconds=max_seconds)
    validate_batch(result, context, pins)
    return result


def assemble_receipts(receipts, context, pins):
    """Validate in-memory batch data; do not trust saved subtotal bounds."""
    rows, seen = [], set()
    if not receipts:
        raise ValueError('At least one batch receipt is required')
    for saved in receipts:
        for row in validate_batch(saved, context, pins):
            ordinal = row['ordinal']
            if ordinal in seen:
                raise ValueError('Completed original-A records overlap across batches')
            seen.add(ordinal)
            rows.append(row)
    rows.sort(key=lambda row: row['ordinal'])
    missing = sorted(set(range(len(SOURCES)))-seen)
    complete = not missing
    bound = context['common_prefix_log2']+GROUPS*float(np.logaddexp2.reduce(
        [row['normalized_group_root_log2'] for row in rows])) if complete else None
    return dict(schema=SCHEMA, kind='assembly', **context, rows=rows,
        complete=complete, full_a4_class_covered=complete, all_a4_bound_log2=bound,
        expected_original_rowsets=len(SOURCES), unique_original_rowsets=len(rows),
        missing_ordinals=missing, missing_original_sources=[list(SOURCES[i]) for i in missing],
        partial_rows_are_not_full_class_bound=not complete,
        labelled_rowset_sum_once_before18th_power=True, reviewed_source_sha256=pins,
        source_pins_unchanged=True, moment_table_builds=0, batch_receipts=len(receipts))


def assemble(paths):
    pins = dependencies()
    bank, measure, prefix = core.frozen.load_measure()
    context = common_context(bank, measure, prefix)
    receipts, records = [], []
    for path in map(Path, paths):
        raw = path.read_bytes()
        receipts.append(json.loads(raw))
        records.append(dict(path=str(path.resolve()), sha256=sha256(raw).hexdigest()))
    result = assemble_receipts(receipts, context, pins)
    if dependencies() != pins or any(sha256(Path(row['path']).read_bytes()).hexdigest() != row['sha256'] for row in records):
        raise ValueError('Source pins or input receipts changed during assembly')
    result['input_receipts'] = records
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    batch = commands.add_parser('batch')
    batch.add_argument('--begin', type=int)
    batch.add_argument('--end', type=int)
    batch.add_argument('--ordinals', nargs='+', type=int)
    batch.add_argument('--max-seconds', type=float, default=120)
    batch.add_argument('--output', type=Path, required=True)
    assembly = commands.add_parser('assemble')
    assembly.add_argument('--inputs', nargs='+', type=Path, required=True)
    assembly.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists(): parser.error('Refuse to overwrite evidence')
    if args.command == 'batch':
        result = run_batch(begin=args.begin, end=args.end, ordinals=args.ordinals,
            max_seconds=args.max_seconds, progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    else:
        result = assemble(args.inputs)
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    keys = ('kind', 'batch_complete', 'completed_sources', 'complete', 'unique_original_rowsets', 'all_a4_bound_log2', 'timeout')
    print(json.dumps({key: result[key] for key in keys if key in result}, allow_nan=False), flush=True)
