"""Complete occupancy unions for unchanged RS16/t64/s16 at explicit lengths.

The setup law is the existing four-row GF16 RS[16,8] outer, independent
uniform GL16 symbol maps, and uniform group-packet and regional shuffles.
The selected t64/s16 inner uses independent uniform GL16 state updates.
Its state starts at zero, continues across regions, and is not flushed.

For K message bits, all occupancies 1..K/128 must be covered. Each component
bounds the expected number of nonzero messages of its occupancy producing
weight at most floor(2K/10). Overlapping tail bounds may be minimized before
the full first-moment sum. A sum below 2^-40 bounds the probability over setup
that minimum distance is at most that cutoff.

assemble() checks saved components and current source identities; it never
claims a fresh certificate. replay() prepares the selected maps once and
recomputes every component. Only a passing, source-stable replay is marked
as a whole-code certificate. Neither path modifies retained proof sources.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
from time import monotonic

from flint import arb, ctx
import packet_rs_k20_whole as retained
import packet_rs_lengths as lengths


q1, sparse, dense = retained.q1, retained.sparse, retained.dense
# The census imports this helper lazily during fresh map preparation.
CORE_FILENAMES = retained.CORE_FILENAMES + ('gf16_packets/feedback_exact.py',)
CONTINUITY = retained.CONTINUITY
SCHEMA = 'rs16-length-complete-occupancy-union-1'
LOW_SCHEMAS = ('finite-packet-q1-diagnostic-1', 'finite-packet-q2-diagnostic-1')
DEFAULT_EXACT = tuple(range(3, 33))


def required_source_paths(schema):
    """Require the audited shared core and each recognized producer's sources."""
    if schema in (dense.SCHEMA, dense.FUGACITY_SCHEMA):
        paths = set(retained.required_source_paths(schema))
    elif schema in (*LOW_SCHEMAS, lengths.TAIL_SCHEMA):
        paths = {str((q1.LOCALITY / name).resolve()) for name in CORE_FILENAMES}
        paths.update(str(Path(module.__file__).resolve()) for module in
                     (q1, retained.uniform, sparse))
        here = Path(__file__).resolve().parent
        paths.add(str((here / 'rs_outer.py').resolve()))
        if schema == LOW_SCHEMAS[1]:
            paths.add(str(Path(retained.q2.__file__).resolve()))
        paths.add(str(q1.kernel_t64.SELECTED_MAP.resolve()))
    else:
        raise ValueError('recognized RS16 length-component schema required')
    paths.update(str((q1.LOCALITY / name).resolve()) for name in CORE_FILENAMES)
    # The new low wrappers fix the length, setup law, and selected-map identity.
    # Dense receipts already carry those premises in their existing producer.
    if schema in (*LOW_SCHEMAS, lengths.TAIL_SCHEMA):
        paths.add(str(Path(lengths.__file__).resolve()))
        paths.add(str(Path(sparse.__file__).resolve()))
    return frozenset(paths)


def dyadic(pair):
    """Validate and normalize a positive finite saved dyadic endpoint exactly."""
    if (not isinstance(pair, list) or len(pair) != 2 or
            any(type(value) is not int for value in pair) or pair[0] <= 0 or
            abs(pair[1]) > 10_000_000):
        raise ValueError('positive bounded-exponent integer dyadic required')
    mantissa, exponent = pair
    trailing = (mantissa & -mantissa).bit_length() - 1
    return mantissa >> trailing, exponent + trailing


def dyadic_less(left, right):
    """Compare normalized positive dyadics without floating-point conversion."""
    lm, le = left
    rm, re = right
    ltop, rtop = lm.bit_length() + le, rm.bit_length() + re
    if ltop != rtop:
        return ltop < rtop
    if le >= re:
        return lm << (le - re) < rm
    return lm < rm << (re - le)


def dyadic_sum(values):
    """Add normalized positive dyadics using exact integer arithmetic."""
    buckets = {}
    for mantissa, exponent in values:
        buckets[exponent] = buckets.get(exponent, 0) + mantissa
    if not buckets:
        return 0, 0
    exponents = sorted(buckets, reverse=True)
    exponent, mantissa = exponents[0], buckets[exponents[0]]
    for next_exponent in exponents[1:]:
        mantissa = (mantissa << (exponent - next_exponent)) + buckets[next_exponent]
        exponent = next_exponent
    trailing = (mantissa & -mantissa).bit_length() - 1
    return mantissa >> trailing, exponent + trailing


def rounded_endpoint(value, precision):
    """Round an exact sum upward once, keeping its JSON mantissa compact."""
    mantissa, exponent = value
    shift = max(0, mantissa.bit_length() - precision)
    if shift:
        mantissa = (mantissa >> shift) + bool(mantissa & ((1 << shift) - 1))
        exponent += shift
    return list(dyadic([mantissa, exponent]))


def _margin(pair):
    # This string is display metadata. Acceptance compares integer dyadics.
    value = arb(pair[0]) * arb(2) ** pair[1]
    return str(-value.log() / arb(2).log())


def _read_record(item):
    if isinstance(item, dict):
        record = item
        raw = json.dumps(record, sort_keys=True, separators=(',', ':')).encode()
        reference = dict(in_memory=True, sha256=hashlib.sha256(raw).hexdigest())
    elif isinstance(item, (str, Path)):
        path = Path(item).resolve()
        raw = path.read_bytes()
        record = json.loads(raw)
        reference = dict(path=str(path), sha256=hashlib.sha256(raw).hexdigest())
    else:
        raise ValueError('component dictionaries or receipt paths required')
    if not isinstance(record, dict):
        raise ValueError('component receipt must be an object')
    return record, reference


def _check_sources(record):
    hashes = record.get('source_sha256')
    if not isinstance(hashes, dict) or not hashes:
        raise ValueError('component mathematical source pins required')
    missing = required_source_paths(record.get('schema')) - hashes.keys()
    if missing:
        raise ValueError('component omits mathematical source pins: ' + ', '.join(sorted(missing)))
    for filename, expected in hashes.items():
        if not isinstance(filename, str) or not isinstance(expected, str):
            raise ValueError('source paths and SHA256 digests must be strings')
        try:
            actual = hashlib.sha256(Path(filename).read_bytes()).hexdigest()
        except OSError as error:
            raise ValueError(f'mathematical source unavailable: {filename}') from error
        if actual != expected:
            raise ValueError(f'mathematical source changed: {filename}')


def _common(record, geometry):
    _check_sources(record)
    expected = dict(K=geometry.K, N=geometry.N, geometry=asdict(geometry),
        threshold=geometry.N // 10, distance='1/10', outer='rs16',
        groups=geometry.group_count, regions=64, group_dimension=128,
        group_output_bits=256, physical_t=64, state_bits=16, macro_t=128,
        physical_steps=geometry.N // 64, macro_steps=geometry.N // 128,
        physical_steps_per_macro=2, macros_per_region=geometry.macros_per_region)
    if record.get('schema') == lengths.TAIL_SCHEMA:
        # This producer stores these facts in its mandatory geometry object,
        # without redundant top-level aliases. Check any aliases if supplied.
        for key in ('groups', 'regions', 'group_dimension', 'group_output_bits'):
            if key not in record:
                del expected[key]
    if any(record.get(key) != value or
           (type(value) is int and type(record.get(key)) is not int)
           for key, value in expected.items()):
        raise ValueError('matching half-rate RS16/t64/s16 finite geometry and cutoff required')
    if (record.get('zero_initial_state') is not True or record.get('final_flush') is not False or
            record.get('whole_code_certificate') is not False or
            record.get('fresh_computation', record.get('schema') == dense.SCHEMA) is not True or
            type(record.get('precision')) is not int or record['precision'] < 128):
        raise ValueError('complete noncertificate component with exact continuous-state scope required')
    lengths.authenticate_record(record.get('map_record'))


def _tail_occupancies(record, groups):
    schema = record['schema']
    if schema == lengths.TAIL_SCHEMA:
        occupancies = record.get('occupancy_covered')
        if (not isinstance(occupancies, list) or not occupancies or
                any(type(q) is not int or not 3 <= q <= groups for q in occupancies) or
                occupancies != sorted(set(occupancies))):
            raise ValueError('explicit sorted distinct exact-tail occupancies required')
    else:
        first, last = record.get('q_min'), record.get('q_max')
        if (type(first) is not int or type(last) is not int or not 3 <= first <= last <= groups or
                record.get('occupancy_covered') != [first, last] or
                record.get('evaluated_every_integer_occupancy') is not True):
            raise ValueError('explicit complete dense-tail interval required')
        occupancies = list(range(first, last + 1))
    if (not isinstance(record.get('occupancy_uppers'), dict) or
            set(record['occupancy_uppers']) != {str(q) for q in occupancies}):
        raise ValueError('exactly one endpoint for every claimed occupancy required')
    return occupancies


def assemble(records, K, *, precision=256):
    """Validate and sum complete saved coverage, minimizing overlapping tails.

    Source hashes authenticate current producer identities, not the numerical
    truth of a saved JSON file. This function therefore never marks a replay.
    """
    geometry = lengths.geometry(K)
    if type(precision) is not int or precision < 128:
        raise ValueError('integer precision >=128 required')
    if isinstance(records, (str, bytes, dict)):
        raise ValueError('a sequence of component records or paths required')
    records = list(records)
    if not records:
        raise ValueError('nonempty component sequence required')
    ctx.prec = precision
    sources = q1.source_snapshot()
    beta, counts = sparse.exact_outer()
    count_hash = sparse.count_hash(counts)
    values, selections, references, summaries, precisions = {}, {}, [], [], []
    map_record = None
    for index, item in enumerate(records):
        record, reference = _read_record(item)
        _common(record, geometry)
        if map_record is None:
            map_record = record['map_record']
        elif record['map_record'] != map_record:
            raise ValueError('identical authenticated selected-map records required')
        schema = record['schema']
        if schema in LOW_SCHEMAS:
            q = LOW_SCHEMAS.index(schema) + 1
            if (q in values or record.get('occupancy_covered') != [q] or
                    record.get('length_component') != f'q{q}' or
                    record.get('count_kind') != ('shells' if q == 1 else 'exact_expected_shells') or
                    record.get('count_sha256') != sparse.count_hash(counts, cumulative=q == 1) or
                    (q == 2 and record.get('all_two_group_support_pairs_covered') is not True)):
                raise ValueError('one complete exact-shell q1/q2 length component required')
            endpoints = {q: dyadic(record.get(f'q{q}_upper'))}
        else:
            if (record.get('every_shell_checked') is not True or
                    record.get('state_continuity') != CONTINUITY or
                    record.get('count_sha256') != count_hash):
                raise ValueError('complete continuous-state outer-envelope tail required')
            try:
                matching_beta = Q(record.get('beta')) == beta
            except (TypeError, ValueError, ZeroDivisionError, OverflowError):
                matching_beta = False
            if not matching_beta:
                raise ValueError('current pointwise RS16 outer comparison factor required')
            endpoints = {q: dyadic(record['occupancy_uppers'][str(q)])
                         for q in _tail_occupancies(record, geometry.group_count)}
        for q, value in endpoints.items():
            if q not in values or dyadic_less(value, values[q]):
                values[q], selections[q] = value, index
        references.append(reference)
        precisions.append(record['precision'])
        summaries.append(dict(schema=schema, occupancy_covered=record['occupancy_covered'],
                              endpoint_count=len(endpoints)))
    if set(values) != set(range(1, geometry.group_count + 1)):
        missing = sorted(set(range(1, geometry.group_count + 1)) - values.keys())
        raise ValueError('incomplete occupancy union; first missing occupancies: ' + str(missing[:16]))
    exact_total = dyadic_sum(values.values())
    exact_tail = dyadic_sum(values[q] for q in range(3, geometry.group_count + 1))
    upper, tail = rounded_endpoint(exact_total, precision), rounded_endpoint(exact_tail, precision)
    if sources != q1.source_snapshot():
        raise RuntimeError('mathematical source changed during union assembly')
    target_met = dyadic_less(dyadic(upper), (1, -40))
    return dict(schema=SCHEMA, outer='rs16', K=geometry.K, N=geometry.N,
        geometry=asdict(geometry), threshold=geometry.N // 10, distance='1/10',
        target_margin_bits=40, target_minimum_distance=geometry.N // 10 + 1,
        state_bits=16, physical_t=64, zero_initial_state=True, final_flush=False,
        state_continuity=CONTINUITY, occupancy_covered=[1, geometry.group_count],
        all_occupancies_covered=True, fresh_replay=False, whole_code_certificate=False,
        target_met=target_met, precision=precision, minimum_component_precision=min(precisions),
        map_record=map_record, input_receipts=references, components=summaries,
        occupancy_uppers={str(q): list(values[q]) for q in sorted(values)},
        selected_components={str(q): selections[q] for q in sorted(values)},
        q1_upper=list(values[1]), q2_upper=list(values[2]), tail_upper=tail, union_upper=upper,
        margin_bits=_margin(upper), tail_margin_bits=_margin(tail),
        arithmetic='Exact integer sum of selected dyadic endpoints; one final upward rounding.',
        source_sha256=sources,
        scope='First-moment union over all nonzero messages in the unchanged ideal RS16/t64/s16 '
              'ensemble. A checked saved-component sum is not a fresh whole-code certificate.')


def recipe(K, *, exact_occupancies=DEFAULT_EXACT, dense_min=3,
           q1_tilts=retained.DEFAULT_Q1, q2_tilts=retained.DEFAULT_Q2,
           exact_tilts=retained.DEFAULT_SPARSE, dense_tilts=retained.DEFAULT_DENSE,
           marker_probabilities=dense.DEFAULT_MARKER_PROBABILITIES, precision=256):
    """Validate a complete replay plan before preparing maps or creating files."""
    geometry = lengths.geometry(K)
    if type(precision) is not int or precision < 256:
        raise ValueError('fresh whole replay requires integer precision >=256')
    if isinstance(exact_occupancies, (str, bytes)):
        raise ValueError('explicit exact-tail occupancy sequence required')
    occupancies = list(exact_occupancies)
    if (any(type(q) is not int or not 3 <= q <= geometry.group_count for q in occupancies) or
            occupancies != sorted(set(occupancies))):
        raise ValueError('sorted distinct valid exact-tail occupancies required')
    if type(dense_min) is not int or not 3 <= dense_min <= geometry.group_count:
        raise ValueError('dense interval must begin in 3..L')
    if set(occupancies) | set(range(dense_min, geometry.group_count + 1)) != set(range(3, geometry.group_count + 1)):
        raise ValueError('replay plan must cover every tail occupancy')
    tilts = {name: list(sparse.checked_options(values, precision, None)[0])
             for name, values in (('q1', q1_tilts), ('q2', q2_tilts),
                                  ('exact', exact_tilts), ('dense', dense_tilts))}
    if isinstance(marker_probabilities, (str, bytes)):
        raise ValueError('a sequence of marker probabilities required')
    try:
        markers = tuple(map(Q, marker_probabilities))
    except (TypeError, ValueError, ZeroDivisionError, OverflowError) as error:
        raise ValueError('finite rational markers required') from error
    if (not markers or len(set(markers)) != len(markers) or any(not 0 < p <= 1 for p in markers) or
            (dense_min < geometry.group_count and not any(p < 1 for p in markers))):
        raise ValueError('distinct valid markers covering the dense interval required')
    return dict(K=K, N=geometry.N, geometry=asdict(geometry), threshold=geometry.N // 10,
        precision=precision, exact_occupancies=occupancies,
        dense_interval=[dense_min, geometry.group_count], tilts=tilts,
        marker_probabilities=list(map(str, markers)))


def _write_new(output, record):
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('x', encoding='utf-8') as stream:
        stream.write(json.dumps(record, indent=2) + '\n')


def replay(output, *, K, **options):
    """Freshly recompute all witnesses from one shared selected-map preparation."""
    start = monotonic()
    plan = recipe(K, **options)
    output = Path(output)
    names = ['q1', 'q2'] + (['exact'] if plan['exact_occupancies'] else []) + ['dense']
    paths = {name: output.with_name(output.stem + '-' + name + '.json') for name in names}
    if output.exists() or any(path.exists() for path in paths.values()):
        raise ValueError('whole replay and every component output must be fresh')
    precision = plan['precision']
    sources_before = q1.source_snapshot()
    data, map_record, sources = lengths.prepare(precision=precision)
    # Fresh preparation can import additional proof helpers lazily. Preserve
    # every earlier pin, then freeze the complete post-preparation snapshot.
    if (any(sources.get(path) != digest for path, digest in sources_before.items()) or
            sources != q1.source_snapshot()):
        raise RuntimeError('mathematical source changed during shared preparation')
    shared = dict(precision=precision, data=data, map_record=map_record)
    records = []
    for q in (1, 2):
        records.append(lengths.run_low(q, K=K, tilts=plan['tilts'][f'q{q}'],
                                      output=paths[f'q{q}'], **shared))
    if plan['exact_occupancies']:
        records.append(lengths.run_tail(K=K, occupancies=plan['exact_occupancies'],
            tilts=plan['tilts']['exact'], output=paths['exact'], **shared))
    records.append(dense.run_fugacity(group_count=plan['geometry']['group_count'],
        q_min=plan['dense_interval'][0], q_max=plan['dense_interval'][1],
        tilts=plan['tilts']['dense'], marker_probabilities=plan['marker_probabilities'],
        output=paths['dense'], **shared))
    result = assemble(records, K, precision=precision)
    if sources != q1.source_snapshot() or result['map_record'] != map_record:
        raise RuntimeError('mathematical source or selected map changed during fresh replay')
    result.update(fresh_replay=True, whole_code_certificate=result['target_met'],
                  recipe=plan, elapsed_seconds=monotonic() - start)
    _write_new(output, result)
    print(f'RS16 K={K} COMPLETE UNION margin_bits={result["margin_bits"]} '
          f'whole_code_certificate={result["whole_code_certificate"]}', flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--K', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--exact-occupancies', type=int, nargs='*', default=list(DEFAULT_EXACT))
    parser.add_argument('--dense-min', type=int, default=3)
    for name, default in (('q1', retained.DEFAULT_Q1), ('q2', retained.DEFAULT_Q2),
                          ('exact', retained.DEFAULT_SPARSE), ('dense', retained.DEFAULT_DENSE)):
        parser.add_argument(f'--{name}-tilts', nargs='+', default=list(default))
    parser.add_argument('--marker-probabilities', nargs='+', default=list(dense.DEFAULT_MARKER_PROBABILITIES))
    parser.add_argument('--receipts', type=Path, nargs='+')
    parser.add_argument('--dry-run', action='store_true')
    args = vars(parser.parse_args())
    output, receipts, dry = args.pop('output'), args.pop('receipts'), args.pop('dry_run')
    if receipts is not None:
        if dry:
            parser.error('--dry-run cannot be combined with --receipts')
        if output.exists():
            raise ValueError('assembly output must be fresh')
        _write_new(output, assemble(receipts, args['K'], precision=args['precision']))
    elif dry:
        print(json.dumps(recipe(**args), indent=2))
    else:
        replay(output, **args)


if __name__ == '__main__':
    main()
