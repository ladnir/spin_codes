"""Validate, repair, and aggregate an explicitly complete sparse receipt.

Input hashes establish provenance and agreement with the retained sources;
they do not authenticate a receipt or replace its original numerical replay.
Each repair uses the frozen selected-weight verifier with its own witness.
The minimum of compatible directed bounds is valid for the same event.
The final union includes each requested weight exactly once. No result here
covers weights outside that interval or constitutes a whole-code certificate.
"""

import ast
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from math import isfinite
from pathlib import Path
from time import monotonic

from flint import arb, ctx

import hd_field128_lift_certificate as cert
import hd_field128_lift_sparse_bands as bands


DEFAULT_REPAIRS = (3, 6, 7, 11, 12, 13, 14, 41)
MAX_INPUT_BYTES = 20 << 20
FROZEN_ENTRYPOINTS = {
    'hd_field128_lift_certificate.py': 'd23a3247fcc1162107821adff91b076b9eb8d0ce8f795719bf5e5ea907cd5473',
    'hd_field128_lift_sparse_bands.py': 'e2f52461a2d5f827ff8a9458918cf05d112cbb460eaaccc9b618d6b87ebd7972',
}


def source_pins():
    directory = Path(__file__).resolve().parent
    pending, seen = [Path(__file__).name], {}
    while pending:
        name = pending.pop()
        if name in seen:
            continue
        path = directory/name
        raw = path.read_bytes()
        seen[name] = sha256(raw).hexdigest()
        for node in ast.walk(ast.parse(raw, filename=str(path))):
            modules = ([item.name for item in node.names] if isinstance(node, ast.Import)
                       else [node.module] if isinstance(node, ast.ImportFrom) and node.module else [])
            for module in modules:
                target = module.split('.')[0]+'.py'
                if (directory/target).is_file() and target not in seen:
                    pending.append(target)
    return dict(sorted(seen.items()))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def same(record, expected, label):
    require(isinstance(record, dict), label+' must be an object')
    for key, value in expected.items():
        require(key in record and type(record[key]) is type(value) and record[key] == value,
                label+': incompatible '+key)


def rational_string(value):
    require(isinstance(value, str) and 0 < len(value) <= 4096, 'Expected a bounded rational string')
    return F(value)


def precision_parameters(precision, bits):
    require(type(precision) is int and 128 <= precision <= 4096, 'Invalid Arb precision')
    require(type(bits) is int and 32 <= bits <= 512, 'Invalid packet mantissa budget')


def common_semantics(q, distance, outer_measure):
    return dict(construction=cert.CONSTRUCTION, field_order=q, field_bits=cert.checked_q(q),
                target_distance=str(distance), forbidden_weight=int(distance*cert.OUTPUT_LENGTH),
                outer_measure=outer_measure, directed=True, whole_code_certificate=False)


def full_semantics(q, distance, outer_measure):
    return dict(common_semantics(q, distance, outer_measure), source_dimension=cert.DIMENSION,
                output_length=cert.OUTPUT_LENGTH, outer_groups=cert.GROUPS,
                projective_divisor=str(q-1), pair_swaps=True,
                state_persists_across_regions=True, full_field_random_coefficients=True)


@lru_cache(64)
def profile_count(r):
    return sum(bool(cert.profile_count_upper(cert.PAIRS, r, a, b))
               for a in range(r//2+1) for b in range(r//2+1))


def validate_weight_row(row, r, maximum, witness, precision, bits):
    same(row, dict(active_groups=r, profiles=profile_count(r),
                   covers_all_supports_at_weight=True, covers_all_even_odd_splits=True,
                   directed=True, whole_code_certificate=False), 'weight row')
    same(row.get('witness'), dict(witness, coefficient_maximum=maximum,
                                 coefficient_marker=str(F(maximum, cert.PAIRS))), 'weight witness')
    z, alpha = map(rational_string, (witness['z'], witness['alpha']))
    require(0 < z < 1 and 0 < alpha <= 1, 'Invalid rational witness')
    precision_parameters(precision, bits)
    bound = rational_string(row.get('log2_upper_rational'))
    enclosure = row.get('log2_upper_enclosure')
    require(isinstance(enclosure, str) and len(enclosure) <= 16384, 'Missing bounded scalar enclosure')
    with ctx.workprec(precision):
        value = arb(enclosure)
        require(value.is_finite() and cert.rational(bound) >= value.upper(),
                'Stored scalar ceiling does not dominate its enclosure')
    return bound


def union_ceiling(rows, precision):
    require(bool(rows), 'Cannot certify an empty sparse union')
    with ctx.workprec(precision):
        mass = sum((arb(2)**cert.rational(rational_string(row['log2_upper_rational']))
                    for row in rows), arb(0))
        return cert.rational_ceiling(mass.log()/arb(2).log())


def validate_schedule(receipt, *, first=1, last=64, q=cert.FIELD_ORDER, distance=F(7, 20),
                      outer_measure='invertible-toeplitz'):
    """Check provenance, compatible semantics, and every indexed input row.

    Scalar enclosure and union checks are consistency checks, not a replay
    of the regional coefficient calculation that produced the input.
    """
    distance = cert.exact(distance)
    require(type(first) is int and type(last) is int and 1 <= first <= last <= 64,
            'Invalid requested sparse interval')
    require(0 < distance < F(1, 2) and outer_measure in cert.OUTER_MEASURES, 'Invalid target semantics')
    pins = bands.source_pins()
    require(all(pins.get(name) == digest for name, digest in FROZEN_ENTRYPOINTS.items()),
            'A frozen verifier entrypoint changed')
    same(receipt, dict(common_semantics(q, distance, outer_measure),
                       requested_weight_range=[first, last], completed_weight_range=[first, last],
                       complete_requested_range=True, pending_weight_ranges=[], failure=None,
                       source_pins_unchanged_at_exit=True, source_sha256=pins), 'input schedule')
    precision, bits = receipt.get('precision'), receipt.get('packet_mantissa_bits')
    precision_parameters(precision, bits)
    completed, planned = receipt.get('completed_bands'), receipt.get('planned_bands')
    require(isinstance(completed, list) and isinstance(planned, list) and bool(completed), 'Missing bands')
    same(receipt, dict(completed_band_count=len(completed), regional_tables_in_completed_bands=len(completed),
                       regional_table_builds_attempted=len(completed)), 'input schedule counts')
    require(len(planned) == len(completed), 'Plan/completion band count differs')
    rows, next_weight = [], first
    for index, (band, plan) in enumerate(zip(completed, planned)):
        require(isinstance(band, dict) and isinstance(plan, dict), 'Invalid band object')
        interval = band.get('weight_range')
        require(isinstance(interval, list) and len(interval) == 2
                and all(type(x) is int for x in interval), 'Invalid band range')
        lo, hi = interval
        require(lo == next_weight and lo <= hi <= last, 'Band gap, overlap, or out-of-range coverage')
        maximum, center = band.get('coefficient_maximum'), band.get('witness_center')
        require(type(maximum) is int and hi <= maximum <= 64, 'Invalid coefficient maximum')
        require(type(center) is int and 1 <= center <= 64, 'Invalid witness center')
        same(plan, dict(center=center, weight_range=interval, coefficient_maximum=maximum), 'planned band')
        same(band, dict(full_semantics(q, distance, outer_measure), source_sha256=pins,
                        coefficient_marker=str(F(maximum, cert.PAIRS)), regional_tables_built=1,
                        weights_evaluated=hi-lo+1, every_weight_evaluated=True,
                        covers_all_supports_in_range=True, covers_all_even_odd_splits=True,
                        coefficient_method='shared-positive-carry-free-Arb-truncated-bivariate'), 'input band')
        witness = band.get('witness')
        same(witness, dict(precision=precision, packet_mantissa_bits=bits, paired=True), 'band witness')
        children = band.get('weight_receipts')
        require(isinstance(children, list) and len(children) == hi-lo+1, 'Incomplete band rows')
        for offset, row in enumerate(children):
            r = lo+offset
            validate_weight_row(row, r, maximum, witness, precision, bits)
            rows.append(dict(row, origin=dict(kind='input-band-row', band_index=index, row_index=offset)))
        require(rational_string(band.get('log2_upper_rational')) >= union_ceiling(children, precision),
                'Stored band union does not dominate its rows')
        next_weight = hi+1
    require(next_weight == last+1 and [row['active_groups'] for row in rows] == list(range(first, last+1)),
            'Input does not cover each requested weight exactly once')
    require(rational_string(receipt.get('completed_union_log2_upper_rational')) >= union_ceiling(completed, precision),
            'Stored schedule union does not dominate its bands')
    return rows


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, 'Duplicate JSON object key: '+key)
        result[key] = value
    return result


def load_receipt(path):
    path = Path(path).resolve()
    require(path.is_file() and path.stat().st_size <= MAX_INPUT_BYTES, 'Missing or oversized input receipt')
    raw = path.read_bytes()
    require(len(raw) <= MAX_INPUT_BYTES, 'Oversized input receipt')
    reject_constant = lambda value: (_ for _ in ()).throw(ValueError('Nonfinite JSON number: '+value))
    receipt = json.loads(raw.decode('utf-8'), object_pairs_hook=_unique_object, parse_constant=reject_constant)
    return receipt, dict(path=str(path), sha256=sha256(raw).hexdigest(), bytes=len(raw))


def repair_provenance(receipts):
    result = []
    for index, receipt in enumerate(receipts):
        raw = json.dumps(receipt, sort_keys=True, separators=(',', ':'), ensure_ascii=True, allow_nan=False).encode('utf-8')
        result.append(dict(repair_index=index, active_groups=receipt['active_groups'], sha256=sha256(raw).hexdigest(),
                           encoding='UTF-8 canonical JSON: sorted keys, comma/colon separators, ensure_ascii=True'))
    return result


def validate_repair(row, r, *, q, distance, outer_measure, precision, bits, source_hashes):
    same(row, dict(full_semantics(q, distance, outer_measure), source_sha256=source_hashes,
                   coefficient_method='positive-carry-free-Arb-truncated-bivariate'), 'repair receipt')
    z, alpha = bands.default_witness(r, q)
    witness = dict(z=str(z), alpha=str(alpha), precision=precision, packet_mantissa_bits=bits, paired=True)
    validate_weight_row(row, r, r, witness, precision, bits)


def repair_schedule(input_path, *, repairs=DEFAULT_REPAIRS, first=1, last=64, q=cert.FIELD_ORDER,
                    distance=F(7, 20), precision=256, bits=80, outer_measure='invertible-toeplitz',
                    max_seconds=3600., max_weight_seconds=600., failure_bits=40, checkpoint=None):
    """Serial repairs; snapshots preserve all original coverage and finished replays."""
    distance = cert.exact(distance)
    precision_parameters(precision, bits)
    require(isfinite(max_weight_seconds) and 0 < max_weight_seconds <= 3600, 'Invalid per-weight deadline')
    require(type(failure_bits) is int and 0 <= failure_bits <= 4096, 'Invalid failure-bit target')
    repairs = tuple(repairs)
    require(all(type(r) is int and first <= r <= last for r in repairs)
            and len(set(repairs)) == len(repairs), 'Repairs must be distinct requested weights')
    budget = cert.Budget(max_seconds)
    receipt, provenance = load_receipt(input_path)
    original = validate_schedule(receipt, first=first, last=last, q=q, distance=distance, outer_measure=outer_measure)
    pins, replay_pins = source_pins(), cert.source_pins()
    selected = {row['active_groups']: row for row in original}
    completed, failure = [], None

    def snapshot():
        rows = [selected[r] for r in range(first, last+1)]
        ceiling = union_ceiling(rows, precision)
        return dict(full_semantics(q, distance, outer_measure),
                    requested_weight_range=[first, last], selected_weight_receipts=rows,
                    every_weight_evaluated=True, complete_requested_range=True,
                    covers_all_supports_in_range=True, covers_all_even_odd_splits=True,
                    selected_weight_count=len(rows), duplicate_weight_count=0,
                    log2_upper_rational=str(ceiling), requested_range_meets_target=ceiling <= -failure_bits,
                    failure_bits=failure_bits, requested_repairs=list(repairs),
                    completed_repairs=completed.copy(), pending_repairs=list(repairs[len(completed):]),
                    completed_repair_provenance=repair_provenance(completed),
                    all_requested_repairs_completed=len(completed) == len(repairs), failure=failure,
                    input_provenance=provenance, input_receipt=receipt, source_sha256=pins,
                    repair_source_sha256=replay_pins, precision=precision, packet_mantissa_bits=bits,
                    max_seconds=max_seconds, max_weight_seconds=max_weight_seconds,
                    elapsed_seconds=monotonic()-budget.started,
                    selection_rule='minimum compatible directed rational ceiling for each weight',
                    input_validation='source pins, semantics, coverage, scalar enclosures, directed union consistency; no full input replay',
                    scope='all supports at each listed sparse weight, exactly once; no other weights claimed')

    if checkpoint:
        checkpoint(snapshot())
    for r in repairs:
        try:
            budget.check()
            remaining = min(max_weight_seconds, max_seconds-(monotonic()-budget.started))
            if remaining <= 0:
                raise TimeoutError('Sparse repair schedule exhausted its deadline')
            z, alpha = bands.default_witness(r, q)
            row = cert.weight_receipt(r, z, alpha, q=q, distance=distance, precision=precision,
                                      bits=bits, paired=True, outer_measure=outer_measure,
                                      max_seconds=remaining, _pins=replay_pins)
            validate_repair(row, r, q=q, distance=distance, outer_measure=outer_measure,
                            precision=precision, bits=bits, source_hashes=replay_pins)
            candidate = dict(row, origin=dict(kind='own-weight-repair', repair_index=len(completed)))
            if rational_string(candidate['log2_upper_rational']) < rational_string(selected[r]['log2_upper_rational']):
                selected[r] = candidate
            completed.append(row)
            if checkpoint:
                checkpoint(snapshot())
        except (TimeoutError, ArithmeticError, RuntimeError) as error:
            failure = dict(active_groups=r, error_type=type(error).__name__, message=str(error))
            break
    require(source_pins() == pins and cert.source_pins() == replay_pins,
            'Pinned source changed during repair/aggregation')
    result = snapshot()
    result['source_pins_unchanged_at_exit'] = True
    result['input_file_unchanged_at_exit'] = sha256(Path(input_path).read_bytes()).hexdigest() == provenance['sha256']
    if checkpoint:
        checkpoint(result)
    return result


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True, type=Path)
    parser.add_argument('--repair', nargs='+', type=int, default=list(DEFAULT_REPAIRS))
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--bits', type=int, default=80)
    parser.add_argument('--max-seconds', type=float, default=3600.)
    parser.add_argument('--max-weight-seconds', type=float, default=600.)
    parser.add_argument('--failure-bits', type=int, default=40)
    parser.add_argument('--checkpoint-dir', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists() or (args.checkpoint_dir and args.checkpoint_dir.exists()):
        parser.error('Output/checkpoint paths must be unused')
    if args.checkpoint_dir:
        args.checkpoint_dir.mkdir(parents=False, exist_ok=False)
    sequence = 0

    def save(result):
        nonlocal sequence
        if args.checkpoint_dir:
            with (args.checkpoint_dir/f'snapshot_{sequence:03d}.json').open('x', encoding='utf-8') as stream:
                json.dump(result, stream, indent=2, allow_nan=False)
                stream.write('\n')
            sequence += 1
        print(json.dumps(dict(completed_repairs=len(result['completed_repairs']), pending_repairs=result['pending_repairs'],
                              log2_upper_rational=result['log2_upper_rational'], failure=result['failure'])), flush=True)

    result = repair_schedule(args.input, repairs=args.repair, precision=args.precision, bits=args.bits,
                             max_seconds=args.max_seconds, max_weight_seconds=args.max_weight_seconds,
                             failure_bits=args.failure_bits, checkpoint=save)
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(output=str(args.output), log2_upper_rational=result['log2_upper_rational'],
                          requested_range_meets_target=result['requested_range_meets_target'],
                          all_requested_repairs_completed=result['all_requested_repairs_completed'])), flush=True)


if __name__ == '__main__':
    main()
