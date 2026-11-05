"""Complete literal-GF128/35% cover from reviewed directed receipt producers.

The input is an explicit partial-index file. Its pinned receipt list is
loaded again, and the frozen index producer reconstructs every selected row.
This wrapper independently checks a disjoint cover of weights 1..8192 and
recomputes the probability union at 256 and 384 bits. The terminal range is
one event, never one probability term per covered weight.

The resulting whole-code statement is conditional on the validity of the
receipts from the listed reviewed producers. This is a provenance, scope,
coverage, and scalar-union replay, not authentication or a full replay of
regional coefficients. No existing partial index is modified or promoted.
"""

from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from importlib.metadata import version
import json
from pathlib import Path

from flint import arb, ctx

import hd_field128_lift_partial_index as index


DISTANCE, FAILURE_BITS = F(7, 20), 40
PRECISIONS = (256, 384)
INDEX = 'hd_field128_lift_partial_index.py'
# Each pair pins the entrypoint bytes and its entire static local-import
# closure. The latter digest is SHA256 of canonical JSON {filename: digest}.
# These are the reviewed literal-h4 producers, not the full-RS candidates.
REVIEWED = {
    INDEX: ('1d63c6a2fce35bfbf3ba596b9620c309fe76fd96397e8eca3b16fadb4c19b166',
            'ccc5e8ed6baa096739e518c62f111021f3c82a6460b7d0393e0d0decaafed079'),
    'hd_field128_lift_certificate.py': ('d23a3247fcc1162107821adff91b076b9eb8d0ce8f795719bf5e5ea907cd5473',
                                        '46fb32a52bebce9abec93f55e1ba56fa0dcd361e0059a7fd89c9ab22cd7c925b'),
    'hd_field128_lift_sparse_repair.py': ('e26796cb12ea060ac9c8bc62bcb77db2ae825a01e9b87b265564bce4f6efcecc',
                                          '5422c6636d2a8550a4b913bcf78435e2436fa2b540188d91efe626dbcdbe9dad'),
    'hd_field128_lift_sparse_wide.py': ('ddddbba5d379d7369a2bd0237087da619f9bf74b0b6638e600e24bb8997c4372',
                                        'b625701a9360be50f810701964833af8e48411c2d61da681509591d5b0ba79cb'),
    'hd_field128_lift_sparse_bands.py': ('e2f52461a2d5f827ff8a9458918cf05d112cbb460eaaccc9b618d6b87ebd7972',
                                         '158eccc43f6eeebe4aef3548c4bad1eca487eed3e5b581b4496651a4014518bb'),
    'hd_field128_lift_alias_cover.py': ('92f5ad1d178d5a88ab05fa1de6c2e18af4a288debf241d9cd4846cbd87665ee2',
                                        '1419378dca5ab27e820cf1d75a6d28276e3e7b4b43b7dd7ef41e56a6b3d51e4c'),
    'hd_field128_lift_alias_batch.py': ('1d4e97fb73e16e76cd5e464d606b19ae8a61b43f108c407e23048e338c650cbd',
                                        'b5aa061e8ffa9732c7e9d7bd348004fb56e960e0b2916f86f095f699cd0e3279'),
    'hd_field128_lift_tail_certificate.py': ('ad5ebdd41456ac940e9321d6ad5dae28906f1e16efbfab4fac848992af30cc61',
                                             'fa5fa4ee0174c542c5addc7c61600816ca43b9f41dffba0ab59d054532c65923'),
    'hd_field128_lift_alias_certificate.py': ('b6dadc8898ee77f52c1ade0d93c98aa46b069a8360a6264355c1fc7ad9829676',
                                              '58a9e752042d0538abd4131e31d9754b58a677f518ab9ffc15aa5ebea436ea42'),
    'hd_field128_lift_acb_alias_certificate.py': ('5883234a792bab9f850d437869ed837a8eb4fd1be7029ce308e72445dc26abf9',
                                                  '37726fffb4633fe1cf258d3d302f5bbb36e6dfa6316b8d97a9e787fdbc017002'),
}


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode('utf-8')


def digest(value):
    return sha256(canonical(value)).hexdigest()


def reviewed_closures():
    """Check fixed reviewed pins before executing the receipt reconstruction."""
    closures, sources = {}, {}
    for producer, (root_pin, closure_pin) in REVIEWED.items():
        closure = index.source_closure(producer)
        index.require(closure.get(producer) == root_pin and digest(closure) == closure_pin,
                      'Reviewed producer or transitive helper changed: ' + producer)
        for name, pin in closure.items():
            index.require(name not in sources or sources[name] == pin, 'Conflicting reviewed source: ' + name)
            sources[name] = pin
        closures[producer] = closure
    return closures, sources


def validate_reviewed_sources(record, closures, sources):
    index.same(record, dict(source_sha256=closures[INDEX]), 'index source closure')
    used = record.get('checked_producer_closures')
    index.require(isinstance(used, dict), 'Missing checked producer closures')
    merged = dict(closures[INDEX])
    for name, pins in used.items():
        index.require(name in closures and pins == closures[name], 'Unreviewed or conflicting producer: ' + name)
        for source, pin in pins.items():
            index.require(source not in merged or merged[source] == pin, 'Conflicting producer source: ' + source)
            merged[source] = pin
    index.require(record.get('checked_source_sha256') == dict(sorted(merged.items())),
                  'Indexed source union differs from the reviewed closures')
    index.require(all(sources.get(name) == pin for name, pin in merged.items()), 'Source-pin conflict')


def event_core(row):
    return {key: row[key] for key in ('active_groups', 'log2_upper_rational', 'provenance')}


def validate_cover(record):
    """Check every weight independently of the index's completion flags."""
    index.same(record, dict(index.Validator(DISTANCE).full(), schema='field128-lift-partial-index-v1'),
               'literal35 index semantics')
    index.same(record, dict(indexed_weight_domain=[1, index.GROUPS], all_weights_covered=True,
                           passing_weight_count=index.GROUPS, passing_weight_ranges=[[1, index.GROUPS]],
                           unresolved_weight_ranges=[], unresolved_weights=[], missing_weight_ranges=[],
                           missing_weights=[], known_above_target_weight_ranges=[], known_above_target_weights=[],
                           whole_code_threshold_conditions_met=True, selected_union_meets_aggregate_target=True),
               'complete index coverage')
    rows, events = record.get('selected_per_weight'), record.get('union_events')
    index.require(isinstance(rows, list) and isinstance(events, list), 'Missing indexed rows or events')
    target = record.get('individual_event_failure_bits')
    index.require(type(target) is int and 0 <= target <= 4096, 'Invalid individual event target')
    selected = {}
    for row in rows:
        weight = row.get('active_groups')
        index.require(type(weight) is int and 1 <= weight <= index.GROUPS and weight not in selected,
                      'Duplicate or invalid selected weight')
        selected[weight] = row
    covered, singles, bounds, terminal_count = set(), set(), [], 0
    terminal = record.get('terminal_event')
    for event in events:
        weight = event.get('active_groups')
        bound = index.fraction(event.get('log2_upper_rational'))
        index.require(bound <= -target, 'Union event does not meet its individual target')
        if weight is None:
            terminal_count += 1
            index.require(terminal_count == 1, 'The terminal union event was counted twice')
            index.require(isinstance(terminal, dict) and event_core(event) == event_core(terminal),
                          'Terminal event differs from the selected terminal receipt')
            index.same(event, dict(weight_range=index.TERMINAL_RANGE), 'terminal event range')
            index.same(terminal, dict(weight_range=index.TERMINAL_RANGE, used_in_union=True,
                                     passes_individual_target=True), 'selected terminal')
            domain = set(range(index.TERMINAL_RANGE[0], index.TERMINAL_RANGE[1] + 1))
        else:
            index.require(type(weight) is int and weight in selected, 'Union event lacks a selected weight')
            index.require(event_core(event) == event_core(selected[weight]), 'Union event differs from its selected row')
            index.same(selected[weight], dict(used_in_union=True, passes_individual_target=True), 'selected event')
            domain = {weight}
            singles.add(weight)
        index.require(not covered.intersection(domain), 'Probability union contains overlapping weight events')
        covered.update(domain)
        bounds.append(bound)
    index.require(covered == set(range(1, index.GROUPS + 1)), 'The disjoint union does not cover all8192weights')
    for weight, row in selected.items():
        index.same(row, dict(used_in_union=weight in singles), 'selected row union use')
    if terminal is not None:
        index.same(terminal, dict(used_in_union=bool(terminal_count)), 'terminal union use')
    index.same(record, dict(union_event_count=len(events), terminal_union_event_count=terminal_count), 'union counts')
    return bounds


def replay_union(bounds, precision):
    """Directed log2 union from exact row ceilings; no binary64 arithmetic."""
    index.require(precision in PRECISIONS and bool(bounds), 'Invalid union replay')
    # Repeated ceilings can be combined as exact integer multiplicities.
    counts = Counter(bounds)
    with ctx.workprec(precision):
        total = sum((count * arb(2)**index.rational_arb(bound) for bound, count in sorted(counts.items())), arb(0))
        value = total.log()/arb(2).log()
        index.require(value.is_finite(), 'Nonfinite complete union')
        ceiling = F(int((value.upper()*10**12).upper().ceil().fmpz()), 10**12)
        return dict(precision=precision, event_count=len(bounds), distinct_ceilings=len(counts),
                    log2_upper_rational=str(ceiling), log2_upper_enclosure=str(value))


def build_manifest(index_path):
    """Create a new complete-cover record or fail; never mutate any receipt."""
    closures, reviewed_sources = reviewed_closures()
    own_path = Path(__file__).resolve()
    own_pin = sha256(own_path.read_bytes()).hexdigest()
    loader = index.Validator(DISTANCE)
    stored, index_pin = loader.load(index_path, role='explicit-partial-index')
    validate_reviewed_sources(stored, closures, reviewed_sources)
    entries = stored.get('explicit_inputs')
    index.require(isinstance(entries, list) and 0 < len(entries) <= index.GROUPS, 'Missing bounded explicit receipt list')
    paths = []
    for entry in entries:
        receipt, pin = loader.dependency(entry, index_pin, 'explicit-reconstructed-input')
        paths.append(pin['path'])
    index.require(paths == sorted(set(paths)), 'Index input list is not canonical and duplicate-free')
    target = stored.get('individual_event_failure_bits')
    aggregate = stored.get('aggregate_failure_bits')
    precision = stored.get('aggregation_precision')
    index.require(type(aggregate) is int and aggregate >= FAILURE_BITS, 'Index aggregate target is weaker than40bits')
    rebuilt = index.build_index(paths, distance=DISTANCE, target_bits=target,
                                aggregate_bits=aggregate, precision=precision)
    # The only accepted extra CLI field is descriptive input_selection; it
    # does not replace the explicit, hashed input list used above.
    index.require(set(stored) <= set(rebuilt) | {'input_selection'}, 'Unexpected partial-index fields')
    index.require(all(key in stored and canonical(stored[key]) == canonical(value)
                      for key, value in rebuilt.items()), 'Partial index differs from reconstruction')
    validate_reviewed_sources(rebuilt, closures, reviewed_sources)
    bounds = validate_cover(rebuilt)
    replays = [replay_union(bounds, bits) for bits in PRECISIONS]
    ceilings = [F(row['log2_upper_rational']) for row in replays]
    index.require(all(value <= -FAILURE_BITS for value in ceilings), 'Complete union fails the40bit target')
    # Keep the weaker of the two valid ceilings; both precision runs must pass.
    ceiling = max(ceilings)
    loader.unchanged()
    for pin in rebuilt['receipt_file_pins']:
        index.require(sha256(Path(pin['path']).read_bytes()).hexdigest() == pin['sha256'],
                      'A transitive receipt changed during verification: ' + pin['path'])
    index.require(reviewed_closures() == (closures, reviewed_sources), 'Reviewed source tree changed during verification')
    index.require(sha256(own_path.read_bytes()).hexdigest() == own_pin, 'Manifest verifier changed during execution')
    return dict(schema='field128-literal35-complete-cover-manifest-v1',
                construction=index.CONSTRUCTION, field_order=index.FIELD_ORDER, field_bits=128,
                source_dimension=index.DIMENSION, output_length=index.OUTPUT_LENGTH, outer_groups=index.GROUPS,
                inner_width=64, state_dimension=16, inner_emission='four-RS[16,4]', feedback='A transpose',
                initial_state='zero', final_state='unrestricted', state_persists_across_regions=True,
                outer_emission='four-RS[64,32] per group',
                outer_randomizer='independent uniform invertible Toeplitz4 per outer group and output column',
                refresh='independent full-field nonzero diagonals D2 H16 D1 at every packet',
                routing='256 independent regional pair permutations and fair pair swaps; alternating cyclic pairings',
                fixed_maps='GF256 maps embedded in GF(2^128)', pair_swaps=True,
                full_field_random_coefficients=True, outer_measure='invertible-toeplitz',
                projective_divisor=str(index.FIELD_ORDER-1), target_distance=str(DISTANCE),
                forbidden_weight=int(DISTANCE*index.OUTPUT_LENGTH), distance_metric='symbol Hamming weight',
                probability_space='ideal independent setup; not a pseudorandom-seed guarantee',
                indexed_weight_domain=[1, index.GROUPS], complete_weight_count=index.GROUPS,
                disjoint_weight_cover=True, all_weights_covered=True, unresolved_weights=[],
                union_events=rebuilt['union_events'], union_event_count=len(bounds),
                terminal_union_event_count=rebuilt['terminal_union_event_count'],
                selected_union_log2_upper_rational=str(ceiling), failure_bits=FAILURE_BITS,
                directed_union_replays=replays, directed=True, whole_code_certificate=True,
                certificate_condition='validity of the pinned receipts from the listed reviewed directed producers',
                verification=dict(construction_and_formula_review='required reviewed-producer premise',
                                  producer_source_pins_checked=True, receipt_file_pins_checked=True,
                                  partial_index_reconstructed=True, scalar_enclosures_checked=True,
                                  disjoint_complete_cover_checked=True, union_replay_precisions=list(PRECISIONS),
                                  full_coefficient_replay=False, receipt_authentication=False),
                full_coefficient_replay=False, original_index_whole_code_certificate=False,
                input_index=index_pin, explicit_inputs=rebuilt['explicit_inputs'],
                receipt_file_pins=rebuilt['receipt_file_pins'], reviewed_producer_pins=REVIEWED,
                checked_producer_closures=rebuilt['checked_producer_closures'],
                checked_source_sha256=rebuilt['checked_source_sha256'],
                source_sha256=dict(closures[INDEX], **{own_path.name: own_pin}),
                source_pins_unchanged_at_exit=True, input_files_unchanged_at_exit=True,
                runtime=dict(python_flint=version('python-flint')),
                scope='complete literal construction at35percent, conditional on reviewed producer receipts; no claim for full-RS variants or timing')


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--index', type=Path, required=True, help='Explicit pinned partial-index JSON')
    parser.add_argument('--output', type=Path, required=True, help='New manifest path; exclusive create')
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output already exists')
    result = build_manifest(args.index)
    with args.output.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({key: result[key] for key in ('construction', 'target_distance', 'complete_weight_count',
                                                  'union_event_count', 'selected_union_log2_upper_rational',
                                                  'whole_code_certificate', 'full_coefficient_replay')}))


if __name__ == '__main__':
    main()
