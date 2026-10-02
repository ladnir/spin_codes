"""Source-checked complete occupancy union, with optional fresh replay.

Without --replay this only combines pinned numerical receipts and explicitly
does not mark a whole-code certificate. --replay recomputes every component
from its proposed tilts before summing. Claims concern the stated ideal
independent-randomness ensemble, never a deterministic seed or every setup.

The certified claim, when a fresh replay meets its target, is the following.
There are 512 outer groups, each mapping 128 message bits to 256 bits using
four parallel RS[16,8] words over GF16 (rs16), or RS[8,4] over GF256 (rs8).
Each symbol is independently mixed by uniform GL16 or GL32, respectively.
Each group independently permutes its 64 four-bit packets. In each of the
64 regions, an independent uniform permutation routes the 512 group packets.
For successive 64-bit inputs x, the fixed selected 64-by-16 expansion A emits
y=x+A*s and updates s=M*s+A^T*x, with independent uniform M in GL(16,2).
The state starts at zero, persists across all steps, and is discarded at end.
This defines a binary [131072,65536] code: the outer is injective, and the
inner map is lower triangular with identity diagonal. Over this ideal setup
law, the probability that its minimum distance is at most 13107 is bounded
by the complete first-moment union computed here. A successful 40-bit replay
therefore gives minimum distance at least 13108 except with probability
at most 2^-40. The union includes every nonzero-message occupancy q=1..512.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
from time import monotonic

from flint import arb, ctx
import packet_q1 as q1
import packet_q2 as q2
import packet_uniform_tail as tail


def dyadic(pair):
    if (not isinstance(pair, list) or len(pair) != 2 or
            any(type(v) is not int for v in pair) or pair[0] <= 0 or abs(pair[1]) > 1000000):
        raise ValueError('positive bounded-exponent dyadic endpoint required')
    return Q(pair[0]) * Q(2)**pair[1]


def counts_for(outer):
    if outer not in ('rs8', 'rs16'):
        raise ValueError('explicit RS8 or RS16 outer required')
    return tail.uniform_envelope(*((8, 4, 8) if outer == 'rs8' else (16, 8, 4)))


def count_hash(counts, cumulative=False):
    if cumulative:
        total, result = Q(0), []
        for count in counts:
            total += count
            result.append(total)
        counts = result
    return hashlib.sha256(json.dumps(list(map(str, counts))).encode()).hexdigest()


def load_receipt(path):
    path = Path(path).resolve()
    raw = path.read_bytes()
    record = json.loads(raw)
    hashes = record.get('source_sha256')
    if not isinstance(hashes, dict) or not hashes:
        raise ValueError('receipt must pin its mathematical sources')
    for filename, expected in hashes.items():
        if hashlib.sha256(Path(filename).read_bytes()).hexdigest() != expected:
            raise ValueError(f'receipt mathematical source changed: {filename}')
    if (record.get('K') != 65536 or record.get('N') != 131072 or record.get('threshold') != 13107
            or record.get('distance') != '1/10'
            or record.get('geometry') != dict(group_count=512, regions=64, group_dimension=128, macro_windows=32, packet_bits=4)
            or record.get('zero_initial_state') is not True or record.get('final_flush') is not False
            or record.get('precision', 0) < 128):
        raise ValueError('matching K16 half-rate continuous unflushed geometry required')
    return record, dict(path=str(path), sha256=hashlib.sha256(raw).hexdigest())


def check_coverage(ranges):
    covered = set()
    for first, last in ranges:
        if type(first) is not int or type(last) is not int or not 1 <= first <= last <= 512:
            raise ValueError('valid explicit occupancy intervals required')
        interval = set(range(first, last + 1))
        if interval & covered:
            raise ValueError('overlapping occupancy intervals')
        covered |= interval
    if covered != set(range(1, 513)):
        raise ValueError('incomplete occupancy union')


def combine(outer, q1_path, q2_path, tail_paths, *, precision=256,
            fresh_replay=False, output=None):
    if type(precision) is not int or precision < 128 or (output is not None and output.exists()):
        raise ValueError('precision >=128 and fresh output required')
    ctx.prec = precision
    beta, counts = counts_for(outer)
    one, one_ref = load_receipt(q1_path)
    two, two_ref = load_receipt(q2_path)
    if (one.get('schema') != 'finite-packet-q1-diagnostic-1' or one.get('occupancy_covered') != [1]
            or one.get('count_kind') != 'shells' or one.get('count_sha256') != count_hash(counts, True)):
        raise ValueError('q1 exact-shell receipt for the selected outer required')
    if (two.get('schema') != 'finite-packet-q2-diagnostic-1' or two.get('occupancy_covered') != [2]
            or two.get('all_two_group_support_pairs_covered') is not True
            or two.get('count_kind') != 'exact_expected_shells' or two.get('count_sha256') != count_hash(counts)):
        raise ValueError('q2 complete exact-shell receipt for the selected outer required')
    map_record = one.get('map_record')
    if not map_record or two.get('map_record') != map_record:
        raise ValueError('identical freshly enumerated physical-map records required')
    values = {1: dyadic(one['q1_upper']), 2: dyadic(two['q2_upper'])}
    references, conversions, ranges = [one_ref, two_ref], [], [(1, 1), (2, 2)]
    for path in tail_paths:
        record, reference = load_receipt(path)
        if (record.get('schema') != 'rs-uniform-regional-tail-diagnostic-1'
                or record.get('map_record') != map_record or record.get('every_shell_checked') is not True
                or record.get('state_continuity') != 'retained_between_every_physical_step_and_region'):
            raise ValueError('matching continuous-state uniform-envelope tail receipt required')
        source_beta, source_counts = counts_for(record.get('outer'))
        if Q(record.get('beta')) != source_beta or record.get('count_sha256') != count_hash(source_counts):
            raise ValueError('freshly checked source outer and exact beta required')
        first, last = record['q_min'], record['q_max']
        if not 3 <= first <= last <= 512 or record.get('occupancy_covered') != [first, last]:
            raise ValueError('valid tail occupancy interval required')
        if set(record['occupancy_uppers']) != {str(q) for q in range(first, last + 1)}:
            raise ValueError('every claimed tail occupancy must have an endpoint')
        ratio = beta / source_beta
        for q in range(first, last + 1):
            if q in values:
                raise ValueError('duplicate occupancy endpoint')
            values[q] = dyadic(record['occupancy_uppers'][str(q)]) * ratio**q
        ranges.append((first, last))
        references.append(reference)
        conversions.append(dict(source_outer=record['outer'], target_outer=outer,
            beta_ratio=str(ratio), operation='multiply each q endpoint by beta_ratio**q before summing',
            occupancy=[first, last]))
    check_coverage(ranges)
    one_upper, two_upper = q1.kernel_t64.aq(values[1]), q1.kernel_t64.aq(values[2])
    tail_upper = q1.kernel_t64.up(sum((q1.kernel_t64.aq(values[q]) for q in range(3, 513)), arb(0)))
    upper = q1.kernel_t64.up(one_upper + two_upper + tail_upper)
    meets = bool(upper < arb(2)**-40)
    record = dict(schema='k16-rs-complete-occupancy-union-1', outer=outer,
        K=65536, N=131072, threshold=13107, distance='1/10', target_margin_bits=40,
        occupancy_covered=[1, 512], all_occupancies_covered=True,
        fresh_replay=fresh_replay, whole_code_certificate=bool(fresh_replay and meets),
        target_met=meets, zero_initial_state=True, final_flush=False,
        geometry=one['geometry'], map_record=map_record, precision=precision,
        input_receipts=references, tail_conversions=conversions,
        q1_upper=q1.endpoint(one_upper), q2_upper=q1.endpoint(two_upper),
        q3_through_512_upper=q1.endpoint(tail_upper), union_upper=q1.endpoint(upper),
        margin_bits=str(-upper.log() / arb(2).log()),
        tail_margin_bits=str(-tail_upper.log() / arb(2).log()),
        source_sha256=q1.source_snapshot(),
        scope='Probability over independent ideal uniform outer symbol GL maps, '
              'uniform per-group column shuffles and independent regional shuffles, '
              'and ideal independent uniform GL16 physical updates, with the selected '
              'fixed maps. Full first-moment union over all nonzero messages. '
              'Not a guarantee for every setup, deterministic seed, or pseudorandom generator. '
              'Without fresh_replay this is a source-checked receipt combination, not a certificate replay.')
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(record, indent=2) + '\n')
    print(f'COMPLETE UNION {outer} margin_bits={record["margin_bits"]} '
          f'tail_margin_bits={record["tail_margin_bits"]} fresh_replay={fresh_replay} '
          f'whole_code_certificate={record["whole_code_certificate"]}', flush=True)
    return record


def replay(outer, q1_path, q2_path, tail_paths, *, precision=256, output):
    """Receipts supply tilt proposals only; no saved numerical bound is reused."""
    start = monotonic()
    old_one, _ = load_receipt(q1_path)
    old_two, _ = load_receipt(q2_path)
    old_tails = [load_receipt(path)[0] for path in tail_paths]
    check_coverage([(1, 1), (2, 2), *[(r['q_min'], r['q_max']) for r in old_tails]])
    paths = [output.with_name(output.stem + suffix + '.json') for suffix in
             ('-q1', '-q2', *(f'-tail-{i}' for i in range(len(old_tails))))]
    if output.exists() or any(path.exists() for path in paths):
        raise ValueError('fresh replay output and component paths required')
    _, counts = counts_for(outer)
    premises = dict(outer='four parallel GF256 RS[8,4]' if outer == 'rs8' else 'four parallel GF16 RS[16,8]',
        label_mixing='independent uniform GL32 per symbol' if outer == 'rs8' else 'independent uniform GL16 per symbol',
        independent_setups_between_groups=True, exact_expected_shell_counts=True)
    # Keep each selected witness spelling, without replaying unused tilts.
    one_used = set(old_one['support_choices']) - {None}
    two_used = {value for column in old_two['support_pair_choices'] for value in column} - {None}
    one_tilts = [tilt for tilt in old_one['tilts'] if tilt in one_used]
    two_tilts = [tilt for tilt in old_two['tilts'] if tilt in two_used]
    if not one_tilts or not two_tilts:
        raise ValueError('nonempty selected q1 and q2 tilt witnesses required')
    q1.evaluate_q1(counts, group_count=512, regions=64, epochs_per_region=16,
        group_dimension=128, count_kind='shells', tilts=one_tilts, precision=precision,
        metadata=premises, output=paths[0])
    q2.evaluate_q2(counts, geometry=q1.Geometry(512, 64, 128), tilts=two_tilts,
        precision=precision, metadata=premises, output=paths[1])
    for path, recipe in zip(paths[2:], old_tails):
        used = set(recipe['occupancy_choices'].values()) - {None}
        selected = [tilt for tilt in recipe['tilts'] if tilt in used]
        if not selected:
            raise ValueError('nonempty selected tail tilt witnesses required')
        tail.run(outer=outer, q_min=recipe['q_min'], q_max=recipe['q_max'],
            precision=precision, tilts=selected, output=path)
    result = combine(outer, paths[0], paths[1], paths[2:], precision=precision,
                     fresh_replay=True, output=output)
    print(f'Fresh replay elapsed={monotonic() - start:.2f}s', flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outer', choices=('rs8', 'rs16'), required=True)
    parser.add_argument('--q1', type=Path, required=True)
    parser.add_argument('--q2', type=Path, required=True)
    parser.add_argument('--tail', type=Path, nargs='+', required=True)
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--replay', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    action = replay if args.replay else combine
    action(args.outer, args.q1, args.q2, args.tail, precision=args.precision, output=args.output)


if __name__ == '__main__':
    main()
