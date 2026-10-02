"""Fresh complete-occupancy replay for the K=2^20 RS16 packet ensemble.

The ensemble has 8192 groups, each encoding 128 bits into 256 bits with four
parallel GF16 RS[16,8] words. Independent uniform GL16 maps mix the sixteen
aligned symbols of each group. Independent uniform group-column shuffles
and regional shuffles route the resulting four-bit packets in 64 regions.
The selected t64/s16 maps and independent uniform GL16 physical updates
define the continuous-state inner, initially zero and without a final flush.

This binary [2097152,1048576] encoder is injective. Its expected number of
nonzero messages producing weight at most 209715 bounds the probability,
over setup, that its minimum distance is at most 209715. A complete sum
below 2^-40 therefore proves minimum distance at least 209716 with the
requested setup-failure margin. Every occupancy q=1..8192 is included.

Saved component receipts may be checked and summed, but that operation is
not a certificate replay. Only replay() freshly enumerates the fixed maps,
computes every component, and marks a passing whole-code certificate. The
existing K16 proof sources and receipts are not modified or relabeled.

Current status: the K20 10%/40-bit target has not closed. The default witness
grid is exploratory, not a passing recipe. A failed replay retains its failed
result and never sets whole_code_certificate=True.
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
import packet_q1 as q1
import packet_q2 as q2
import packet_uniform_tail as uniform
import packet_rs_k20_sparse as sparse
import packet_rs_k20_dense as dense
# Preparation otherwise imports this lazily. Importing the module performs no
# census and makes the full capped-density dependency visible to source pins.
import kernel_birth_density


GEOMETRY = sparse.GEOMETRY
DEFAULT_Q1 = ('.00016', '.00024', '.00032')
DEFAULT_Q2 = ('.00032', '.00048', '.00064')
DEFAULT_SPARSE = ('.0008', '.0016', '.0032', '.0064', '.0128')
DEFAULT_DENSE = ('.0032', '.0064', '.0128', '.0256', '.0512', '.1024',
                 '.2048', '.4096', '.8192', '1.6384', '2.1972246')
CONTINUITY = 'retained_between_every_physical_step_and_region'

# Conservative audited import graph of the shared proof engine. This includes
# imported proposal helpers as well as the active outward arithmetic: either
# can affect which complete valid envelope is selected. Keep this explicit so
# callers, tests, and unrelated modules cannot change the required receipt pins.
# Paths are relative to q1.LOCALITY; the fixed-map JSON is required separately.
CORE_FILENAMES = tuple('''
basis_lattice.py
bch_joint_support.py
cancellation_joint.py
feedback_character_census.py
feedback_convolution.py
fresh_collision.py
full_feedback_census.py
full_feedback_refinement.py
gf16_packets/birth_classes.py
gf16_packets/birth_classes_probe.py
gf16_packets/birth_refresh.py
gf16_packets/birth_refresh_probe.py
gf16_packets/conditioned_return.py
gf16_packets/gf_refresh.py
gf16_packets/occupancy_birth_classes.py
gf16_packets/occupancy_kernel.py
gf16_packets/occupancy_rank.py
gf16_packets/profile_return.py
gf16_packets/rank_return.py
gf16_packets/refresh_kernel.py
gf16_packets/scalar_cover.py
gf16_packets/shape_return.py
gf16_packets/single_group.py
gf16_packets/sparse_cover.py
gf16_packets/trimmed_return.py
gf16_packets/weighted_return.py
group_moment.py
group_rank_one_verify.py
independent_rows/candidates/density_extend.py
independent_rows/candidates/exact_feedback.py
independent_rows/candidates/mass_density.py
independent_rows/candidates/mass_density_screen.py
independent_rows/candidates/overlap_mean/mean.py
independent_rows/candidates/overlap_second/density_second.py
independent_rows/candidates/overlap_second/quadratic.py
independent_rows/candidates/overlap_second/second.py
independent_rows/column_density.py
independent_rows/dense_closure/kernel.py
independent_rows/dense_closure/mixture.py
independent_rows/feedback_density.py
independent_rows/support.py
independent_rows/universal_density.py
independent_rows/verify.py
joint_oa_caps.py
joint_support.py
mature_tail.py
mixing_rounds.py
occupancy_adaptive.py
occupancy_allones.py
occupancy_cdf_cover.py
occupancy_count_gap.py
occupancy_fresh_moment.py
occupancy_memory.py
occupancy_memory_verify.py
occupancy_model.py
occupancy_multi_average.py
occupancy_screen.py
occupancy_sensitivity.py
occupancy_weighted.py
occupancy_window_average.py
packed_mixing/s16_closure/kernel_birth_density.py
packed_mixing/s16_closure/kernel_maps.py
packed_mixing/s16_closure/kernel_t64.py
packed_mixing/s16_closure/sparse_kernel.py
packed_mixing/s16_maps.py
pair_tail.py
random_group_verify.py
random_rank_one.py
rank_two_moment.py
rank_two_types.py
shortened_bound.py
two_bit/heterogeneous.py
two_bit/iid_kernel.py
two_bit/mixture_cover.py
two_bit/model.py
two_bit/probe.py
two_bit/row_mixture.py
two_column_moment.py
two_group_moment.py
two_group_screen.py
two_group_verify.py
window_feedback.py
window_histogram.py
zero_moment.py
'''.split())


def required_source_paths(schema):
    """Return the complete audited core plus the schema's actual producer."""
    producers = {
        'finite-packet-q1-diagnostic-1': ('packet_q1.py',),
        'finite-packet-q2-diagnostic-1': ('packet_q2.py',),
        sparse.TAIL_SCHEMA: ('packet_rs_k20_sparse.py', 'packet_uniform_tail.py'),
        dense.SCHEMA: ('packet_rs_k20_dense.py',),
        dense.FUGACITY_SCHEMA: ('packet_rs_k20_dense.py',),
    }
    if schema not in producers:
        raise ValueError('recognized K20 component schema required')
    here = Path(__file__).resolve().parent
    paths = {q1.LOCALITY / name for name in CORE_FILENAMES}
    paths.add(q1.kernel_t64.SELECTED_MAP)
    paths.update(here / name for name in
                 ('packet_q1.py', 'rs_outer.py', 'rs_uniform_envelope.py', *producers[schema]))
    return frozenset(str(path.resolve()) for path in paths)


def endpoint_value(pair):
    """Read a positive finite dyadic directly into Arb, without giant fractions."""
    if (not isinstance(pair, list) or len(pair) != 2 or
            any(type(v) is not int for v in pair) or pair[0] <= 0 or
            abs(pair[1]) > 10_000_000):
        raise ValueError('positive bounded-exponent dyadic endpoint required')
    value = arb(pair[0]) * arb(2) ** pair[1]
    if not value.is_finite() or not value > 0:
        raise ValueError('positive finite endpoint required')
    return value


def check_coverage(ranges):
    covered = set()
    for first, last in ranges:
        if (type(first) is not int or type(last) is not int or
                not 1 <= first <= last <= GEOMETRY.group_count):
            raise ValueError('explicit valid occupancy intervals required')
        interval = set(range(first, last + 1))
        if interval & covered:
            raise ValueError('overlapping occupancy intervals')
        covered.update(interval)
    if covered != set(range(1, GEOMETRY.group_count + 1)):
        raise ValueError('incomplete occupancy union')


def load_receipt(path):
    path = Path(path).resolve()
    raw = path.read_bytes()
    record = json.loads(raw)
    hashes = record.get('source_sha256')
    if not isinstance(hashes, dict) or not hashes:
        raise ValueError('receipt must pin its mathematical sources')
    missing = required_source_paths(record.get('schema')) - hashes.keys()
    if missing:
        raise ValueError('receipt omits mathematical source pins: ' + ', '.join(sorted(missing)))
    for filename, expected in hashes.items():
        if hashlib.sha256(Path(filename).read_bytes()).hexdigest() != expected:
            raise ValueError(f'mathematical source changed: {filename}')
    if (record.get('geometry') != asdict(GEOMETRY) or
            record.get('K') != GEOMETRY.K or record.get('N') != GEOMETRY.N or
            record.get('threshold') != sparse.THRESHOLD or record.get('distance') != '1/10' or
            record.get('zero_initial_state') is not True or record.get('final_flush') is not False or
            type(record.get('precision')) is not int or record['precision'] < 128 or
            record.get('whole_code_certificate') is not False or record.get('outer') != 'rs16'):
        raise ValueError('matching K20 half-rate continuous-state RS16 component required')
    return record, dict(path=str(path), sha256=hashlib.sha256(raw).hexdigest())


def combine(one_path, two_path, tail_paths, *, precision=256):
    """Check pinned components and sum their endpoints; never mark a replay."""
    if type(precision) is not int or precision < 128:
        raise ValueError('integer precision >=128 required')
    ctx.prec = precision
    beta, counts = sparse.exact_outer()
    one, one_ref = load_receipt(one_path)
    two, two_ref = load_receipt(two_path)
    if (one.get('schema') != 'finite-packet-q1-diagnostic-1' or
            one.get('occupancy_covered') != [1] or one.get('count_kind') != 'shells' or
            one.get('count_sha256') != sparse.count_hash(counts, cumulative=True)):
        raise ValueError('complete q1 exact-shell component required')
    if (two.get('schema') != 'finite-packet-q2-diagnostic-1' or
            two.get('occupancy_covered') != [2] or
            two.get('all_two_group_support_pairs_covered') is not True or
            two.get('count_kind') != 'exact_expected_shells' or
            two.get('count_sha256') != sparse.count_hash(counts)):
        raise ValueError('complete q2 exact-shell component required')
    map_record = one.get('map_record')
    if not map_record or two.get('map_record') != map_record:
        raise ValueError('identical enumerated physical maps required')
    values = {1: endpoint_value(one['q1_upper']), 2: endpoint_value(two['q2_upper'])}
    references, ranges = [one_ref, two_ref], [(1, 1), (2, 2)]
    component_precisions = [one['precision'], two['precision']]
    summaries = []
    for path in tail_paths:
        record, reference = load_receipt(path)
        if (record.get('schema') not in (sparse.TAIL_SCHEMA, dense.SCHEMA, dense.FUGACITY_SCHEMA) or
                record.get('map_record') != map_record or record.get('every_shell_checked') is not True or
                record.get('state_continuity') != CONTINUITY or
                record.get('count_sha256') != sparse.count_hash(counts) or Q(record.get('beta')) != beta):
            raise ValueError('matching continuous-state RS16 uniform-envelope component required')
        first, last = record['q_min'], record['q_max']
        if (type(first) is not int or type(last) is not int or
                not 3 <= first <= last <= GEOMETRY.group_count or
                record.get('occupancy_covered') != [first, last]):
            raise ValueError('valid explicit tail occupancy interval required')
        if set(record['occupancy_uppers']) != {str(q) for q in range(first, last + 1)}:
            raise ValueError('every claimed occupancy needs an endpoint')
        subtotal = arb(0)
        for q in range(first, last + 1):
            if q in values:
                raise ValueError('duplicate occupancy endpoint')
            values[q] = endpoint_value(record['occupancy_uppers'][str(q)])
            subtotal += values[q]
        ranges.append((first, last))
        references.append(reference)
        component_precisions.append(record['precision'])
        subtotal = q1.kernel_t64.up(subtotal)
        summaries.append(dict(occupancy=[first, last], schema=record['schema'],
            upper=q1.endpoint(subtotal), margin_bits=str(-subtotal.log() / arb(2).log())))
    check_coverage(ranges)
    tail_upper = q1.kernel_t64.up(sum((values[q] for q in range(3, GEOMETRY.group_count + 1)), arb(0)))
    upper = q1.kernel_t64.up(values[1] + values[2] + tail_upper)
    meets = bool(upper < arb(2) ** -40)
    return dict(schema='k20-rs16-complete-occupancy-union-1', outer='rs16',
        K=GEOMETRY.K, N=GEOMETRY.N, geometry=asdict(GEOMETRY),
        threshold=sparse.THRESHOLD, distance='1/10', target_margin_bits=40,
        guaranteed_minimum_distance=sparse.THRESHOLD + 1,
        occupancy_covered=[1, GEOMETRY.group_count], all_occupancies_covered=True,
        fresh_replay=False, whole_code_certificate=False, target_met=meets,
        zero_initial_state=True, final_flush=False, state_continuity=CONTINUITY,
        map_record=map_record, precision=precision, minimum_component_precision=min(component_precisions),
        input_receipts=references, tail_components=summaries,
        q1_upper=q1.endpoint(values[1]), q2_upper=q1.endpoint(values[2]),
        q3_through_8192_upper=q1.endpoint(tail_upper), union_upper=q1.endpoint(upper),
        margin_bits=str(-upper.log() / arb(2).log()),
        tail_margin_bits=str(-tail_upper.log() / arb(2).log()),
        source_sha256=q1.source_snapshot(),
        scope='Probability over independent uniform GL16 outer symbol maps, '
              'uniform per-group packet permutations and uniform regional permutations, '
              'and independent uniform GL16 physical state updates with the selected fixed maps. '
              'Full first-moment union over all nonzero messages. A saved-component sum '
              'is not marked as a fresh whole-code certificate.')


def recipe(*, sparse_max=32, q1_tilts=DEFAULT_Q1, q2_tilts=DEFAULT_Q2,
           sparse_tilts=DEFAULT_SPARSE, dense_tilts=DEFAULT_DENSE,
           marker_probabilities=dense.DEFAULT_MARKER_PROBABILITIES, precision=256):
    if type(sparse_max) is not int or not 3 <= sparse_max < GEOMETRY.group_count:
        raise ValueError('sparse endpoint must lie in 3..8191')
    if type(precision) is not int or precision < 256:
        raise ValueError('fresh whole replay requires at least 256-bit precision')
    tilts = {}
    for name, values in [('q1', q1_tilts), ('q2', q2_tilts),
                         ('sparse', sparse_tilts), ('dense', dense_tilts)]:
        tilts[name] = list(sparse.checked_options(values, precision, None)[0])
    if isinstance(marker_probabilities, (str, bytes)):
        raise ValueError('a sequence of marker probabilities is required')
    markers = tuple(map(Q, marker_probabilities))
    if (not markers or len(set(markers)) != len(markers) or
            any(not 0 < p <= 1 for p in markers) or not any(p < 1 for p in markers)):
        raise ValueError('distinct valid markers covering the dense interval required')
    return dict(K=GEOMETRY.K, N=GEOMETRY.N, geometry=asdict(GEOMETRY),
        threshold=sparse.THRESHOLD, precision=precision, sparse_max=sparse_max,
        intervals=[[1, 1], [2, 2], [3, sparse_max], [sparse_max + 1, GEOMETRY.group_count]],
        tilts=tilts, marker_probabilities=list(map(str, markers)))


def replay(output, **options):
    """Recompute every component; the only inputs are geometry and tilt witnesses."""
    start = monotonic()
    plan = recipe(**options)
    output = Path(output)
    paths = [output.with_name(output.stem + '-' + name + '.json')
             for name in ('q1', 'q2', 'sparse', 'dense')]
    if output.exists() or any(path.exists() for path in paths):
        raise ValueError('whole replay and all component paths must be fresh')
    precision, split, tilts = plan['precision'], plan['sparse_max'], plan['tilts']
    ctx.prec = precision
    data, map_record = q1.kernel_t64.prepare(birth_density='capped')
    sources = q1.source_snapshot()
    shared = dict(precision=precision, data=data, map_record=map_record)
    sparse.run_q1(tilts=tilts['q1'], output=paths[0], **shared)
    sparse.run_q2(tilts=tilts['q2'], output=paths[1], **shared)
    sparse.run_tail(q_min=3, q_max=split, tilts=tilts['sparse'], output=paths[2], **shared)
    dense.run_fugacity(q_min=split + 1, q_max=GEOMETRY.group_count,
              tilts=tilts['dense'], marker_probabilities=plan['marker_probabilities'],
              output=paths[3], **shared)
    result = combine(paths[0], paths[1], paths[2:], precision=precision)
    if sources != q1.source_snapshot() or result['map_record'] != map_record:
        raise RuntimeError('mathematical source or map changed during whole replay')
    result.update(fresh_replay=True, whole_code_certificate=result['target_met'],
                  recipe=plan, elapsed_seconds=monotonic() - start)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('x', encoding='utf-8') as stream:
        stream.write(json.dumps(result, indent=2) + '\n')
    print(f'K20 RS16 COMPLETE UNION margin_bits={result["margin_bits"]} '
          f'tail_margin_bits={result["tail_margin_bits"]} '
          f'whole_code_certificate={result["whole_code_certificate"]} '
          f'elapsed={result["elapsed_seconds"]:.2f}s', flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--sparse-max', type=int, default=32)
    parser.add_argument('--q1-tilts', nargs='+', default=list(DEFAULT_Q1))
    parser.add_argument('--q2-tilts', nargs='+', default=list(DEFAULT_Q2))
    parser.add_argument('--sparse-tilts', nargs='+', default=list(DEFAULT_SPARSE))
    parser.add_argument('--dense-tilts', nargs='+', default=list(DEFAULT_DENSE))
    parser.add_argument('--marker-probabilities', nargs='+', default=list(dense.DEFAULT_MARKER_PROBABILITIES))
    parser.add_argument('--dry-run', action='store_true')
    args = vars(parser.parse_args())
    output, dry = args.pop('output'), args.pop('dry_run')
    if dry:
        print(json.dumps(recipe(**args), indent=2))
    else:
        replay(output, **args)


if __name__ == '__main__':
    main()
