"""Assemble fresh larger-outer q1/q2 with the authenticated derived tail.

This is a derived certificate, not a fresh replay of all occupancies. The
tail is regenerated from the pinned K18 receipt by the coset transfer lemma.
No old exact-shell q1/q2 endpoint is accepted.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import json
from pathlib import Path

from flint import arb, ctx
import packet_outer_block_transfer as transfer

SCHEMA = 'rs-larger-outer-derived-whole-1'
REQUIRED_LOCAL_SOURCES = (
    'packet_larger_outer_sparse.py', 'packet_inner_quadratic_extension.py',
    'packet_inner_small_extension.py', 'packet_inner_s20.py', 'packet_q1.py',
    'packet_regional_power.py', 'packet_uniform_tail.py', 'rs_uniform_envelope.py',
    'rs_outer.py',
)


def authenticate_sparse(path, *, multiplier, source_whole):
    """Validate the fresh partial producer, source pins, geometry and map."""
    transfer.beta_ratio(multiplier)
    path = Path(path).resolve()
    raw = path.read_bytes()
    record = json.loads(raw)
    expected = dict(schema='larger-rs-outer-uniform-sparse-1', multiplier=multiplier,
        K=transfer.BASE_K*multiplier, N=2*transfer.BASE_K*multiplier,
        groups=2048, regions=64*multiplier, outer_n=8*multiplier, outer_k=4*multiplier,
        parallel_rows=4, base_field_size=256, symbol_mixer='independent uniform GL32',
        outer_dimension=128*multiplier, outer_length=256*multiplier,
        beta=str(transfer.BASE_BETA**multiplier*transfer.beta_ratio(multiplier)),
        distance='1/10', cutoff=2*transfer.BASE_K*multiplier//10,
        occupancy_covered=[1, 2], whole_code_certificate=False,
        fresh_computation=True, physical_t=64, state_bits=22,
        zero_initial_state=True, final_flush=False, continuous_state=True)
    if any(record.get(key) != value for key, value in expected.items()):
        raise ValueError('fresh sparse scope, outer geometry or beta mismatch')
    if type(record.get('precision')) is not int or record['precision'] < 192:
        raise ValueError('fresh sparse outward precision must be at least 192 bits')
    if record.get('map_record') != source_whole['map_record']:
        raise ValueError('fresh sparse must use the exact retained S22 map and census scope')
    sources = record.get('source_sha256')
    if not isinstance(sources, dict) or not sources:
        raise ValueError('fresh sparse mathematical source pins are required')
    here = Path(__file__).resolve().parent
    required = {str((here / name).resolve()) for name in REQUIRED_LOCAL_SOURCES}
    if not required <= sources.keys():
        raise ValueError('required sparse producer/core source pins are missing')
    for filename, digest in sources.items():
        if transfer._digest(Path(filename).read_bytes()) != digest:
            raise ValueError(f'fresh sparse source changed: {filename}')
    endpoints, choices = record.get('occupancy_uppers'), record.get('occupancy_choices')
    if (not isinstance(endpoints, dict) or set(endpoints) != {'1', '2'} or
            not isinstance(choices, dict) or set(choices) != {'1', '2'}):
        raise ValueError('exactly fresh q1 and q2 endpoints and witnesses are required')
    tilts = record.get('tilts')
    if not isinstance(tilts, list) or not tilts or any(Q(t) <= 0 for t in tilts):
        raise ValueError('positive sparse rational trial tilts required')
    for q in ('1', '2'):
        transfer._dyadic(endpoints[q])
        if choices[q] not in tilts or Q(choices[q]) <= 0:
            raise ValueError('selected sparse tilt was not a computed trial')
    return record, dict(path=str(path), sha256=transfer._digest(raw),
                       canonical_sha256=transfer._digest(transfer._canonical(record)))


def assemble(source, sparse, *, multiplier, precision=256, output=None):
    """Cover all q with an exact dyadic sum and one final upward rounding."""
    if type(precision) is not int or precision < 256:
        raise ValueError('whole-certificate assembly requires at least 256-bit arithmetic')
    if output is not None and Path(output).exists():
        raise ValueError('output must be fresh')
    old, _, _ = transfer.authenticate(source)
    sparse_record, sparse_ref = authenticate_sparse(sparse, multiplier=multiplier, source_whole=old)
    tail = transfer.run(source, multiplier=multiplier, precision=precision)
    values = {int(q): transfer._dyadic(v) for q, v in tail['occupancy_uppers'].items()}
    for q, endpoint in sparse_record['occupancy_uppers'].items():
        if int(q) in values:
            raise ValueError('sparse and tail occupancies overlap')
        values[int(q)] = transfer._dyadic(endpoint)
    if set(values) != set(range(1, 2049)):
        raise ValueError('every occupancy from 1 through 2048 must be covered exactly once')
    previous = ctx.prec
    try:
        ctx.prec = precision
        # This addition is exact; only its final conversion is rounded upward.
        mantissa, exponent = transfer._sum_dyadics(values.values())
        upper = transfer._up(arb(mantissa)*arb(2)**exponent)
        passes = transfer._fraction(upper) < Q(1, 1 << 40)
        result = dict(schema=SCHEMA, K=tail['K'], N=tail['N'], groups=2048,
            regions=tail['regions'], outer=tail['outer'],
            group_dimension=tail['group_dimension'], group_output_bits=tail['group_output_bits'],
            physical_t=64, state_bits=22, map_record=old['map_record'],
            inner_distribution='uniform_gl', independent_physical_updates=True,
            independent_outer_symbol_maps=True, independent_regional_shuffles=True,
            zero_initial_state=True, final_flush=False,
            state_continuity='retained_between_every_physical_step_and_region',
            distance='1/10', threshold=tail['threshold'], target_minimum_distance=tail['threshold']+1,
            target_margin_bits=40, precision=precision, occupancy_covered=[1, 2048],
            all_occupancies_covered=True,
            occupancy_uppers={str(q): list(values[q]) for q in sorted(values)},
            union_upper=upper, margin_bits=transfer._margin(upper),
            q1_upper=list(values[1]), q2_upper=list(values[2]),
            tail_upper=tail['tail_upper'], tail_margin_bits=tail['tail_margin_bits'],
            sparse_source=sparse_ref,
            tail_source=dict(canonical_sha256=transfer._digest(transfer._canonical(tail)),
                source=tail['source'], helper_sha256=tail['helper_sha256'],
                proof_note_sha256=tail['proof_note_sha256']),
            source_sha256=dict(sparse_record['source_sha256']),
            assembly_source_sha256=transfer._digest(Path(__file__).read_bytes()),
            arithmetic='Exact sum of disjoint dyadic endpoints; one final upward Arb rounding.',
            target_met=passes, whole_code_certificate=passes,
            derived_certificate=passes, fresh_replay=False,
            scope='Full first-moment bound for the specified ideal larger-outer ensemble. '
                  'Fresh direct q1/q2 plus the proved transfer of authenticated K18 tail bounds; '
                  'not an all-occupancy fresh replay, implementation, or seeded guarantee.')
        # Reject a source change during assembly, including the fresh sparse source.
        transfer.authenticate(source)
        authenticate_sparse(sparse, multiplier=multiplier, source_whole=old)
    finally:
        ctx.prec = previous
    if output is not None:
        with Path(output).open('x') as handle:
            json.dump(result, handle, indent=2)
            handle.write('\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source')
    parser.add_argument('sparse')
    parser.add_argument('--multiplier', type=int, choices=(2, 4), required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    result = assemble(args.source, args.sparse, multiplier=args.multiplier, output=args.output)
    print(json.dumps({k: result[k] for k in ('K', 'margin_bits', 'derived_certificate')}))
