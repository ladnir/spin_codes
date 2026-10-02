"""Reproduce the complete K20, rate-one-half, t64/s20 RS16 packet bound.

From the repository root, run:

    python -B research/workstreams/k16_design/reproduce_rs16_k20.py \
        --output tmp/rs16-k20-reproduction/whole-p256.json

Python, python-flint, NumPy and SciPy must already be installed. No saved
numerical receipt, native encoder, remote host or benchmark is required.
The fresh 256-bit replay takes about eight minutes on the development CPU.
It counts every nonzero occupancy q=1..4096 at distance cutoff 209715.

The outer has eight parallel GF16 RS[16,8] rows and independent symbol
randomizers with uniform nonzero fixed-input images; adjoints of uniform
nonzero GF(2^32) multipliers are one supported family. Packet routing is
independently uniform within groups and regions. The recorded t64/s20 inner
starts at zero and keeps its state across every step and region, without a
final flush. The claim concerns the ideal setup ensemble, not every setup.

Only rational tilt and marker choices are encoded here, not saved endpoints.
The frozen whole replay freshly prepares the fixed inner maps, computes the
exact-placement prefix and conditioned-iid suffix, and sums their disjoint
outward endpoints exactly. Its source snapshot also pins this entry point
when it is run as a script. Thus source hashes differ from earlier invocations
while the numerical union endpoint is reproducible.

Use --dry-run to inspect the plan without importing numerical dependencies,
creating directories or files, or performing a bound calculation. All generated
filenames must be fresh; existing files are never overwritten.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import json
from pathlib import Path
from time import monotonic


PRECISION = 256
K, N, GROUPS, REGIONS, STATE_BITS = 1 << 20, 1 << 21, 4096, 128, 20
PREFIX_END = 442
PREFIX_RECIPES = (
    (1, 1, '1/2500'), (2, 4, '3/2500'), (5, 10, '2/625'),
    (11, 21, '4/625'), (22, 40, '7/500'), (41, 63, '3/125'),
    (64, 84, '7/200'), (85, 100, '17/400'), (101, 116, '1/20'),
    (117, 132, '23/400'), (133, 148, '13/200'), (149, 169, '29/400'),
    (170, 197, '17/200'), (198, 238, '1/10'), (239, 298, '1/8'),
    (299, 442, '4/25'),
)
DENSE_TILTS = ('1/64', '1/32', '3/64', '1/16', '5/64', '3/32',
    '1/8', '3/16', '1/4', '3/8', '1/2', '3/4', '1', '3/2',
    '10986123/5000000')
MARKERS = tuple(str(Q(j, 64)) for j in range(1, 65))


def prefix_choices():
    """Return parameter-selection hints without computed support bounds."""
    result = {}
    for first, last, tilt in PREFIX_RECIPES:
        if type(first) is not int or type(last) is not int or first > last or Q(tilt) <= 0:
            raise ValueError('valid positive rational prefix recipes required')
        for q in range(first, last+1):
            if str(q) in result:
                raise ValueError('prefix occupancy recipes overlap')
            result[str(q)] = dict(tilt=str(Q(tilt)))
    if set(result) != {str(q) for q in range(1, PREFIX_END+1)}:
        raise ValueError('prefix recipes must cover every occupancy exactly once')
    return result


def output_paths(output):
    output = Path(output).resolve()
    if output.suffix.lower() != '.json':
        raise ValueError('the output filename must end in .json')
    proposal = output.with_name(output.stem+'-proposal-prefix.json')
    components = [Path(str(output)+'.exact.json'), Path(str(output)+'.dense.json')]
    return output, proposal, components


def checked_plan(output):
    output, proposal, components = output_paths(output)
    prefix_choices()
    if any(Q(tilt) <= 0 for tilt in DENSE_TILTS) or len(set(DENSE_TILTS)) != len(DENSE_TILTS):
        raise ValueError('distinct positive dense tilts required')
    if tuple(map(Q, MARKERS)) != tuple(Q(j, 64) for j in range(1, 65)):
        raise ValueError('the fixed marker grid must be j/64 for j=1..64')
    if not 0 < PREFIX_END < GROUPS:
        raise ValueError('proper exact prefix and disjoint nonempty dense suffix required')
    paths = [output, proposal, *components]
    if len(set(paths)) != len(paths):
        raise ValueError('all generated output paths must be distinct')
    existing = [str(path) for path in paths if path.exists()]
    if existing:
        raise FileExistsError('fresh output paths required: '+', '.join(existing))
    return dict(K=K, N=N, groups=GROUPS, regions=REGIONS, state_bits=STATE_BITS,
        physical_t=64, outer_dimension=256, precision=PRECISION,
        cutoff=N//10, minimum_distance=N//10+1, distance='1/10',
        target_margin_bits=40, occupancy_covered=[1, GROUPS],
        exact_prefix=[1, PREFIX_END], dense_suffix=[PREFIX_END+1, GROUPS],
        prefix_recipes=[dict(q_min=a, q_max=b, tilt=t) for a, b, t in PREFIX_RECIPES],
        dense_tilts=list(DENSE_TILTS), marker_probabilities=list(MARKERS),
        whole_code_certificate=False, dry_run_is_not_a_certificate=True,
        numerical_endpoints_present=False, output=str(output),
        proposal_path=str(proposal), fresh_component_paths=list(map(str, components)))


def proposal_record(envelope_metadata):
    return dict(schema='rs16-k20-standalone-prefix-proposal-1', proposal_only=True,
        whole_code_certificate=False, numerical_endpoints_present=False,
        K=K, state_bits=STATE_BITS, envelope=envelope_metadata, best=prefix_choices(),
        scope='Rational tilt choices only; the replay recomputes every numerical endpoint.')


def reproduce(output):
    plan = checked_plan(output)
    # Import only after the read-only preflight, keeping --dry-run dependency-free.
    import packet_outer_cost_whole as whole

    start = monotonic()
    reference, _ = whole.outer_options(256)
    recipe = proposal_record(reference.metadata())
    output, proposal = Path(plan['output']), Path(plan['proposal_path'])
    output.parent.mkdir(parents=True, exist_ok=True)
    with proposal.open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(recipe, handle, indent=2)
        handle.write('\n')
    print('Fresh K20 RS16/s20 replay: q=1..4096, precision=256, no saved endpoints.', flush=True)
    result = whole.replay(proposal, output, K=K, bits=STATE_BITS,
        group_dimension=256, prefix_end=PREFIX_END, precision=PRECISION,
        dense_tilts=tuple(map(Q, DENSE_TILTS)), markers=tuple(map(Q, MARKERS)),
        progress=lambda value: print(json.dumps(value), flush=True))
    if (result.get('fresh_replay') is not True or result.get('all_occupancies_covered') is not True
            or result.get('whole_code_certificate') is not True or result.get('target_met') is not True):
        raise ArithmeticError('fresh reproduction did not certify the complete requested 40-bit claim')
    print(json.dumps(dict(output=str(output), margin_bits=result['margin_bits'],
        whole_code_certificate=True, elapsed_seconds=monotonic()-start)), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    if args.dry_run:
        print(json.dumps(checked_plan(args.output), indent=2))
    else:
        reproduce(args.output)


if __name__ == '__main__':
    main()
