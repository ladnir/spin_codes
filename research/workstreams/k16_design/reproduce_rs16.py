"""Reproduce the complete K16 RS16 bound from repository sources and fixed maps.

From the repository root, run:

    python -B research/workstreams/k16_design/reproduce_rs16.py \
        --output tmp/rs16-reproduction/whole-p256.json

Python, python-flint, NumPy, and SciPy must already be installed, and the
retained t64/s16 selected-map declaration must be present. No saved
numerical receipt, remote host, native encoder, or benchmark is required.
The run uses 256-bit Arb arithmetic and takes a few minutes on a typical CPU.
Every generated filename must be fresh; existing files are never overwritten.

This entry point writes four explicitly proposal-only recipes containing
fixed tilt lists, geometry, and source hashes, but no numerical endpoints.
The recipe choice fields merely select tilts for the existing replay API;
they do not claim optimality or encode previously computed support bounds.
packet_rs_whole.replay then freshly computes q1, q2, q3..32, and q33..512,
checks their source/scope/coverage, and emits the full union receipt.

Use --dry-run to inspect the plan without importing numerical dependencies,
creating files, or running any bound computation. For example:

    python -B research/workstreams/k16_design/reproduce_rs16.py \
        --output tmp/rs16-reproduction/whole-p256.json --dry-run
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from time import monotonic


PRECISION = 256
Q1_TILTS = ('.00512', '.01024')
Q2_TILTS = ('.00512', '.01024')
TAIL_RECIPES = (
    (3, 32, ('.01024', '.0256', '.0512', '.1024')),
    (33, 512, ('.1024', '.1536', '.2048', '.3072', '.4096', '.6144',
               '.8192', '1.2288', '1.6384', '2.1972246')),
)


def output_paths(output):
    output = Path(output).resolve()
    if output.suffix.lower() != '.json':
        raise ValueError('the output filename must end in .json')
    proposals = [output.with_name(output.stem + suffix + '.json') for suffix in
                 ('-proposal-q1', '-proposal-q2', '-proposal-tail-0', '-proposal-tail-1')]
    replayed = [output.with_name(output.stem + suffix + '.json') for suffix in
                ('-q1', '-q2', '-tail-0', '-tail-1')]
    return output, proposals, replayed


def checked_plan(output):
    output, proposals, replayed = output_paths(output)
    coverage = [1, 2]
    for first, last, tilts in TAIL_RECIPES:
        if not tilts or len(set(tilts)) != len(tilts):
            raise ValueError('distinct nonempty tilt lists required')
        coverage.extend(range(first, last + 1))
    if sorted(coverage) != list(range(1, 513)):
        raise ValueError('the reproduction plan must cover each occupancy exactly once')
    paths = [output, *proposals, *replayed]
    if len(set(paths)) != len(paths):
        raise ValueError('all output paths must be distinct')
    existing = [str(path) for path in paths if path.exists()]
    if existing:
        raise FileExistsError('fresh output paths required: ' + ', '.join(existing))
    return dict(outer='rs16', precision=PRECISION, K=65536, N=131072,
        threshold=13107, distance='1/10', occupancy_covered=[1, 512],
        whole_code_certificate=False, dry_run_is_not_a_certificate=True,
        q1_tilts=list(Q1_TILTS), q2_tilts=list(Q2_TILTS),
        tail_recipes=[dict(q_min=first, q_max=last, tilts=list(tilts))
                      for first, last, tilts in TAIL_RECIPES],
        output=str(output), proposal_paths=list(map(str, proposals)),
        fresh_component_paths=list(map(str, replayed)))


def proposal_records(source_hashes):
    """Parameter-selection recipes, deliberately without numerical bounds."""
    common = dict(schema='rs16-standalone-tilt-proposal-1', proposal_only=True,
        choice_fields_are_tilt_selection_hints_only=True, numerical_endpoints_present=False,
        K=65536, N=131072, threshold=13107, distance='1/10', precision=PRECISION,
        geometry=dict(group_count=512, regions=64, group_dimension=128,
                      macro_windows=32, packet_bits=4),
        zero_initial_state=True, final_flush=False, whole_code_certificate=False,
        source_sha256=source_hashes)
    records = [dict(common, target_occupancy=[1], tilts=list(Q1_TILTS),
                    support_choices=list(Q1_TILTS)),
               dict(common, target_occupancy=[2], tilts=list(Q2_TILTS),
                    support_pair_choices=[list(Q2_TILTS)])]
    for first, last, tilts in TAIL_RECIPES:
        records.append(dict(common, q_min=first, q_max=last, tilts=list(tilts),
            occupancy_choices={f'tilt-selection-{i}': value for i, value in enumerate(tilts)}))
    return records


def reproduce(output):
    plan = checked_plan(output)
    # Import only after a read-only path preflight. The real replay, not these
    # recipe files, enumerates maps, rebuilds outer counts, and computes bounds.
    import packet_rs_whole as whole

    start = monotonic()
    output = Path(plan['output'])
    proposals = list(map(Path, plan['proposal_paths']))
    sources = whole.q1.source_snapshot()
    records = proposal_records(sources)
    output.parent.mkdir(parents=True, exist_ok=True)
    for path, record in zip(proposals, records):
        with path.open('x', encoding='utf-8', newline='\n') as handle:
            json.dump(record, handle, indent=2)
            handle.write('\n')
    print('Fresh RS16 reproduction: q=1..512, precision=256, no saved numerical bounds.', flush=True)
    result = whole.replay('rs16', proposals[0], proposals[1], proposals[2:],
                          precision=PRECISION, output=output)
    if (result.get('fresh_replay') is not True or result.get('all_occupancies_covered') is not True
            or result.get('whole_code_certificate') is not True or result.get('target_met') is not True):
        raise ArithmeticError('fresh reproduction did not certify the complete requested 40-bit claim')
    print(f'Reproduction complete in {monotonic() - start:.2f}s: {output}', flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--output', type=Path, required=True, help='fresh final JSON filename; component files are created alongside it')
    parser.add_argument('--dry-run', action='store_true', help='validate and display the plan without files or numerical work')
    args = parser.parse_args()
    if args.dry_run:
        print(json.dumps(checked_plan(args.output), indent=2))
    else:
        reproduce(args.output)


if __name__ == '__main__':
    main()
