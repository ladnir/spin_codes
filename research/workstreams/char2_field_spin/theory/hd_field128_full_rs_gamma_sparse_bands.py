"""Directed gamma sparse bands, using one frozen regional down-set per witness.

Every integer weight receives its own prefix and complete matching-profile
sum. The frozen gamma provider's GG-only refinement is a valid weaker cap
for the same construction; this driver does not claim all-pair coefficients.
Immutable completed rows survive timeouts. Partial bands have no full bound.
"""

from fractions import Fraction as F
from hashlib import sha256
from importlib.metadata import version
import json
from pathlib import Path
import sys
from time import monotonic

from flint import ctx

import hd_field128_full_rs_gamma_certificate as certificate
import hd_field128_full_rs_state_sparse_bands as frozen_driver


base, local = certificate.base, certificate.local
_weight_from_regions = frozen_driver._weight_from_regions


def source_pins():
    pins = dict(certificate.source_pins(), **frozen_driver.source_pins())
    path = Path(__file__).resolve()
    pins[path.name] = sha256(path.read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def check_pins(pins):
    directory = Path(__file__).resolve().parent
    if any(sha256((directory/name).read_bytes()).hexdigest() != digest for name, digest in pins.items()):
        raise RuntimeError('A pinned gamma sparse-band helper changed during replay')


def write_new(path, value):
    return frozen_driver.write_new(path, value)


def batch_receipt(first, last, z, alpha, *, maximum=None, state=30,
                  distance=F(2, 5), precision=256, bits=80, paired=True,
                  failure_bits=40, max_seconds=600., witness_source=None,
                  checkpoint=None, progress=None):
    z, alpha, distance = map(base.exact, (z, alpha, distance))
    name = local.construction(state)
    maximum = last if maximum is None else maximum
    if (any(type(value) is not int for value in (first, last, maximum, precision, bits, failure_bits))
            or not 1 <= first <= last <= maximum <= certificate.MAXIMUM
            or not 128 <= precision <= 4096 or not 32 <= bits <= 512
            or not 0 <= failure_bits <= 4096 or not 0 < z < 1
            or not 0 < alpha <= 1 or not 0 < distance < F(1, 2) or type(paired) is not bool):
        raise ValueError('Invalid exact gamma sparse band, witness, or directed precision')
    pins, budget = source_pins(), base.Budget(max_seconds)
    witness = dict(z=str(z), alpha=str(alpha), precision=precision, packet_mantissa_bits=bits,
                   paired=paired, coefficient_maximum=maximum, coefficient_marker=str(F(maximum, base.PAIRS)))
    metadata = dict(construction=name, state_dimension=state, inner_width=64,
                    inner_emission=f'single-RS[64,{state}]', emission_minimum=65-state,
                    feedback='gamma A transpose', CA_zero=True,
                    gamma_extension_degree=16, gamma_base_field_order=256,
                    gamma_choice='any fixed element of degree16 over GF256',
                    principal_nonsingularity='deterministic symmetric unit-diagonal determinant degree-halving lemma',
                    fixed_map_field_order=256, fixed_nodes='GF256 additive subspace0..63',
                    state_nodes=f'128..{127+state}',
                    refresh=f'D2 H_{state} D1; fixed MDS Lagrange from0..{state-1} to{state}..{2*state-1}',
                    field_order=certificate.FIELD_ORDER, field_bits=128, outer_groups=base.GROUPS,
                    source_dimension=base.DIMENSION, output_length=base.OUTPUT_LENGTH,
                    forbidden_weight=int(distance*base.OUTPUT_LENGTH), target_distance=str(distance),
                    outer_measure='invertible-toeplitz',
                    outer_randomizer='independent full-field invertible4x4 Toeplitz per outer group and output column',
                    projective_divisor=str(certificate.FIELD_ORDER-1), pair_swaps=True,
                    state_persists_across_regions=True, initial_state='zero', final_state='unrestricted',
                    full_field_random_coefficients=True, joint_rank_return_cap=True,
                    composition_steps=2 if paired and alpha != 1 else 1,
                    all_pair_before_alpha=False, witness=witness, witness_source=witness_source,
                    coefficient_method='shared-gamma-positive-carry-free-Arb-truncated-bivariate',
                    local_caps='exact rational joint active/inactive return count; upward dyadic packet rounding',
                    old_numerical_bounds_reused=False, requested_failure_bits=failure_bits,
                    source_sha256=pins, pin_scope='complete recursive static local-import closure, including lazy imports',
                    runtime=dict(python=sys.version, python_flint=version('python-flint'), numpy=version('numpy')),
                    arithmetic='directed Arb with exact rational local caps', directed=True, whole_code_certificate=False)
    rows, failure, tables = [], None, 0
    phase, active = 'regional-table', None
    with ctx.workprec(precision):
        try:
            budget.check()
            if progress:
                progress(dict(stage=phase, coefficient_maximum=maximum, state=state))
            regions = certificate.regional_directed(z, alpha, maximum, state=state, bits=bits,
                                                    paired=paired, budget=budget)
            tables = 1
            phase = 'weight'
            for active in range(first, last+1):
                result = _weight_from_regions(active, regions, z, alpha, q=certificate.FIELD_ORDER,
                                              distance=distance, outer_measure='invertible-toeplitz', budget=budget)
                check_pins(pins)
                row = dict(metadata, **result,
                           meets_requested_selected_weight_target=F(result['log2_upper_rational']) <= -failure_bits,
                           elapsed_seconds=monotonic()-budget.started, max_seconds=max_seconds,
                           scope='all supports at this exact selected weight; shared down-set, no interpolation')
                rows.append(row)
                if checkpoint:
                    checkpoint(row)
                if progress:
                    progress(dict(stage='completed-weight', active_groups=active,
                                  log2_upper_rational=row['log2_upper_rational'], weights_evaluated=len(rows)))
        except (TimeoutError, ArithmeticError) as error:
            failure = dict(phase=phase, active_groups=active, error_type=type(error).__name__, message=str(error))
    check_pins(pins)
    complete = failure is None and len(rows) == last-first+1
    union = frozen_driver._union_ceiling(rows, precision)
    return dict(metadata, status='complete' if complete else 'incomplete',
                requested_weight_range=[first, last], weight_range=[first, last],
                coefficient_maximum=maximum, coefficient_marker=str(F(maximum, base.PAIRS)),
                regional_tables_built=tables, weight_receipts=rows, weights_evaluated=len(rows),
                every_integer_weight_evaluated=complete, complete_requested_range=complete,
                completed_weight_range=[first, first+len(rows)-1] if rows else None,
                pending_weight_ranges=[] if len(rows) == last-first+1 else [[first+len(rows), last]],
                completed_union_log2_upper_rational=union, log2_upper_rational=union if complete else None,
                meets_requested_band_target=complete and F(union) <= -failure_bits,
                covers_all_supports_in_range=complete, covers_all_even_odd_splits=complete,
                no_weight_interpolation=True, failure=failure, source_pins_unchanged_at_exit=True,
                elapsed_seconds=monotonic()-budget.started, max_seconds=max_seconds,
                scope='only completed integer weights in this sparse band; not a whole-code certificate')


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--first', type=int, required=True)
    parser.add_argument('--last', type=int, required=True)
    parser.add_argument('--maximum', type=int)
    parser.add_argument('--z', required=True)
    parser.add_argument('--alpha', required=True)
    parser.add_argument('--state', type=int, choices=certificate.STATES, default=30)
    parser.add_argument('--distance', default='2/5')
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--bits', type=int, default=80)
    parser.add_argument('--failure-bits', type=int, default=40)
    parser.add_argument('--max-seconds', type=float, default=600.)
    parser.add_argument('--unpaired', action='store_true')
    parser.add_argument('--witness-source', type=Path)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        parser.error('Output directory must be unused')
    raw = args.witness_source.read_bytes() if args.witness_source else None
    provenance = (dict(path=str(args.witness_source.resolve()), sha256=sha256(raw).hexdigest(),
                       role='witness provenance only; no numerical bound reused') if raw is not None else None)
    args.output_dir.mkdir(exist_ok=False)

    def checkpoint(row):
        if raw is not None and args.witness_source.read_bytes() != raw:
            raise RuntimeError('Pinned witness source changed')
        write_new(args.output_dir/f"weight_{row['active_groups']:04d}.json", row)

    result = batch_receipt(args.first, args.last, args.z, args.alpha, maximum=args.maximum,
                           state=args.state, distance=args.distance, precision=args.precision,
                           bits=args.bits, paired=not args.unpaired, failure_bits=args.failure_bits,
                           max_seconds=args.max_seconds, witness_source=provenance, checkpoint=checkpoint,
                           progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    if raw is not None and args.witness_source.read_bytes() != raw:
        raise RuntimeError('Pinned witness source changed')
    write_new(args.output_dir/'summary.json', result)
    print(json.dumps({key: result[key] for key in ('status', 'weights_evaluated', 'pending_weight_ranges',
                                                  'completed_union_log2_upper_rational', 'failure')}, allow_nan=False))
    if not result['complete_requested_range']:
        raise SystemExit(2)


if __name__ == '__main__':
    main()
