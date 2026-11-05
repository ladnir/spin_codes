"""Bounded directed full-RS alias bands with immutable complete row receipts.

One previously built state alias bank and one frozen RegionalCache serve
every explicitly requested integer weight. The frozen normalization,
support sum, exact potential candidates, and union arithmetic are reused.
The loop only adds complete row checkpoints and preserves pending weights
after a timeout. No endpoints are interpolated and no old bound is reused.
"""

from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from time import monotonic

from flint import ctx

import hd_field128_full_rs_state_alias_certificate as provider
import hd_field128_lift_alias_cover as cover


base, tail = provider.base, provider.tail


def source_pins():
    pins = dict(provider.source_pins(), **cover.source_pins())
    path = Path(__file__).resolve()
    pins[path.name] = sha256(path.read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def write_new(path, value):
    raw = (json.dumps(value, indent=2, allow_nan=False)+'\n').encode('utf-8')
    with Path(path).open('xb') as stream:
        stream.write(raw)
    return dict(path=str(Path(path).resolve()), sha256=sha256(raw).hexdigest())


def band_receipt(bank, first, last, ratios, *, method='runs', max_weights=256,
                 max_seconds=600., failure_bits=40, witness_source=None, checkpoint=None):
    cover.checked_range(first, last, max_weights)
    ratios = tuple(dict.fromkeys(map(base.exact, ratios)))
    if (not ratios or len(ratios) > 32 or min(ratios) <= 0
            or type(failure_bits) is not int or not 0 <= failure_bits <= 4096
            or method not in ('runs', 'marginal')):
        raise ValueError('Expected1..32 positive exact potentials and a valid support method/target')
    if not isinstance(bank, provider.AliasBank):
        raise ValueError('This driver requires the explicit full-RS state directed alias provider')
    bank.check_pins()
    pins = source_pins()
    budget = base.Budget(max_seconds)
    metadata = dict(bank.metadata(), source_sha256=pins, coverage_driver_sha256=pins,
                    witness_source=witness_source, potential_candidates=list(map(str, ratios)),
                    support_method=method, requested_failure_bits=failure_bits,
                    old_numerical_bounds_reused=False)
    cache, rows, failure = cover.RegionalCache(bank), [], None
    with ctx.workprec(bank.precision):
        for weight in range(first, last+1):
            try:
                budget.check()
                regions = cache.regional(weight, budget)
                candidates = []
                for ratio in ratios:
                    support = tail.support_sum(base.PAIRS, weight, regions, ratio, method=method, budget=budget)
                    value = bank.constant+bank.slope*weight+tail.log2arb(support)
                    candidates.append((base.rational_ceiling(value), ratio, str(value)))
                ceiling, ratio, enclosure = min(candidates, key=lambda item: item[0])
                budget.check()
                bank.check_pins()
                if source_pins() != pins:
                    raise RuntimeError('A pinned alias-band helper changed before the row checkpoint')
                row = dict(metadata, active_groups=weight, regional_histograms=len(regions),
                           log2_upper_rational=str(ceiling), log2_upper_enclosure=enclosure,
                           state_potential_ratio=str(ratio), covers_all_supports_at_weight=True,
                           covers_all_even_odd_splits=True, meets_requested_selected_weight_target=ceiling <= -failure_bits,
                           meets_40_bit_selected_weight_target=ceiling <= -40,
                           elapsed_seconds=monotonic()-budget.started, max_seconds=max_seconds,
                           scope='all supports at this exact selected weight; no interpolation or pending-weight claim')
                rows.append(row)
                if checkpoint:
                    checkpoint(row)
            except (TimeoutError, ArithmeticError) as error:
                failure = dict(active_groups=weight, error_type=type(error).__name__, message=str(error))
                break
    bank.check_pins()
    if source_pins() != pins:
        raise RuntimeError('A pinned alias-band helper changed during the band replay')
    ceiling = cover.union_ceiling(rows, bank.precision)
    next_weight, complete = first+len(rows), len(rows) == last-first+1
    return dict(metadata, requested_weight_range=[first, last], weight_range=[first, last],
                weight_receipts=rows, weights_evaluated=len(rows), every_integer_weight_evaluated=complete,
                completed_weight_range=[first, next_weight-1] if rows else None,
                pending_weight_ranges=[] if complete else [[next_weight, last]],
                complete_requested_range=complete, covers_all_supports_in_range=complete,
                covers_all_even_odd_splits=True, no_weight_interpolation=True,
                completed_union_log2_upper_rational=ceiling, log2_upper_rational=ceiling if complete else None,
                requested_range_meets_target=complete and F(ceiling) <= -failure_bits,
                weights_above_target=[row['active_groups'] for row in rows
                                      if F(row['log2_upper_rational']) > -failure_bits],
                regional_caches_built=1, source_pins_unchanged_at_exit=True, failure=failure,
                max_weights=max_weights, elapsed_seconds=monotonic()-budget.started, max_seconds=max_seconds,
                scope='only completed integer weights in this band; not a whole-code certificate')


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--first', type=int, required=True)
    parser.add_argument('--last', type=int, required=True)
    parser.add_argument('--state', type=int, choices=provider.STATES, default=24)
    parser.add_argument('--z', required=True)
    parser.add_argument('--alpha', required=True)
    parser.add_argument('--marker', nargs=2, action='append', required=True)
    parser.add_argument('--ratio', action='append', required=True)
    parser.add_argument('--distance', default='2/5')
    parser.add_argument('--field-bits', type=int, default=128)
    parser.add_argument('--grid', type=int, default=128)
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--bits', type=int, default=80)
    parser.add_argument('--method', choices=('runs', 'marginal'), default='runs')
    parser.add_argument('--max-weights', type=int, default=256)
    parser.add_argument('--failure-bits', type=int, default=40)
    parser.add_argument('--build-seconds', type=float, default=600.)
    parser.add_argument('--max-seconds', type=float, default=600.)
    parser.add_argument('--max-table-bytes', type=int, default=384 << 20)
    parser.add_argument('--witness-source', type=Path)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        parser.error('Output directory must be unused')
    if not 8 <= args.field_bits <= 128 or args.field_bits%8:
        parser.error('Field bits must be a multiple of8in8..128')
    cover.checked_range(args.first, args.last, args.max_weights)
    provenance = None
    if args.witness_source:
        raw = args.witness_source.read_bytes()
        provenance = dict(path=str(args.witness_source.resolve()), sha256=sha256(raw).hexdigest(),
                          role='exact witness provenance only; no numerical bound reused')
    args.output_dir.mkdir(exist_ok=False)
    bank = provider.AliasBank(args.z, args.alpha, args.marker, state=args.state,
                              q=1 << args.field_bits, distance=args.distance, grid=args.grid,
                              precision=args.precision, bits=args.bits, max_seconds=args.build_seconds,
                              max_table_bytes=args.max_table_bytes,
                              progress=lambda row: print(json.dumps(row), flush=True))

    def checkpoint(row):
        write_new(args.output_dir/f"weight_{row['active_groups']:04d}.json", row)
        print(json.dumps({key: row[key] for key in ('active_groups', 'log2_upper_rational',
                                                   'state_potential_ratio', 'meets_requested_selected_weight_target')},
                         allow_nan=False), flush=True)

    result = band_receipt(bank, args.first, args.last, args.ratio, method=args.method,
                          max_weights=args.max_weights, failure_bits=args.failure_bits,
                          max_seconds=args.max_seconds, witness_source=provenance, checkpoint=checkpoint)
    result.update(build_seconds=args.build_seconds, max_table_bytes=args.max_table_bytes)
    write_new(args.output_dir/'summary.json', result)
    print(json.dumps({key: result[key] for key in ('weights_evaluated', 'pending_weight_ranges',
                                                  'completed_union_log2_upper_rational', 'failure')}, allow_nan=False))


if __name__ == '__main__':
    main()
