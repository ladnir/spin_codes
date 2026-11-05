"""Directed gamma alias bank with every packet pair composed before alpha.

Cache all65**2 chronological raw occupancy products raised to alpha once.
Aggregate all561**2 exact histogram multiplicities outside that power. The
full degree128 pair polynomial replaces GG/GD/DG/DD; both single components
are zero. The frozen interval DFT,64-pair regional power, and4096-slot
normalization are unchanged. This is a cap refinement, not a construction
change, and no FLOAT coefficient or prior numerical bound is imported.
"""

import ast
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from math import comb
from pathlib import Path

from flint import arb, ctx

import hd_field128_full_rs_gamma_alias_certificate as frozen


base, tail, local, interval = frozen.base, frozen.tail, frozen.local, frozen.interval
STATES, FIELD_ORDER = frozen.STATES, frozen.FIELD_ORDER


@lru_cache(4)
def histogram_geometry(slots=32):
    if type(slots) is not int or not 1 <= slots <= 32:
        raise ValueError('Expected1..32 pair slots')
    return tuple((s, f, comb(slots, s)*comb(slots-s, f))
                 for s in range(slots+1) for f in range(slots-s+1))


@lru_cache(4)
def grouped_pair_geometry(slots=32):
    """Group all ordered histogram pairs by occupancies and total histogram.

    Only exact integer multiplicities are combined here. Their common raw
    matrix product is powered once later, before these counts multiply it.
    """
    geometry = histogram_geometry(slots)
    by_occupancy = [[] for unused in range(2*slots+1)]
    for s, f, count in geometry:
        by_occupancy[s+2*f].append((f, count))
    result = []
    for j, left in enumerate(by_occupancy):
        for k, right in enumerate(by_occupancy):
            counts = [0]*(2*slots+1)
            for f, count in left:
                for g, other in right:
                    counts[f+g] += count*other
            result.append((j, k, tuple((j+k-2*f, f, count)
                                      for f, count in enumerate(counts) if count)))
    return tuple(result)


def all_pair_coefficients(table, alpha, *, slots=32, budget=None):
    """Generic positive-Arb coefficient interface for independent tiny tests."""
    histogram_geometry(slots)
    alpha = base.exact(alpha)
    side = 2*slots+1
    if (not 0 < alpha <= 1 or set(table) != set(range(side))
            or any(len(row) != 4 for row in table.values())):
        raise ValueError('Expected a complete occupancy table and exact0<alpha<=1')
    budget = budget or base.Budget()
    budget.check()
    rows = [tuple(base.upper(value) for value in table[j]) for j in range(side)]
    exponent = base.rational(alpha)
    powered = []
    for left in rows:
        budget.check()
        for right in rows:
            raw = base.matrix_product(left, right)
            powered.append(tuple(base.upper(value**exponent) if not value.is_zero() else arb(0)
                                 for value in raw))
    budget.check()
    grouped = grouped_pair_geometry(slots)
    budget.check()
    composed = {(s, f): [arb(0) for unused in range(4)]
                for s in range(side) for f in range(side-s)}
    for index, (j, k, histograms) in enumerate(grouped):
        if not index % side:
            budget.check()
        a, b, c, d = powered[j*side+k]
        for s, f, count in histograms:
            target = composed[s, f]
            target[0] += count*a
            target[1] += count*b
            target[2] += count*c
            target[3] += count*d
    budget.check()
    return {key: tuple(base.upper(value) for value in row) for key, row in composed.items()}


def directed_parts(z, alpha, *, state=24, bits=80, budget=None):
    budget = budget or base.Budget()
    z, alpha = base.exact(z), base.exact(alpha)
    if not 0 < alpha <= 1:
        raise ValueError('Expected exact0<alpha<=1')
    table = frozen.sparse.packet_table(z, 64, state=state, bits=bits, budget=budget)
    # Empty single components are identically zero; all routes are in the pair.
    return {}, all_pair_coefficients(table, alpha, budget=budget)


def source_pins():
    directory, pending, seen = Path(__file__).resolve().parent, [Path(__file__).name], {}
    while pending:
        name = pending.pop()
        if name in seen:
            continue
        raw = (directory/name).read_bytes()
        seen[name] = sha256(raw).hexdigest()
        for node in ast.walk(ast.parse(raw)):
            modules = ([item.name for item in node.names] if isinstance(node, ast.Import)
                       else [node.module] if isinstance(node, ast.ImportFrom) and node.module else [])
            for module in modules:
                target = module.split('.')[0]+'.py'
                if (directory/target).is_file() and target not in seen:
                    pending.append(target)
    return dict(sorted(seen.items()))


class AliasBank(frozen.AliasBank):
    """Fresh full-pair coefficients; inherited directed normalization and receipts."""

    def __init__(self, z, alpha, markers, *, state=24, q=FIELD_ORDER,
                 distance=F(2, 5), grid=128, precision=256, bits=80,
                 outer_measure='invertible-toeplitz', max_seconds=180.,
                 max_table_bytes=384 << 20, progress=None):
        self.construction = local.construction(state, q)
        self.state = state
        self.z, self.alpha, self.distance = map(base.exact, (z, alpha, distance))
        self.q, self.field_bits = q, base.checked_q(q)
        self.markers = tuple(dict.fromkeys(tuple(map(base.exact, row)) for row in markers))
        if (not 0 < self.z < 1 or not 0 < self.alpha <= 1 or not 0 < self.distance < F(1, 2)
                or not 1 <= len(self.markers) <= 128
                or any(len(row) != 2 or min(row) <= 0 for row in self.markers)
                or type(precision) is not int or not 128 <= precision <= 4096
                or type(bits) is not int or not 32 <= bits <= 512
                or type(grid) is not int or grid not in (64, 128, 256)
                or outer_measure not in base.OUTER_MEASURES
                or type(max_table_bytes) is not int or not 1 <= max_table_bytes <= 512 << 20
                or len(self.markers)*4*grid*grid*8 > max_table_bytes):
            raise ValueError('Invalid exact all-pair gamma alias witness or retained-table budget')
        self.precision, self.bits, self.outer_measure, self.grid = precision, bits, outer_measure, grid
        self.forbidden = int(self.distance*base.OUTPUT_LENGTH)
        self.pins = source_pins()
        budget = base.Budget(max_seconds)
        with ctx.workprec(precision):
            if progress:
                progress(dict(stage='exact-gamma-all-pair-coefficients', state=state))
            single, composed = directed_parts(self.z, self.alpha, state=state, bits=bits, budget=budget)
            self.atlases, self.numerators, self.diagnostics = [], [], []
            for index, row in enumerate(self.markers):
                budget.check()
                if progress:
                    progress(dict(stage='directed-interval-all-pair-alias', marker=index, grid=grid, state=state))
                atlas, fallback, diagnostic = interval.component_alias(
                    *interval.components(single, composed, row), grid, budget)
                self.atlases.append(atlas)
                self.numerators.append(fallback)
                self.diagnostics.append(diagnostic)
            self.density_log = tail.log2arb(base.rational(base.outer_density(q, outer_measure)))
            self.slope = base.rational(self.alpha)*self.density_log
            self.constant = -base.rational(self.alpha)*(tail.log2arb(arb(q-1))
                             +self.forbidden*tail.log2arb(base.rational(self.z)))
            self.log_markers = tuple(tuple(tail.log2arb(base.rational(x)) for x in row) for row in self.markers)
        budget.check()
        self.check_pins()

    def check_pins(self):
        if source_pins() != self.pins:
            raise RuntimeError('A pinned all-pair gamma alias helper changed during replay')

    def metadata(self):
        result = super().metadata()
        result.pop('good_geometry_types')
        result.pop('raw_good_pair_products')
        result.update(coefficient_method='gamma-all-pair-before-alpha-directed-interval-radix2-DFT',
                      pair_refinement='every ordered pair; all intermediate-state paths summed before alpha',
                      raw_occupancy_pair_products=65**2, raw_histogram_pair_products=561**2,
                      composed_histogram_terms=2145, single_components='identically zero',
                      counts_outside_alpha=True, all_pair_before_alpha=True,
                      construction_unchanged=True, no_coefficient_truncation=True,
                      bank_source_sha256=self.pins,
                      bank_pin_scope='complete all-pair subclass closure, retained through outer receipt wrappers')
        result['witness'].update(all_pair_before_alpha=True)
        return result


def main():
    import argparse
    import hd_field128_gamma_profile_certificate as profiles
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--r', type=int, required=True)
    parser.add_argument('--state', type=int, choices=STATES, default=24)
    parser.add_argument('--z', required=True)
    parser.add_argument('--alpha', required=True)
    parser.add_argument('--marker', nargs=2, action='append', required=True)
    parser.add_argument('--distance', default='2/5')
    parser.add_argument('--grid', type=int, default=128)
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--bits', type=int, default=80)
    parser.add_argument('--failure-bits', type=int, default=40)
    parser.add_argument('--build-seconds', type=float, default=600.)
    parser.add_argument('--max-seconds', type=float, default=600.)
    parser.add_argument('--max-table-bytes', type=int, default=384 << 20)
    parser.add_argument('--witness-source', type=Path)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists() or not 1 <= args.r <= base.GROUPS:
        parser.error('Expected an unused output directory and selected weight in1..8192')
    base.Budget(args.build_seconds)
    base.Budget(args.max_seconds)
    raw = args.witness_source.read_bytes() if args.witness_source else None
    provenance = (dict(path=str(args.witness_source.resolve()), sha256=sha256(raw).hexdigest(),
                       role='witness provenance only; no numerical bounds reused') if raw is not None else None)
    args.output_dir.mkdir(exist_ok=False)
    progress = lambda row: print(json.dumps(row, allow_nan=False), flush=True)
    try:
        bank = AliasBank(args.z, args.alpha, args.marker, state=args.state, distance=args.distance,
                         grid=args.grid, precision=args.precision, bits=args.bits,
                         max_seconds=args.build_seconds, max_table_bytes=args.max_table_bytes,
                         progress=progress)
    except (TimeoutError, ArithmeticError) as error:
        profiles.write_new(args.output_dir/'summary.json', dict(status='bank-build-incomplete',
                  construction=local.construction(args.state), active_groups=args.r,
                  complete_selected_weight=False, covers_all_supports_at_weight=False,
                  log2_upper_rational=None, directed=False, whole_code_certificate=False,
                  all_pair_before_alpha=True, source_sha256=source_pins(), witness_source=provenance,
                  failure=dict(error_type=type(error).__name__, message=str(error))))
        raise SystemExit(2)

    def checkpoint(row):
        if raw is not None and args.witness_source.read_bytes() != raw:
            raise RuntimeError('Pinned witness source changed')
        profiles.write_new(args.output_dir/f"even_full_{row['even_full']:04d}.json", row)

    result = profiles.weight_receipt(bank, args.r, max_seconds=args.max_seconds,
                                    failure_bits=args.failure_bits, witness_source=provenance,
                                    checkpoint=checkpoint, progress=progress)
    if raw is not None and args.witness_source.read_bytes() != raw:
        raise RuntimeError('Pinned witness source changed')
    result.update(build_seconds=args.build_seconds, max_table_bytes=args.max_table_bytes)
    profiles.write_new(args.output_dir/'summary.json', result)
    progress({key: result[key] for key in ('status', 'completed_even_rows', 'pending_even_full_range',
                                          'log2_upper_rational', 'failure')})
    if not result['complete_selected_weight']:
        raise SystemExit(2)


if __name__ == '__main__':
    main()
