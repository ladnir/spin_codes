"""Directed alias banks for full-RS[64,s], with s in {16,24,32}.

The fixed maps use A from RS[64,s], C=A^T, and an MDS H_s between
independent nonzero full-field diagonals. The outer Toeplitz setup and
paired regional route are unchanged. These are distinct constructions
from the literal four-RS[16,4] code, including the single-RS state16 control.

The frozen exact local formulas supply all 65 occupancy caps. Upward
dyadic rounding precedes the 4225 chronological good-good matrix products;
each raw product is raised entrywise to alpha before its arrangement count.
No single or paired histogram is truncated. Reviewed interval DFT and
common-potential support helpers cover all matching profiles at a selected
weight in 1..8192. No FLOAT coefficient or old probability ceiling is used.
Each receipt is selected-weight only, not a whole-code certificate.
"""

import ast
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
from time import monotonic

from flint import arb, ctx

import hd_field128_full_rs_state_certificate as sparse
import hd_field128_full_rs_state_screen as local
import hd_field128_lift_alias_certificate as interval


base, tail = interval.base, interval.tail
STATES = local.STATES
DEFAULT_STATE = 24


def directed_parts(z, alpha, *, q=base.FIELD_ORDER, state=DEFAULT_STATE,
                   bits=80, joint=True, budget=None):
    """All 561 single histograms and 4225 ordered raw good-good products."""
    z, alpha = base.exact(z), base.exact(alpha)
    if not 0 < alpha <= 1:
        raise ValueError('Expected exact 0<alpha<=1')
    budget = budget or base.Budget()
    table = sparse.packet_table(z, 64, q=q, state=state, bits=bits, joint=joint, budget=budget)
    exponent = base.rational(alpha)
    powered = {j: tuple(base.upper(value**exponent) if not value.is_zero() else arb(0)
                        for value in entries) for j, entries in table.items()}
    single = {(s, f): tuple(count * value for value in powered[s + 2*f])
              for s, f, count in local.fixed16.geometry(128)}
    good = local.fixed16.geometry(128, good=True)
    composed = {}
    for s, f, count in good:
        budget.check()
        for t, g, other in good:
            values = composed.setdefault((s + t, f + g), [arb(0) for unused in range(4)])
            raw = base.matrix_product(table[s + 2*f], table[t + 2*g])
            for entry, value in enumerate(raw):
                if not value.is_zero():
                    values[entry] += count * other * value**exponent
    budget.check()
    return single, composed


def source_pins():
    """Complete recursive static local-import closure, including lazy imports."""
    directory = Path(__file__).resolve().parent
    pending, seen = [Path(__file__).name], {}
    while pending:
        name = pending.pop()
        if name in seen:
            continue
        raw = (directory / name).read_bytes()
        seen[name] = sha256(raw).hexdigest()
        for node in ast.walk(ast.parse(raw)):
            modules = ([item.name for item in node.names] if isinstance(node, ast.Import)
                       else [node.module] if isinstance(node, ast.ImportFrom) and node.module else [])
            for module in modules:
                target = module.split('.')[0] + '.py'
                if (directory / target).is_file() and target not in seen:
                    pending.append(target)
    return dict(sorted(seen.items()))


class AliasBank(interval.AliasBank):
    """New exact coefficient producer; frozen directed transform/support code."""

    def __init__(self, z, alpha, markers, *, state=DEFAULT_STATE, joint=True,
                 q=base.FIELD_ORDER, distance=F(2, 5), grid=128, precision=256,
                 bits=80, outer_measure='invertible-toeplitz', max_seconds=180.,
                 max_table_bytes=384 << 20, progress=None):
        self.construction = local.construction(state)
        self.state, self.joint = state, joint
        self.z, self.alpha, self.distance = map(base.exact, (z, alpha, distance))
        self.q, self.field_bits = q, base.checked_q(q)
        self.markers = tuple(dict.fromkeys(tuple(map(base.exact, row)) for row in markers))
        if (not 0 < self.z < 1 or not 0 < self.alpha <= 1 or not 0 < self.distance < F(1, 2)
                or not 1 <= len(self.markers) <= 128 or any(len(row) != 2 or min(row) <= 0 for row in self.markers)
                or type(precision) is not int or not 128 <= precision <= 4096
                or type(bits) is not int or not 32 <= bits <= 512 or type(joint) is not bool
                or type(grid) is not int or grid not in (64, 128, 256)
                or outer_measure not in base.OUTER_MEASURES
                or type(max_table_bytes) is not int or not 1 <= max_table_bytes <= 512 << 20
                or len(self.markers)*4*grid*grid*8 > max_table_bytes):
            raise ValueError('Invalid exact state alias-bank witness or retained-table budget')
        self.precision, self.bits, self.outer_measure, self.grid = precision, bits, outer_measure, grid
        self.forbidden = int(self.distance * base.OUTPUT_LENGTH)
        self.pins = source_pins()
        budget = base.Budget(max_seconds)
        with ctx.workprec(precision):
            if progress:
                progress(dict(stage='exact-full-RS-state-coefficients', state=state))
            single, composed = directed_parts(self.z, self.alpha, q=q, state=state,
                                               bits=bits, joint=joint, budget=budget)
            self.atlases, self.numerators, self.diagnostics = [], [], []
            for index, row in enumerate(self.markers):
                budget.check()
                if progress:
                    progress(dict(stage='directed-interval-alias', marker=index, grid=grid, state=state))
                atlas, fallback, diagnostic = interval.component_alias(
                    *interval.components(single, composed, row), grid, budget)
                self.atlases.append(atlas)
                self.numerators.append(fallback)
                self.diagnostics.append(diagnostic)
            self.density_log = tail.log2arb(base.rational(base.outer_density(q, outer_measure)))
            self.slope = base.rational(self.alpha)*self.density_log
            self.constant = -base.rational(self.alpha)*(tail.log2arb(arb(q - 1))
                             +self.forbidden*tail.log2arb(base.rational(self.z)))
            self.log_markers = tuple(tuple(tail.log2arb(base.rational(x)) for x in row) for row in self.markers)
        budget.check()
        self.check_pins()

    def check_pins(self):
        if source_pins() != self.pins:
            raise RuntimeError('A pinned transitive full-RS alias helper changed during replay')

    def metadata(self):
        result = super().metadata()
        result.update(construction=self.construction, inner_width=64, state_dimension=self.state,
                      inner_emission=f'single-RS[64,{self.state}]', emission_minimum=65 - self.state,
                      feedback='A transpose', CA_zero=True, initial_state='zero', final_state='unrestricted',
                      fixed_map_field_order=256, fixed_nodes='GF256 additive subspace0..63',
                      state_nodes=f'128..{127 + self.state}',
                      refresh=f'D2 H_{self.state} D1; H fixed MDS Lagrange from0..{self.state - 1} to{self.state}..{2*self.state - 1}',
                      refresh_diagonal_fullfield_products_per_packet=2*self.state,
                      single_histogram_terms=561, good_geometry_types=65, raw_good_pair_products=4225,
                      maximum_packet_occupancy=64, maximum_paired_occupancy=128,
                      pair_slots_per_packet=32, pair_slots_per_region=4096,
                      packets_per_region=128, composed_pairs_per_region=64, regions=256,
                      joint_rank_return_cap=self.joint,
                      coefficient_method='full-RS-state-positive-residue-aliases-directed-interval-radix2-DFT',
                      local_caps='exact rational generic MDS rank/fiber return bound; upward dyadic packet rounding',
                      no_coefficient_truncation=True, old_numerical_ceilings_reused=False,
                      scope='new full-RS variable-state construction; requested supports only; no whole-code or timing claim')
        result['witness'].update(state_dimension=self.state, inner_width=64, joint_rank_return_cap=self.joint)
        return result

    def weight_receipt(self, weight, ratio, *, method='runs', max_seconds=60.):
        result = super().weight_receipt(weight, ratio, method=method, max_seconds=max_seconds)
        result['scope'] = 'new full-RS variable-state construction; every support at this selected weight only; no encoder or timing claim'
        return result


def main():
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--r', type=int, required=True)
    parser.add_argument('--state', type=int, choices=STATES, required=True)
    parser.add_argument('--z', required=True, help='Exact rational or decimal witness')
    parser.add_argument('--alpha', required=True, help='Exact rational or decimal witness')
    parser.add_argument('--marker', nargs=2, action='append', required=True, help='Exact positive SINGLE FULL markers')
    parser.add_argument('--ratio', required=True, help='Exact positive common state potential')
    parser.add_argument('--method', choices=('marginal', 'runs'), default='runs')
    parser.add_argument('--distance', required=True)
    parser.add_argument('--field-bits', type=int, default=128)
    parser.add_argument('--grid', type=int, default=128)
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--bits', type=int, default=80)
    parser.add_argument('--simple-return', action='store_true')
    parser.add_argument('--build-seconds', type=float, default=180.)
    parser.add_argument('--max-seconds', type=float, default=60.)
    parser.add_argument('--max-table-bytes', type=int, default=384 << 20)
    parser.add_argument('--output', type=Path, required=True, help='New output path; exclusive create')
    args = parser.parse_args()
    if not 8 <= args.field_bits <= 128 or args.field_bits % 8:
        parser.error('Field bits must be a multiple of 8 in 8..128')
    if not 1 <= args.r <= base.GROUPS:
        parser.error('Selected weight must lie in 1..8192')
    if args.output.exists():
        parser.error('Output already exists')
    started = monotonic()
    bank = AliasBank(args.z, args.alpha, args.marker, state=args.state, joint=not args.simple_return,
                     q=1 << args.field_bits, distance=args.distance, grid=args.grid, precision=args.precision,
                     bits=args.bits, max_seconds=args.build_seconds, max_table_bytes=args.max_table_bytes,
                     progress=lambda row: print(json.dumps(row), flush=True))
    result = bank.weight_receipt(args.r, args.ratio, method=args.method, max_seconds=args.max_seconds)
    result.update(elapsed_seconds=monotonic() - started, build_seconds=args.build_seconds,
                  selected_weight_seconds=args.max_seconds, max_table_bytes=args.max_table_bytes)
    with args.output.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({key: result[key] for key in ('construction', 'active_groups', 'target_distance',
                                               'log2_upper_rational', 'directed', 'whole_code_certificate')}))


if __name__ == '__main__':
    main()
