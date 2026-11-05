"""FLOAT diagnostic: compose every gamma packet pair before taking alpha.

For a packet histogram h=(s,f), let P_h=P_(s+2f) be the positive local
2-by-2 cap before the fractional power. The paired polynomial is
sum_(h,k) count(h) count(k) u^(s_h+s_k) v^(f_h+f_k) (P_h P_k)^circ-alpha.
The matrix product is chronological; counts remain outside the power.
For 0<alpha<=1, (AB)^circ-alpha <= A^circ-alpha B^circ-alpha entrywise.
Thus this replaces the frozen GG/GD/DG/DD cap by a no-larger pair cap.

All 561 packet histograms and 314721 ordered histogram pairs are retained.
Only 4225 distinct raw occupancy products are needed. The complete paired
degree128 polynomial feeds the frozen alias engine, with its single-packet
components set to zero. The engine still raises the pair to power64 and
normalizes by the same4096-slot histogram. No construction or route changes.

This is only a FLOAT proposal diagnostic. It does not produce directed
receipts, whole-code certificates, or encoder performance measurements.
"""

import ast
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from math import comb, isfinite, log2
from pathlib import Path

import numpy as np

import hd_field128_full_rs_gamma_screen as gamma
import hd_field128_full_rs_gamma_alias_search as search
import hd_field128_gamma_profile_diagnostic as profiles


old, alias = gamma.old, gamma.frozen_alias.alias
FIELD_ORDER, NEG = gamma.FIELD_ORDER, gamma.NEG


@lru_cache(4)
def histogram_geometry(slots=32):
    """Ordered slot assignments; fair singleton orientations are averaged out."""
    if type(slots) is not int or not 1 <= slots <= 32:
        raise ValueError('Expected1..32 pair slots per packet')
    return tuple((s, f, comb(slots, s)*comb(slots-s, f))
                 for s in range(slots+1) for f in range(slots-s+1))


def all_pair_coefficients(logs, alpha, *, slots=32, _budget=None):
    """Full paired log polynomial from an occupancy-indexed local cap table.

    The generic slots parameter exists for tiny independent algebra tests.
    Public construction banks always use slots32 and gamma's q128 caps.
    """
    geometry = histogram_geometry(slots)
    logs = np.asarray(logs, dtype=float)
    if (logs.shape != (2*slots+1, 2, 2) or np.any(np.isnan(logs))
            or np.any(np.isposinf(logs)) or not isfinite(alpha) or not 0 < alpha <= 1):
        raise ValueError('Invalid occupancy table or fractional exponent')
    budget = _budget or old.Budget()
    budget.check()
    s = np.array([row[0] for row in geometry], dtype=int)
    f = np.array([row[1] for row in geometry], dtype=int)
    weights = np.array([log2(row[2]) for row in geometry])
    occupancy = s+2*f
    side = 2*slots+1
    # One state-path sum per ordered occupancy pair, before alpha.
    raw = np.logaddexp2(logs[:, None, :, 0, None]+logs[None, :, None, 0, :],
                       logs[:, None, :, 1, None]+logs[None, :, None, 1, :])
    raw = raw.reshape((side, side, 4))
    degrees = ((s[:, None]+s[None, :])*side+f[:, None]+f[None, :]).ravel()
    count_logs = weights[:, None]+weights[None, :]
    result = np.full((4, side*side), NEG)
    for entry in range(4):
        budget.check()
        values = count_logs+alpha*raw[occupancy[:, None], occupancy[None, :], entry]
        np.logaddexp2.at(result[entry], degrees, values.ravel())
    budget.check()
    return result.reshape((4, side, side))


class AliasCoefficients(gamma.AliasCoefficients):
    """Isolated adapter: the full pair polynomial replaces all four partitions."""

    def __init__(self, z, alpha, *, q=FIELD_ORDER, state=24, max_seconds=60., _budget=None):
        gamma.construction(state, q)
        if not isfinite(z) or not isfinite(alpha) or not 0 < z < 1 or not 0 < alpha <= 1:
            raise ValueError('Expected finite valid mark and fractional power')
        budget = _budget or old.Budget(max_seconds)
        self.z, self.alpha, self.q, self.state, self.joint = float(z), float(alpha), q, state, True
        self.gamma_degree, self.all_pairs = 16, True
        logs = gamma.local_logs(self.z, q, state=state)
        self.composed = all_pair_coefficients(logs, self.alpha, _budget=budget)
        # The frozen transform adds GD+DG+DD. Zeroing both single partitions
        # makes it transform exactly the complete paired polynomial above.
        self.single = np.full((4, 33, 33), NEG)
        self.packet = self
        budget.check()

    def pair_logs(self, odds, *, max_seconds=60., _budget=None):
        odds = alias.checked_markers((odds,))[0]
        budget = _budget or old.Budget(max_seconds)
        budget.check()
        degree = np.arange(65)
        tilt = degree[:, None]*odds[0]+degree[None, :]*odds[1]
        result = np.logaddexp2.reduce((self.composed+tilt[None, :, :]).reshape((4, -1)),
                                     axis=1).reshape((2, 2))
        if np.any(np.isnan(result)) or np.any(np.isposinf(result)):
            raise ArithmeticError('Invalid all-pair positive logarithm')
        budget.check()
        return result


class AliasBankScreen(gamma.AliasBankScreen):
    """Same gamma construction and frozen support engine, stronger pair caps."""

    def __init__(self, z, alpha, marker_log_odds, log_ratio, *, state=24, q=FIELD_ORDER,
                 max_seconds=60., _coefficients=None, _budget=None, **kwargs):
        budget = _budget or old.Budget(max_seconds)
        coefficients = _coefficients or AliasCoefficients(z, alpha, q=q, state=state, _budget=budget)
        if not isinstance(coefficients, AliasCoefficients) or not coefficients.all_pairs:
            raise ValueError('An all-pair gamma coefficient producer is required')
        super().__init__(z, alpha, marker_log_odds, log_ratio, state=state, q=q,
                         _coefficients=coefficients, _budget=budget, **kwargs)

    def metadata(self):
        result = super().metadata()
        result.pop('good_geometry_types')
        result.pop('raw_good_pair_products')
        result.update(coefficient_method='all-pair-before-alpha-normalized-FFT-alias-FLOAT',
                      pair_refinement='every ordered pair; all intermediate-state paths summed before alpha',
                      raw_occupancy_pair_products=4225, raw_histogram_pair_products=314721,
                      composed_histogram_terms=2145, maximum_packet_occupancy=64,
                      maximum_paired_occupancy=128, pair_slots_per_packet=32,
                      pair_slots_per_region=4096, packets_per_region=128,
                      composed_pairs_per_region=64, regions=256,
                      no_coefficient_truncation=True, construction_unchanged=True,
                      old_numerical_bounds_reused=False)
        return result


def selected_screen(bank, weight, *, support_method='runs', optimize_potential=False,
                    max_evaluations=40, _budget=None):
    """Choose frozen run/potential or sharpened-count chronological support."""
    budget = _budget or old.Budget()
    if support_method == 'runs':
        return dict(bank.evaluate(weight, method='runs', optimize_potential=optimize_potential,
                                  max_evaluations=max_evaluations, _budget=budget),
                    support_method='runs')
    if support_method != 'exact-profiles' or optimize_potential:
        raise ValueError('Use runs or exact-profiles; exact-profiles has no potential to optimize')
    regions = bank.regional(weight)
    result = profiles.support_sum(old.PAIRS, weight, regions, budget=budget)
    prefix = bank.alpha*bank.density*weight+bank.constant
    result.update(bank.metadata(), active_groups=weight, target_distance=str(bank.distance),
                  support_method='exact-profiles', log2_outer_prefix=prefix,
                  log2_upper_screen=prefix+result['log2_support_sum_screen'],
                  log2_potential_ratio=bank.log_ratio,
                  potential_role='witness retained only; chronological support does not use a common potential',
                  regional_alias_diagnostics=bank.regional_diagnostics,
                  scope='one selected weight; all-pair caps, sharpened profile upper counts, chronological transfer; FLOAT only')
    return result


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


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True, help='Saved FLOAT row; witness literals only')
    parser.add_argument('--r', type=int, required=True)
    parser.add_argument('--state', type=int, choices=gamma.STATES, default=24)
    parser.add_argument('--distance', default='2/5')
    parser.add_argument('--z')
    parser.add_argument('--alpha')
    parser.add_argument('--marker-weight', type=int)
    parser.add_argument('--grid', type=int)
    parser.add_argument('--log-ratio', type=float)
    parser.add_argument('--support-method', choices=('runs', 'exact-profiles'), default='runs')
    parser.add_argument('--optimize-potential', action='store_true')
    parser.add_argument('--max-evaluations', type=int, default=40)
    parser.add_argument('--max-seconds', type=float, default=120.)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output must be a new file')
    if args.support_method == 'exact-profiles' and args.optimize_potential:
        parser.error('exact-profiles does not use a common potential; omit --optimize-potential')
    raw, pins = args.source.read_bytes(), source_pins()
    proposal = search.prepare_proposal(json.loads(raw),
                                      dict(path=str(args.source.resolve()), sha256=sha256(raw).hexdigest()),
                                      r=args.r, state=args.state, distance=F(args.distance),
                                      z=args.z, alpha=args.alpha, marker_weight=args.marker_weight)
    witness = proposal['witness']
    odds = tuple(tuple(float(F(x)) for x in row) for row in witness['marker_log_odds_exact'])
    ratio = float(F(witness['log2_potential_ratio_exact'])) if args.log_ratio is None else args.log_ratio
    budget = old.Budget(args.max_seconds)
    bank = AliasBankScreen(float(F(witness['z'])), float(F(witness['alpha'])), odds, ratio,
                           state=args.state, distance=F(args.distance),
                           grid=witness['grid'] if args.grid is None else args.grid,
                           floor_policy='positive-fallback', _budget=budget)
    result = selected_screen(bank, args.r, support_method=args.support_method,
                             optimize_potential=args.optimize_potential,
                             max_evaluations=args.max_evaluations, _budget=budget)
    export = bank.witness()
    fitted_ratio = result['log2_potential_ratio']
    export.update(log2_potential_ratio=fitted_ratio, log2_potential_ratio_exact=str(F(fitted_ratio)))
    result.update(witness=export, proposal=proposal, source_sha256=pins)
    if source_pins() != pins or args.source.read_bytes() != raw:
        raise RuntimeError('A pinned all-pair source changed during the screen')
    search.write_new(args.output, result)
    print(json.dumps({key: result[key] for key in ('construction', 'active_groups', 'log2_upper_screen',
                                                  'directed', 'whole_code_certificate')}))


if __name__ == '__main__':
    main()
