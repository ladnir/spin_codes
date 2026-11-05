"""FLOAT alias bank for full-RS[64,s] states16/24/32, never a certificate.

Only the coefficient producer is new. The frozen FFT normalization, residue
aliasing, positive fallback, histogram normalization, and all-profile common
potential code are reused without mutation. Source artifacts supply witness
literals only; their construction or numerical bounds are never inherited.
"""

import ast
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isfinite, log2
from pathlib import Path

import numpy as np

import hd_field128_full_rs_state_screen as local
import hd_field128_lift_alias_screen as alias
import hd_field128_lift_alias_cover as cover


old = local.old


class AliasCoefficients:
    """561 single histograms and65 good types; raw GG is composed before alpha."""

    def __init__(self, z, alpha, *, q=old.FIELD_ORDER, state=16, joint=True,
                 max_seconds=60., _budget=None):
        local.construction(state)
        old.checked_q(q, construction=True)
        if not isfinite(z) or not isfinite(alpha) or not 0 < z < 1 or not 0 < alpha <= 1:
            raise ValueError('Expected finite valid mark and fractional power')
        budget = _budget or old.Budget(max_seconds)
        budget.check()
        self.z, self.alpha, self.q, self.state, self.joint = float(z), float(alpha), q, state, joint
        logs = local.local_logs(self.z, q, state=state, joint=joint)
        self.single = np.full((4, 33, 33), old.NEG)
        for s, f, count in local.fixed16.geometry(128):
            self.single[:, s, f] = log2(count)+self.alpha*logs[s+2*f].ravel()
        good = local.fixed16.geometry(128, good=True)
        s = np.array([row[0] for row in good], dtype=int)
        f = np.array([row[1] for row in good], dtype=int)
        weights = np.log2([row[2] for row in good])
        raw = logs[s+2*f]
        composed = np.logaddexp2(raw[:, None, :, 0, None]+raw[None, :, None, 0, :],
                                 raw[:, None, :, 1, None]+raw[None, :, None, 1, :])
        degrees = ((s[:, None]+s[None, :])*65+f[:, None]+f[None, :]).ravel()
        values = weights[:, None, None]+weights[None, :, None]+self.alpha*composed.reshape((65, 65, 4))
        result = np.full((4, 3*65), old.NEG)
        for entry in range(4):
            budget.check()
            np.logaddexp2.at(result[entry], degrees, values[:, :, entry].ravel())
        self.composed = result.reshape((4, 3, 65))
        # The frozen alias interface only needs packet.pair_logs.
        self.packet = self
        budget.check()

    def pair_logs(self, odds, *, max_seconds=60., _budget=None):
        odds = alias.checked_markers((odds,))[0]
        budget = _budget or old.Budget(max_seconds)
        budget.check()
        tilt = np.arange(33)[:, None]*odds[0]+np.arange(33)[None, :]*odds[1]
        marked = self.single+tilt[None, :, :]
        good = np.logaddexp2.reduce(marked[:, :2, :].reshape((4, -1)), axis=1).reshape((2, 2))
        dirty = np.logaddexp2.reduce(marked[:, 2:, :].reshape((4, -1)), axis=1).reshape((2, 2))
        tilt2 = np.arange(3)[:, None]*odds[0]+np.arange(65)[None, :]*odds[1]
        result = np.logaddexp2.reduce((self.composed+tilt2[None, :, :]).reshape((4, -1)), axis=1).reshape((2, 2))
        for left, right in ((good, dirty), (dirty, good), (dirty, dirty)):
            result = np.logaddexp2(result, np.asarray(old.multiply(left, right)))
        if np.any(np.isnan(result)) or np.any(np.isposinf(result)):
            raise ArithmeticError('Invalid positive full-RS paired logarithm')
        budget.check()
        return result


class AliasBankScreen(alias.AliasBankScreen):
    """Frozen numerical alias/support engine with state-specific local caps."""

    def __init__(self, z, alpha, marker_log_odds, log_ratio, *, state=16, joint=True,
                 q=old.FIELD_ORDER, distance=F(7, 20), outer_measure='invertible-toeplitz',
                 grid=128, floor_policy='positive-fallback', max_seconds=60.,
                 max_table_bytes=384 << 20, _coefficients=None, _budget=None):
        self.construction = local.construction(state)
        self.state, self.joint = state, joint
        budget = _budget or old.Budget(max_seconds)
        coefficients = _coefficients or AliasCoefficients(z, alpha, q=q, state=state, joint=joint, _budget=budget)
        if (coefficients.state, coefficients.joint) != (state, joint):
            raise ValueError('Cached coefficients use a different state or return cap')
        super().__init__(z, alpha, marker_log_odds, log_ratio, q=q, distance=distance,
                         outer_measure=outer_measure, grid=grid, floor_policy=floor_policy,
                         max_table_bytes=max_table_bytes, _coefficients=coefficients, _budget=budget)

    def metadata(self):
        return dict(construction=self.construction, state_dimension=self.state, inner_width=64,
                    emission_minimum=65-self.state, joint_rank_return_cap=self.joint,
                    pair_swaps=True, state_persists_across_regions=True, full_field_random_coefficients=True,
                    single_histogram_terms=561, good_geometry_types=65, raw_good_pair_products=4225,
                    directed=False, whole_code_certificate=False, arithmetic='FLOAT')

    def evaluate(self, weight, **kwargs):
        return dict(super().evaluate(weight, **kwargs), **self.metadata(),
                    scope='selected active-group weight only; no native implementation or global-distance claim')

    def witness(self):
        return dict(super().witness(), **self.metadata())

    def box(self, ranges, **kwargs):
        return dict(super().box(ranges, **kwargs), **self.metadata())


def proposal_from_artifact(artifact, provenance, *, state, distance, old_source=False, z=None, alpha=None):
    """Copy only exact proposal literals; never import an old numerical ceiling."""
    construction = local.construction(state)
    witness = cover.witness_from_artifact(artifact, old=old_source)
    # Short rational proposals avoid carrying enormous exact z^16 denominators.
    mark = F(str(float(witness['z']))) if z is None else cover.base.exact(z)
    exponent = F(str(float(witness['alpha']))) if alpha is None else cover.base.exact(alpha)
    distance = F(distance)
    if not 0 < mark < 1 or not 0 < exponent <= 1 or not 0 < distance < F(1, 2):
        raise ValueError('Inadmissible new-construction proposal')
    return dict(schema='field128-full-rs-state-alias-witness-only-v1', construction=construction,
                state_dimension=state, inner_width=64, field_order=old.FIELD_ORDER, field_bits=128,
                target_distance=str(distance), pair_swaps=True, full_field_random_coefficients=True,
                source=provenance, source_role='exact proposal literals only', old_numerical_bounds_reused=False,
                status='new-construction untested witness proposal', directed=False, whole_code_certificate=False,
                state_potential_ratio=str(witness['ratio']),
                witness=dict(z=str(mark), alpha=str(exponent), grid=witness['grid'],
                             markers=[list(map(str, pair)) for pair in witness['markers']]))


def log_rational(value):
    value = F(value)
    if value <= 0:
        raise ValueError('Expected a positive exact marker or potential')
    return log2(value.numerator)-log2(value.denominator)


def source_pins():
    directory = Path(__file__).resolve().parent
    pending, seen = [Path(__file__).name], {}
    while pending:
        name = pending.pop()
        if name in seen:
            continue
        raw = (directory/name).read_bytes()
        seen[name] = sha256(raw).hexdigest()
        for node in ast.walk(ast.parse(raw)):
            modules = ([alias.name for alias in node.names] if isinstance(node, ast.Import)
                       else [node.module] if isinstance(node, ast.ImportFrom) and node.module else [])
            for module in modules:
                target = module.split('.')[0]+'.py'
                if (directory/target).is_file() and target not in seen:
                    pending.append(target)
    return dict(sorted(seen.items()))


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--old-source', action='store_true')
    parser.add_argument('--state', type=int, choices=local.STATES, required=True)
    parser.add_argument('--distance', required=True)
    parser.add_argument('--r', type=int, required=True)
    parser.add_argument('--z')
    parser.add_argument('--alpha')
    parser.add_argument('--grid', type=int)
    parser.add_argument('--log-ratio', type=float)
    parser.add_argument('--floor-policy', choices=('positive-fallback', 'legacy-cushion'), default='positive-fallback')
    parser.add_argument('--optimize-potential', action='store_true')
    parser.add_argument('--max-evaluations', type=int, default=40)
    parser.add_argument('--max-seconds', type=float, default=120.)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output must be a new file')
    budget = old.Budget(args.max_seconds)
    artifact, provenance = cover.read_json(args.source)
    proposal = proposal_from_artifact(artifact, provenance, state=args.state, distance=F(args.distance),
                                      old_source=args.old_source, z=args.z, alpha=args.alpha)
    w = proposal['witness']
    markers = tuple(tuple(log_rational(x) for x in pair) for pair in w['markers'])
    ratio = log_rational(proposal['state_potential_ratio']) if args.log_ratio is None else args.log_ratio
    bank = AliasBankScreen(float(F(w['z'])), float(F(w['alpha'])), markers, ratio,
                           state=args.state, distance=F(args.distance), grid=args.grid or w['grid'],
                           floor_policy=args.floor_policy, _budget=budget)
    result = bank.evaluate(args.r, method='runs', optimize_potential=args.optimize_potential,
                            max_evaluations=args.max_evaluations, _budget=budget)
    result.update(proposal=proposal, witness=bank.witness(), source_sha256=source_pins())
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('proposal', 'witness', 'source_sha256')}))


if __name__ == '__main__':
    main()
