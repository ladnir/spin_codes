"""Isolated FLOAT full-RS[64,s], s in16/24/28/30, with C'=gamma A^T.

Gamma has degree16 over GF256. Symmetry and the unit diagonal of A_D A_D^T
halve the degree of every principal determinant's square factor, so these
states have nonsingular active principal matrices for every coordinate set.
Only the returning-state cap changes. The new cap counts active and inactive
output-zero constraints jointly. No frozen provider is modified.

Public construction screens require GF(2^128). The generic exact local
oracle is conditional on the same principal-nonsingularity property; small
field tests supply a separate map/gamma satisfying it. No FLOAT result is
a directed certificate, whole-code distance claim, or encoder benchmark.
"""

import ast
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from math import comb, isfinite, log2
from pathlib import Path

import numpy as np

import hd_field128_full_rs_state_screen as frozen
import hd_field128_full_rs_state_alias_screen as frozen_alias


old = frozen.old
FIELD_ORDER, NEG = old.FIELD_ORDER, old.NEG
STATES = (16, 24, 28, 30)


def construction(state=24, q=FIELD_ORDER):
    if type(state) is not int or state not in STATES or q != FIELD_ORDER:
        raise ValueError('This construction uses state16/24/28/30 and GF(2^128) with gamma of degree16 over GF256')
    return f'toeplitz-fullfield-full-rs64-{state}-gamma16-feedback-alternating-swapped-pairs-k20-rate-half-v1'


@lru_cache(4)
def fixed_maps(state=24):
    construction(state)
    a = frozen.lagrange(range(64), range(128, 128+state))
    h = frozen.lagrange(range(state, 2*state), range(state))
    return a, tuple(zip(*a)), h


def joint_probability(j, c, e, q, state):
    """Nonzero syndrome, c active zeros and e inactive zeros, under uniform X_S.

    The hypothesis is that I+gamma C_D A_D is invertible for every D.
    The complement S\\D has rank u=min(state,j-c) under C. The syndrome
    intersection has dimension at most min(u,state-e); every syndrome has
    q**(j-c-u) preimages on this complement. The active entries are fixed.
    """
    if not 0 <= c <= j or not 0 <= e:
        raise ValueError('Invalid active/inactive zero counts')
    if e >= state:
        return F(0)
    u = min(state, j-c)
    dimension = min(u, state-e)
    return F((q**dimension-1)*q**(j-c-u), q**j)


@lru_cache(256)
def exact_moments(j, z, q=FIELD_ORDER, n=64, state=24):
    """Exact K,birth,return-cap,V; conditional local oracle, no field embedding claim."""
    quiet, birth, basic_return, live = frozen.exact_moments(j, z, q, n, state, False)
    z = F(z)
    inactive = n-j
    bound = sum((comb(j, c)*comb(inactive, e)*z**(n-c-e)*(1-z)**(c+e)
                 * joint_probability(j, c, e, q, state)
                 for c in range(j) for e in range(min(state, inactive+1))), F(0))
    return quiet, birth, min(basic_return, bound), live


def exact_packet(j, z, q=FIELD_ORDER, n=64, state=24):
    quiet, birth, ret, live = exact_moments(j, z, q, n, state)
    d = q**state-1
    density = F(d, (q-1)**state)
    return quiet, density*birth, ret/d, min(density*live, F(z)**(n-state+1)) if not j else density*live


@lru_cache(64)
def local_logs(z, q=FIELD_ORDER, n=64, state=24):
    """Positive log expansion: c=j terms are mathematically zero, not clipped."""
    bits = frozen.parameters(q, n, state)
    if not isfinite(z) or not 0 < z < 1:
        raise ValueError('Expected a finite mark strictly between zero and one')
    result = frozen.local_logs(z, q, n, state, False).copy()
    lz, lw = log2(z), log2(1-z)
    log_d = state*bits+old.log1mexp2(-state*bits)
    for j in range(1, n+1):
        inactive, terms = n-j, []
        for c in range(j):
            u = min(state, j-c)
            for e in range(min(state, inactive+1)):
                dimension = min(u, state-e)
                probability = (dimension-c-u)*bits+old.log1mexp2(-dimension*bits)
                terms.append(log2(comb(j, c)*comb(inactive, e))+(n-c-e)*lz+(c+e)*lw+probability)
        bound = float(np.logaddexp2.reduce(terms))
        result[j, 1, 0] = min(result[j, 1, 0], bound-log_d)
    if np.any(np.isnan(result)) or np.any(np.isposinf(result)):
        raise ArithmeticError('Invalid gamma local logarithm')
    result.setflags(write=False)
    return result


class AliasCoefficients:
    """Frozen width64 coefficient layout, new gamma-local caps only."""

    def __init__(self, z, alpha, *, q=FIELD_ORDER, state=24, max_seconds=60., _budget=None):
        construction(state, q)
        if not isfinite(z) or not isfinite(alpha) or not 0 < z < 1 or not 0 < alpha <= 1:
            raise ValueError('Expected finite valid mark and fractional power')
        budget = _budget or old.Budget(max_seconds)
        budget.check()
        self.z, self.alpha, self.q, self.state, self.joint = float(z), float(alpha), q, state, True
        self.gamma_degree = 16
        logs = local_logs(self.z, q, state=state)
        self.single = np.full((4, 33, 33), NEG)
        for s, f, count in frozen.fixed16.geometry(128):
            self.single[:, s, f] = log2(count)+self.alpha*logs[s+2*f].ravel()
        good = frozen.fixed16.geometry(128, good=True)
        s = np.array([row[0] for row in good], dtype=int)
        f = np.array([row[1] for row in good], dtype=int)
        weights = np.array([log2(row[2]) for row in good])
        raw = logs[s+2*f]
        composed = np.logaddexp2(raw[:, None, :, 0, None]+raw[None, :, None, 0, :],
                                 raw[:, None, :, 1, None]+raw[None, :, None, 1, :])
        degrees = ((s[:, None]+s[None, :])*65+f[:, None]+f[None, :]).ravel()
        values = weights[:, None, None]+weights[None, :, None]+self.alpha*composed.reshape((65, 65, 4))
        result = np.full((4, 3*65), NEG)
        for entry in range(4):
            budget.check()
            np.logaddexp2.at(result[entry], degrees, values[:, :, entry].ravel())
        self.composed = result.reshape((4, 3, 65))
        self.packet = self
        budget.check()

    pair_logs = frozen_alias.AliasCoefficients.pair_logs


class AliasBankScreen(frozen_alias.alias.AliasBankScreen):
    """The frozen alias/support computation, with explicit new construction metadata."""

    def __init__(self, z, alpha, marker_log_odds, log_ratio, *, state=24, q=FIELD_ORDER,
                 distance=F(9, 20), max_seconds=60., _coefficients=None, _budget=None, **kwargs):
        name = construction(state, q)
        budget = _budget or old.Budget(max_seconds)
        coefficients = _coefficients or AliasCoefficients(z, alpha, q=q, state=state, _budget=budget)
        if (not isinstance(coefficients, AliasCoefficients) or coefficients.gamma_degree != 16
                or coefficients.state != state):
            raise ValueError('A gamma-specific coefficient producer is required')
        self.state, self.joint = state, True
        super().__init__(z, alpha, marker_log_odds, log_ratio, q=q, distance=distance,
                         _coefficients=coefficients, _budget=budget, **kwargs)
        self.construction = name

    def metadata(self):
        return dict(construction=self.construction, state_dimension=self.state, inner_width=64,
                    emission_minimum=65-self.state, joint_rank_return_cap=True,
                    pair_swaps=True, state_persists_across_regions=True, full_field_random_coefficients=True,
                    single_histogram_terms=561, good_geometry_types=65, raw_good_pair_products=4225,
                    directed=False, whole_code_certificate=False, arithmetic='FLOAT',
                    feedback='gamma A transpose', gamma_extension_degree=16,
                    gamma_base_field_order=256, principal_nonsingularity='deterministic degree-halving lemma',
                    return_cap='joint active/inactive zero constraints with invertible active feedback block',
                    CA_zero=True, native_implementation=False)

    def evaluate(self, weight, **kwargs):
        return dict(super().evaluate(weight, **kwargs), **self.metadata(),
                    scope='new gamma construction, selected weight only; not a directed or whole-code certificate')

    def witness(self):
        return dict(super().witness(), **self.metadata())

    def box(self, ranges, **kwargs):
        return dict(super().box(ranges, **kwargs), **self.metadata())


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
    parser.add_argument('--source', required=True, type=Path, help='Saved FLOAT row or search summary; witness literals only')
    parser.add_argument('--r', required=True, type=int)
    parser.add_argument('--state', type=int, choices=STATES, default=24)
    parser.add_argument('--distance', default='9/20')
    parser.add_argument('--max-seconds', type=float, default=120.)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output must be a new file')
    pins, raw = source_pins(), args.source.read_bytes()
    artifact = json.loads(raw)
    row = artifact.get('best', artifact)
    witness = row['witness']
    z, alpha = float(F(witness['z'])), float(F(witness['alpha']))
    ratio = row.get('log2_potential_ratio', witness['log2_potential_ratio'])
    budget = old.Budget(args.max_seconds)
    bank = AliasBankScreen(z, alpha, witness['marker_log_odds'], ratio, state=args.state,
                           distance=F(args.distance), grid=witness['grid'], _budget=budget)
    result = bank.evaluate(args.r, method='runs', optimize_potential=True, _budget=budget)
    export = bank.witness()
    export.update(log2_potential_ratio=result['log2_potential_ratio'],
                  log2_potential_ratio_exact=str(F(result['log2_potential_ratio'])))
    result.update(witness=export, source=dict(path=str(args.source.resolve()), sha256=sha256(raw).hexdigest()),
                  old_numerical_bounds_reused=False, source_sha256=pins)
    if source_pins() != pins or args.source.read_bytes() != raw:
        raise RuntimeError('A pinned source changed during the screen')
    with args.output.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({key: result[key] for key in ('construction', 'active_groups', 'log2_upper_screen',
                                                  'log2_potential_ratio', 'directed', 'whole_code_certificate')}))


if __name__ == '__main__':
    main()
