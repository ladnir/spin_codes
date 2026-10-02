"""Geometry-explicit RS occupancy-tail comparison screen.

Positive Bernoulli components dominate every expected outer support shell.
A component belongs to a group for all regions; it is not resampled between
regions. A tilted, shuffled IID comparison bounds each interval of component
means. Point intervals do not cover the tail. Even a complete tail cover here
does not certify q=1, q=2, or the whole code.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
from fractions import Fraction as Q
import hashlib
import json
from math import comb, log
from pathlib import Path
from time import monotonic

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import logsumexp
from flint import arb, ctx

import packet_q1 as q1
import rs_outer
from heterogeneous import capped_density_loss


def logq(value):
    return log(value.numerator) - log(value.denominator)


def exact_mixture(n=8, k=4, packets_per_symbol=8):
    """One rational coefficient majorant per exact symbol-support layer."""
    regions = n * packets_per_symbol
    if regions != 64 or n != 2 * k:
        raise ValueError('this candidate requires 64 packet regions and half rate')
    layers = rs_outer.group_support_layers(n=n, k=k, packet_bits=4,
                                          packets_per_symbol=packets_per_symbol)
    counts = rs_outer.expected_group_support_counts(n=n, k=k, packet_bits=4,
                                                    packets_per_symbol=packets_per_symbol)
    components = []
    for layer in layers:
        p = Q(15 * layer.symbol_weight, 16 * n)
        mass = max(Q(c, layer.denominator) /
                   (comb(regions, u) * p**u * (1 - p)**(regions - u))
                   for u, c in enumerate(layer.numerators) if c)
        components.append((mass, p, 1))
    for u, count in enumerate(counts):
        upper = comb(regions, u) * sum((mass * p**u * (1 - p)**(regions - u)
                                        for mass, p, _ in components), Q(0))
        if upper < count:
            raise ArithmeticError(f'exact shell majorant fails at support {u}')
    # The actual zero message remains distinct from artificial all-zero
    # draws of a positive comparison component, which retains active label 1.
    return [(Q(1), Q(0), 0), *components], counts


def probabilities(p):
    return [1 - p, *(p * Q(comb(4, w), 15) for w in range(1, 5))]


def log_power(matrix, exponent):
    """Rescaled floating proposal; its exponent is always supplied explicitly."""
    row = np.zeros(len(matrix)); row[0] = 1
    row_log, power_log = 0., 0.
    power = matrix.copy()
    while exponent:
        if exponent & 1:
            row = row @ power
            scale = row.max()
            if not scale > 0:
                return -np.inf
            row /= scale
            row_log += power_log + np.log(scale)
        exponent >>= 1
        if exponent:
            power = power @ power
            scale = power.max()
            if not scale > 0:
                return -np.inf
            power /= scale
            power_log = 2 * power_log + np.log(scale)
    return row_log + np.log(row.sum())


def outer_dual(logs, features, active, target, occupancy, positive):
    """Floating proposal only; exact rational witnesses are replayed below."""
    features, active = np.asarray(features), np.asarray(active, dtype=bool)
    lo, hi = (0., 20000.) if positive else (-20000., 0.)

    def evaluate(eta):
        scores = logs + eta * features
        inactive_log = logsumexp(scores[~active])
        active_log = logsumexp(scores[active])
        mu = (20000. if occupancy >= 1 else
              np.clip(np.log(occupancy) - np.log1p(-occupancy) + inactive_log - active_log, 0., 20000.))
        scores = scores + mu * active
        norm = logsumexp(scores)
        return norm - eta * target - mu * occupancy, float(np.exp(scores - norm) @ features - target), float(mu)

    if evaluate(lo)[1] >= 0:
        eta = lo
    elif evaluate(hi)[1] <= 0:
        eta = hi
    else:
        eta = brentq(lambda e: evaluate(e)[1], lo, hi, xtol=1e-10)
    _, _, mu = evaluate(eta)
    return Q(round(eta * 10**8), 10**8), Q(round(mu * 10**8), 10**8)


class TailModel:
    def __init__(self, geometry, components, data, *, q_min=3, tilt=Q(3, 16)):
        if (not isinstance(geometry, q1.Geometry) or type(q_min) is not int
                or not 1 <= q_min <= geometry.group_count or not 0 < Q(tilt) <= 1):
            raise ValueError('checked geometry, occupancy, and positive base tilt required')
        self.geometry, self.components, self.data = geometry, tuple(components), data
        self.q_min, self.tilt = q_min, Q(tilt)
        if (not self.components or self.components[0] != (1, 0, 0)
                or any(c <= 0 or not 0 < p < 1 or active != 1 for c, p, active in self.components[1:])):
            raise ValueError('one genuine zero component and positive active Bernoulli components required')
        self.features = tuple(p * self.tilt / (1 - p + p * self.tilt) for _, p, _ in self.components)
        self.active = tuple(active for _, _, active in self.components)
        self.normalized_masses = tuple(c * (1 - p + p * self.tilt)**geometry.regions
                                       for c, p, _ in self.components)
        self.logs = np.array([logq(c) for c in self.normalized_masses])
        self.feature_array = np.array(list(map(float, self.features)))
        self.active_array = np.array(self.active)
        self.min_positive = (min(1 - f for f in self.features if f < 1),
                             min(f for f in self.features if f > 0))
        self.root = (Q(q_min, geometry.group_count) * min(self.features[1:]), max(self.features))
        self.threshold = geometry.N // 10
        self.packets = geometry.group_count * geometry.regions
        self.epochs = geometry.N // 128
        assert self.packets * 4 == geometry.N and self.epochs == geometry.regions * geometry.macros_per_region

    def check_cell(self, cell):
        if len(cell) != 2 or not self.root[0] <= cell[0] <= cell[1] <= self.root[1]:
            raise ValueError('ordered mean cell inside this model domain required')

    def loss(self, cell):
        self.check_cell(cell)
        groups = self.geometry.group_count
        caps = tuple(max(0, min(groups, int(groups * x / p)))
                     for x, p in zip((1 - cell[0], cell[1]), self.min_positive))
        return capped_density_loss(groups, caps)

    def weights(self, cell):
        self.check_cell(cell)
        return 1 - cell[0], cell[1] / self.tilt

    def propose(self, cell):
        self.check_cell(cell)
        groups = self.geometry.group_count
        best = None
        for positive in (False, True):
            target = cell[0] if positive else cell[1]
            eta, mu = outer_dual(self.logs, self.feature_array, self.active,
                                 float(target), self.q_min / groups, positive)
            outside = groups * (logsumexp(self.logs + float(eta) * self.feature_array +
                float(mu) * self.active_array) - float(eta * target)) - float(mu * self.q_min)
            if best is None or outside < best[0]:
                best = outside, eta, mu
        outside, eta, mu = best
        weights = self.weights(cell)
        scale = sum(weights)
        law = probabilities(weights[1] / scale)

        def objective(lam):
            return log_power(q1.kernel_t64.floating(self.data, law, lam), self.epochs) + lam * self.threshold

        grid = [0.000001, .0001, .001, .004, .016, .064, .256, 1., 4.]
        scores = [objective(lam) for lam in grid]
        index = int(np.argmin(scores))
        fit = minimize_scalar(objective, bounds=(grid[max(0, index - 1)], grid[min(len(grid) - 1, index + 1)]), method='bounded')
        lam = Q(max(1, round(float(fit.x) * 10**10)), 10**10)
        inner = objective(float(lam))
        value = (outside + inner + self.packets * logq(scale) +
                 self.geometry.regions * logq(self.loss(cell))) / np.log(2)
        return value, dict(tilt=str(self.tilt), parameters=list(map(str, (lam, eta, mu))))

    def outward(self, cell, witness):
        self.check_cell(cell)
        if Q(witness['tilt']) != self.tilt or len(witness['parameters']) != 3:
            raise ValueError('matching tilt and three rational parameters required')
        lam, eta, mu = map(Q, witness['parameters'])
        if lam <= 0 or mu < 0:
            raise ValueError('positive output tilt and nonnegative occupancy dual required')
        weights = self.weights(cell)
        scale = sum(weights)
        matrix = q1.kernel_t64.outward(self.data, probabilities(weights[1] / scale), lam)**self.epochs
        inner = sum((matrix[0, j] for j in range(matrix.ncols())), arb(0))
        total = sum((q1.kernel_t64.aq(c) * q1.kernel_t64.aq(eta * f + mu * active).exp()
                     for c, f, active in zip(self.normalized_masses, self.features, self.active)), arb(0))
        log_upper = (inner.log() + self.packets * q1.kernel_t64.aq(scale).log() +
            q1.kernel_t64.aq(lam) * self.threshold +
            self.geometry.regions * q1.kernel_t64.aq(self.loss(cell)).log() +
            self.geometry.group_count * total.log() -
            q1.kernel_t64.aq(self.geometry.group_count * min(eta * cell[0], eta * cell[1]) + mu * self.q_min))
        result = q1.kernel_t64.up(log_upper.exp())
        if not result.is_finite() or not result > 0:
            raise ArithmeticError('positive finite outward endpoint required')
        return result


def run(*, outer='rs8', precision=192, tilt=Q(3, 16), means=None,
        max_cells=0, target_bits=55, output=None):
    if precision < 128 or max_cells < 0 or target_bits < 1 or (output is not None and output.exists()):
        raise ValueError('valid precision, cell budget, target, and fresh output required')
    start = monotonic()
    ctx.prec = precision
    n, k, packets = (8, 4, 8) if outer == 'rs8' else (16, 8, 4) if outer == 'rs16' else (0, 0, 0)
    components, counts = exact_mixture(n, k, packets)
    geometry = q1.Geometry(512, 64, 128)
    data, map_record = q1.kernel_t64.prepare(birth_density='capped')
    model = TailModel(geometry, components, data, tilt=Q(tilt))
    sources = q1.source_snapshot()
    record = dict(schema='finite-packet-tail-screen-1', outer=outer, geometry=asdict(geometry),
        K=geometry.K, N=geometry.N, threshold=model.threshold, distance='1/10',
        groups=geometry.group_count, regions=geometry.regions, packets=model.packets,
        macros=model.epochs, q_min=model.q_min, q_max=geometry.group_count,
        proposed_occupancy_scope=[model.q_min, geometry.group_count],
        occupancy_excluded=[1, 2], whole_code_certificate=False, tail_cover_complete=False,
        zero_initial_state=True, final_flush=False, component_fixed_across_regions=True,
        precision=precision, base_tilt=str(tilt), root=list(map(str, model.root)),
        shell_majorant_verified=True, shell_counts_sha256=hashlib.sha256(json.dumps(list(map(str, counts))).encode()).hexdigest(),
        mixture=[dict(mass=str(c), activity=str(p), active=a) for c, p, a in components],
        map_record=map_record, source_sha256=sources, points=[], accepted_cells=[], unresolved_cells=[],
        scope='Positive shell comparison for q>=3 only, with independent ideal GL setup, '
              'uniform shuffles, selected t64/s16 maps, continuous unflushed state. '
              'Point bounds cover only their exact mean fibers. No numerical bound from '
              'an older code or geometry is reused; no whole-code claim is made.')

    def save():
        if sources != q1.source_snapshot():
            raise RuntimeError('loaded local mathematical source changed during screen')
        record['elapsed_seconds'] = monotonic() - start
        if output is not None:
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text(json.dumps(record, indent=2) + '\n')

    print(f'TAIL SCREEN {outer} q=3..512 K={geometry.K} N={geometry.N} '
          f'root={model.root} theta={tilt}', flush=True)
    chosen = ([model.root[0], Q('.0015'), Q('.002'), Q('.003'), Q('.005'), Q('.01'),
               Q('.02'), Q('.04'), Q('.08'), Q('.125'), Q('.2'), Q('.35'), Q('.5'),
               Q('.7'), model.root[1]] if means is None else list(map(Q, means)))
    for mean in chosen:
        if not model.root[0] <= mean <= model.root[1]:
            continue
        trial = monotonic()
        proposed, witness = model.propose((mean, mean))
        bound = model.outward((mean, mean), witness)
        bits = str(-bound.log() / arb(2).log())
        record['points'].append(dict(mean=str(mean), proposal_log2=proposed,
            upper=q1.endpoint(bound), margin_bits=bits, witness=witness,
            elapsed_seconds=monotonic() - trial))
        save()
        print(f'mean={mean} point_margin_bits={bits} '
              f'proposal_log2={proposed:.6f} elapsed={monotonic() - trial:.2f}s', flush=True)
    if max_cells:
        pending = [model.root]
        total = arb(0)
        visited = 0
        while pending and visited < max_cells:
            cell = pending.pop()
            proposed, witness = model.propose(cell)
            bound = model.outward(cell, witness)
            visited += 1
            if bound < arb(2)**(-target_bits):
                record['accepted_cells'].append(dict(cell=list(map(str, cell)),
                    witness=witness, upper=q1.endpoint(bound)))
                total = q1.kernel_t64.up(total + bound)
            else:
                middle = (cell[0] + cell[1]) / 2
                pending.extend([(middle, cell[1]), (cell[0], middle)])
            record.update(visited_cells=visited, unresolved_cells=[list(map(str, c)) for c in pending],
                accepted_cell_union_upper=q1.endpoint(total), tail_cover_complete=not pending)
            save()
            print(f'cover visited={visited} accepted={len(record["accepted_cells"])} '
                  f'pending={len(pending)} cell_margin={-bound.log() / arb(2).log()}', flush=True)
        if not pending:
            record['tail_upper'] = q1.endpoint(total)
            record['tail_margin_bits'] = str(-total.log() / arb(2).log())
            save()
    print('Whole-code certificate: false. q=1 and q=2 are excluded; point tests are not a cover.', flush=True)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outer', choices=('rs8', 'rs16'), default='rs8')
    parser.add_argument('--precision', type=int, default=192)
    parser.add_argument('--tilt', type=Q, default=Q(3, 16))
    parser.add_argument('--means', nargs='*', type=Q)
    parser.add_argument('--max-cells', type=int, default=0)
    parser.add_argument('--target-bits', type=int, default=55)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    run(outer=args.outer, precision=args.precision, tilt=args.tilt, means=args.means,
        max_cells=args.max_cells, target_bits=args.target_bits, output=args.output)


if __name__ == '__main__':
    main()
