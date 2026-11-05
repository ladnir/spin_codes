"""FLOAT alias proposals for the literal field lift, never certificates.

Packet entries and positive coefficients stay in logs until alpha, marker
tilts, and state similarity scaling have all been applied. A normalized FFT
then computes regional residue sums. Small or negative inverse-transform
values are never claimed to be zero. The default bank uses the positive
regional bound for entries below its cancellation floor. Legacy-cushion
mode reproduces the older FLOAT proposal convention and labels that choice.
"""

from dataclasses import dataclass
from fractions import Fraction as F
from hashlib import sha256
from math import isfinite, log2
from pathlib import Path
from time import monotonic

import numpy as np

from hd_field128_lift_screen import (Budget, FIELD_ORDER, GROUPS, PAIRS, OUTPUT_LENGTH,
    NEG, checked_q, local_logs, log1mexp2, outer_log_density)
from hd_field128_lift_tail_screen import (PositiveBankScreen, PositivePacket,
    LOG_FACTORIAL, full_counts, source_pins as tail_source_pins)
from hd_randomized_outer_state_profile import multiply, power
from hd_shifted_pairs_screen import three_type_geometry


def checked_markers(rows):
    rows = tuple(dict.fromkeys(tuple(map(float, row)) for row in rows))
    if (not rows or len(rows) > 128 or any(len(row) != 2 or not all(map(isfinite, row))
                                         or any(abs(x) > 4096 for x in row) for row in rows)):
        raise ValueError('Expected 1..128 finite marker pairs with log magnitudes <=4096')
    return rows


def perron_logs(matrix):
    """Cancellation-safe log state ratio and positive row-cap scale."""
    a, b, c, d = np.asarray(matrix, dtype=float).ravel()
    if not all(map(isfinite, (a, b, c, d))):
        raise ArithmeticError('A strictly positive paired matrix is required for FFT scaling')
    high, low = max(a, d), min(a, d)
    gap = NEG if high == low else high+log1mexp2(low-high)
    discriminant = .5*float(np.logaddexp2(2*gap, 2+b+c))
    if a >= d:
        lambda_minus_a = 1+b+c-float(np.logaddexp2(discriminant, gap))
    else:
        lambda_minus_a = float(np.logaddexp2(discriminant, gap))-1
    ratio = lambda_minus_a-b
    rho = max(float(np.logaddexp2(a, b+ratio)), float(np.logaddexp2(c-ratio, d)))
    if not isfinite(ratio) or not isfinite(rho):
        raise ArithmeticError('Nonfinite log similarity witness')
    return rho, ratio


class AliasCoefficients:
    """Full unmarked single and good-good coefficients in log2 arithmetic."""

    def __init__(self, z, alpha, *, q=FIELD_ORDER, max_seconds=60., _budget=None):
        budget = _budget or Budget(max_seconds)
        self.packet = PositivePacket(z, alpha, q=q, _budget=budget)
        self.z, self.alpha, self.q = self.packet.z, self.packet.alpha, q
        counts, weights, degrees, indices = three_type_geometry()
        budget.check()
        half_weights = weights-counts[:, 0]-counts[:, 1]
        weights = (half_weights[:, None]+half_weights[None, :]).ravel()
        degrees = (degrees//33**2+(degrees//33) % 33)*33+degrees % 33
        packet_logs = local_logs(self.z, q).reshape((-1, 4))
        single = np.full((4, 33**2), NEG)
        for entry in range(4):
            budget.check()
            np.logaddexp2.at(single[entry], degrees, weights+self.alpha*packet_logs[indices, entry])
        self.single = single.reshape((4, 33, 33))
        s, f = self.packet.good_s, self.packet.good_f
        degrees = ((s[:, None]+s[None, :])*65+f[:, None]+f[None, :]).ravel()
        weights = self.packet.good_weights[:, None]+self.packet.good_weights[None, :]
        composed = np.full((4, 3*65), NEG)
        for entry in range(4):
            budget.check()
            values = weights+self.packet.composed_logs[:, :, entry]
            np.logaddexp2.at(composed[entry], degrees, values.ravel())
        self.composed = composed.reshape((4, 3, 65))
        budget.check()


@dataclass
class AliasTable:
    logs: np.ndarray
    reliable: np.ndarray
    positive_numerator: np.ndarray
    diagnostic: dict


def _frequency(log_coefficients, odds, log_factors, grid):
    degree = np.arange(log_coefficients.shape[1])[:, None]*odds[0]
    degree = degree+np.arange(log_coefficients.shape[2])[None, :]*odds[1]
    scaled = log_coefficients+degree[None, :, :]+np.asarray(log_factors)[:, None, None]
    finite = np.isfinite(scaled)
    if np.any(scaled[finite] > 1000):
        raise ArithmeticError('Normalized FFT component still exceeds its FLOAT range')
    with np.errstate(under='ignore'):
        coefficients = np.exp2(scaled)
    lost = finite & (coefficients == 0)
    polynomial = np.zeros((4, grid, grid))
    rows, columns = np.indices(log_coefficients.shape[1:])
    for entry in range(4):
        np.add.at(polynomial[entry], (rows.ravel() % grid, columns.ravel() % grid),
                  coefficients[entry].ravel())
    return np.fft.fftn(polynomial, axes=(1, 2)), int(lost.sum())


def _matrix_product(left, right):
    a, b, c, d = left
    e, f, g, h = right
    return np.stack((a*e+b*g, a*f+b*h, c*e+d*g, c*f+d*h))


def alias_table(z, alpha, odds, *, q=FIELD_ORDER, grid=128, cushion=1e-12,
                cancellation_floor=1e-10, max_seconds=60., _coefficients=None, _budget=None):
    """Full unnormalized log alias atlas plus reliability diagnostics."""
    checked_q(q, construction=True)
    odds = checked_markers((odds,))[0]
    if (type(grid) is not int or not 64 <= grid <= 256 or grid & (grid-1)
            or not isfinite(cushion) or not 0 < cushion <= 1e-8
            or not isfinite(cancellation_floor) or not 100*cushion <= cancellation_floor <= 1e-4):
        raise ValueError('Use grid64/128/256 and positive cushion/floor with floor>=100*cushion')
    budget = _budget or Budget(max_seconds)
    coefficients = _coefficients or AliasCoefficients(z, alpha, q=q, _budget=budget)
    if (coefficients.z, coefficients.alpha, coefficients.q) != (float(z), float(alpha), q):
        raise ValueError('Cached coefficient parameters differ from this witness')
    pair = coefficients.packet.pair_logs(odds, _budget=budget)
    log_rho, log_ratio = perron_logs(pair)
    log_vector = np.array((0., log_ratio))
    similarity = (log_vector[None, :]-log_vector[:, None]).ravel()
    good, dirty = coefficients.single.copy(), coefficients.single.copy()
    good[:, 2:, :] = NEG
    dirty[:, :2, :] = NEG
    transforms, underflows = [], 0
    for polynomial, scale in ((good, log_rho/2), (dirty, log_rho/2),
                              (coefficients.composed, log_rho)):
        budget.check()
        frequency, count = _frequency(polynomial, odds, similarity-scale, grid)
        transforms.append(frequency)
        underflows += count
    good, dirty, composed = transforms
    values = composed+_matrix_product(good, dirty)+_matrix_product(dirty, good)+_matrix_product(dirty, dirty)
    del transforms, good, dirty, composed
    row_sums = values[:, 0, 0].real.reshape((2, 2)).sum(axis=1)
    row_error = float(np.max(abs(row_sums-1)))
    if not np.all(np.isfinite(values)) or row_error > 1e-9:
        raise ArithmeticError('Normalized positive row sums failed their FLOAT check')
    for _ in range(6):
        budget.check()
        a, b, c, d = values
        trace, bc = a+d, b*c
        values = np.stack((a*a+bc, b*trace, c*trace, d*d+bc))
    aliases = np.fft.ifftn(values, axes=(1, 2))
    if not np.all(np.isfinite(aliases)):
        raise ArithmeticError('Nonfinite inverse FFT')
    imaginary = float(np.max(abs(aliases.imag)))
    minimum = float(np.min(aliases.real))
    if imaginary > 1e-10 or minimum < -1e-10:
        raise ArithmeticError('FLOAT alias cancellation diagnostic failed')
    reliable = aliases.real >= cancellation_floor
    # Every low entry receives a positive diagnostic cushion, never an exact
    # zero claim. The default bank ignores these entries in favor of a
    # positive-marker bound; legacy mode exposes the old proposal convention.
    logs = np.log2(np.maximum(aliases.real, 0)+cushion)
    logs += 64*log_rho-similarity[:, None, None]
    numerator = np.array(power(pair, 64)).ravel()
    budget.check()
    return AliasTable(logs, reliable, numerator, dict(
        grid=grid, log2_pair_scale=log_rho, log2_similarity_ratio=log_ratio,
        normalized_row_sum_error=row_error, maximum_imaginary_alias=imaginary,
        minimum_real_alias=minimum, cancellation_floor=cancellation_floor,
        diagnostic_cushion=cushion, reliable_entries=int(reliable.sum()),
        below_floor_entries=int(reliable.size-reliable.sum()),
        normalized_input_underflow_terms=underflows,
        underflow_interpretation='FLOAT terms below representation; not exact zeros',
        directed=False, whole_code_certificate=False))


class AliasBankScreen(PositiveBankScreen):
    """Same all-profile potential API, with explicit FFT fallback policy."""

    def __init__(self, z, alpha, marker_log_odds, log_ratio, *, q=FIELD_ORDER,
                 distance=F(7, 20), outer_measure='invertible-toeplitz', grid=128,
                 floor_policy='positive-fallback', max_seconds=60., max_table_bytes=384 << 20,
                 _coefficients=None, _budget=None):
        self.bits = checked_q(q, construction=True)
        self.q, self.distance = q, F(distance)
        if not 0 < self.distance < F(1, 2):
            raise ValueError('Expected a distance strictly between zero and one half')
        self.z, self.alpha, self.log_ratio = float(z), float(alpha), float(log_ratio)
        self.odds = checked_markers(marker_log_odds)
        if floor_policy not in ('positive-fallback', 'legacy-cushion'):
            raise ValueError('Unknown cancellation-floor policy')
        if (type(grid) is not int or grid not in (64, 128, 256)
                or not 1 <= max_table_bytes <= 512 << 20
                or len(self.odds)*4*grid*grid*9 > max_table_bytes):
            raise ValueError('Alias bank exceeds its explicit table memory budget')
        self.grid, self.floor_policy = grid, floor_policy
        self.outer_measure = outer_measure
        self.density = outer_log_density(q, outer_measure)
        self.forbidden = int(self.distance*OUTPUT_LENGTH)
        self.constant = -self.alpha*(self.bits+log1mexp2(-self.bits)+self.forbidden*log2(self.z))
        budget = _budget or Budget(max_seconds)
        coefficients = _coefficients or AliasCoefficients(z, alpha, q=q, _budget=budget)
        self.tables = tuple(alias_table(self.z, self.alpha, odds, q=q, grid=grid,
                                       _coefficients=coefficients, _budget=budget) for odds in self.odds)
        self.numerators = np.array([table.positive_numerator.reshape((2, 2)) for table in self.tables])
        self.regional_diagnostics = None
        budget.check()

    def regional(self, weight):
        if type(weight) is not int or not 1 <= weight <= GROUPS:
            raise ValueError('Selected weight must be in 1..8192')
        full = full_counts(PAIRS, weight)
        single, empty = weight-2*full, PAIRS-weight+full
        normalizer = LOG_FACTORIAL[PAIRS]-LOG_FACTORIAL[empty]-LOG_FACTORIAL[single]-LOG_FACTORIAL[full]
        best = np.full((4, len(full)), np.inf)
        selected_reliable = np.zeros((4, len(full)), dtype=bool)
        any_reliable = np.zeros((4, len(full)), dtype=bool)
        for odds, table in zip(self.odds, self.tables):
            numerator = table.logs[:, single % self.grid, full % self.grid]
            reliable = table.reliable[:, single % self.grid, full % self.grid]
            any_reliable |= reliable
            if self.floor_policy == 'positive-fallback':
                numerator = np.where(reliable, numerator, table.positive_numerator[:, None])
            candidate = numerator-normalizer-single*odds[0]-full*odds[1]
            improve = candidate < best
            selected_reliable[improve] = reliable[improve]
            np.minimum(best, candidate, out=best)
        self.regional_diagnostics = dict(
            floor_policy=self.floor_policy, selected_below_floor_entries=int((~selected_reliable).sum()),
            selected_positive_fallback_entries=(int((~selected_reliable).sum())
                                                if self.floor_policy == 'positive-fallback' else 0),
            selected_cushioned_entries=(int((~selected_reliable).sum())
                                       if self.floor_policy == 'legacy-cushion' else 0),
            entries_without_reliable_alias=int((~any_reliable).sum()), entries=int(best.size),
            no_zero_coefficient_claims=True)
        return best.T.reshape((-1, 2, 2))

    def evaluate(self, weight, **kwargs):
        result = super().evaluate(weight, **kwargs)
        result['regional_alias_diagnostics'] = dict(self.regional_diagnostics)
        return result

    def witness(self):
        return dict(super().witness(), grid=self.grid, floor_policy=self.floor_policy,
                    coefficient_method='normalized-FFT-alias-FLOAT',
                    atlas_diagnostics=[table.diagnostic for table in self.tables])


def source_pins():
    result = tail_source_pins()
    path = Path(__file__).resolve()
    result[path.name] = sha256(path.read_bytes()).hexdigest()
    return result


def main():
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--field-bits', type=int, default=128)
    parser.add_argument('--distance', default='7/20')
    parser.add_argument('--r', type=int, required=True)
    parser.add_argument('--z', required=True)
    parser.add_argument('--alpha', required=True)
    parser.add_argument('--marker', type=float, nargs=2, action='append', required=True)
    parser.add_argument('--log-ratio', type=float, default=-8.)
    parser.add_argument('--grid', type=int, default=128)
    parser.add_argument('--method', choices=('marginal', 'runs'), default='runs')
    parser.add_argument('--floor-policy', choices=('positive-fallback', 'legacy-cushion'), default='positive-fallback')
    parser.add_argument('--optimize-potential', action='store_true')
    parser.add_argument('--max-seconds', type=float, default=60.)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if not 8 <= args.field_bits <= 128 or args.field_bits % 8:
        parser.error('--field-bits must be a multiple of eight in 8..128')
    if args.output and args.output.exists():
        parser.error('Output exists; choose an unused filename')
    budget = Budget(args.max_seconds)
    screen = AliasBankScreen(float(F(args.z)), float(F(args.alpha)), args.marker, args.log_ratio,
                             q=1 << args.field_bits, distance=F(args.distance), grid=args.grid,
                             floor_policy=args.floor_policy, _budget=budget)
    result = screen.evaluate(args.r, method=args.method, optimize_potential=args.optimize_potential, _budget=budget)
    result.update(witness=screen.witness(), source_sha256=source_pins(),
                  elapsed_seconds=monotonic()-budget.started, max_seconds=args.max_seconds)
    encoded = json.dumps(result, indent=2, allow_nan=False)
    if args.output:
        with args.output.open('x', encoding='utf-8') as target:
            target.write(encoded+'\n')
    print(encoded)


if __name__ == '__main__':
    main()
