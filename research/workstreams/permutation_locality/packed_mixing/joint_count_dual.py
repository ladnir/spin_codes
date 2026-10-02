"""Checked bounds on E[f(J)] from one joint count distribution.

The caller supplies atom caps for a normalized law on J=0,...,n. Each raw
MGF pair (t,L) asserts E[exp(t*J)] <= exp(L). L must bound the LOG of the
entire MGF; an atom-corrected quantity is not a valid substitute. Optional
moment pairs (g,B) assert E[g(J)] <= B and permit signed feature vectors.

prepare() builds fresh reusable constraints. Prepared.bound() proposes a
dual with floating LP, then verifies every coordinate using Arb. Prepared
objects are precision-specific; witnesses contain only rational duals and
a power-of-two objective scale. Saved numerical bounds are never inputs.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from math import isfinite

import numpy as np
from scipy.optimize import linprog
from flint import arb, ctx

SCHEMA = 'joint-count-dual-1'


def _arb(value):
    if isinstance(value, arb):
        result = value
    else:
        if isinstance(value, (bool, np.bool_)):
            raise ValueError('booleans are not numeric premises')
        value = Q(value)
        result = arb(value.numerator)/value.denominator
    if not result.is_finite():
        raise ValueError('finite numerical premises required')
    return result


def _up(value):
    if not value.is_finite():
        raise ArithmeticError('nonfinite outward arithmetic')
    mantissa, exponent = value.upper().man_exp()
    return arb(int(mantissa))*arb(2)**int(exponent)


def _down(value):
    if not value.is_finite():
        raise ArithmeticError('nonfinite outward arithmetic')
    mantissa, exponent = value.lower().man_exp()
    return arb(int(mantissa))*arb(2)**int(exponent)


def _power_scale(value):
    """A nearby exact power of two, strictly above a positive upper bound."""
    value = _up(value)
    if value < 0:
        raise ValueError('nonnegative scale required')
    if value == 0:
        return 0
    mantissa, exponent = value.man_exp()
    return int(exponent)+int(mantissa).bit_length()


def _points(values, *, nonnegative=False):
    result = tuple(_up(_arb(value)) for value in values)
    if not result or nonnegative and any(value < 0 for value in result):
        raise ValueError('nonempty finite nonnegative vector required')
    return result


def _rational(value):
    if isinstance(value, (bool, np.bool_, float, np.floating)):
        raise ValueError('witness coefficients must be exact rationals')
    return Q(value)


@dataclass(frozen=True)
class Prepared:
    """Fresh immutable constraints; construct with prepare(), not directly.

    Let u_j be the effective atom cap and x_j=p_j/u_j when u_j>0. Then
    0<=x_j<=1, sum u_j*x_j=1, and each stored row satisfies A*x<=B.
    Zero-cap coordinates have p_j=0 and contribute zero everywhere.
    """
    _caps: tuple
    _rows: tuple
    _rhs: tuple
    _float_caps: np.ndarray
    _float_rows: np.ndarray
    _float_rhs: np.ndarray
    precision: int

    @property
    def size(self):
        return len(self._caps)

    @property
    def constraint_count(self):
        return len(self._rows)

    def _check_precision(self):
        if ctx.prec != self.precision:
            raise ArithmeticError('prepared constraints belong to a different precision')

    def _values(self, values):
        self._check_precision()
        values = _points(values, nonnegative=True)
        if len(values) != self.size:
            raise ValueError('one objective value per count required')
        return values

    def _objective(self, values, exponent):
        scale = arb(2)**exponent
        return tuple(_up(cap*value/scale) for cap, value in zip(self._caps, values)), scale

    def verify(self, values, witness):
        """Verify a proposed dual afresh; no solver or saved bound is trusted.

        For each coordinate derive gamma_j>=0 so that
        u_j*f_j/S <= alpha*u_j + sum beta_k*A_kj + gamma_j.
        Taking expectations in x gives S*(alpha+sum beta_k*B_k+sum gamma_j).
        The repair is charged to atom caps, never divided by tiny u_j.
        """
        values = self._values(values)
        if (not isinstance(witness, dict)
                or set(witness) != {'schema', 'size', 'scale_exponent', 'alpha', 'beta'}
                or witness['schema'] != SCHEMA or type(witness['size']) is not int
                or witness['size'] != self.size
                or type(witness['scale_exponent']) is not int
                or not -(1 << 30) <= witness['scale_exponent'] <= (1 << 30)
                or not isinstance(witness['beta'], (list, tuple))
                or len(witness['beta']) != self.constraint_count):
            raise ValueError('matching finite rational joint-count witness required')
        alpha = _rational(witness['alpha'])
        beta = tuple(map(_rational, witness['beta']))
        if any(value < 0 for value in beta):
            raise ValueError('inequality dual coefficients must be nonnegative')
        a, scale = self._objective(values, witness['scale_exponent'])
        alpha_ball, beta_balls = _arb(alpha), tuple(map(_arb, beta))
        active = [(coefficient, row, rhs) for coefficient, row, rhs in
                  zip(beta_balls, self._rows, self._rhs) if coefficient != 0]
        result = alpha_ball+sum((coefficient*rhs for coefficient, _, rhs in active), arb(0))
        for j, (value, cap) in enumerate(zip(a, self._caps)):
            deficit = _up(value-alpha_ball*cap-sum(
                (coefficient*row[j] for coefficient, row, _ in active), arb(0)))
            result += max(arb(0), deficit)
        upper = _up(scale*result)
        if upper < 0:
            raise ArithmeticError('incompatible probability premises or negative dual bound')
        self._check_precision()
        return upper

    def bound(self, values):
        """Return (checked upper, rational witness), with a verified fallback."""
        values = self._values(values)
        baseline = _up(sum((cap*value for cap, value in zip(self._caps, values)), arb(0)))
        exponent = _power_scale(baseline)
        witness = dict(schema=SCHEMA, size=self.size, scale_exponent=exponent,
                       alpha='0', beta=['0']*self.constraint_count)
        best = self.verify(values, witness)
        if baseline == 0:
            return best, witness
        objective, _ = self._objective(values, exponent)
        fit = linprog(-np.array([float(value) for value in objective]),
            A_ub=self._float_rows if self.constraint_count else None,
            b_ub=self._float_rhs if self.constraint_count else None,
            A_eq=self._float_caps[None, :], b_eq=[1.], bounds=(0., 1.), method='highs')
        if fit.success:
            coefficients = [-float(fit.eqlin.marginals[0]),
                            *(-float(x) for x in fit.ineqlin.marginals)]
            if all(isfinite(value) for value in coefficients):
                candidate = dict(witness, alpha=str(Q(str(coefficients[0]))),
                    beta=[str(Q(str(max(0., value)))) for value in coefficients[1:]])
                try:
                    checked = self.verify(values, candidate)
                except ArithmeticError:
                    checked = None
                if checked is not None and checked < best:
                    best, witness = checked, candidate
        return best, witness


def prepare(atom_caps, raw_log_mgfs=(), *, moment_bounds=()):
    """Build reusable fresh constraints for a normalized count law.

    raw_log_mgfs contains (exact tilt, outward LOG upper) pairs. A valid
    source is regional_count.mgf_upper before any tilted-atom correction.
    moment_bounds contains (signed feature vector, outward upper RHS) pairs.
    For example, (j,mean_hi), (-j,-mean_lo), and ((j-a)^2,second_hi).
    All vector positions correspond to the complete domain j=0,...,n.
    """
    precision = ctx.prec
    caps = tuple(min(arb(1), value) for value in _points(atom_caps, nonnegative=True))
    raw = []
    for item in raw_log_mgfs:
        if not isinstance(item, (list, tuple)) or len(item) != 2:
            raise ValueError('raw MGF premises are (exact tilt, log upper) pairs')
        tilt, logarithm = _rational(item[0]), _up(_arb(item[1]))
        raw.append((tilt, logarithm))
        caps = tuple(min(cap, _up((logarithm-_arb(tilt*j)).exp())) for j, cap in enumerate(caps))
    if _up(sum(caps, arb(0))) < 1:
        raise ValueError('atom/MGF caps cannot contain a normalized law')
    rows, rhs = [], []
    for tilt, logarithm in raw:
        # Lower coefficients and upper RHS define a valid relaxed constraint.
        rows.append(tuple(max(arb(0), _down(cap*(_arb(tilt*j)-logarithm).exp()))
                          for j, cap in enumerate(caps)))
        rhs.append(arb(1))
    for item in moment_bounds:
        if not isinstance(item, (list, tuple)) or len(item) != 2:
            raise ValueError('moment premises are (feature vector, upper) pairs')
        features, upper = tuple(map(_arb, item[0])), _up(_arb(item[1]))
        if len(features) != len(caps):
            raise ValueError('moment feature vector must cover every count')
        magnitude = max([_up(abs(upper)), *(_up(abs(value)) for value in features)])
        scale = arb(2)**_power_scale(magnitude)
        rows.append(tuple(_down(cap*feature/scale) for cap, feature in zip(caps, features)))
        rhs.append(_up(upper/scale))
    arrays = [np.array([float(value) for value in caps]),
              np.array([[float(value) for value in row] for row in rows]).reshape(len(rows), len(caps)),
              np.array([float(value) for value in rhs])]
    if any(not np.all(np.isfinite(array)) for array in arrays):
        raise ArithmeticError('constraint scaling failed to produce finite proposals')
    for array in arrays:
        array.flags.writeable = False
    if ctx.prec != precision:
        raise ArithmeticError('precision changed while preparing constraints')
    return Prepared(caps, tuple(rows), tuple(rhs), *arrays, precision)


def bound(values, atom_caps, raw_log_mgfs=(), *, moment_bounds=()):
    return prepare(atom_caps, raw_log_mgfs, moment_bounds=moment_bounds).bound(values)


def verify(values, atom_caps, raw_log_mgfs, witness, *, moment_bounds=()):
    return prepare(atom_caps, raw_log_mgfs, moment_bounds=moment_bounds).verify(values, witness)
