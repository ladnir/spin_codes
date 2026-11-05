"""Directed positive atoms for the frozen bit-syndrome A8/F1 witness.

This layer receives fresh Arb regional upper bounds; it never builds a
provider or reads saved FLOAT bounds. Only exact logarithmic witness
parameters are read from the pinned proposal. Three slopes remain fixed.
Their atoms are rounded upward to exact dyadics, and the four intercepts
are raised outward until every regional entry is dominated. The 00 masses
are one exact dyadic p and its exact complement 1-p.

Primitive matrices and terminal weights use the same positive measure.
Information-cover powers are enclosed by Arb and checked by exact integer
power comparisons. Binary64 export uses explicit common power-of-two
scales and exact rational comparisons, including subnormal values. This
module does not run a column DP or claim a distance certificate.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isfinite, nextafter
from pathlib import Path

import numpy as np
from flint import arb, ctx

from hd_post_refresh_quiet_rows import transition_geometry


DIRECTORY = Path(__file__).resolve().parent
PROPOSAL = 'hd_bit_syndrome_orbit_993_x0_float.json'
PROPOSAL_SHA256 = '4d5d2374edca033663a18de558b1eef56751e9e4fc6b8081df5ec13c0f89e369'
ENDPOINTS = ((0, 0), (0, 0), (0, 1), (1, 0), (1, 1))
PARAMETERS = {
    'probability_log': F('-4264401345954727/140737488355328'),
    'atom00_log': F('-8719986785324751/1125899906842624'),
    'atom01_log': F('-2928190749607293/1125899906842624'),
    'atom10_log': F('-365228954560611/140737488355328'),
    'scale01_log': F('-93571663151647/1125899906842624'),
    'scale10_log': F('-282371531011411/2199023255552'),
    'scale11_log': F('-2499131197676909/35184372088832'),
}


def rational(value):
    value = F(value)
    return arb(value.numerator) / value.denominator


def upper_fraction(value):
    """Exact rational value of an Arb upper endpoint, not its midpoint."""
    if not isinstance(value, arb) or not value.is_finite():
        raise ValueError('A finite Arb enclosure is required')
    # Exact dyadic endpoints must not be rounded again at a possibly lower
    # ambient precision after their producing workprec block has ended.
    endpoint = value if value.is_exact() else value.upper()
    q = endpoint.fmpq()
    answer = F(int(q.p), int(q.q))
    if answer < 0:
        raise ValueError('A nonnegative upper endpoint is required')
    return answer


def _power_two(exponent):
    return F(1 << exponent) if exponent >= 0 else F(1, 1 << -exponent)


def _floor_log2(value):
    value = F(value)
    if value <= 0:
        raise ValueError('Positive rational required')
    e = value.numerator.bit_length() - value.denominator.bit_length()
    return e - int(value < _power_two(e))


def dyadic_up(value, bits):
    """Ceiling to a relative dyadic grid, proved by integer arithmetic."""
    value = F(value)
    if value < 0 or type(bits) is not int or bits < 2:
        raise ValueError('Nonnegative rational and at least two bits required')
    if not value:
        return value
    quantum = _power_two(_floor_log2(value) - bits + 1)
    scaled = value / quantum
    result = ((scaled.numerator + scaled.denominator - 1) // scaled.denominator) * quantum
    if result < value:
        raise ArithmeticError('Upward dyadic rounding failed')
    return result


def _precision(precision):
    if type(precision) is not int or not 128 <= precision <= 4096:
        raise ValueError('Explicit Arb precision128..4096 required')
    return precision


def frozen_parameters(proposal_path=None):
    """Authenticate literals and scope; no saved curve or estimate is used."""
    path = DIRECTORY / PROPOSAL if proposal_path is None else Path(proposal_path)
    raw = path.read_bytes()
    if sha256(raw).hexdigest() != PROPOSAL_SHA256:
        raise ValueError('Frozen atom-witness proposal changed')
    data = json.loads(raw)
    required = dict(mark='993/1000', regional_alpha=1, active_groups=18, arity=8,
        columns=16, log2_activity_marker=0.0, initial_state=0, terminal_vector=[1, 1],
        regional_full_pair_restriction='f<=min(1,floor(m/2))',
        matrix_entry_order=['00', '01', '10', '11'], complete_requested_control=True,
        directed=False, whole_code_certificate=False)
    if any(data.get(key) != value for key, value in required.items()):
        raise ValueError('Wrong frozen atom-witness scope')
    actual = {key: F(value) for key, value in data['state_measure']['exact_float_parameters'].items()}
    if actual != PARAMETERS:
        raise ValueError('The seven frozen exact logarithmic literals changed')
    return dict(actual)


def _curve_fractions(curves):
    values = np.asarray(curves, dtype=object)
    if values.shape != (4, 19):
        raise ValueError('Four complete regional curves at m0..18 required')
    answer = tuple(tuple(upper_fraction(value) for value in row) for row in values)
    if answer[0][0] > 1:
        raise ValueError('The 00 zeroth moment cannot exceed one')
    return answer


def build_atoms(curves, *, precision=256, proposal_path=None):
    """Choose a fixed exact dyadic measure dominating all fresh upper bounds.

    The input layout is (entry,m), entries00,01,10,11. Values must be Arb
    upper bounds, not logs. Choosing larger intercepts is a deterministic
    upward verification step, not an optimization or a saved-bound reuse.
    """
    precision = _precision(precision)
    parameters = frozen_parameters(proposal_path)
    caps = _curve_fractions(curves)
    with ctx.workprec(precision):
        # Reserve room for the complement of a small dyadic probability.
        bits = precision - 64
        candidates = {key: dyadic_up(upper_fraction(arb(2) ** rational(value)), bits)
                      for key, value in parameters.items()}
        b00, b01, b10 = (candidates[key] for key in ('atom00_log', 'atom01_log', 'atom10_log'))
        p = dyadic_up(max(candidates['probability_log'],
                         *(caps[0][m] / b00**m for m in range(1, 19))), bits)
        if p > 1:
            raise ValueError('Frozen 00 shape requires probability greater than one')
        c01 = dyadic_up(max(candidates['scale01_log'],
                           *(caps[1][m] / b01**m for m in range(19))), bits)
        c10 = dyadic_up(max(candidates['scale10_log'],
                           *(caps[2][m] / b10**m for m in range(19))), bits)
        c11 = dyadic_up(max(candidates['scale11_log'], *caps[3]), bits)
    atoms = (F(0), b00, b01, b10, F(1))
    masses = (1-p, p, c01, c10, c11)
    for m in range(19):
        moments = (F(1) if m == 0 else p*b00**m, c01*b01**m, c10*b10**m, c11)
        if any(moments[entry] < caps[entry][m] for entry in range(4)):
            raise ArithmeticError('Exact all-occupancy domination failed')
    if masses[0] + masses[1] != 1:
        raise ArithmeticError('The two 00 masses must sum exactly to one')
    return dict(atoms=atoms, masses=masses, endpoints=ENDPOINTS, precision=precision,
        maximum=18, frozen_log_parameters=parameters, proposal_sha256=PROPOSAL_SHA256,
        all_76_moments_verified=True, probability_complement_exact=True,
        saved_float_bounds_used=False, exact_dyadic_inputs={
            'atoms': list(map(str, atoms)), 'masses': list(map(str, masses))},
        whole_code_certificate=False)


def power_upper(base, exponent):
    """Dyadic upper bound on a nonnegative rational power, integer-checked.

    The caller controls workprec. For exponent n/d>0, u^d>=base^n is
    precisely the required direction because all quantities are positive.
    """
    base, exponent = F(base), F(exponent)
    if base < 0 or exponent < 0:
        raise ValueError('Nonnegative base and exponent required')
    if not exponent:
        return arb(1)
    if not base:
        return arb(0)
    upper = (rational(base) ** rational(exponent)).upper()
    exact = upper_fraction(upper)
    if exact**exponent.denominator < base**exponent.numerator:
        raise ArithmeticError('Rational-power upper endpoint check failed')
    return upper


def _measure(values):
    atoms, masses = tuple(map(F, values['atoms'])), tuple(map(F, values['masses']))
    if (len(atoms) != 5 or len(masses) != 5 or any(x < 0 for x in atoms + masses)
            or atoms[0] != 0 or atoms[4] != 1 or masses[0] + masses[1] != 1
            or tuple(values.get('endpoints', ENDPOINTS)) != ENDPOINTS):
        raise ValueError('Five nonnegative atoms with one exact 00 probability pair required')
    return atoms, masses


def _caps(table, arity=None):
    table = tuple(tuple(row) for row in table)
    if (not table or any(not row for row in table) or len({len(row) for row in table}) != 1
            or not 2 <= len(table[0]) <= 9
            or (arity is not None and len(table[0]) != arity+1)
            or any(type(n) is not int or n < 0 for row in table for n in row)):
        raise ValueError('Exact nonnegative integer coefficient rows required')
    return table


def _elementary(markers):
    values = [arb(1)]
    for marker in markers:
        values.append(arb(0))
        for j in range(len(values)-1, 0, -1):
            values[j] += marker*values[j-1]
    return values


def primitive_weights(measure, cover_values, quiet_caps, outside_caps, *, groups=18):
    """Arb upper primitive weights under ONE common five-atom measure.

    Production uses four cap rows, eight marked positions and groups18.
    Smaller cap tables/groups support independent exact test oracles.
    Terminal count ordering is transition_geometry(arity). All terminals
    use X=1 and omit the edge masses, which the ordered DP charges once.
    """
    atoms, masses = _measure(measure)
    precision = _precision(measure['precision'])
    quiet_caps = _caps(quiet_caps)
    arity = len(quiet_caps[0])-1
    outside_caps = _caps(outside_caps, arity)
    covers = tuple(map(F, cover_values))
    if (len(outside_caps) != len(quiet_caps) or any(row[0] for row in outside_caps)
            or not covers or any(c < 0 or c > 1 for c in covers)
            or type(groups) is not int or groups < 1):
        raise ValueError('Matched cap classes, cover exponents[0,1] and positive group count required')
    counts, _ = transition_geometry(arity)
    with ctx.workprec(precision):
        av = tuple(rational(value) for value in atoms)
        mv = tuple(rational(value) for value in masses)
        empty = np.array([[arb(0) for _ in range(2)] for _ in range(2)], dtype=object)
        for (i, j), mass in zip(ENDPOINTS, mv):
            empty[i, j] += mass
        # This is an exact identity of the chosen measure, not a numerical
        # cancellation between independently rounded intervals.
        empty[0, 0] = arb(1)
        cover_matrices = np.empty((len(covers), 2, 2), dtype=object)
        for ci, c in enumerate(covers):
            matrix = [[arb(0) for _ in range(2)] for _ in range(2)]
            for atom, mass, (i, j) in zip(atoms, mv, ENDPOINTS):
                if c:
                    inner = upper_fraction(power_upper(atom, 1/c))
                    weight = power_upper(1+255*inner, groups*c)
                else:
                    weight = rational(max(F(1), atom)**groups)
                matrix[i][j] += mass*weight
            for i in range(2):
                for j in range(2):
                    cover_matrices[ci, i, j] = matrix[i][j].upper()
        terminal_q = np.empty((len(quiet_caps), len(counts)), dtype=object)
        terminal_n = np.empty_like(terminal_q)
        for ci, count in enumerate(counts):
            markers = tuple(atom for atom, number in zip(av, count) for _ in range(number))
            elementary = _elementary(markers)
            for target, table in ((terminal_q, quiet_caps), (terminal_n, outside_caps)):
                for cls, coefficients in enumerate(table):
                    total = sum((coefficient*value for coefficient, value in zip(coefficients, elementary)), arb(0))
                    target[cls, ci] = (total**groups).upper()
        empty = np.array([[x.upper() for x in row] for row in empty], dtype=object)
        edges = np.array([x.upper() for x in mv], dtype=object)
    return dict(edge_masses=edges, empty_matrix=empty, cover_matrices=cover_matrices,
        terminal_quiet=terminal_q, terminal_nonquiet=terminal_n, count_order=counts,
        cover_values=covers, groups=groups, arity=arity, activity_marker=F(1),
        precision=precision, exact_dyadic_inputs={'atoms': list(map(str, atoms)),
            'masses': list(map(str, masses))}, common_measure=True,
        edge_masses_in_terminals=False, whole_code_certificate=False)


def export_scaled(values, *, exponent=None):
    """Export exact Arb uppers as binary64 mantissas times one power of two.

    Every returned mantissa is checked as an exact rational. A positive
    value below binary64's subnormal range becomes the minimum subnormal,
    never zero. Zero is preserved only for an exact zero upper endpoint.
    """
    array = np.asarray(values, dtype=object)
    exact = tuple(upper_fraction(value) for value in array.flat)
    if not exact:
        raise ValueError('A nonempty array of Arb upper endpoints required')
    if exponent is None:
        largest = max(exact)
        exponent = _floor_log2(largest)+1 if largest else 0
    if type(exponent) is not int:
        raise ValueError('The common binary scale must be an integer')
    scale = _power_two(exponent)
    converted = []
    for value in exact:
        target = value/scale
        if not target:
            converted.append(0.0)
            continue
        try:
            answer = float(target)
        except OverflowError as error:
            raise OverflowError('Requested export scale overflows binary64') from error
        if not isfinite(answer):
            raise OverflowError('Requested export scale overflows binary64')
        if F(answer) < target:
            answer = nextafter(answer, float('inf'))
        if not isfinite(answer) or F(answer) < target or answer <= 0:
            raise ArithmeticError('Exact binary64 upper-export verification failed')
        converted.append(answer)
    return dict(mantissas=np.array(converted, dtype=np.float64).reshape(array.shape), exponent=exponent)
