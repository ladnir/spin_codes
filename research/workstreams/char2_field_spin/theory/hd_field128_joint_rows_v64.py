"""Exact positive active-row refinement of V60; same construction/event.

Only local_logs converts exact caps to nondirected FLOAT proposals.
No searches, coefficient tables, artifacts or jobs are produced here.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from math import comb, log2
from pathlib import Path

import hd_field128_wide_state_v60 as old

HERE = Path(__file__).resolve().parent
NAME = Path(__file__).name
FIELD_ORDER, WIDTH, STATES = old.FIELD_ORDER, old.WIDTH, old.STATES
FROZEN = {
    'hd_field128_wide_state_v60.py':
        '939ba095547900a6bc234f01426a48073f6097d0cfbfe16e398ddaf7a23813a0',
    'hd_field128_return_tail_moments_v57.py':
        '8c1ed22d4ccd3a7033bbd384f23e46291785c0c1bafc6f332a277b08aec61bef',
    '../HD_FIELD128_JOINT_ROWS.md':
        'edb927c905c49ea5ada8a8655340132816ac71f7055820a79ec1bab311ab06f1',
}


def verify_frozen():
    for module in (old, old.tail):
        if Path(module.__file__).resolve().parent != HERE:
            raise RuntimeError('Nonlocal mathematical provider')
    for name, digest in FROZEN.items():
        if sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise RuntimeError('Changed joint-row dependency: '+name)
    old.verify_frozen()


def source_pins():
    verify_frozen()
    pins = old.source_pins()
    for name, digest in {**FROZEN, NAME: sha256(Path(__file__).read_bytes()).hexdigest()}.items():
        if name in pins and pins[name] != digest:
            raise RuntimeError('Conflicting joint-row source: '+name)
        pins[name] = digest
    verify_frozen()
    return dict(sorted(pins.items()))


def check_source_pins(expected):
    if source_pins() != expected:
        raise RuntimeError('Joint-row source closure changed')


@dataclass(frozen=True)
class Row:
    active_zeros: int
    bounds: tuple
    clipped: tuple
    tails: tuple


@dataclass(frozen=True)
class RowTable:
    active: int
    inactive: int
    masses: tuple
    active_moments: tuple
    rows: tuple


def _marginal(j, q, n, state):
    old._geometry(j, q, n, state)
    moments = tuple(F(comb(j, c), q**c)*(1-F(1, q**min(state, j-c)))
                    for c in range(j+1))
    masses = tuple(sum(((-1)**(c-a)*comb(c, a)*moments[c]
                        for c in range(a, j+1)), F(0)) for a in range(j+1))
    if (any(value < 0 for value in masses) or masses[j] != 0
            or sum(masses, F(0)) != 1-F(1, q**min(state, j))):
        raise ArithmeticError('Invalid exact active-zero marginal')
    if any(sum((comb(a, c)*masses[a] for a in range(c, j+1)), F(0)) != moments[c]
           for c in range(j+1)):
        raise ArithmeticError('Active factorial inversion mismatch')
    return masses, moments


@lru_cache(maxsize=65)
def _row_table(j, q, n, state):
    old._geometry(j, q, n, state)
    masses, moments = _marginal(j, q, n, state)
    inactive = n-j
    choose_j = tuple(comb(j, c) for c in range(j+1))
    choose_r = tuple(comb(inactive, e) for e in range(inactive+1))
    joint = tuple(tuple(choose_j[c]*choose_r[e]
                        *old.joint_probability(j, c, e, q, n, state)
                        for e in range(inactive+1)) for c in range(j+1))
    if tuple(row[0] for row in joint) != moments:
        raise ArithmeticError('V60 active marginal mismatch')
    rows = []
    for a, mass in enumerate(masses):
        if not mass:
            continue
        choose_a = tuple(comb(a, c) for c in range(a+1))
        bounds = (mass,) + tuple(min(choose_r[e]*mass,
            min(joint[c][e]/choose_a[c] for c in range(a+1)))
            for e in range(1, inactive+1))
        maximum, clipped, tails = old.tail.factorial_envelope(bounds)
        if (len(bounds) != inactive+1 or clipped[0] != mass
                or len(tails) != maximum+1 or any(value < 0 for value in bounds)):
            raise ArithmeticError('Invalid exact row envelope')
        rows.append(Row(a, bounds, clipped, tails))
    return RowTable(j, inactive, masses, moments, tuple(rows))


def _powers(z, n):
    zp, wp, w = [F(1)], [F(1)], 1-z
    for unused in range(n):
        zp.append(zp[-1]*z)
        wp.append(wp[-1]*w)
    return zp, wp


def _row_cap(table, z, zp, wp):
    r, result = table.inactive, F(0)
    for row in table.rows:
        bernstein = sum((value*zp[r-e]*wp[e] for e, value in enumerate(row.clipped)
                         if value), F(0))
        tail = row.bounds[0]*zp[r]+(1-z)*sum(
            (value*zp[r-e] for e, value in enumerate(row.tails) if e and value), F(0))
        selected = min(bernstein, tail)
        if not row.bounds[0]*zp[r] <= selected <= row.bounds[0]:
            raise ArithmeticError('Exact inactive-row moment invariant failed')
        result += zp[table.active-row.active_zeros]*selected
    return result


@lru_cache(maxsize=128)
def _components(j, z, q, n, state):
    old._geometry(j, q, n, state)
    z = old._mark(z)
    base = old._moments(j, z, q, n, state, True)
    table = _row_table(j, q, n, state)
    zp, wp = _powers(z, n)
    rows = _row_cap(table, z, zp, wp)
    selected = min(base[2], rows)
    if not zp[n]*table.active_moments[0] <= selected <= base[2]:
        raise ArithmeticError('Exact refined return invariant failed')
    if rows > base[1]:
        raise ArithmeticError('Row marginal exceeds exact birth cap')
    return (base[0], base[1], selected, base[3]), base[2], rows


def _parameters(j, z, state):
    state, z = old._state(state), old._mark(z)
    old._geometry(j, FIELD_ORDER, WIDTH, state)
    return z, state


def factorial_rows(j, state):
    state = old._state(state)
    old._geometry(j, FIELD_ORDER, WIDTH, state)
    verify_frozen()
    try:
        return _row_table(j, FIELD_ORDER, WIDTH, state)
    finally:
        verify_frozen()


def exact_moments(j, z, state):
    """Quiet,birth,return,live before dividing return by d."""
    z, state = _parameters(j, z, state)
    verify_frozen()
    try:
        return _components(j, z, FIELD_ORDER, WIDTH, state)[0]
    finally:
        verify_frozen()


def cap_components(j, z, state):
    z, state = _parameters(j, z, state)
    verify_frozen()
    try:
        refined, base, rows = _components(j, z, FIELD_ORDER, WIDTH, state)
        return dict(v60_return=base, row_return=rows, selected=refined[2])
    finally:
        verify_frozen()


def _packet(moments, j, z, state):
    quiet, birth, ret, live = moments
    d = FIELD_ORDER**state-1
    density = F(d, (FIELD_ORDER-1)**state)
    return (quiet, density*birth, ret/d,
            min(density*live, z**(WIDTH-state+1)) if j == 0 else density*live)


def exact_packet(j, z, state):
    z, state = _parameters(j, z, state)
    return _packet(exact_moments(j, z, state), j, z, state)


def local_logs(z, state):
    """All65 packets to FLOAT; phase checks, no per-row I/O."""
    import numpy as np
    z, state = _parameters(0, z, state)
    verify_frozen()
    try:
        rows = []
        for j in range(WIDTH+1):
            moments = _components(j, z, FIELD_ORDER, WIDTH, state)[0]
            packet = _packet(moments, j, z, state)
            rows.append([log2(v.numerator)-log2(v.denominator) if v else -np.inf
                         for v in packet])
        result = np.asarray(rows, dtype=float).reshape((WIDTH+1, 2, 2))
        if np.any(np.isnan(result)) or np.any(np.isposinf(result)):
            raise ArithmeticError('Invalid exact-sourced local log')
        result.setflags(write=False)
        return result
    finally:
        verify_frozen()


def principal_bad_gamma_upper(state):
    return old.principal_bad_gamma_upper(state)


def construction(state):
    return old.construction(state)


def metadata(state):
    meta = old.metadata(state)
    meta.update(schema='wide-state-exact-joint-rows-v64', source_sha256=source_pins(),
        mathematical_provider_sha256=FROZEN[Path(old.__file__).name],
        return_refinement='min(V60, exact-active-row V57 clipped/tail sum)',
        exact_active_zero_marginal=True, nonnegative_exact_row_bounds=True,
        row_tables_independent_of_z=True, full_joint_LP_solved=False,
        construction_changed=False, new_setup_event=False,
        numerical_improvement_measured=False, all_profiles_replayed=False)
    return meta
