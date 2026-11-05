"""Raw eight-cone caps for post-refresh coupled states.

Y=X+AU; (Unew,Vnew)=(V,M(U+CX)+V), with CA=0 and fresh MDS M.
State order Z,L,R,B,D,P,F,G. D is a flat nonzero diagonal; P is an
arbitrary nonzero diagonal point mixture; F has arbitrary nonzero first
coordinate and conditional flat nonzero second coordinate; G is an
arbitrary nonzero first coordinate with zero second coordinate. These
are positive comparison measures, not claims about actual independence.

The generic exact helper reuses frozen occupancy-only component bounds.
The production FLOAT adapter additionally uses the old T10/T11 row where
the post-refresh map permits it. No bank, state potential, or gate runs
on import. Neither helper is a directed distance certificate.
"""
from fractions import Fraction as Rat
from math import isfinite

import numpy as np

from hd_coupled_state_packet import _components, _validate


STATE_ORDER = ('Z', 'L', 'R', 'B', 'D', 'P', 'F', 'G')
Z, L, R, B, D, P, F, G = range(8)


def _assemble(occupancy, z, q, length, rank, *, exact, old_row=None):
    number = Rat if exact else float
    values = _components(occupancy, z, q, length, rank, exact=exact)
    beta, size = values['beta'], values['dimension']
    quiet, nz = values['quiet'], values['nonzero_feedback']
    emitted = z**values['weight']
    t10 = nz*values['return_mark']/size
    t11 = beta*values['uniform_live_emission']
    if old_row is not None:
        # The old density of W=M(U+CX) is independent of the reference V.
        # This permits the same T11 at B->B and T11/size at B->R.
        t10, t11 = min(t10, old_row[0]), min(t11, old_row[1])
    point_return = values['sup_point_return']
    point_live = beta*values['sup_point_live']
    matrix = [number(0)]*64
    matrix[8*Z+Z] = quiet*emitted
    matrix[8*Z+L] = beta*nz*emitted
    matrix[8*L+D] = quiet*emitted
    matrix[8*L+R] = beta*nz*emitted/size
    matrix[8*L+B] = beta*nz*emitted
    matrix[8*R+Z] = t10
    matrix[8*R+L] = t11
    matrix[8*B+D] = t10
    matrix[8*B+R] = t11/size
    matrix[8*B+B] = t11
    matrix[8*D+P] = t10
    # This is the fixed root-selected D decomposition. Do not mix its
    # entries with the alternative averaged-J decomposition into G/F.
    matrix[8*D+R] = point_live/size
    matrix[8*D+B] = point_live
    matrix[8*P+P] = point_return
    matrix[8*P+G] = point_live/size
    matrix[8*P+F] = point_live
    matrix[8*F+D] = point_return
    matrix[8*F+R] = point_live/size
    matrix[8*F+B] = point_live
    matrix[8*G+Z] = point_return
    matrix[8*G+L] = point_live
    return tuple(matrix), dict(values, old_return_cap=t10, old_live_cap=t11)


def packet_caps(occupancy, z, *, q=256, length=16, rank=4):
    """Exact rational flat64 cap under the frozen block-MDS assumptions.

    For each fixed nonzero refresh input, every output atom must be at
    most (q-1)^(-blocks*rank). P/G coefficients bound mass. An F coefficient
    bounds a mixture of conditional flat-second references, not just mass
    of an unrestricted pair distribution. Every cone has terminal mass one.
    """
    z = Rat(z)
    occupancy = _validate(occupancy, z, q, length, rank)
    return _assemble(occupancy, z, q, length, rank, exact=True)[0]


def packet_details(occupancy, z, *, q=256, length=16, rank=4):
    z = Rat(z)
    occupancy = _validate(occupancy, z, q, length, rank)
    return _assemble(occupancy, z, q, length, rank, exact=True)[1]


def raw_packet_table(z):
    """Production FLOAT8x8 table, with no geometry average or state scaling.

    Fresh nonzero X is independent of all past weighted histories after
    conditioning on the complete outer support word. Old T10/T11 only
    bound the unchanged marked map U->M(U+CX) from a flat nonzero U.
    In B, the flat V is independent of that whole marked map; translating
    W by V gives the stated D/R/B row, including the old empty sharpening.
    """
    if isinstance(z, bool) or not isfinite(z) or not 0 < z <= 1:
        raise ValueError('Finite mark in (0,1] required')
    from hd_xor_multigroup import nonzero_matrices, profile_arrays
    occupancies, _ = profile_arrays()
    _, _, old_caps = nonzero_matrices(float(z))
    table = {}
    for occupancy, old in zip(occupancies, old_caps):
        key = tuple(map(int, occupancy))
        values, _ = _assemble(key, float(z), 256, 16, 4, exact=False,
                              old_row=(old[1, 0], old[1, 1]))
        matrix = np.array(values).reshape(8, 8)
        if np.any(~np.isfinite(matrix)) or np.any(matrix < 0):
            raise ArithmeticError('Invalid post-refresh eight-cone cap')
        table[key] = matrix
    return table


def metadata():
    return dict(state_order=list(STATE_ORDER),
        recurrence='Y=X+AU; (Unew,Vnew)=(V,M(U+CX)+V)',
        field_order=256, packet_width=64, state_symbols=32, refresh_symbols=16,
        regional_alpha=1, diagonal_point_cone_is_arbitrary_mass=True,
        first_point_cone_is_arbitrary_mass=True,
        F_is_arbitrary_first_with_conditional_flat_nonzero_second=True,
        product_cones_assert_actual_independence=False, state_scaling_applied=False,
        D_decomposition='P:T10; R:beta*Lsup/Ddim; B:beta*Lsup',
        old_packet_row_reused_at=['R->Z,L', 'B->D,R,B', 'D->P'],
        inverse_refresh_empty_cap_used_for_B=True,
        terminal_vector=[1]*8, parameter_fits=0, directed=False,
        whole_code_certificate=False)
