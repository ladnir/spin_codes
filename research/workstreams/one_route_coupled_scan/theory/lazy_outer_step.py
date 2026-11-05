"""Active envelopes when each input has its own fresh outer coefficient.

Requires the lazy Toeplitz exposure and retained nonzero pivot in
LAZY_OUTER_COEFFICIENTS.md. It is not valid for already-exposed outer inputs.
No files or computations at import.
"""

from fractions import Fraction as F

from check_two_groups import fresh_terms


def terminals(vector, mask, z, tau):
    psi = min(tau, F(256))
    k = mask.bit_count()
    reset = [[F(0) for _ in range(16)] for _ in range(5)]
    live = [[F(0) for _ in range(16)] for _ in range(5)]
    live[0][0] = z**k * max(vector[r][mask] for r in range(1, k + 1))
    for h in range(1, 16):
        # With both prices psi, this is z^|J| psi^-|(I union H)\(J intersect H')|.
        terms = fresh_terms(mask, h, z, psi)
        for r in range(1, h.bit_count() + 1):
            reset[r][h] = sum((weight for following, weight, _ in terms if not following), F(0))
            live[r][h] = sum((weight * max(vector[s][following]
                                          for s in range(1, min(r + k, following.bit_count()) + 1))
                              for following, weight, _ in terms if following), F(0))
    return reset, live


def pivot_envelope(vector, tau):
    """First-nonzero RS-row partition before ONE group's lazy column."""
    quiet = sum((tau**(-i) for i in range(4)), F(0))
    live = max(sum((tau**(-i) * max(vector[s][h] for s in range(max(1, r - i), r + 1))
                    for i in range(4)), F(0)) / vector[r][h]
               for h in range(1, 16) for r in range(1, h.bit_count() + 1))
    return [[quiet, F(0)], [F(0), live]]
