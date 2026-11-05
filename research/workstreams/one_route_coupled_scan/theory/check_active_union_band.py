#!/usr/bin/env python3
"""Exact partial certificate for one complete active-union size band.

This is NOT a full-code certificate. It covers all exact support masks,
amplitudes, weights <=734003, and marked excursion placements with
8096<=v<=8192 and 1<=c<=32 at the fixed K20 packet-four parameters.
ACTIVE_UNION_GATE.md supplies the count; the final comparison is integer.
The one-diagonal code and its setup distribution are unchanged.

Root owns execution. Use python -B, timeout60s and a256MiB memory cap.
No files or raw large integers are emitted. No parameter search.
"""

from fractions import Fraction
from math import comb
from time import monotonic


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def main():
    started = monotonic()
    n, b, B, p, field_bits = 16384, 128, 64, 4, 32
    N = n*b
    e = (35*N)//100
    vl, vh, cl, ch = 8096, 8192, 1, 32
    G = N//(p*B)
    tau = Fraction(33, 32)
    t = tau**p
    require(e == 734003 and vh*2 == n, "fixed rectangle geometry")
    require(1 < tau <= (1 << field_bits), "fractional marker")

    # C <= 2^L <= q^-m, with the entire volume/mask count inside the
    # fractional power. nu=log_q(tau) need never be approximated.
    u_min = b*vl
    exponent = 16*e - 14*u_min + 15*b*ch + 2048
    require(exponent < 0, "rectangle does not have negative amplitude exponent")
    m = (-exponent)//field_bits
    require(exponent <= -field_bits*m, "integer exponent rounding direction")

    # Bound J(v,c) uniformly. Both lower indices are below their modes,
    # so increasing the top and bottom arguments here only enlarges it.
    top_first, top_second = vh-cl, n-vl+ch
    require(ch-1 <= top_first//2 and ch <= top_second//2, "binomial monotonicity")
    max_layouts = comb(top_first, ch-1)*comb(top_second, ch)
    classes = (vh-vl+1)*(ch-cl+1)

    # rho<=1/2. Compute H at rho=1/2 by exact rational arithmetic.
    # For U made of whole packets, positive D(U) is at least p.
    H = sum((Fraction(comb(B,j), 1 << B)*t**max(0,j-B//2)
             for j in range(B+1)), Fraction())
    require(H > 1, "positive guarded moment")
    numerator_power = H.numerator**G
    denominator_power = H.denominator**G
    require(numerator_power > denominator_power, "exact excess sign")

    # classes*max_layouts*tau^-m * t/(t-1) * (H^G-1).
    # Avoid Fraction's unnecessary large gcd: compare numerator/denominator
    # directly, retaining the exact subtraction of the constant one.
    numerator = (classes*max_layouts*(numerator_power-denominator_power)
                 *tau.denominator**m*t.numerator)
    denominator = (denominator_power*tau.numerator**m
                   *(t.numerator-t.denominator))
    margin = denominator.bit_length()-numerator.bit_length()-1
    require(margin >= 40 and (numerator << margin) < denominator,
            "partial rectangle does not close40 bits")
    elapsed = monotonic()-started
    require(elapsed < 60, "60-second cap")
    print(dict(kind="exact_complete_active_union_band_not_full_code",
               K=N//2, N=N, B=B, packet=p, batch=b, field_bits=field_bits,
               bad_weight_threshold=e, v_range=(vl,vh), c_range=(cl,ch),
               classes=classes, amplitude_log2_upper=exponent,
               integer_q_penalty=m, tau=str(tau),
               failure_margin_at_least_bits=margin,
               exact_numerator_bits=numerator.bit_length(),
               exact_denominator_bits=denominator.bit_length(),
               seconds=round(elapsed,3)))
    print("All masks/weights/placements in this band are covered; all other (v,c) remain open here.")


if __name__ == "__main__":
    main()
