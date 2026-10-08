"""Fixed affine-placement spectral receipt; NOT a distance certificate.

The64 odd generators define a bipartite Cayley graph on32768 source rows.
The directed-Arb test certifies its largest nonconstant positive eigenvalue.
It does not model triple coordinates as independent random samples.
"""

import argparse
from fractions import Fraction as F
from math import comb


ROWS = 32768
STRIDES = (
    12623, 13781, 1327, 8485, 15923, 13269, 9939, 15617,
    11733, 7157, 4563, 9235, 4579, 3107, 8209, 4815,
    10163, 3237, 2417, 10819, 15471, 3299, 11593, 14227,
    10361, 6701, 15631, 14507, 8535, 2041, 461, 3057,
    13069, 37, 16173, 10917, 7993, 10657, 2063, 6261,
    7265, 7819, 4669, 14679, 2989, 2637, 10487, 16033,
    3573, 9877, 9539, 4089, 10903, 6659, 9425, 14581,
    3003, 12613, 10389, 7933, 9513, 6025, 6205, 6119,
)
LAMBDA_UPPER = F(323, 1000)


def spectral_screen():
    import numpy as np
    kernel = np.zeros(ROWS)
    for stride in STRIDES:
        kernel[stride] += 1 / (2 * len(STRIDES))
        kernel[-stride] += 1 / (2 * len(STRIDES))
    eigenvalues = np.fft.fft(kernel).real
    mode = int(np.argmax(eigenvalues[1:])) + 1
    return dict(maximum=float(eigenvalues[mode]), mode=mode,
                bipartite_eigenvalue=float(eigenvalues[ROWS // 2]), rigorous=False)


def certify_spectrum(precision=96):
    from flint import arb, ctx
    if len(STRIDES) != 64 or len(set(STRIDES)) != 64 or any(a % 2 != 1 for a in STRIDES):
        raise ArithmeticError("Invalid fixed stride family")
    previous = ctx.prec
    ctx.prec = precision
    try:
        angle = 2 * arb.pi() / ROWS
        cosine = [(angle * j).cos() for j in range(ROWS // 2 + 1)]
        bound = arb(LAMBDA_UPPER.numerator) / LAMBDA_UPPER.denominator
        maximum_upper = arb(-1)
        maximum_mode = None
        for mode in range(1, ROWS // 2 + 1):
            value = arb(0)
            for stride in STRIDES:
                residue = (mode * stride) % ROWS
                value += cosine[min(residue, ROWS - residue)]
            value /= len(STRIDES)
            if not value < bound:
                raise ArithmeticError(("Spectral inequality not verified", mode, str(value)))
            if value.upper() > maximum_upper:
                maximum_upper, maximum_mode = value.upper(), mode
        return dict(rows=ROWS, generators=len(STRIDES), modes_checked=ROWS // 2,
                    lambda_upper=str(LAMBDA_UPPER), maximum_upper=str(maximum_upper),
                    maximum_mode=maximum_mode, precision=precision,
                    rigorous_spectral_inequality=True, distance_certificate=False)
    finally:
        ctx.prec = previous


def safe_column_bound(a=8192, cancelled_middle=0, pre_order=2):
    """All row subsets of sizea; one or two scans; explicit boundary allowance."""
    if pre_order not in (1, 2) or (pre_order == 1 and cancelled_middle):
        raise ValueError("Middle cancellations apply only to the two-prefix case")
    b, n = len(STRIDES), len(STRIDES) * ROWS
    alpha = F(a, ROWS)
    structural = n * (alpha * alpha + LAMBDA_UPPER * alpha * (1 - alpha))
    upper = structural.numerator // structural.denominator + pre_order * b + cancelled_middle
    return dict(a=a, pre_order=pre_order, check_capacity=(b // 2) * a, structural_upper=structural,
                boundary_allowance=pre_order * b, cancelled_middle_allowance=cancelled_middle,
                safe_column_upper=upper, fewer_safe_columns_than_checks=upper < (b // 2) * a)


def cancellation_tail_receipt(cutoff=1024, q=65536):
    """Factorial-moment bound for Binomial(N-1,1/(q-1)).

    Adjacent-equality indicators in an iid alphabet sequence are independent:
    any prescribed path-edge set imposes one independent equality per edge.
    """
    trials = len(STRIDES) * ROWS - 1
    numerator = comb(trials, cutoff)
    denominator = (q - 1) ** cutoff
    # Compare exact integers; display only the certified binary exponent.
    bits = denominator.bit_length() - numerator.bit_length() - 1
    while (numerator << (bits + 1)) < denominator:
        bits += 1
    while (numerator << bits) >= denominator:
        bits -= 1
    return dict(trials=trials, alphabet=q - 1, event=f"cancelled_middle >= {cutoff}",
                failure_upper_strictly_below_power_of_two=-bits,
                rigorous_integer_comparison=True, distance_certificate=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certify", action="store_true")
    args = parser.parse_args()
    print(certify_spectrum() if args.certify else spectral_screen())
    print(safe_column_bound(cancelled_middle=1023))
    print(safe_column_bound(pre_order=1))
    print(cancellation_tail_receipt())


if __name__ == "__main__":
    main()
