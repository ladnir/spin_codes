"""Exact integer replay of the dense-state class bounds in PARTIAL_PROOF.md.

No floats, no experiment, and no whole-code certificate.
"""


def main():
    k = 1 << 20
    n, s, b, field = 2 * k, 64, 128, 1 << 128
    assert n % b == 0 and s * (n // b) == k
    # Bernoulli: (1-1/q)^R >= 1-R/q >= 1-K/q > 1/2.
    assert 2 * k < field
    intervals = 256 * 257 // 2
    interval_bits = (intervals - 1).bit_length()
    assert intervals < (1 << interval_bits)
    for numerator, denominator, deficit in [(49, 100, 2293), (491, 1000, 1245)]:
        weight = numerator * n // denominator
        exponent = n + 128 * (weight - k + 2 * deficit)
        margin = -exponent - 1  # density factor <2
        assert margin >= 40
        assert margin - interval_bits >= 40
        print(f"delta={numerator}/{denominator}, W={weight}, M<={deficit}: "
              f"whole-scan class <2^-{margin}; aligned-interval union <2^-{margin-interval_bits}")
    print(f"Region intervals={intervals} <2^{interval_bits}; exact integer checks PASS")
    print("Partial-state trajectories exceeding these budgets remain uncovered.")


if __name__ == '__main__':
    main()
