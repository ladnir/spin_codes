"""Exact arithmetic checks for DENSE_SUPPORT_BOUND.md; no file output.

The analytic row-factor argument covers every width and support. Here we
check all width-four rows and selected width-eight rows against the source
transfer, then verify the finite K20 inequality using integer arithmetic.
"""

from fractions import Fraction as F

from coupled_rank_transfer import local_row, pair_antipodal_cycle, rotate_state_mask


def require(condition, context):
    if not condition:
        raise AssertionError(context)


def check_rows():
    tau, z = 256, F(7, 3341)
    a = F(1, tau)
    f = z + a + z * a
    factors = {
        (0, 0): (F(1), z),
        (0, 1): (F(1), 2 * z),
        (1, 0): (z, z + a),
        (1, 1): (z * (1 + a), f),
    }
    count = 0
    for width in (4, 8):
        full = (1 << width) - 1
        masks = range(1 << width) if width == 4 else (0, 1, 3, 5, 15, 85, 127, 255)
        for cycle in (tuple(range(width)), pair_antipodal_cycle(width)):
            for state in masks:
                rotated = rotate_state_mask(width, state, cycle)
                for inputs in masks:
                    product = F(1)
                    for j in range(width):
                        product *= factors[((state >> j) & 1, (rotated >> j) & 1)][(inputs >> j) & 1]
                    actual = sum(local_row(width, inputs, state, z, tau, cycle))
                    require(actual <= product, (width, cycle, inputs, state, actual, product))
                    require(product <= f**width / z**(width - inputs.bit_count()),
                            ("missing-coordinate bound", width, cycle, inputs, state))
                    if inputs == full:
                        require(actual <= f**width, ("full row", width, state))
                    count += 1
    return count


def check_finite_bound():
    tau, z = 256, F(7, 3341)
    f = z + (1 + z) / tau
    require(f == F(5, 832), "chosen f")
    require(0 < z <= F(1, tau), "row-factor parameter domain")
    b = tau**50 / z**35 * f**100 * F(100, 99)**100 * (99 / z)
    require(b < F(1, 256), "one-percent band bound")
    left = 256**51 * 3341**36 * 500**100
    right = 7**36 * 832**100 * 99**99
    require(left < right, "cleared integer comparison")
    require(F(left, right) == 256 * b, "algebra of integer comparison")
    n = 1 << 21
    # U^100 < 2^(-800-8N), hence U < 2^-167780.
    require(800 + 8 * n > 100 * 167780, "reported integer margin")
    return 167780


if __name__ == "__main__":
    rows = check_rows()
    margin = check_finite_bound()
    print(f"PASS: {rows} exact row checks; high-source-weight band failure < 2^-{margin}.")
    print("Scope: at most floor(N/100) missing outer coordinates; not a whole-code certificate.")
