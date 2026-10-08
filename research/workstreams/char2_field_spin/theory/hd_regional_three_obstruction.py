"""Exact distance obstructions for coordinate-major one-route scan codes.

This is a proof helper, not an encoder benchmark or a distance certificate.
All finite inequalities and probability bounds use integers or Fractions.
"""

from fractions import Fraction
from math import comb, log2


def check_parameters(b, rows, pre_order, post_order):
    if b < 2 or b % 2 or rows < 2:
        raise ValueError("Need even outer length and at least two outer rows")
    if not 0 <= pre_order < rows or not 1 <= post_order < rows:
        raise ValueError("Invalid inverse bandwidth")


def cyclic_ceiling(b, rows, pre_order=2, post_order=1):
    """All routes: select a cyclic interval of outer parity-check groups."""
    check_parameters(b, rows, pre_order, post_order)
    p, q = pre_order, post_order
    a = (2 * q * rows + 2 * (q + 1) * p) // (2 * q + 1) + 1
    if a > rows:
        return None
    capacity = (b // 2) * a
    expected_lower = b * ((q + 1) * (a - p) - q * rows)
    assert expected_lower > capacity
    return dict(a=a, capacity=capacity, confined_lower=expected_lower,
                distance_upper=capacity + 1)


def regional_prefix_ceiling(b, rows, pre_order=2, post_order=1):
    """Every regional route: select a fixed prefix of outer check groups."""
    check_parameters(b, rows, pre_order, post_order)
    p, q = pre_order, post_order
    a = (2 * q * rows + 2 * (q + 1) * p + 2 * q) // (2 * q + 1) + 1
    if a > rows:
        return None
    capacity = (b // 2) * a
    safe = a - p
    confined = b * max(0, (q + 1) * safe - q * rows - q)
    assert confined > capacity
    return dict(a=a, safe=safe, capacity=capacity, confined_lower=confined,
                distance_upper=capacity + 1)


def window_count(bits, order):
    """Number of wholly safe within-region windows of length order+1."""
    return sum(all(bits[i:i + order + 1])
               for i in range(len(bits) - order))


def regional_failure_lower(b, rows, pre_order=2, post_order=1,
                           target=Fraction(2, 5)):
    """Lower bound on failing the target distance under uniform regional P.

    The coefficients may be fixed arbitrarily. Independence between regions
    is not needed for this first-moment LOWER bound; uniform marginals suffice.
    """
    check_parameters(b, rows, pre_order, post_order)
    if not 0 < target <= Fraction(1, 2):
        raise ValueError("Expected a target in (0, 1/2]")
    p, q, rho = pre_order, post_order, b // 2
    n = b * rows
    bad_max = (n * target.numerator - 1) // target.denominator
    a = (bad_max - 1) // rho
    safe = a - p
    if safe < q + 1:
        raise ValueError("Target is too small for this sufficient event")
    capacity = rho * a
    expected = b * (rows - q) * Fraction(comb(safe, q + 1),
                                        comb(rows, q + 1))
    maximum = b * (safe - q)
    lower = (max(Fraction(0), (expected - capacity) / (maximum - capacity))
             if maximum > capacity else Fraction(0))
    return dict(a=a, safe=safe, capacity=capacity, bad_max=bad_max,
                distance_upper_on_event=capacity + 1, expected=expected,
                maximum=maximum, failure_lower=lower)


def source_group_supports(b, rows, pre_order, post_order, inverse_route):
    """Conservative H-column group supports, including global carry.

    inverse_route[j] is the pre-route coordinate reached from routed index j.
    No modular wrap is used for the actual scan indices. Only outer row
    GROUP labels are reduced modulo rows.
    """
    n = b * rows
    if sorted(inverse_route) != list(range(n)):
        raise ValueError("The inverse route must be a permutation")
    result = []
    for column in range(n):
        groups = set()
        for routed in range(column, min(n, column + post_order + 1)):
            middle = inverse_route[routed]
            for source in range(middle, min(n, middle + pre_order + 1)):
                groups.add(source % rows)
        result.append(groups)
    return result


def block_placement_failure_lower(b, rows, block_size, pre_order=2,
                                  target=Fraction(2, 5)):
    """Scalar-post gate when each placement preserves common row blocks.

    Each coordinate may permute the blocks and their internal orders
    arbitrarily, but it may not mix different source blocks within a block.
    The intervening route has uniform regional marginals.
    """
    check_parameters(b, rows, pre_order, 1)
    if rows % block_size or block_size <= pre_order:
        raise ValueError("Need whole blocks wider than the pre bandwidth")
    rho = b // 2
    bad_max = (b * rows * target.numerator - 1) // target.denominator
    blocks = (bad_max - 1) // (rho * block_size)
    a = blocks * block_size
    safe = blocks * (block_size - pre_order)
    capacity = rho * a
    expected = Fraction(b * safe * (safe - 1), rows)
    maximum = b * max(0, safe - 1)
    lower = (max(Fraction(0), (expected - capacity) / (maximum - capacity))
             if maximum > capacity else Fraction(0))
    return dict(blocks=blocks, a=a, safe=safe, capacity=capacity,
                distance_upper_on_event=capacity + 1, expected=expected,
                failure_lower=lower)


def block_degree_ceiling(b, rows, block_size, pre_order=2, post_order=1):
    """Deterministic grouped-check ceiling for common-block placements."""
    check_parameters(b, rows, pre_order, post_order)
    if block_size < 1 or rows % block_size:
        raise ValueError("Need a partition into whole common source blocks")
    groups = rows // block_size
    checks_per_group = (b // 2) * block_size
    pre_groups = (pre_order + block_size - 1) // block_size + 1
    degree = pre_groups * (post_order + 1)
    if degree > groups:
        return None
    denominator = comb(groups - 1, degree - 1)
    low, high = degree, groups
    while low < high:
        middle = (low + high) // 2
        if 2 * comb(middle - 1, degree - 1) > denominator:
            high = middle
        else:
            low = middle + 1
    return dict(groups=groups, checks_per_group=checks_per_group, degree=degree,
                selected_groups=low, distance_upper=checks_per_group * low + 1)


def main():
    b, rows, p = 64, 32768, 2
    print("K20; [64,32]; coordinate-major; global scan carry; pre order 2")
    for q in range(1, 5):
        cyclic = cyclic_ceiling(b, rows, p, q)
        regional = regional_prefix_ceiling(b, rows, p, q)
        failure = regional_failure_lower(b, rows, p, q)
        lower = failure["failure_lower"]
        print(dict(post_order=q, cyclic=cyclic, fixed_prefix=regional,
                   distance_on_random_event=failure["distance_upper_on_event"],
                   expected=float(failure["expected"]),
                   failure_lower=float(lower),
                   reliability_bits=(-log2(float(lower)) if lower else None)))
    for block in (2, 4, 8, 16, 32, 64, 128):
        print(dict(block_size=block, grouped=block_degree_ceiling(b, rows, block, p)))
        if block > p:
            result = block_placement_failure_lower(b, rows, block, p)
            print(dict(block_size=block, a=result["a"], safe=result["safe"],
                       distance_upper_on_event=result["distance_upper_on_event"],
                       failure_lower=float(result["failure_lower"])))


if __name__ == "__main__":
    main()
