"""Exact binary weight-enumerator utilities.  Standard library only."""

from math import comb


def full_from_symmetric_half(half: dict[int, int], n: int = 256) -> list[int]:
    """Expand A_w=A_{n-w} from entries 0<=w<=n/2."""
    A = [0] * (n + 1)
    for w, value in half.items():
        if not (0 <= w <= n // 2):
            raise ValueError(f"half-spectrum weight out of range: {w}")
        A[w] = int(value)
        A[n - w] = int(value)
    return A


def krawtchouk(n: int, j: int, w: int) -> int:
    """Binary Krawtchouk K_j(w), computed exactly."""
    lo = max(0, j - (n - w))
    hi = min(j, w)
    return sum(
        (-1) ** s * comb(w, s) * comb(n - w, j - s)
        for s in range(lo, hi + 1)
    )


def krawtchouk_row(n: int, j: int) -> list[int]:
    return [krawtchouk(n, j, w) for w in range(n + 1)]


def macwilliams(A: list[int], k: int, n: int | None = None) -> list[int]:
    """Exact MacWilliams transform of a binary [n,k] linear code enumerator."""
    if n is None:
        n = len(A) - 1
    if len(A) != n + 1:
        raise ValueError("enumerator length does not match n")
    if sum(A) != 1 << k:
        raise ValueError(f"enumerator sum is {sum(A)}, expected 2^{k}")
    den = 1 << k
    nz = [(w, a) for w, a in enumerate(A) if a]
    B = []
    for j in range(n + 1):
        s = sum(a * krawtchouk(n, j, w) for w, a in nz)
        q, r = divmod(s, den)
        if r:
            raise ArithmeticError(f"nonintegral MacWilliams coefficient j={j}")
        B.append(q)
    return B


def dual_min_distance(B: list[int]) -> int | None:
    for j in range(1, len(B)):
        if B[j]:
            return j
    return None


def oa_moment(A: list[int], t: int) -> int:
    return sum(a * comb(w, t) for w, a in enumerate(A))


def expected_oa_moment(n: int, k: int, t: int) -> int:
    if t > k:
        # Formula remains an integer only when the OA strength hypothesis supports it.
        raise ValueError("use only for t within the claimed OA strength")
    return (1 << (k - t)) * comb(n, t)


def symmetric_half_linear_coeff(n: int, w: int, f) -> int:
    """Coefficient of half-spectrum variable x_w in sum_i x_i f(i)."""
    if w == n - w:
        return f(w)
    return f(w) + f(n - w)


def symmetric_kraw_coeff(n: int, j: int, w: int) -> int:
    return symmetric_half_linear_coeff(n, w, lambda x: krawtchouk(n, j, x))


def symmetric_moment_coeff(n: int, t: int, w: int) -> int:
    return symmetric_half_linear_coeff(n, w, lambda x: comb(x, t))
