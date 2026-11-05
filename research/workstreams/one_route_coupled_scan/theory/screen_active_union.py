#!/usr/bin/env python3
"""Bounded floating diagnostic of complete (v,c) excursion classes.

The formula is in ACTIVE_UNION_GATE.md. Each reported score sums all marked
layouts and all output/state masks at the specified active-batch count v and
excursion count c. It does NOT cover other (v,c), and floating arithmetic is
not a certificate. The regional occupancy vector is relaxed to its mean;
this can destroy useful route information. Positive scores are failed bounds,
not codewords. No output files or third-party modules. Root owns execution.

Hard caps: 30 seconds, 256 MiB externally, at most 16 v values, 10 c values
per v, 13 fractional markers. No optimization or unbounded parameter search.
"""

from math import exp, expm1, isfinite, lgamma, log, log1p
from time import monotonic


B, K_LOCAL, PACKET, BATCH = 64, 32, 4, 128
N_BATCHES, FIELD_BITS = 16384, 32
N = N_BATCHES * BATCH
GROUPS = N // (PACKET * B)
REGION_BATCHES = N_BATCHES // B
CUTOFF = 35 * N // 100
LN2 = log(2)
Q_HALF = 1 << (FIELD_BITS // 2)
LOG_TAUS = (1 / 256, 1 / 64, 1 / 16, 0.125, 0.25, 0.5,
            1, 2, 4, 8, 16, 24, 32)
V_GRID = (33, 64, 128, 256, 512, 1024, 2048, 4096,
          5734, 6000, 8192, 10000, 12000, 14000, 16000, 16384)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def log_choose(n, k):
    if not 0 <= k <= n:
        return -float("inf")
    return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)


def log_add(values):
    values = tuple(values)
    largest = max(values, default=-float("inf"))
    if not isfinite(largest):
        return largest
    return largest + log(sum(exp(x - largest) for x in values))


def log_expm1(x):
    require(x > 0, "positive expm1 argument required")
    if x > 40:
        return x + log1p(-exp(-x))
    return log(expm1(x))


def log_one_plus_exp(x):
    if x > 40:
        return x + log1p(exp(-x))
    return log1p(exp(x))


def log_guarded_moment(v, log_tau):
    """Positive H-1 sum, followed by exact algebra before floating evaluation."""
    rho = v / N_BATCHES
    log_t = PACKET * log_tau * LN2
    if v == N_BATCHES:
        log_h_minus = log_expm1((B - K_LOCAL) * log_t)
    else:
        log_h_minus = log_add(
            log_choose(B, j) + j * log(rho) + (B-j) * log1p(-rho)
            + log_expm1((j-K_LOCAL) * log_t)
            for j in range(K_LOCAL + 1, B + 1))
    log_h = log_one_plus_exp(log_h_minus)
    # H>=1, so subtracting one occurs before any logarithmic evaluation.
    log_power_minus = log_expm1(GROUPS * log_h)
    return log_t - log_expm1(log_t) + log_power_minus


def log_weight_volume(u, cutoff):
    """Upper bound on sum_{w<=cutoff} binom(u,w) Q_HALF^w.

    Include weights below the necessary reset minimum: safe overcount.
    Use the full binomial sum above its mode; otherwise use the positive
    Chernoff bound with the minimizing output marker.
    """
    e = min(u, cutoff)
    if e == 0:
        return 0.0
    if e * (1 + Q_HALF) >= Q_HALF * u:
        return u * log1p(Q_HALF)
    fraction = e / u
    entropy = -fraction * log(fraction) - (1-fraction) * log1p(-fraction)
    return u * entropy + e * log(Q_HALF)


def touched_cap(v, c, region_batches=REGION_BATCHES):
    return min(B, 2*c + (v-2*c) // region_batches)


def score(v, c, log_tau):
    require(1 <= c and 2*c-1 <= v <= N_BATCHES, "invalid layout class")
    if (BATCH+1)*(c-1)+1 > min(CUTOFF, BATCH*v):
        return -float("inf")
    if touched_cap(v, c) <= K_LOCAL:
        return -float("inf")
    active_after_start = BATCH * (v-c)
    log_layouts = log_choose(v-c, c-1) + log_choose(N_BATCHES-v+c, c)
    log_amplitudes = (active_after_start * LN2
        + (BATCH-active_after_start) * FIELD_BITS * LN2 / 2
        + log_weight_volume(BATCH*v, CUTOFF)
        - log((1 << FIELD_BITS)-1))
    nu = log_tau / FIELD_BITS
    return (log_layouts + nu*log_amplitudes + log_guarded_moment(v, log_tau)) / LN2


def exact_small_layout_checks():
    """Enumerate interval boundaries, not field elements; arithmetic is integer."""
    from math import comb
    checked = 0
    for n in range(1, 9):
        actual = {}
        # Ordered, non-overlapping intervals may abut. All but the last
        # have length >=2; the last may have length1 (free terminal).
        def visit(end, lengths):
            if lengths:
                v, c = sum(lengths), len(lengths)
                actual[v, c] = actual.get((v, c), 0) + 1
            if lengths and lengths[-1] < 2:
                return
            for start in range(end, n):
                for stop in range(start+1, n+1):
                    visit(stop, lengths + (stop-start,))
        visit(0, ())
        for v in range(1, n+1):
            for c in range(1, (v+1)//2+1):
                expected = comb(v-c, c-1) * comb(n-v+c, c)
                require(actual.get((v, c), 0) == expected, "layout formula")
                checked += 1
    # Every interval's exact region count is bounded, including singletons
    # and starts immediately before a region boundary.
    for length in range(1, 17):
        for region_length in range(1, 9):
            cap = 2 + (length-2) // region_length
            for start in range(region_length):
                actual = (start+length-1)//region_length - start//region_length + 1
                require(actual <= cap, "interval touched-region cap")
                checked += 1
    return checked


def main():
    deadline = monotonic() + 30
    checks = exact_small_layout_checks()
    print(f"Exact local checks: {checks}; floating diagnostic follows, NOT a certificate.")
    print("v,c,touched_cap,best_log2_bound,log2_tau")
    classes = excluded = 0
    for v in V_GRID:
        largest_c = min((v+1)//2, 1 + (CUTOFF-1)//(BATCH+1))
        c_grid = sorted({1, 2, 8, 16, 17, max(1, v//100),
                         max(1, v//10), max(1, v//3), largest_c})
        for c in c_grid:
            if c > largest_c:
                continue
            if touched_cap(v, c) <= K_LOCAL:
                excluded += 1
                continue
            candidates = [(score(v, c, marker), marker) for marker in LOG_TAUS]
            best, marker = min(candidates)
            print(f"{v},{c},{touched_cap(v,c)},{best:.6f},{marker:g}")
            classes += 1
            if monotonic() > deadline:
                raise TimeoutError("30-second cap; diagnostic incomplete")
    print(f"Screened {classes} complete (v,c) classes; {excluded} others excluded by touched regions.")
    print("Other classes are not covered. Positive bounds do not exhibit a failing codeword.")


if __name__ == "__main__":
    main()
