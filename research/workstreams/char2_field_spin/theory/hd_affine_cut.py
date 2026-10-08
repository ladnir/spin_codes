"""Bounded row-set cut diagnostics for the frozen affine placement.

Actual safe counts are exact integers. The interval-pressure screen uses
stationary floating arithmetic and is NOT a distance certificate or an
impossibility result. Only the separately requested exact-moment path can
produce a finite failure lower bound.
"""
import argparse
from functools import lru_cache
from fractions import Fraction as F
import json

import numpy as np

from hd_two_state_regional_strides import ROWS, STRIDES, LAMBDA_UPPER


def safe_counts(selected, strides=STRIDES, offsets=None):
    """Actual pre1 inverse-column support, including global carry/end.

    pi_c(j)=stride_c*j+offset_c. Interior supports are adjacent pairs.
    At a regional boundary the next endpoint is pi_(c+1)(0); the last
    global inverse column has only its first endpoint.
    """
    selected = np.asarray(selected, dtype=bool)
    rows, regions = len(selected), len(strides)
    if offsets is None:
        offsets = (0,)*regions
    counts = []
    for c, (stride, offset) in enumerate(zip(strides, offsets)):
        cyclic = int(np.count_nonzero(selected & np.roll(selected, -stride)))
        last, first = (offset-stride) % rows, offset % rows
        cyclic -= int(selected[last] and selected[first])
        if c+1 < regions:
            cyclic += int(selected[last] and selected[offsets[c+1] % rows])
        else:
            cyclic += int(selected[last])
        counts.append(cyclic)
    return tuple(counts)


def spectral_safe_upper(a):
    """Directed spectral receipt, with only one all-zero-offset terminal."""
    alpha = F(a, ROWS)
    bound = len(STRIDES)*ROWS*(alpha*alpha+LAMBDA_UPPER*alpha*(1-alpha))
    return bound.numerator//bound.denominator+1


def grouped_safe_counts(selected, group, strides=STRIDES):
    """Zero-offset affine permutations of common groups; partial groups allowed."""
    selected = np.asarray(selected, dtype=bool)
    if len(selected) % group:
        raise ValueError("group must divide row count")
    blocks = selected.reshape((-1, group))
    within = int(np.count_nonzero(blocks[:, :-1] & blocks[:, 1:]))
    first, last = blocks[:, 0], blocks[:, -1]
    counts = [within+int(np.count_nonzero(last & np.roll(first, -stride)))
              for stride in strides]
    counts[-1] += int(last[-strides[-1] % len(blocks)] and not first[0])
    return tuple(counts)


@lru_cache(maxsize=None)
def stationary_pressure(safe, numerator, denominator):
    """Mean optimal score per position for G-(numerator/denominator)*C.

    The two-tail interval DP has clipped score-difference state. For one
    Bernoulli region this routine solves its finite stationary chain.
    The safe-count/R density is used only for this labeled numerical screen.
    """
    if safe == 0:
        return 0.0
    if safe == ROWS:
        return 1.0
    p = safe/ROWS
    opening = 2*(numerator+denominator)
    values, indices, pending = [], {}, [-opening]
    while pending:
        d = pending.pop()
        if d in indices:
            continue
        indices[d] = len(values)
        values.append(d)
        clipped = min(0, max(-opening, d))
        pending.extend((clipped+denominator, clipped-numerator))
    n = len(values)
    matrix = -np.eye(n)
    for i, d in enumerate(values):
        clipped = min(0, max(-opening, d))
        matrix[indices[clipped+denominator], i] += p
        matrix[indices[clipped-numerator], i] += 1-p
    matrix[-1] = 1
    right = np.zeros(n)
    right[-1] = 1
    stationary = np.linalg.solve(matrix, right)
    return float(np.dot(stationary, np.maximum(0, values))/denominator)


MARKS = ((1, 8), (1, 6), (1, 4), (1, 3), (1, 2), (2, 3),
         (3, 4), (1, 1), (4, 3), (3, 2), (2, 1), (3, 1), (4, 1))


def screen_profile(a, counts, distance=F(2, 5), extra_marks=()):
    regions, n = len(counts), len(counts)*ROWS
    h = regions*a/2
    total = sum(counts)
    record = dict(a=a, density=a/ROWS, safe_total=total,
                  safe_min=min(counts), safe_max=max(counts), check_capacity=int(h),
                  safe_surplus=total-int(h), counts=list(counts),
                  finite_failure_certificate=False)
    if total <= h:
        return dict(record, status="safe-count prerequisite fails exactly")
    cap = float(distance)-h/n
    if cap <= 0:
        return dict(record, status="check-capacity prerequisite fails exactly")
    candidates = []
    for num, den in (*MARKS, *extra_marks):
        pressure = sum(stationary_pressure(safe, num, den) for safe in counts)/regions
        candidates.append((pressure+num/den*cap-h/n, num, den))
    gap, num, den = min(candidates)
    return dict(record, status="stationary diagnostic only", dual_gap=gap,
                mark=f"{num}/{den}", mean_safe_density=total/n)


def hillclimb(selected, max_swaps=512):
    """Bounded deterministic swap ascent of exact internal-edge count."""
    selected = np.array(selected, dtype=bool, copy=True)
    neighbors = np.array([*STRIDES, *(-s for s in STRIDES)], dtype=np.int64)
    degree = sum((np.roll(selected, int(shift)).astype(np.int16)
                  for shift in neighbors), np.zeros(ROWS, dtype=np.int16))
    moves = 0
    for _ in range(max_swaps):
        remove = int(np.argmin(np.where(selected, degree, 1000)))
        add = int(np.argmax(np.where(~selected, degree, -1000)))
        adjacent = int(np.any((remove+neighbors) % ROWS == add))
        if int(degree[add])-int(degree[remove])-adjacent <= 0:
            break
        selected[remove], selected[add] = False, True
        degree[(remove+neighbors) % ROWS] -= 1
        degree[(add+neighbors) % ROWS] += 1
        moves += 1
    return selected, moves


def candidate_sets():
    """Fixed bounded menu, not an exhaustive search over source-row sets."""
    x = np.arange(ROWS)
    sizes = [int(ROWS*v/40) for v in range(11, 32)]
    for a in sizes:
        yield f"cyclic-{a}", x < a
    for modulus in (2, 4, 8, 16, 32, 64):
        for take in range(1, modulus):
            if modulus > 8 and take not in (modulus//4, 3*modulus//8, modulus//2,
                                            5*modulus//8, 3*modulus//4):
                continue
            yield f"lowbit-{modulus}-{take}", x % modulus < take
    # Largest positive modes from the frozen spectral receipt, plus simple
    # low frequencies. Stable sorting resolves every cosine tie explicitly.
    kernel = np.zeros(ROWS)
    for stride in STRIDES:
        kernel[stride] += 1
        kernel[-stride] += 1
    eigenvalues = np.fft.rfft(kernel).real/(2*len(STRIDES))
    modes = [1, 2, 4, *map(int, np.argsort(-eigenvalues[1:])[:5]+1)]
    for mode in dict.fromkeys(modes):
        order = np.argsort(-np.cos(2*np.pi*mode*x/ROWS), kind="stable")
        for fraction in (0.3, 0.4, 0.5, 0.6, 0.7):
            a = int(fraction*ROWS)
            selected = np.zeros(ROWS, dtype=bool)
            selected[order[:a]] = True
            yield f"fourier-{mode}-{a}", selected
    for fraction in (0.3, 0.4, 0.5, 0.6, 0.7):
        a = int(fraction*ROWS)
        order = np.argsort(-np.cos(2*np.pi*1082*x/ROWS), kind="stable")
        selected = np.zeros(ROWS, dtype=bool)
        selected[order[:a]] = True
        selected, moves = hillclimb(selected)
        yield f"hillclimb-1082-{a}-{moves}", selected


def grouped_candidate_sets(group):
    """Bounded full-block and partial-block row sets for a locality variant."""
    size = ROWS//group
    x = np.arange(size)
    fractions = (.025, .05, .075, .1, .125, .15, .175, .2, .25,
                 .3, .35, .4, .45, .5, .55, .6, .65, .7, .75)
    for fraction in fractions:
        selected = x < int(fraction*size)
        yield f"group{group}-cyclic-{int(selected.sum())}", np.repeat(selected, group)
    kernel = np.zeros(size)
    for stride in STRIDES:
        kernel[stride % size] += 1
        kernel[-stride % size] += 1
    eigenvalues = np.fft.rfft(kernel).real/(2*len(STRIDES))
    modes = np.argsort(-eigenvalues[1:])[:5]+1
    for mode in modes:
        order = np.argsort(-np.cos(2*np.pi*mode*x/size), kind="stable")
        for fraction in (.05, .1, .2, .3, .4, .5, .6, .7):
            selected = np.zeros(size, dtype=bool)
            selected[order[:int(fraction*size)]] = True
            yield f"group{group}-fourier-{mode}-{int(selected.sum())}", np.repeat(selected, group)
    for pattern in range(1, (1 << group)-1):
        lane = np.array([(pattern >> i) & 1 for i in range(group)], dtype=bool)
        for fraction in (.1, .25, .5, .75, 1):
            selected = (x < int(fraction*size))[:, None] & lane[None, :]
            yield f"group{group}-partial-{pattern}-{int(fraction*size)}", selected.ravel()


def exact_profile_moments(a, counts, block=256):
    """Optional exact finite G/C cut using the fixed 3G-C block policy.

    Counts may differ across regions. The imported exact pair moments
    retain dependence between fixed-cardinality blocks in each region.
    This path is deliberately separate from the floating stationary screen.
    """
    from hd_two_state_cut import regional_moments, word_moments
    sums = word_moments(block)
    cached = {}
    mean_g = var_g = mean_c = var_c = F(0)
    for safe in counts:
        if safe not in cached:
            cached[safe] = (regional_moments(ROWS, safe, block, sums, 1),
                            regional_moments(ROWS, safe, block, sums, 2))
        (g, vg), (c, vc) = cached[safe]
        mean_g += g
        var_g += vg
        mean_c += c
        var_c += vc
    n, h = len(counts)*ROWS, len(counts)*a//2
    cutoff = (2*n+4)//5-1
    gap_g, gap_c = mean_g-h, cutoff-h-1-mean_c
    lower = (max(F(0), 1-var_g/gap_g**2-var_c/gap_c**2)
             if gap_g > 0 and gap_c > 0 else F(0))
    return dict(mean_gain=mean_g, variance_gain=var_g, mean_cost=mean_c,
                variance_cost=var_c, gap_gain=gap_g, gap_cost=gap_c,
                failure_lower=lower, finite_failure_certificate=lower > 0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--details", action="store_true")
    parser.add_argument("--group", type=int, choices=(1, 2, 4), default=1)
    args = parser.parse_args()
    records = []
    candidates = candidate_sets() if args.group == 1 else grouped_candidate_sets(args.group)
    for name, selected in candidates:
        a = int(selected.sum())
        counts = safe_counts(selected) if args.group == 1 else grouped_safe_counts(selected, args.group)
        # The 1/12 mark is the retained refinement of the three initially
        # positive coarse-grid grouped diagnostics, not an optimizer run.
        extra_marks = () if args.group == 1 else ((1, max(1, round(ROWS/(0.6*a)))), (1, 12))
        record = dict(name=name, **screen_profile(a, counts, extra_marks=extra_marks))
        records.append(record)
        if args.details:
            print(json.dumps(record), flush=True)
    screened = [r for r in records if "dual_gap" in r]
    print(json.dumps(dict(group=args.group, candidates=len(records),
                         safe_prerequisite_failures=sum(r["safe_surplus"] <= 0 for r in records),
                         largest_dual_gaps=sorted(screened, key=lambda r:r["dual_gap"], reverse=True)[:10],
                         rigorous=False, finite_failure_certificate=False), indent=2))
