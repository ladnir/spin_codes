"""Exact finite obstruction for a rational scalar postscan and regional routing.

The written rank and Doob arguments are in HD_RATIONAL_POST_OBSTRUCTION.md.
All certificate arithmetic here uses integers/Fraction; floating displays
are diagnostic only. No Monte Carlo, field-size approximation, or optimizer.
"""
import argparse
from fractions import Fraction as F
import json


def selected_runs(bits):
    """Maximal half-open one-runs of length at least two, in one region."""
    if any(x not in (0, 1) for x in bits):
        raise ValueError("bits must be zero or one")
    result = []
    start = 0
    while start < len(bits):
        if not bits[start]:
            start += 1
            continue
        end = start + 1
        while end < len(bits) and bits[end]:
            end += 1
        if end-start >= 2:
            result.append((start, end))
        start = end
    return result


def run_statistics(bits):
    """Return T=sum(length-1), r=number of selected runs, and |S|."""
    runs = selected_runs(bits)
    total = sum(end-start for start, end in runs)
    return total-len(runs), len(runs), total


def exact_means(regions, rows, safe):
    """Independent uniform safe-subsets of fixed size in each region."""
    if regions < 1 or rows < 2 or not 0 <= safe <= rows:
        raise ValueError("invalid regional geometry")
    adjacent = F(regions*safe*(safe-1), rows)
    runs = F(regions*(rows-safe+1)*safe*(safe-1), rows*(rows-1))
    return adjacent, runs


def finite_bound(b=64, k=32, rows=32768, a=None, pre_order=2,
                 distance=F(2, 5)):
    """A whole-code FAILURE lower bound, not a distance certificate.

    The precoder inverse has lower bandwidth pre_order. The first a outer
    row-groups give h=k*a check rows and a-pre_order safe columns/region.
    The two-statistic Doob proof supplies Var(T),Var(r) <=4N independently
    of coefficient choices. The routing law must be independent uniform
    scalar permutations inside the b regions.
    """
    if a is None:
        a = 13*rows//25
    if b < 2 or b != 2*k or rows < 2 or pre_order < 0:
        raise ValueError("invalid code geometry")
    if not pre_order <= a <= rows or not 0 < distance < 1:
        raise ValueError("invalid cut")
    n, h, safe = b*rows, k*a, a-pre_order
    # Forbidden weights are through ceil(distance*N)-1, including when
    # distance*N is an integer.
    target = distance*n
    cutoff = (target.numerator+target.denominator-1)//target.denominator-1
    run_cap = cutoff-h-1
    mean_t, mean_r = exact_means(b, rows, safe)
    gap_t, gap_r = mean_t-h, run_cap-mean_r
    if gap_t <= 0 or gap_r <= 0:
        raise ValueError("this cut has no positive two-sided expectation gaps")
    variance = 4*n
    bad_t, bad_r = F(variance)/(gap_t*gap_t), F(variance)/(gap_r*gap_r)
    lower = max(F(0), 1-bad_t-bad_r)
    return dict(b=b, k=k, rows=rows, n=n, message_dimension=k*rows,
                pre_order=pre_order, a=a, safe_per_region=safe,
                check_capacity=h, forbidden_cutoff=cutoff, run_cap=run_cap,
                mean_t=mean_t, mean_runs=mean_r, gap_t=gap_t, gap_runs=gap_r,
                variance_upper=variance, bad_t_upper=bad_t, bad_runs_upper=bad_r,
                failure_lower=lower, exceeds_19_over_20=lower > F(19, 20),
                route="independent uniform scalar permutations within regions",
                coefficient_law="arbitrary valid nonzero coefficients, even route-dependent",
                rigorous_arithmetic=True, bound_type="whole-code failure lower bound")


def json_record(result):
    return {name: ({"exact": str(value), "display": float(value)}
                   if isinstance(value, F) else value)
            for name, value in result.items()}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rows", type=int, default=32768)
    parser.add_argument("--outer", type=int, nargs=2, default=(64, 32))
    parser.add_argument("--a", type=int)
    parser.add_argument("--pre-order", type=int, default=2)
    parser.add_argument("--distance", type=F, default=F(2, 5))
    args = parser.parse_args()
    print(json.dumps(json_record(finite_bound(*args.outer, args.rows, args.a,
                     args.pre_order, args.distance)), indent=2))
