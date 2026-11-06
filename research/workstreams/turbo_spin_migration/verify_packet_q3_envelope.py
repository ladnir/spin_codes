"""Outward interval check for fixed occupation three, all packet bands.

The three-packet kernel is evaluated by enumerating its four-state paths.
For a path with k live spacings, the live-time fraction has a Beta(k,4-k)
law. This gives explicit positive Laplace integrals and avoids interval
matrix exponentiation. Witness optimization is only a search heuristic.
"""

import argparse
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path

from mpmath import iv

from q3_packet_screen import screen as numeric_screen
from verify_packet_q1_interval import DELTA, M, P0, SCHEDULE, perron_2x2
from research.workstreams.paper_architecture.certificates.single_sampled_ba_rm2sub.linear_time_audit.verify_outer_interval import (
    HALF, add_upper, as_interval, interval, negative_entropy_upper,
    upper_float,
)


def kernels(gamma):
    g = as_interval(gamma)
    e = iv.exp(-g)
    i0 = (1-e)/g
    i1 = (1-(1+g)*e)/g**2
    i2 = (2-(2+2*g+g**2)*e)/g**3
    h = 2*(i0-i1)
    j = 2*i1
    r = iv.mpf(1)/M
    p = ((iv.mpf(0), iv.mpf(1)), (r, 1-r))
    weights = (iv.mpf(1), 3*(i0-2*i1+i2),
               6*(i1-i2), 3*i2, e)
    k3 = [[iv.mpf(0), iv.mpf(0)], [iv.mpf(0), iv.mpf(0)]]
    for path in itertools.product((0, 1), repeat=4):
        weight = weights[sum(path)]
        for i in range(3):
            weight *= p[path[i]][path[i+1]]
        k3[path[0]][path[-1]] += weight
    return (((iv.mpf(1), iv.mpf(0)), (iv.mpf(0), e)),
            ((iv.mpf(0), i0), (r*i0, (1-r)*e)),
            ((r*h, (1-r)*j), (r*(1-r)*j, r*j+(1-r)**2*e)),
            k3)


def bound(rows, witness):
    u = [as_interval(Fraction(format(math.exp(x), ".17g")))
         for x in witness[:3]]
    gamma = Fraction(format(math.exp(witness[3]), ".17g"))
    k = kernels(gamma)
    s1 = sum(u)
    s2 = u[0]*u[1]+u[0]*u[2]+u[1]*u[2]
    s3 = u[0]*u[1]*u[2]
    t = [[k[0][i][j]+s1*k[1][i][j]+s2*k[2][i][j]+s3*k[3][i][j]
          for j in range(2)] for i in range(2)]
    radius = perron_2x2(t[0][0], t[0][1], t[1][0], t[1][1])
    value = (iv.log(radius)
             +as_interval(gamma)*as_interval(DELTA)/as_interval(P0)) / 4
    other = [upper_float(3*iv.log(iv.mpf(2))/as_interval(SCHEDULE))]
    for row, ui in zip(rows, u):
        lo, hi = Fraction(row["rho_lo"]), Fraction(row["rho_hi"])
        value -= interval(lo, hi)*iv.log(ui)/4
        other.extend((row["upper_nats_per_bit"],
                      negative_entropy_upper(lo, hi)/4))
    return add_upper(upper_float(value), *other)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("envelope")
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--progress-every", type=int, default=500)
    args = parser.parse_args()
    with open(args.envelope, encoding="utf-8") as source:
        bands = json.load(source)["bands"]
    selected = list(range(0, len(bands), args.stride))
    if selected[-1] != len(bands)-1:
        selected.append(len(bands)-1)
    worst = -math.inf
    count = 0
    for indices in itertools.combinations_with_replacement(selected, 3):
        rows = [bands[i] for i in indices]
        numeric = numeric_screen(indices, bands)
        upper = bound(rows, numeric[3])
        worst = max(worst, upper)
        count += 1
        if upper >= 0:
            print(json.dumps({"result": "open", "indices": indices,
                              "upper": upper, "checked": count}))
            raise SystemExit(1)
        if args.progress_every and count % args.progress_every == 0:
            print(json.dumps({"checked": count, "worst_upper": worst}), flush=True)
    print(json.dumps({"result": "closed", "triples": count,
                      "worst_upper": worst}))


if __name__ == "__main__":
    main()
