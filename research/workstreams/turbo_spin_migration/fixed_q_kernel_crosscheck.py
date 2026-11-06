"""Independent high-precision fixed-Q kernel computation via Beta moments.

This is a diagnostic for the large-Q saddle. It enumerates the zero/live
path polynomial rather than exponentiating a block matrix.
"""

import argparse

import mpmath as mp

from fixed_q_transfer_screen import kernels as block_exponential_kernels


def beta_kernels(q, gamma, dps=70):
    mp.mp.dps = dps
    r = mp.mpf(1)/((1 << 22)-1)
    p = ((mp.mpf(0), mp.mpf(1)), (r, 1-r))
    out = []
    for a in range(q+1):
        matrix = mp.matrix(2)
        for start in (0, 1):
            paths = {(start, start): mp.mpf(1)}
            for _ in range(a):
                next_paths = {}
                for (state, k), weight in paths.items():
                    for target in (0, 1):
                        key = (target, k+target)
                        next_paths[key] = next_paths.get(key, 0)+weight*p[state][target]
                paths = next_paths
            for (end, k), weight in paths.items():
                moment = mp.hyper([k], [a+1], -gamma)
                matrix[start, end] += weight*moment
        out.append(matrix)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--q", type=int, default=16)
    parser.add_argument("--gamma", type=float, default=58.85484378041179)
    args = parser.parse_args()
    high = beta_kernels(args.q, args.gamma)
    fast = block_exponential_kernels(args.q, args.gamma)
    absolute = max(abs(float(high[a][i, j])-fast[a][i, j])
                   for a in range(args.q+1) for i in range(2) for j in range(2))
    relative = max(abs(float(high[a][i, j])-fast[a][i, j])
                   /max(1e-300, abs(float(high[a][i, j])))
                   for a in range(args.q+1) for i in range(2) for j in range(2))
    print({"q": args.q, "gamma": args.gamma,
           "max_absolute_difference": absolute,
           "max_relative_difference": relative,
           "Kq": [[str(high[-1][i, j]) for j in range(2)] for i in range(2)]})


if __name__ == "__main__":
    main()
