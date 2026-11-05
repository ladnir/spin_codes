"""Bounded floating-point diagnostic for SPARSE_GATE.md, not a certificate.

Width four, exact single-source insertion distribution and free terminal
state. Evaluates the relaxed lookahead operator for the full-support
one-group gate. Positive matrix powers are scaled to avoid overflow.
Underflow and rounding are NOT enclosed, so no success claim is certified.

Root owns numerical runs. Uses one BLAS thread, <128 MiB working arrays,
at most 16 parameter pairs per invocation, and writes no files. Use -B.
"""

import os

for _variable in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[_variable] = "1"

import argparse
import math
import time

import numpy as np

from coupled_rank_transfer import (
    lookahead_transition_terms, pair_antipodal_cycle, symbolic_transitions,
)


W, MASKS, STATES = 4, 16, 256


def symbolic_operators():
    cycle = pair_antipodal_cycle(W)
    operators, terminals = [], []
    for inputs in (0, 1, 2, 4, 8):
        entries, ending = [], []
        for h in range(MASKS):
            ordinary = symbolic_transitions(W, inputs, h, cycle)
            for h1 in range(MASKS):
                row = h * MASKS + h1
                for weight, charge, multiplicity in ordinary.get(h1, ()):
                    ending.append((row, weight, charge, multiplicity))
                for h2 in range(MASKS):
                    for weight, charge, multiplicity in lookahead_transition_terms(W, inputs, h, h1, h2, cycle):
                        entries.append((row, h1 * MASKS + h2, weight, charge, multiplicity))
        operators.append(np.asarray(entries, dtype=np.int64).T)
        terminals.append(np.asarray(ending, dtype=np.int64).T)
    return operators, terminals


def evaluate(symbols, z, tau):
    operators, terminals = symbols
    matrices, vectors = [], []
    for entries, ending in zip(operators, terminals):
        matrix, vector = np.zeros((STATES, STATES)), np.zeros(STATES)
        row, column, weight, charge, count = entries
        np.add.at(matrix, (row, column), count * np.power(z, weight) * np.power(tau, -charge.astype(float)))
        row, weight, charge, count = ending
        np.add.at(vector, row, count * np.power(z, weight) * np.power(tau, -charge.astype(float)))
        matrices.append(matrix)
        vectors.append(vector)
    return matrices[0], sum(matrices[1:]) / 4, vectors[0], sum(vectors[1:]) / 4


def normalize_pair(p, s, logscale):
    scale = max(float(np.max(p)), float(np.max(s)))
    if not math.isfinite(scale) or scale <= 0:
        raise ArithmeticError("invalid positive matrix scale")
    return p / scale, s / scale, logscale + math.log2(scale)


def compose(left, right):
    """(E^a, insertion_a) times (E^b, insertion_b), common scaling."""
    p, s, logscale = left
    q, t, other_scale = right
    return normalize_pair(p @ q, s @ q + p @ t, logscale + other_scale)


def insertion_power(e, a, exponent):
    result = (np.eye(STATES), np.zeros_like(e), 0.0)
    factor = normalize_pair(e, a, 0.0)
    while exponent:
        if exponent & 1:
            result = compose(result, factor)
        exponent >>= 1
        if exponent:
            factor = compose(factor, factor)
    return result


def screen(symbols, z, tau, steps=2048, regions=256, groups=8192, dimension=128):
    e, a, f0, f1 = evaluate(symbols, z, tau)
    p, s, logscale = insertion_power(e, a, steps - 1)
    # B1=(S_(m-1) E + E^(m-1) A)/m. Last region omits lookahead
    # in the last step: b1=(S_(m-1) f0+E^(m-1) f1)/m.
    b = (s @ e + p @ a) / steps
    terminal = (s @ f0 + p @ f1) / steps
    vector = np.zeros(STATES)
    vector[:MASKS] = 1
    total_scale = 0.0
    for _ in range(regions - 1):
        vector = vector @ b
        scale = float(np.max(vector))
        if scale == 0:
            return float("-inf"), float("nan")
        vector /= scale
        total_scale += logscale + math.log2(scale)
    final = float(vector @ terminal)
    if final <= 0:
        return float("-inf"), float("nan")
    n = steps * W * regions
    threshold = (7 * n) // 20
    logbound = (math.log2(groups) + (dimension - 1) * math.log2(tau)
                - threshold * math.log2(z) + total_scale + logscale + math.log2(final))
    live = [h * MASKS + h1 for h in range(1, MASKS) for h1 in range(1, MASKS)]
    # A diagnostic eigenvalue only; not a componentwise verified potential.
    eigenvalue = float(np.max(np.abs(np.linalg.eigvals(e[np.ix_(live, live)]))))
    return logbound, eigenvalue


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--z", type=float, nargs="+", default=[0.99, 0.994, 0.996, 0.998, 0.999])
    parser.add_argument("--tau", type=float, nargs="+", default=[256.0])
    args = parser.parse_args()
    if len(args.z) * len(args.tau) > 16 or not all(0 < z < 1 for z in args.z) or not all(1 <= t <= 256 for t in args.tau):
        parser.error("need at most 16 pairs, 0<z<1, and 1<=tau<=256")
    symbols = symbolic_operators()
    print("FLOATING DIAGNOSTIC ONLY: full-support one-group gate, width4 antipodal, free terminal", flush=True)
    for tau in args.tau:
        for z in args.z:
            started = time.monotonic()
            value, eigenvalue = screen(symbols, z, tau)
            print(f"tau={tau:g} z={z:.9g} log2_U={value:.6f} live_rho={eigenvalue:.9g} seconds={time.monotonic()-started:.3f}", flush=True)


if __name__ == "__main__":
    main()
