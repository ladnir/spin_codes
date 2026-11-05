"""Propose and EXACTLY check a positive live empty-gap envelope at width4.

Float iteration only proposes a vector. A successful result requires every
inequality T(v)<=lambda*v to pass using Fraction arithmetic after rounding
the proposed vector to rationals. No file output or full-code certificate.
Bounded to16384 iterations and small state tables. Run serially with -B.
"""

from fractions import Fraction as F
from math import ceil, log2

from coupled_rank_transfer import pair_antipodal_cycle
from empty_span_rank import empty_step, terminal_ones


def check(z=F(999, 1000), contraction=F(499, 500), *, tau=F(256), verbose=True):
    width = 4
    tau = F(tau)
    # z need not be below lambda: a step emits up to four symbols.
    if not 0 < z < 1 or not 0 < contraction < 1 or tau < 1:
        raise ValueError("need output marker and contraction in (0,1)")
    cycle = pair_antipodal_cycle(width)
    valid = [(rank, mask) for rank in range(1, width + 1)
             for mask in range(1, 1 << width) if rank <= mask.bit_count()]
    vector = terminal_ones(width, 1.0)
    vector[0][0] = 0.0
    for iteration in range(1, 16385):
        transformed = empty_step(vector, width, float(z), float(tau), cycle=cycle)
        if iteration % 16 == 0 and all(
                transformed[r][h] <= float(contraction) * vector[r][h] - 1e-4
                for r, h in valid):
            break
        for r, h in valid:
            vector[r][h] = 1 + transformed[r][h] / float(contraction)
    else:
        raise RuntimeError("bounded potential proposal did not converge; no claim")
    denominator = 1 << 20
    rational = [[F(0) for _ in range(1 << width)] for _ in range(width + 1)]
    for r, h in valid:
        rational[r][h] = F(ceil(vector[r][h] * denominator), denominator)
    evaluated = empty_step(rational, width, z, tau, cycle=cycle)
    slacks = [contraction * rational[r][h] - evaluated[r][h] for r, h in valid]
    if not all(slack >= 0 for slack in slacks):
        raise AssertionError("exact positive-vector verification failed")
    minimum = min(rational[r][h] for r, h in valid)
    maximum = max(rational[r][h] for r, h in valid)
    ratio = maximum / minimum
    if verbose:
        print(f"PASS: {len(valid)} EXACT live-state inequalities; z={z}, tau={tau}, lambda={contraction}.")
        print(f"Proposal iterations={iteration}; verified min slack={min(slacks)}.")
        print(f"For free live terminal: v_L(r,H) <= ({contraction})^L v(r,H)/min(v).")
        print(f"Diagnostic only: log2(max(v)/min(v))={log2(float(ratio)):.6f}; "
              f"uniform log2 envelope at L2048={2048 * log2(float(contraction)) + log2(float(ratio)):.6f}.")
    return rational, z, tau, contraction


if __name__ == "__main__":
    check()
