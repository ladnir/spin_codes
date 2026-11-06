"""Outward sparse and fixed-row continuum gates for the shared Turbo inner.

These are local inequalities. The finite-region limit and first-moment
assembly still require a proof for the new inner distribution.
"""

from __future__ import annotations

from fractions import Fraction as F
import json
from pathlib import Path
import sys

from flint import arb, ctx

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from research.workstreams.turbo_spin_migration.shared_inner_dense_certificate import (
    exact_segments, num, upper,
)
from research.workstreams.turbo_spin_migration.shared_inner_dense_screen import (
    M, sparse_gate,
)


ROOT = Path(__file__).resolve().parents[1]
PRIOR_FIXED = ROOT / "inner_design/imt_asymptotic/FIXED_LIMIT.json"


def sparse_outward():
    # The exact rational Collatz residuals are checked by sparse_gate().
    sparse = sparse_gate()
    amax = F(1, 10_000)
    row_values = []
    for lo, hi, slope, intercept in exact_segments():
        for x in (lo, hi):
            xx = num(x)
            cost = (num(slope * x + intercept) + num(2).log()
                    + xx * xx.log() + (1 - xx) * (1 - xx).log())
            row_values.append((upper(cost), x))
    row_bound, worst_x = max(row_values)
    expression = (
        num(row_bound) - num(F(8, 5)).log()
        + num(F(3, 5) / (1 - F(8, 5) * amax))
        - num(F(3, 4))
        + num(F(11, 100) * F(8, 5) / (1 - F(8, 5) * amax))
        + num(2).log() / num(F(39, 4))
    )
    bound = upper(expression)
    assert bound < 0
    return {"worst_outer_endpoint": str(worst_x),
            "sparse_exponent_per_alpha_upper": str(bound),
            "sparse_exponent_per_alpha_decimal": float(bound),
            "collatz_residuals_exact": sparse["collatz_residual_upper_over_alpha"]}


def fixed_one_two_outward():
    r, p0 = num(F(1, M)), num(F(1 << 21, M))
    ln2 = num(2).log()
    rows = []
    for q in (1, 2):
        gamma = q * ln2
        theta = gamma / p0
        a = num(F(1, 1 << q))
        f = (1 - a) / gamma
        gap = 2 * (gamma - 1 + a) / gamma**2
        edge = 2 * (1 - (1 + gamma) * a) / gamma**2
        k0 = [[arb(1), arb(0)], [arb(0), a]]
        k1 = [[arb(0), f], [r * f, a]]
        k2 = [[r * gap, edge], [r * edge, r * edge + a]]
        matrix = [[k0[i][j] + q * k1[i][j]
                   + (k2[i][j] if q == 2 else 0)
                   for j in range(2)] for i in range(2)]
        x, y = matrix[0]
        z, w = matrix[1]
        radius = (x + w + ((x - w)**2 + 4 * y * z).sqrt()) / 2
        exponent = (radius.log() + num(F(11, 100)) * theta
                    - q * ln2 / 2 + q * num(F(1281, 100000))
                    + q * ln2 / num(F(39, 4)))
        bound = upper(exponent / q)
        assert bound < 0
        rows.append({"q": q, "per_active_row_bit_upper": str(bound),
                     "decimal": float(bound)})
    return rows


def fixed_q_ge_three():
    saved = json.loads(PRIOR_FIXED.read_text(encoding="utf-8"))
    assert saved["status"] == "OUTWARD_CONDITIONAL_IMT_CONTINUUM_Q_GE_3_INEQUALITIES"
    assert saved["delta"] == "11/100" and saved["block_constant"] == "39/4"
    old = F(1 << 18, (1 << 19) - 1)
    new = F(1 << 21, M)
    increase = F(11, 100) * F(133, 125) * (1 / new - 1 / old)
    old_max = max(F(row["exponent_upper"]) for row in saved["checks"])
    # In fixed_limit.py, the reset probability enters only through a
    # positive term in the norm. The new probability is smaller. The
    # live-density change affects only the delta*theta term above.
    assert F(1, M) < F(1, (1 << 19) - 1)
    bound = old_max + increase
    assert bound < 0
    return {"old_endpoint_upper": str(old_max),
            "new_chernoff_increase": str(increase),
            "new_endpoint_upper_by_monotonicity": str(bound),
            "decimal": float(bound)}


def main():
    if not __debug__:
        raise RuntimeError("Certificate checks require Python assertions")
    ctx.prec = 256
    print(json.dumps({"status": "outward_local_inequalities",
                      "precision_bits": ctx.prec,
                      "sparse": sparse_outward(),
                      "fixed_one_two": fixed_one_two_outward(),
                      "fixed_q_ge_three": fixed_q_ge_three(),
                      "full_asymptotic_theorem": False}))


if __name__ == "__main__":
    main()
