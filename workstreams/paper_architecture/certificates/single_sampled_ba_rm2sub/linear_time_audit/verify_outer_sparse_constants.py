"""Outward-rounded constants for the sparse Golay--BA-3 endpoint proof."""

from __future__ import annotations

import json

from mpmath import iv


iv.dps = 60
delta = iv.mpf(13) / 125
b0 = iv.mpf(1) / 1000
c = 4 * delta / (1 - delta)
q = iv.sqrt(c)
kappa = 4 * q / (1 - q)
x_max = 12 * b0

log_derivative = (
    iv.log(kappa * x_max / (1 - 2 * x_max))
    + 1
    + 2 * x_max / (1 - 2 * x_max)
)
constant = (
    iv.e
    * (2**12 - 1)
    / 24
    * (4 * kappa / (1 - 8 * b0)) ** 4
)
rho = constant * b0**3
odd_factor = 2 * delta / (1 - delta)


def upper(value) -> float:
    return float(value.b)


result = {
    "status": "proved",
    "claims": {
        "accumulator_tail_base_below_one": upper(c) < 1.0,
        "odd_tail_factor_below_one": upper(odd_factor) < 1.0,
        "first_moment_decreases_with_r": upper(log_derivative) < 0.0,
        "sparse_geometric_base_below_0_657": upper(rho) < 0.657,
    },
    "upper_endpoints": {
        "c": upper(c),
        "odd_factor": upper(odd_factor),
        "log_derivative_at_x_0_012": upper(log_derivative),
        "K": upper(kappa),
        "C_star": upper(constant),
        "C_star_times_b0_cubed": upper(rho),
    },
    "interval_dps": iv.dps,
}

if not all(result["claims"].values()):
    result["status"] = "failed"
print(json.dumps(result, indent=2, sort_keys=True))
if result["status"] != "proved":
    raise SystemExit(1)
