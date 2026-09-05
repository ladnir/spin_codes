#!/usr/bin/env python3
"""Fit and upper-bound a uniform per-early-start tau on the isolated slice."""

from __future__ import annotations

import math
from pathlib import Path

from check_q1_law import parse_enum
from check_isolated_earlycount import exact_slice_prob, isolated_early_model


def main() -> int:
    configs = [
        (Path(r"C:\Users\peter\repo\libOTe\enum_WC12_WC120.120.12.txt"), 120, 0.115217),
        (Path(r"C:\Users\peter\repo\libOTe\enum_WC16_WC128.128.16.txt"), 128, 0.125),
    ]

    datasets: list[tuple[str, list[list[float]], int, int, int]] = []
    for path, n, delta in configs:
        rows = parse_enum(path)
        cut = math.floor(delta * n)
        l = math.ceil(2.0 * delta * n)
        datasets.append((path.name, rows, n, cut, l))

    best_fit: tuple[float, float] | None = None
    for tau_i in range(1, 301):
        tau = tau_i / 1000.0
        err = 0.0
        for _, rows, n, cut, l in datasets:
            for h in range(1, 11):
                exact = exact_slice_prob(rows, n, h, cut)
                pred = isolated_early_model(n, h, l, tau)
                gap = math.log2(pred / exact)
                err += gap * gap
        cand = (err, tau)
        if best_fit is None or cand < best_fit:
            best_fit = cand

    safe_tau = None
    for tau_i in range(1, 301):
        tau = tau_i / 1000.0
        ok = True
        worst_over = 0.0
        for _, rows, n, cut, l in datasets:
            for h in range(1, 11):
                exact = exact_slice_prob(rows, n, h, cut)
                pred = isolated_early_model(n, h, l, tau)
                gap = math.log2(pred / exact)
                if gap < 0.0:
                    ok = False
                    break
                worst_over = max(worst_over, gap)
            if not ok:
                break
        if ok:
            safe_tau = (tau, worst_over)
            break

    out = Path(__file__).with_name("isolated_early_tau_report.txt")
    with out.open("w") as f:
        f.write("isolated early-count tau fit\n")
        f.write(f"  best-fit tau over h=1..10: {best_fit[1]:.3f}\n")
        f.write(f"  best-fit objective      : {best_fit[0]:.6f}\n")
        if safe_tau is not None:
            f.write(f"  first safe tau          : {safe_tau[0]:.3f}\n")
            f.write(f"  worst over-bound bits   : {safe_tau[1]:.6f}\n")
        f.write("\n")
        for name, rows, n, cut, l in datasets:
            tau = safe_tau[0] if safe_tau is not None else best_fit[1]
            f.write(f"{name}: n={n} cut={cut} L={l} tau={tau:.3f}\n")
            f.write("h, exact, model, gap_bits\n")
            for h in range(1, 11):
                exact = exact_slice_prob(rows, n, h, cut)
                pred = isolated_early_model(n, h, l, tau)
                gap = math.log2(pred / exact)
                f.write(f"{h}, {exact:.6e}, {pred:.6e}, {gap:.6f}\n")
            f.write("\n")

    print(out.read_text(), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
