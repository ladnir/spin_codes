#!/usr/bin/env python3
"""Create a block-outer finite certificate manifest from existing artifacts.

This is the candidate-code workflow bridge:

1. run block_outer_upgrade_probe.py with --model spectrum-csv;
2. run certify_geometric_tail_from_csv.py on the prefix CSV;
3. run this script to collect paths, hashes, prefix constants, tail constants,
   and residual constants into one manifest;
4. run verify_block_certificate_manifest.py on the emitted manifest.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

from verify_block_certificate_manifest import (
    large_r_low_slot_max_block,
    product_spectrum_log2_w,
    sha256_file,
)
from verify_checkpoint_residual import bounded_r_log2_ratio, high_slot_exponent
from verify_dense_claims import ETA_CRIT, scan_gap


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def read_one(path: Path) -> dict[str, str]:
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    if len(rows) != 1:
        raise SystemExit(f"{path} must contain exactly one row, found {len(rows)}")
    return rows[0]


def repo_rel(path: Path) -> str:
    path = path.resolve()
    try:
        return path.relative_to(repo_root()).as_posix()
    except ValueError:
        return str(path)


def local_spectrum_stats(path: Path) -> tuple[int, int, int, int]:
    total = 0
    min_nonzero = None
    max_weight = 0
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            weight = int(row["weight"])
            count = int(row["count"])
            total += count
            max_weight = max(max_weight, weight)
            if weight > 0 and count > 0:
                min_nonzero = weight if min_nonzero is None else min(min_nonzero, weight)
    if min_nonzero is None:
        raise SystemExit(f"{path} has no nonzero local codewords")
    dim = total.bit_length() - 1
    if total != 1 << dim:
        raise SystemExit(f"{path} spectrum total is {total}, not a power of two")
    return total, dim, min_nonzero, max_weight


def as_float(row: dict[str, str], key: str) -> float:
    return float(row[key])


def as_int(row: dict[str, str], key: str) -> int:
    return int(row[key])


def residual_constants(
    *,
    spectrum_path: Path,
    k: int,
    n: int,
    blocks: int,
    sigma: int,
    delta: float,
    h_min: int,
    r_cap: int,
    theta_slot: float,
    z: float,
    eta: float,
    linear_eta_hi: float,
    run_theta: float,
    xi: float,
    gap_step: float,
) -> dict[str, float]:
    offset_s = sigma - math.ceil(math.log2(k))
    eta_lo = h_min / n
    log2_w = product_spectrum_log2_w(spectrum_path=spectrum_path, blocks=blocks, z=z)
    log2_q = bounded_r_log2_ratio(n=n, h=h_min, r_cap=r_cap, theta=theta_slot, z=z)
    high_slot, _ = high_slot_exponent(delta=delta, theta=theta_slot, eta=eta, z=z)
    low_value, _, _, _, _ = large_r_low_slot_max_block(
        log2_w=log2_w,
        n=n,
        offset_s=offset_s,
        h_min=h_min,
        r_cap=r_cap,
        eta=eta,
        theta=theta_slot,
        z=z,
    )
    low_union = math.log2(3.0) + 4.0 * math.log2(n) + low_value
    gap = scan_gap(eta_lo=eta_lo, eta_hi=linear_eta_hi, theta=run_theta, delta=delta, xi=xi, step=gap_step)
    return {
        "log2_w_out": log2_w,
        "bounded_r_low_slot_log2_q": log2_q,
        "bounded_r_high_slot_exponent": high_slot,
        "large_r_low_slot_ledger_log2": low_value,
        "large_r_low_slot_union_log2": low_union,
        "linear_handoff_gap": gap.worst_gap,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--name", required=True)
    parser.add_argument("--description", default="")
    parser.add_argument("--local-code-name", required=True)
    parser.add_argument("--local-spectrum-csv", type=Path, required=True)
    parser.add_argument("--summary-csv", type=Path, required=True)
    parser.add_argument("--prefix-csv", type=Path, required=True)
    parser.add_argument("--tailcert-csv", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--spectrum-source", default="")
    parser.add_argument("--certified-total-log2-max", type=float, required=True)
    parser.add_argument("--max-observed-log2-ratio", type=float, default=-0.9)
    parser.add_argument("--max-tail-residual-log2", type=float, default=-500.0)
    parser.add_argument("--analytic-ratio-log2", type=float, required=True)
    parser.add_argument("--max-analytic-tail-residual-log2", type=float, default=-500.0)
    parser.add_argument("--r-cap", type=int, default=12)
    parser.add_argument("--theta-slot", type=float, default=0.35)
    parser.add_argument("--z", type=float, default=0.4)
    parser.add_argument("--eta", type=float, default=0.25)
    parser.add_argument("--linear-eta-hi", type=float, default=ETA_CRIT)
    parser.add_argument("--run-theta", type=float, default=0.001)
    parser.add_argument("--xi", type=float, default=8.0)
    parser.add_argument("--gap-step", type=float, default=1e-4)
    parser.add_argument("--log2-w-out-max", type=float, default=40.0)
    parser.add_argument("--bounded-r-low-slot-log2-q-max", type=float, default=-0.1)
    parser.add_argument("--bounded-r-high-slot-exponent-min", type=float, default=0.05)
    parser.add_argument("--large-r-low-slot-ledger-log2-max", type=float, default=-250.0)
    parser.add_argument("--large-r-low-slot-union-log2-max", type=float, default=-160.0)
    parser.add_argument("--linear-handoff-gap-max", type=float, default=-0.002)
    args = parser.parse_args()

    spectrum_path = args.local_spectrum_csv.resolve()
    summary_path = args.summary_csv.resolve()
    prefix_path = args.prefix_csv.resolve()
    tailcert_path = args.tailcert_csv.resolve()
    summary = read_one(summary_path)
    tail = read_one(tailcert_path)
    _, spectrum_dim, spectrum_min, spectrum_max = local_spectrum_stats(spectrum_path)

    k = as_int(summary, "k")
    n = as_int(summary, "N")
    blocks = as_int(summary, "blocks")
    sigma = as_int(summary, "sigma")
    delta = as_float(summary, "delta")
    h_max = as_int(summary, "h_max")
    block_bits = as_int(summary, "block_bits")
    d0 = as_int(summary, "d0")
    if d0 != spectrum_min:
        raise SystemExit(f"summary d0={d0} but spectrum minimum distance is {spectrum_min}")
    if block_bits != spectrum_dim:
        raise SystemExit(f"summary block_bits={block_bits} but spectrum dimension is {spectrum_dim}")

    residual = residual_constants(
        spectrum_path=spectrum_path,
        k=k,
        n=n,
        blocks=blocks,
        sigma=sigma,
        delta=delta,
        h_min=h_max,
        r_cap=args.r_cap,
        theta_slot=args.theta_slot,
        z=args.z,
        eta=args.eta,
        linear_eta_hi=args.linear_eta_hi,
        run_theta=args.run_theta,
        xi=args.xi,
        gap_step=args.gap_step,
    )
    manifest: dict[str, Any] = {
        "schema": "permute-conv.block-certificate.v1",
        "name": args.name,
        "description": args.description,
        "parameters": {
            "k": k,
            "n": n,
            "delta": delta,
            "d": as_int(summary, "d"),
            "sigma": sigma,
            "block_bits": block_bits,
            "blocks": blocks,
            "h_max": h_max,
            "certified_total_log2_max": args.certified_total_log2_max,
        },
        "local_code": {
            "name": args.local_code_name,
            "length": spectrum_max,
            "dimension": spectrum_dim,
            "distance": spectrum_min,
            "spectrum_csv": repo_rel(spectrum_path),
            "spectrum_sha256": sha256_file(spectrum_path),
            "spectrum_source": args.spectrum_source,
        },
        "artifacts": {
            "summary_csv": {"path": repo_rel(summary_path), "sha256": sha256_file(summary_path)},
            "prefix_csv": {"path": repo_rel(prefix_path), "sha256": sha256_file(prefix_path)},
            "tailcert_csv": {"path": repo_rel(tailcert_path), "sha256": sha256_file(tailcert_path)},
        },
        "prefix": {
            "total_log2": as_float(summary, "total_log2"),
            "peak_h": as_int(summary, "peak_h"),
            "peak_term_log2": as_float(summary, "peak_term_log2"),
            "peak_outer_log2": as_float(summary, "peak_outer_log2"),
            "peak_inner_log2": as_float(summary, "peak_inner_log2"),
        },
        "tail": {
            "check_from": as_int(tail, "check_from"),
            "tail_after": as_int(tail, "tail_after"),
            "observed_worst_log2_ratio": as_float(tail, "observed_worst_log2_ratio"),
            "max_observed_log2_ratio": args.max_observed_log2_ratio,
            "tail_residual_log2": as_float(tail, "tail_residual_log2"),
            "max_tail_residual_log2": args.max_tail_residual_log2,
            "analytic_ratio_log2": args.analytic_ratio_log2,
            "max_analytic_tail_residual_log2": args.max_analytic_tail_residual_log2,
        },
        "residual": {
            "h_min": h_max,
            "r_cap": args.r_cap,
            "theta_slot": args.theta_slot,
            "z": args.z,
            "eta": args.eta,
            "linear_eta_hi": args.linear_eta_hi,
            "run_theta": args.run_theta,
            "xi": args.xi,
            "gap_step": args.gap_step,
            **residual,
            "thresholds": {
                "log2_w_out_max": args.log2_w_out_max,
                "bounded_r_low_slot_log2_q_max": args.bounded_r_low_slot_log2_q_max,
                "bounded_r_high_slot_exponent_min": args.bounded_r_high_slot_exponent_min,
                "large_r_low_slot_ledger_log2_max": args.large_r_low_slot_ledger_log2_max,
                "large_r_low_slot_union_log2_max": args.large_r_low_slot_union_log2_max,
                "linear_handoff_gap_max": args.linear_handoff_gap_max,
            },
        },
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as f:
        json.dump(manifest, f, indent=2)
        f.write("\n")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
