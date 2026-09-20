#!/usr/bin/env python3
"""Verify a block-outer finite certificate manifest.

The manifest records the local spectrum source, finite-prefix artifacts, tail
certificate, and residual constants for a spectrum-csv block outer.  This
script checks that the referenced files match their hashes and recomputes the
block-spectrum residual constants from the manifest parameters.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

from block_outer_upgrade_probe import load_local_spectrum, local_spectrum_log2
from verify_checkpoint_residual import (
    bounded_r_log2_ratio,
    high_slot_exponent,
    h_candidates_for_fixed_rc,
    ledger_value,
)
from verify_dense_claims import scan_gap


DEFAULT_MANIFEST = Path(__file__).with_name("certificates") / "rm512_256_sig32_delta009.json"


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def resolve_repo_path(path_text: str) -> Path:
    path = Path(path_text)
    if path.is_absolute():
        return path
    return repo_root() / path


def read_one(path: Path) -> dict[str, str]:
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    if len(rows) != 1:
        raise SystemExit(f"{path} must contain exactly one row, found {len(rows)}")
    return rows[0]


def parse_log2(text: str) -> float:
    text = text.strip()
    if text == "-inf":
        return float("-inf")
    return float(text)


def read_term(path: Path, h: int) -> float:
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            if int(row["h"]) == h:
                return parse_log2(row["term_log2"])
    raise SystemExit(f"{path} has no row h={h}")


def read_spectrum_sum(path: Path) -> tuple[int, int, int]:
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
    return total, min_nonzero, max_weight


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def log2add(x: float, y: float) -> float:
    if x == float("-inf"):
        return y
    if y == float("-inf"):
        return x
    hi = max(x, y)
    lo = min(x, y)
    return hi + math.log2(1.0 + 2.0 ** (lo - hi))


def product_spectrum_log2_w(*, spectrum_path: Path, blocks: int, z: float) -> float:
    local = local_spectrum_log2(load_local_spectrum(str(spectrum_path)), z)
    all_log = blocks * local
    if all_log == float("-inf"):
        return float("-inf")
    if all_log <= 1e-10:
        val = math.expm1(all_log * math.log(2.0))
        return math.log2(val) if val > 0.0 else float("-inf")
    return all_log + math.log2(1.0 - 2.0 ** (-all_log))


def large_r_low_slot_max_block(
    *,
    log2_w: float,
    n: int,
    offset_s: int,
    h_min: int,
    r_cap: int,
    eta: float,
    theta: float,
    z: float,
) -> tuple[float, int, int, int, int]:
    h_high_global = math.floor(eta * n)
    m_slot = math.floor(theta * n)
    best_value = float("-inf")
    best_h = -1
    best_r = -1
    best_c = -1

    for r in range(r_cap + 1, h_high_global // 2 + 1):
        for c in (2 * r, 2 * r + 1):
            h_low = max(h_min, c)
            if h_low > h_high_global:
                continue
            for h in h_candidates_for_fixed_rc(n=n, h_low=h_low, h_high=h_high_global, c=c, m_slot=m_slot, z=z):
                value = ledger_value(
                    log2_w=log2_w,
                    n=n,
                    offset_s=offset_s,
                    z=z,
                    h=h,
                    r=r,
                    c=c,
                    m_slot=m_slot,
                )
                if value > best_value:
                    best_value = value
                    best_h = h
                    best_r = r
                    best_c = c
    return best_value, best_h, best_r, best_c, m_slot


def close(name: str, actual: float, expected: float, *, tol: float) -> bool:
    passed = abs(actual - expected) <= tol
    print(f"{name}: {'PASS' if passed else 'FAIL'} ({actual:.12g} vs {expected:.12g})")
    return passed


def check(name: str, passed: bool, detail: str) -> bool:
    print(f"{name}: {'PASS' if passed else 'FAIL'} ({detail})")
    return passed


def artifact_path(manifest: dict[str, Any], key: str) -> Path:
    return resolve_repo_path(manifest["artifacts"][key]["path"])


def verify_hashes(manifest: dict[str, Any]) -> bool:
    ok = True
    spectrum_path = resolve_repo_path(manifest["local_code"]["spectrum_csv"])
    ok &= check(
        "spectrum sha256",
        sha256_file(spectrum_path) == manifest["local_code"]["spectrum_sha256"],
        str(spectrum_path),
    )
    for key, meta in manifest["artifacts"].items():
        path = resolve_repo_path(meta["path"])
        ok &= check(f"{key} sha256", sha256_file(path) == meta["sha256"], str(path))
    return ok


def verify_prefix_and_tail(manifest: dict[str, Any], *, tol: float) -> bool:
    params = manifest["parameters"]
    local = manifest["local_code"]
    prefix_expected = manifest["prefix"]
    tail_expected = manifest["tail"]
    summary = read_one(artifact_path(manifest, "summary_csv"))
    tail = read_one(artifact_path(manifest, "tailcert_csv"))
    prefix_csv = artifact_path(manifest, "prefix_csv")
    spectrum_path = resolve_repo_path(local["spectrum_csv"])
    spectrum_total, spectrum_min, spectrum_max = read_spectrum_sum(spectrum_path)

    ok = True
    ok &= check("spectrum total", spectrum_total == (1 << local["dimension"]), f"sum={spectrum_total}")
    ok &= check("spectrum min distance", spectrum_min == local["distance"], f"d0={spectrum_min}")
    ok &= check("spectrum length", spectrum_max == local["length"], f"max_weight={spectrum_max}")
    ok &= check("k", int(summary["k"]) == params["k"], summary["k"])
    ok &= check("N", int(summary["N"]) == params["n"], summary["N"])
    ok &= check("d", int(summary["d"]) == params["d"], summary["d"])
    ok &= check("sigma", int(summary["sigma"]) == params["sigma"], summary["sigma"])
    ok &= check("block bits", int(summary["block_bits"]) == params["block_bits"], summary["block_bits"])
    ok &= check("blocks", int(summary["blocks"]) == params["blocks"], summary["blocks"])
    ok &= check("h_max", int(summary["h_max"]) == params["h_max"], summary["h_max"])
    ok &= check("model", summary["model"] == "spectrum-csv", summary["model"])
    ok &= close("delta", float(summary["delta"]), params["delta"], tol=1e-15)
    ok &= close("prefix total log2", float(summary["total_log2"]), prefix_expected["total_log2"], tol=tol)
    ok &= check("peak h", int(summary["peak_h"]) == prefix_expected["peak_h"], summary["peak_h"])
    ok &= close("peak term log2", float(summary["peak_term_log2"]), prefix_expected["peak_term_log2"], tol=tol)
    ok &= close("peak outer log2", float(summary["peak_outer_log2"]), prefix_expected["peak_outer_log2"], tol=tol)
    ok &= close("peak inner log2", float(summary["peak_inner_log2"]), prefix_expected["peak_inner_log2"], tol=tol)
    ok &= check(
        "prefix margin",
        float(summary["total_log2"]) <= params["certified_total_log2_max"],
        f"{summary['total_log2']} <= {params['certified_total_log2_max']}",
    )

    observed_ratio = float(tail["observed_worst_log2_ratio"])
    tail_residual = float(tail["tail_residual_log2"])
    tail_after = int(tail["tail_after"])
    last_term = read_term(prefix_csv, tail_after)
    analytic_ratio = tail_expected["analytic_ratio_log2"]
    q = 2.0 ** analytic_ratio
    analytic_tail = last_term + analytic_ratio - math.log2(1.0 - q)
    total_with_observed = log2add(float(summary["total_log2"]), tail_residual)
    total_with_analytic = log2add(float(summary["total_log2"]), analytic_tail)

    ok &= check("tail after", tail_after == tail_expected["tail_after"], str(tail_after))
    ok &= check("tail check from", int(tail["check_from"]) == tail_expected["check_from"], tail["check_from"])
    ok &= close("observed ratio log2", observed_ratio, tail_expected["observed_worst_log2_ratio"], tol=tol)
    ok &= check(
        "observed ratio threshold",
        observed_ratio <= tail_expected["max_observed_log2_ratio"],
        f"{observed_ratio:.6f} <= {tail_expected['max_observed_log2_ratio']:.6f}",
    )
    ok &= close("tail residual log2", tail_residual, tail_expected["tail_residual_log2"], tol=tol)
    ok &= check(
        "tail residual threshold",
        tail_residual <= tail_expected["max_tail_residual_log2"],
        f"{tail_residual:.6f} <= {tail_expected['max_tail_residual_log2']:.6f}",
    )
    ok &= check(
        "analytic tail threshold",
        analytic_tail <= tail_expected["max_analytic_tail_residual_log2"],
        f"{analytic_tail:.6f} <= {tail_expected['max_analytic_tail_residual_log2']:.6f}",
    )
    ok &= check(
        "certified total threshold",
        total_with_observed <= params["certified_total_log2_max"]
        and total_with_analytic <= params["certified_total_log2_max"],
        f"{total_with_observed:.6f}, {total_with_analytic:.6f} <= {params['certified_total_log2_max']:.6f}",
    )
    return ok


def verify_residual(manifest: dict[str, Any], *, tol: float) -> bool:
    params = manifest["parameters"]
    local = manifest["local_code"]
    residual = manifest["residual"]
    thresholds = residual["thresholds"]
    spectrum_path = resolve_repo_path(local["spectrum_csv"])
    offset_s = params["sigma"] - math.ceil(math.log2(params["k"]))
    eta_lo = residual["h_min"] / params["n"]

    log2_w = product_spectrum_log2_w(spectrum_path=spectrum_path, blocks=params["blocks"], z=residual["z"])
    log2_q = bounded_r_log2_ratio(
        n=params["n"],
        h=residual["h_min"],
        r_cap=residual["r_cap"],
        theta=residual["theta_slot"],
        z=residual["z"],
    )
    high_slot, _ = high_slot_exponent(
        delta=params["delta"],
        theta=residual["theta_slot"],
        eta=residual["eta"],
        z=residual["z"],
    )
    low_value, _, _, _, _ = large_r_low_slot_max_block(
        log2_w=log2_w,
        n=params["n"],
        offset_s=offset_s,
        h_min=residual["h_min"],
        r_cap=residual["r_cap"],
        eta=residual["eta"],
        theta=residual["theta_slot"],
        z=residual["z"],
    )
    low_union = math.log2(3.0) + 4.0 * math.log2(params["n"]) + low_value
    gap = scan_gap(
        eta_lo=eta_lo,
        eta_hi=residual["linear_eta_hi"],
        theta=residual["run_theta"],
        delta=params["delta"],
        xi=residual["xi"],
        step=residual["gap_step"],
    )

    ok = True
    ok &= close("log2 W_out(z)", log2_w, residual["log2_w_out"], tol=tol)
    ok &= close("bounded-r low-slot log2 q", log2_q, residual["bounded_r_low_slot_log2_q"], tol=tol)
    ok &= close("bounded-r high-slot exponent", high_slot, residual["bounded_r_high_slot_exponent"], tol=tol)
    ok &= close("large-r low-slot ledger", low_value, residual["large_r_low_slot_ledger_log2"], tol=tol)
    ok &= close("large-r low-slot union", low_union, residual["large_r_low_slot_union_log2"], tol=tol)
    ok &= close("linear handoff gap", gap.worst_gap, residual["linear_handoff_gap"], tol=tol)
    ok &= check("outer generating threshold", log2_w < thresholds["log2_w_out_max"], f"{log2_w:.6f}")
    ok &= check(
        "bounded-r ratio threshold",
        log2_q <= thresholds["bounded_r_low_slot_log2_q_max"],
        f"{log2_q:.6f}",
    )
    ok &= check(
        "high-slot exponent threshold",
        high_slot >= thresholds["bounded_r_high_slot_exponent_min"],
        f"{high_slot:.6f}",
    )
    ok &= check(
        "large-r ledger threshold",
        low_value <= thresholds["large_r_low_slot_ledger_log2_max"],
        f"{low_value:.6f}",
    )
    ok &= check(
        "large-r union threshold",
        low_union <= thresholds["large_r_low_slot_union_log2_max"],
        f"{low_union:.6f}",
    )
    ok &= check(
        "linear handoff threshold",
        gap.worst_gap <= thresholds["linear_handoff_gap_max"],
        f"{gap.worst_gap:.6f}",
    )
    return ok


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--tolerance", type=float, default=1e-6)
    args = parser.parse_args()

    with args.manifest.open() as f:
        manifest = json.load(f)

    print(f"Block certificate manifest: {manifest['name']}")
    ok = True
    ok &= check("schema", manifest["schema"] == "permute-conv.block-certificate.v1", manifest["schema"])
    ok &= verify_hashes(manifest)
    ok &= verify_prefix_and_tail(manifest, tol=args.tolerance)
    ok &= verify_residual(manifest, tol=args.tolerance)
    print(f"STATUS: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
