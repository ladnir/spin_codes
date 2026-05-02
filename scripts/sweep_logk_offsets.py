#!/usr/bin/env python3
"""Sweep sigma = ceil(log2(k)) + c for the dense large-k evaluator."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

from dense_largek_eval import (
    ETA_CRIT,
    GapResult,
    adaptive_small_h_localtail,
    adaptive_small_h_rigorous,
    exact_smallw_inner_bound,
    geometric_window_bound,
    linear_window_gap,
    outer_generating_bound,
    outer_small_h_prefix,
    run_mixture_model_log2,
    runmix_gamma_from_q,
    single_run_product_log2,
    small_h_heuristic_log2,
    total_length,
    tiny_window_bound,
)


def parse_deltas(text: str) -> list[float]:
    return [float(x) for x in text.split(",") if x.strip()]


def evaluate_point(
    *,
    k_msg: int,
    sigma: int,
    delta: float,
    z: float,
    rho: float,
    kappa: float,
    theta: float,
    xi: float,
    eta_hi: float,
    gap_step: float,
    heuristic_h_cap: int,
    n_mode: str,
    outer_smallh_mode: str,
    bracket_h_cap: int,
    only_bracket: bool,
) -> dict[str, float]:
    n = total_length(k_msg, sigma, n_mode)
    parity_n = n - k_msg
    h_mid_hi = int(math.floor(ETA_CRIT * n))
    bracket_h_hi = min(bracket_h_cap, h_mid_hi)

    w_out = None
    h0 = 0
    s_tiny = 0.0
    s_low = 0.0
    gap = None
    log2_s_lin = float("-inf")
    log2_s_top = float("-inf")
    log2_total = float("-inf")
    h_cap = 0
    log2_h_small = float("-inf")
    log2_h_resid = float("-inf")
    log2_h_total = float("-inf")
    log2_r_small = float("-inf")
    log2_r_total = float("-inf")
    a_stop_h = 0
    a_peak_h = 0
    a_peak_log2 = float("-inf")
    log2_a_small = float("-inf")
    log2_a_resid = float("-inf")
    log2_a_total = float("-inf")
    l_stop_h = 0
    l_peak_h = 0
    l_peak_log2 = float("-inf")
    log2_l_tail = float("-inf")
    log2_l_total = float("-inf")

    if not only_bracket:
        w_out = outer_generating_bound(k_msg, sigma, z)
        h0, s_tiny = tiny_window_bound(n, sigma, delta, z, xi, kappa, w_out)
        s_low = geometric_window_bound(
            n=n,
            h_lo=h0 + 1,
            h_hi=h_mid_hi,
            z=z,
            rho=rho,
            w_out=w_out,
        )

        gap = linear_window_gap(delta, theta, xi, eta_hi, gap_step)
        lin_lo = int(math.ceil(ETA_CRIT * n))
        lin_hi = int(math.floor(min(eta_hi, 1.0 - theta) * n))
        lin_count = max(0, lin_hi - lin_lo + 1)
        log2_s_lin = float("-inf") if lin_count == 0 else math.log2(lin_count) + n * gap.worst_gap
        s_lin = 0.0 if log2_s_lin < -900.0 else 2.0 ** log2_s_lin

        top_lo = int(math.ceil(eta_hi * n))
        top_count = max(0, n - top_lo + 1)
        top_exp = (
            float("-inf")
            if eta_hi >= 1.0
            else (
                -(10**9)
                if eta_hi <= 0.0
                else (
                    -eta_hi * math.log2(eta_hi)
                    - (1.0 - eta_hi) * math.log2(1.0 - eta_hi)
                    - 0.5
                )
            )
        )
        log2_s_top = float("-inf") if top_count == 0 else math.log2(top_count) + n * top_exp
        s_top = 0.0 if log2_s_top < -900.0 else 2.0 ** log2_s_top

        total = s_tiny + s_low + s_lin + s_top
        log2_total = math.log2(total) if total > 0.0 else float("-inf")

        h_cap = min(heuristic_h_cap, h_mid_hi)
        log2_h_small = small_h_heuristic_log2(
            n=n,
            k_msg=k_msg,
            sigma=sigma,
            h_lo=1,
            h_hi=h_cap,
            theta=theta,
            outer_mode=outer_smallh_mode,
            parity_n=parity_n,
        )
        log2_r_small = float("-inf")
        outer_logs = outer_small_h_prefix(k_msg, sigma, h_cap, outer_mode=outer_smallh_mode, parity_n=parity_n)
        for h in range(1, h_cap + 1):
            out = outer_logs[h]
            inn = exact_smallw_inner_bound(n, h, sigma, delta, xi)
            if inn <= 0.0:
                continue
            term = out + math.log2(inn)
            if log2_r_small == float("-inf"):
                log2_r_small = term
            else:
                hi = max(log2_r_small, term)
                lo = min(log2_r_small, term)
                log2_r_small = hi + math.log2(1.0 + 2.0 ** (lo - hi))
        h_resid = geometric_window_bound(
            n=n,
            h_lo=h_cap + 1,
            h_hi=h_mid_hi,
            z=z,
            rho=rho,
            w_out=w_out,
        )
        log2_h_resid = math.log2(h_resid) if h_resid > 0.0 else float("-inf")
        if log2_h_small == float("-inf"):
            log2_h_total = log2_h_resid
        elif log2_h_resid == float("-inf"):
            log2_h_total = log2_h_small
        else:
            hi = max(log2_h_small, log2_h_resid)
            lo = min(log2_h_small, log2_h_resid)
            log2_h_total = hi + math.log2(1.0 + 2.0 ** (lo - hi))

        if log2_r_small == float("-inf"):
            log2_r_total = log2_h_resid
        elif log2_h_resid == float("-inf"):
            log2_r_total = log2_r_small
        else:
            hi = max(log2_r_small, log2_h_resid)
            lo = min(log2_r_small, log2_h_resid)
            log2_r_total = hi + math.log2(1.0 + 2.0 ** (lo - hi))

        log2_a_small, a_stop_h, a_peak_h, a_peak_log2 = adaptive_small_h_rigorous(
            n=n,
            k_msg=k_msg,
            sigma=sigma,
            delta=delta,
            xi=xi,
            h_hi=h_mid_hi,
            outer_mode=outer_smallh_mode,
            parity_n=parity_n,
        )
        a_resid = geometric_window_bound(
            n=n,
            h_lo=a_stop_h + 1,
            h_hi=h_mid_hi,
            z=z,
            rho=rho,
            w_out=w_out,
        )
        log2_a_resid = math.log2(a_resid) if a_resid > 0.0 else float("-inf")
        if log2_a_small == float("-inf"):
            log2_a_total = log2_a_resid
        elif log2_a_resid == float("-inf"):
            log2_a_total = log2_a_small
        else:
            hi = max(log2_a_small, log2_a_resid)
            lo = min(log2_a_small, log2_a_resid)
            log2_a_total = hi + math.log2(1.0 + 2.0 ** (lo - hi))

        log2_l_total, l_stop_h, l_peak_h, l_peak_log2, log2_l_tail = adaptive_small_h_localtail(
            n=n,
            k_msg=k_msg,
            sigma=sigma,
            delta=delta,
            xi=xi,
            h_hi=bracket_h_hi if outer_smallh_mode == "banded" else h_mid_hi,
            outer_mode=outer_smallh_mode,
            parity_n=parity_n,
        )
        if log2_s_lin != float("-inf"):
            hi = max(log2_l_total, log2_s_lin)
            lo = min(log2_l_total, log2_s_lin)
            log2_l_total = hi + math.log2(1.0 + 2.0 ** (lo - hi))
        if log2_s_top != float("-inf"):
            hi = max(log2_l_total, log2_s_top)
            lo = min(log2_l_total, log2_s_top)
            log2_l_total = hi + math.log2(1.0 + 2.0 ** (lo - hi))
    else:
        lin_gap = linear_window_gap(delta, theta, xi, eta_hi, gap_step)
        lin_lo = int(math.ceil(ETA_CRIT * n))
        lin_hi = int(math.floor(min(eta_hi, 1.0 - theta) * n))
        lin_count = max(0, lin_hi - lin_lo + 1)
        log2_s_lin = float("-inf") if lin_count == 0 else math.log2(lin_count) + n * lin_gap.worst_gap
        top_lo = int(math.ceil(eta_hi * n))
        top_count = max(0, n - top_lo + 1)
        top_exp = (
            float("-inf")
            if eta_hi >= 1.0
            else (
                -(10**9)
                if eta_hi <= 0.0
                else (
                    -eta_hi * math.log2(eta_hi)
                    - (1.0 - eta_hi) * math.log2(1.0 - eta_hi)
                    - 0.5
                )
            )
        )
        log2_s_top = float("-inf") if top_count == 0 else math.log2(top_count) + n * top_exp

        log2_l_total, l_stop_h, l_peak_h, l_peak_log2, log2_l_tail = adaptive_small_h_localtail(
            n=n,
            k_msg=k_msg,
            sigma=sigma,
            delta=delta,
            xi=xi,
            h_hi=bracket_h_hi if outer_smallh_mode == "banded" else h_mid_hi,
            outer_mode=outer_smallh_mode,
            parity_n=parity_n,
        )
        if log2_s_lin != float("-inf"):
            hi = max(log2_l_total, log2_s_lin)
            lo = min(log2_l_total, log2_s_lin)
            log2_l_total = hi + math.log2(1.0 + 2.0 ** (lo - hi))
        if log2_s_top != float("-inf"):
            hi = max(log2_l_total, log2_s_top)
            lo = min(log2_l_total, log2_s_top)
            log2_l_total = hi + math.log2(1.0 + 2.0 ** (lo - hi))

    q_single = min(1.0, 2.0 * delta)
    log2_q_total, q_peak_h, q_peak_log2 = single_run_product_log2(
        k_msg=k_msg,
        sigma=sigma,
        h_lo=1,
        h_hi=bracket_h_hi,
        q_single=q_single,
        outer_mode=outer_smallh_mode,
        parity_n=parity_n,
    )
    if log2_s_lin != float("-inf"):
        hi = max(log2_q_total, log2_s_lin)
        lo = min(log2_q_total, log2_s_lin)
        log2_q_total = hi + math.log2(1.0 + 2.0 ** (lo - hi))
    if log2_s_top != float("-inf"):
        hi = max(log2_q_total, log2_s_top)
        lo = min(log2_q_total, log2_s_top)
        log2_q_total = hi + math.log2(1.0 + 2.0 ** (lo - hi))

    gamma_runmix = runmix_gamma_from_q(q_single)
    log2_m_total, m_peak_h, m_peak_log2 = run_mixture_model_log2(
        n=n,
        k_msg=k_msg,
        sigma=sigma,
        h_lo=1,
        h_hi=bracket_h_hi,
        q_single=q_single,
        gamma=gamma_runmix,
        outer_mode=outer_smallh_mode,
        parity_n=parity_n,
    )
    if log2_s_lin != float("-inf"):
        hi = max(log2_m_total, log2_s_lin)
        lo = min(log2_m_total, log2_s_lin)
        log2_m_total = hi + math.log2(1.0 + 2.0 ** (lo - hi))
    if log2_s_top != float("-inf"):
        hi = max(log2_m_total, log2_s_top)
        lo = min(log2_m_total, log2_s_top)
        log2_m_total = hi + math.log2(1.0 + 2.0 ** (lo - hi))

    gamma_runmix_legacy = 0.2 * q_single
    log2_m0_total, m0_peak_h, m0_peak_log2 = run_mixture_model_log2(
        n=n,
        k_msg=k_msg,
        sigma=sigma,
        h_lo=1,
        h_hi=bracket_h_hi,
        q_single=q_single,
        gamma=gamma_runmix_legacy,
        outer_mode=outer_smallh_mode,
        parity_n=parity_n,
    )
    if log2_s_lin != float("-inf"):
        hi = max(log2_m0_total, log2_s_lin)
        lo = min(log2_m0_total, log2_s_lin)
        log2_m0_total = hi + math.log2(1.0 + 2.0 ** (lo - hi))
    if log2_s_top != float("-inf"):
        hi = max(log2_m0_total, log2_s_top)
        lo = min(log2_m0_total, log2_s_top)
        log2_m0_total = hi + math.log2(1.0 + 2.0 ** (lo - hi))

    return {
        "n": float(n),
        "w_out_log2": math.log2(w_out) if (w_out is not None and w_out > 0.0) else float("-inf"),
        "h0": float(h0),
        "s_tiny_log2": math.log2(s_tiny) if s_tiny > 0.0 else float("-inf"),
        "s_low_log2": math.log2(s_low) if s_low > 0.0 else float("-inf"),
        "gap_worst": gap.worst_gap if gap is not None else lin_gap.worst_gap,
        "gap_eta": gap.worst_eta if gap is not None else lin_gap.worst_eta,
        "s_lin_log2": log2_s_lin,
        "s_top_log2": log2_s_top,
        "total_log2": log2_total,
        "h_cap": float(h_cap),
        "h_small_log2": log2_h_small,
        "h_tail_log2": log2_h_resid,
        "h_total_log2": log2_h_total,
        "r_small_log2": log2_r_small,
        "r_total_log2": log2_r_total,
        "a_stop_h": float(a_stop_h),
        "a_peak_h": float(a_peak_h),
        "a_peak_log2": a_peak_log2,
        "a_small_log2": log2_a_small,
        "a_tail_log2": log2_a_resid,
        "a_total_log2": log2_a_total,
        "l_stop_h": float(l_stop_h),
        "l_peak_h": float(l_peak_h),
        "l_peak_log2": l_peak_log2,
        "l_tail_log2": log2_l_tail,
        "l_total_log2": log2_l_total,
        "q_single": q_single,
        "q_peak_h": float(q_peak_h),
        "q_peak_log2": q_peak_log2,
        "q_total_log2": log2_q_total,
        "m_gamma": gamma_runmix,
        "m_peak_h": float(m_peak_h),
        "m_peak_log2": m_peak_log2,
        "m_total_log2": log2_m_total,
        "m0_gamma": gamma_runmix_legacy,
        "m0_peak_h": float(m0_peak_h),
        "m0_peak_log2": m0_peak_log2,
        "m0_total_log2": log2_m0_total,
        "b_opt_log2": log2_l_total,
        "b_ctr_log2": log2_m_total,
        "b_pess_log2": log2_q_total,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--power-min", type=int, default=10)
    parser.add_argument("--power-max", type=int, default=20)
    parser.add_argument("--powers", type=str, default="", help="Optional comma-separated explicit powers.")
    parser.add_argument("--c-min", type=int, default=0)
    parser.add_argument("--c-max", type=int, default=4)
    parser.add_argument("--deltas", type=str, default="0.10,0.11,0.12")
    parser.add_argument("--z", type=float, default=1.0 / 3.0)
    parser.add_argument("--rho", type=float, default=1.0 / 4.0)
    parser.add_argument("--kappa", type=float, default=0.5)
    parser.add_argument("--theta", type=float, default=0.005)
    parser.add_argument("--xi", type=float, default=8.0)
    parser.add_argument("--eta-hi", type=float, default=0.99)
    parser.add_argument("--gap-step", type=float, default=1e-5)
    parser.add_argument("--heuristic-h-cap", type=int, default=64)
    parser.add_argument("--n-mode", choices=("paper", "extended"), default="paper")
    parser.add_argument("--outer-smallh-mode", choices=("conv", "banded"), default="banded")
    parser.add_argument("--bracket-h-cap", type=int, default=24)
    parser.add_argument("--only-bracket", action="store_true")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    if args.powers.strip():
        powers = [int(x) for x in args.powers.split(",") if x.strip()]
    else:
        powers = list(range(args.power_min, args.power_max + 1, 2))
    deltas = parse_deltas(args.deltas)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "power",
                "k",
                "sigma",
                "c",
                "delta",
                "n",
                "w_out_log2",
                "h0",
                "s_tiny_log2",
                "s_low_log2",
                "gap_worst",
                "gap_eta",
                "s_lin_log2",
                "s_top_log2",
                "total_log2",
                "h_cap",
                "h_small_log2",
                "h_tail_log2",
                "h_total_log2",
                "r_small_log2",
                "r_total_log2",
                "a_stop_h",
                "a_peak_h",
                "a_peak_log2",
                "a_small_log2",
                "a_tail_log2",
                "a_total_log2",
                "l_stop_h",
                "l_peak_h",
                "l_peak_log2",
                "l_tail_log2",
                "l_total_log2",
                "q_single",
                "q_peak_h",
                "q_peak_log2",
                "q_total_log2",
                "m_gamma",
                "m_peak_h",
                "m_peak_log2",
                "m_total_log2",
                "m0_gamma",
                "m0_peak_h",
                "m0_peak_log2",
                "m0_total_log2",
                "b_opt_log2",
                "b_ctr_log2",
                "b_pess_log2",
            ],
        )
        writer.writeheader()

        for power in powers:
            k_msg = 1 << power
            base_sigma = power
            for c in range(args.c_min, args.c_max + 1):
                sigma = base_sigma + c
                print(f"running k=2^{power}, sigma={sigma} (c={c})...")
                for delta in deltas:
                    row = evaluate_point(
                        k_msg=k_msg,
                        sigma=sigma,
                        delta=delta,
                        z=args.z,
                        rho=args.rho,
                        kappa=args.kappa,
                        theta=args.theta,
                        xi=args.xi,
                        eta_hi=args.eta_hi,
                        gap_step=args.gap_step,
                        heuristic_h_cap=args.heuristic_h_cap,
                        n_mode=args.n_mode,
                        outer_smallh_mode=args.outer_smallh_mode,
                        bracket_h_cap=args.bracket_h_cap,
                        only_bracket=args.only_bracket,
                    )
                    row.update(
                        {
                            "power": power,
                            "k": k_msg,
                            "sigma": sigma,
                            "c": c,
                            "delta": delta,
                        }
                    )
                    writer.writerow(row)
                    f.flush()

    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
