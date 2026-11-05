"""Directed neighboring-weight cover from one shared hybrid coefficient build.

At the fixed z/alpha of the frozen r512 proof, all local powered coefficients
and seven low-full tables are independent of the requested source weight.
The high exact rectangle is enlarged to F<=floor(rhi/2),S<=rhi-434. While
each table is alive, extract every requested (r,f,s=r-2f), divide by its own
M4096(f,s)*u^f*v^s, and retain only small normalized matrices and terminals.

Each r has its own exact indexed source counts and joint-missing Omega_r.
Keep r fixed throughout its alternating regional chain, including the
output-only terminal. Sum the complete per-r bounds; no interpolation,
monotonicity assumption, cross-r matrix envelope, or new alpha is used.
Odd source weights use the same audited exact counter. The band, not the
whole code, is the scope of any successful bound.

Frozen helpers remain unchanged. At most17 sizes in480..544,900 or1200s,
512MiB process peak,20GiB free space, stdout only, zero generated files.
"""

import argparse
from fractions import Fraction as F
import gc
import hashlib
import importlib.util
from itertools import combinations_with_replacement
import json
from math import comb
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True
HYBRID_PIN = "985c3fe04f055811976d8060a4b5a93063a6cabd81dd6f4c285769f60e096727"
CONTROL_512 = ("-60.8663747947622293153415391337338207056409052227494556026125591895557569",
               "2.47e-71")
CONTROL_223 = ("-82.3387645228211550551233638766999240969165462066241805609379602210861984",
               "5.63e-71")


def load_hybrid():
    path = Path(__file__).with_name("cheap_lane_cauchy_hybrid_directed.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != HYBRID_PIN:
        raise RuntimeError("frozen hybrid helper hash mismatch")
    spec = importlib.util.spec_from_file_location("neighbor_hybrid", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def join_weight(v, counts, anchor, core, r, regular, terminal, alpha, z, guard, started):
    """All ordered profiles at one r, never replacing counts by envelopes."""
    arb = v.arb
    run_bounds, run_counts = [arb(0) for _ in range(r+1)], [0]*(r+1)
    ordered = 0; worst = selected = None
    fullmax = r//2
    for f0 in range(fullmax+1):
        guard()
        for f1 in range(f0, fullmax+1):
            number = counts.profile_count(4096, r, f0, f1)
            if not number:
                assert r % 2 == 0 and f0 == f1 == fullmax
                continue
            runs = r-f0-f1
            r0, r1 = regular[f0], regular[f1]
            shared = v.matpow(v.matmul(r0, r1), 127)
            terms = [(f0, f1, v.positive(v.matmul(v.matmul(shared, r0), terminal[f1])[0][0]))]
            if f0 != f1:
                terms.append((f1, f0, v.positive(v.matmul(v.matmul(r1, shared), terminal[f0])[0][0])))
            run_counts[runs] += len(terms)*number; ordered += len(terms)
            for a, b, value in terms:
                mass = number*value; run_bounds[runs] += mass
                if worst is None or mass.upper() > worst[2].upper():
                    worst = (a, b, mass)
                if a == b == 223:
                    selected = mass
        if f0 % 64 == 0:
            print(json.dumps({"phase": "weight_profile_rows", "weight": r, "through_full0": f0,
                "ordered_profiles": ordered, "seconds": time.monotonic()-started}), flush=True)
    assert ordered == (fullmax+1)**2-int(r % 2 == 0)
    assert sum(run_counts) == comb(8192, r)
    for c in range(1, r+1):
        assert run_counts[c] == counts.run_total(4096, r, c)
    h, c_h, eta = anchor.source_factor(r)
    source = (v.exact(eta)*v.exact(core.BETA)**r*v.exact(z)**(-core.LIMIT))**v.exact(alpha)
    bound = source*sum(run_bounds, arb(0)); assert bound > 0
    logtwo = arb(2).log()
    bits = bound.log()/logtwo
    worst_run = max(range(1, r+1), key=lambda c: run_bounds[c].upper())
    selected_bits = (source*selected).log()/logtwo
    if r == 512:
        assert bits.overlaps(arb(*CONTROL_512))
        assert selected_bits.overlaps(arb(*CONTROL_223))
    report = {"weight": r, "log2_bound_interval": str(bits),
        "below_2_pow_minus_40": bool(bound < arb(2)**(-40)),
        "below_2_pow_minus_50": bool(bound < arb(2)**(-50)),
        "ordered_profiles": ordered, "all_indexed_support_counts_checked": True,
        "all_run_subtotals_checked": r, "source_factor": {"h": h, "c_h": str(c_h), "eta": str(eta)},
        "worst_ordered_profile": {"full": list(worst[:2]),
            "log2_bound_interval": str((source*worst[2]).log()/logtwo)},
        "worst_runs": worst_run, "worst_run_log2_interval": str((source*run_bounds[worst_run]).log()/logtwo),
        "selected_f223_log2_interval": str(selected_bits),
        "frozen512_control_replayed": r == 512,
        "seconds": time.monotonic()-started}
    return bound, report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rlo", type=int, default=508)
    parser.add_argument("--rhi", type=int, default=516)
    parser.add_argument("--seconds", type=int, choices=(900, 1200), default=900)
    parser.add_argument("--precision", type=int, choices=(256, 384), default=256)
    args = parser.parse_args()
    if not (480 <= args.rlo <= args.rhi <= 544 and args.rhi-args.rlo < 17):
        parser.error("require480<=rlo<=rhi<=544 and at most17 sizes")
    fullmax, singlemax = args.rhi//2, args.rhi-434
    if (fullmax+1)*(singlemax+1) > 32768:
        parser.error("high exact rectangle exceeds32768 coefficients per entry")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    started = time.monotonic(); deadline = started+args.seconds
    hybrid = load_hybrid(); old = hybrid.load_directed(); base = old.load_base()
    fine = base.load("cheap_lane_cauchy_fine.py", hybrid.FINE_PIN, "neighbor_exact")
    anchor = base.load("cheap_lane_cone_cauchy_anchors.py", base.ANCHOR_PIN, "neighbor_anchor")
    counts = base.load("cheap_lane_fast_profile_counts.py", base.COUNT_PIN, "neighbor_counts")
    screen = base.load("cheap_lane_cauchy_profile_screen.py", old.SCREEN_PIN, "neighbor_screen")
    counts.DEADLINE = deadline
    peak_memory = anchor.peak_memory_reader()

    def guard():
        if time.monotonic() > deadline:
            raise TimeoutError("neighbor cover time budget exceeded")
        if peak_memory() > 512*(1 << 20):
            raise MemoryError("512MiB neighbor cover budget exceeded")

    adapter = anchor.load_adapter(); core = adapter.load_verifier()
    core.Z, core.ALPHA = hybrid.Z, hybrid.ALPHA
    from flint import arb, ctx
    ctx.prec = args.precision

    class NeighborVerifier(core.Verifier):
        def guard(self):
            guard()

    v = NeighborVerifier(2, deadline)
    tests = {"mixed_exact_toy_coefficients": hybrid.mixed_tests(v, adapter),
        "rectangle_exact_toy_coefficients": fine.exact_tests(v, adapter),
        "cyclic_exact_toy_coefficients": base.exact_tests(v, adapter),
        "orientation_checks": old.orientation_tests(), **counts.checks()}
    h511, c511, _ = anchor.source_factor(511)
    h512, c512, eta512 = anchor.source_factor(512)
    assert (h511, c511) == (1, F(1, 512))
    assert (h512, c512, eta512) == (2, F(1, 3), F(43861, 7225))
    weights = tuple(range(args.rlo, args.rhi+1))
    schedule = hybrid.SCHEDULE[:-1]+((232, 217, fullmax, "exact"),)
    assert [f for _, lo, hi, _ in schedule for f in range(lo, hi+1)] == list(range(fullmax+1))
    regular = {r: [None]*(r//2+1) for r in weights}
    terminal = {r: [None]*(r//2+1) for r in weights}
    caps = v.local_caps(list(combinations_with_replacement(range(17), 4)))
    degrees, normalizers = base.local_degrees(v, caps, hybrid.ALPHA, F(1), F(1))
    del caps
    witnesses = []; cyclic_dc = truncated_dc = extractions = 0
    print(json.dumps({"phase": "shared_local_coefficients_complete", "weights": list(weights),
        "high_rectangle": [fullmax, singlemax], "seconds": time.monotonic()-started}), flush=True)
    for center, lo, hi, kind in schedule:
        guard()
        u, w = screen.marker_curve(center, old.CURVE_FULL, old.CURVE_SINGLE)
        up, wp = [v.exact(u)**k for k in range(33)], [v.exact(w)**k for k in range(33)]
        scaled = {(f, s): [up[f]*wp[s]*x for x in row] for (f, s), row in degrees.items()}
        plain = [sum((row[k] for row in scaled.values()), arb(0)) for k in range(6)]
        plain_p, plain_t = [[plain[0], plain[1]], [plain[2], plain[3]]], [[plain[4]], [plain[5]]]
        plain_r, plain_v = v.matpow(plain_p, 128), v.matmul(v.matpow(plain_p, 127), plain_t)
        if kind == "cyclic":
            full, single, stride = 127, 127, 255
            rr, tt = base.fold(v, scaled, 128, 128)
            product = lambda a, b: base.multiply(v, a, b, 128, 128)
        elif kind == "mixed":
            full, single, stride = hi, 127, 2*hi+1
            rr, tt = hybrid.mixed_pack(v, scaled, full, 128)
            product = lambda a, b: hybrid.mixed_product(v, a, b, full, 128)
        else:
            full, single, stride = fullmax, singlemax, 2*fullmax+1
            rr, tt = fine.rectangle_pack(v, scaled, full, single)
            product = lambda a, b: fine.rectangle_product(v, a, b, full, single)
        for _ in range(7):
            tt = product(rr, tt)
            rr = product(rr, rr)
        for table, evaluated in ((rr, plain_r), (tt, plain_v)):
            for i in range(2):
                for j in range(len(table[0])):
                    dc = sum((table[i][j][f+stride*s]
                              for s in range(single+1) for f in range(full+1)), arb(0))
                    if kind == "cyclic":
                        assert dc.overlaps(evaluated[i][j]); cyclic_dc += 1
                    else:
                        assert dc <= evaluated[i][j] or dc.overlaps(evaluated[i][j])
                        truncated_dc += 1
        band_targets = 0
        for r in weights:
            for f in range(max(lo, r-4096, 0), min(hi, r//2)+1):
                s = r-2*f
                ff = f % 128 if kind == "cyclic" else f
                ss = s if kind == "exact" else s % 128
                assert 0 <= ff <= full and 0 <= ss <= single
                index = ff+stride*ss
                denominator = arb(comb(4096, f)*comb(4096-f, s))*v.exact(u)**f*v.exact(w)**s
                assert denominator > 0 and regular[r][f] is None
                regular[r][f] = [[v.positive(item[index])/denominator for item in row] for row in rr]
                terminal[r][f] = [[v.positive(item[index])/denominator for item in row] for row in tt]
                band_targets += 1
        extractions += band_targets
        witnesses.append({"center": center, "first_full": lo, "last_full": hi, "method": kind,
            "full_limit_or_residue": full+1 if kind == "cyclic" else full,
            "single_limit_or_residue": single if kind == "exact" else single+1,
            "u": str(u), "v": str(w), "targets_extracted": band_targets})
        del rr, tt, table, evaluated, scaled, plain, plain_p, plain_t, plain_r, plain_v, product
        gc.collect(); guard()
        print(json.dumps({"phase": "shared_coefficient_band_complete", **witnesses[-1],
            "seconds": time.monotonic()-started, "peak_process_bytes": peak_memory()}), flush=True)
    assert extractions == sum(r//2+1 for r in weights)
    assert all(x is not None for r in weights for x in regular[r]+terminal[r])
    reports = []; total = new_total = arb(0)
    for r in weights:
        bound, report = join_weight(v, counts, anchor, core, r, regular[r], terminal[r],
                                    hybrid.ALPHA, hybrid.Z, guard, started)
        total += bound
        if r != 512:
            new_total += bound
        reports.append(report)
        print(json.dumps({"phase": "weight_complete", **report}), flush=True)
        del regular[r], terminal[r]
        guard()
    logtwo = arb(2).log()
    parameters = {"rlo": args.rlo, "rhi": args.rhi, "z": str(hybrid.Z),
        "alpha": str(hybrid.ALPHA), "distance_threshold": core.LIMIT, "shared_bands": witnesses}
    witness_hash = hashlib.sha256(json.dumps(parameters, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    print(json.dumps({"scope": "all actual supports and state paths at every explicitly listed weight",
        "directed": True, "all_profiles_each_weight_covered": True, "whole_code_certificate": False,
        "log2_band_union_interval": str(total.log()/logtwo),
        "band_below_2_pow_minus_40": bool(total < arb(2)**(-40)),
        "band_below_2_pow_minus_50": bool(total < arb(2)**(-50)),
        "log2_new_weights_excluding512_interval": str(new_total.log()/logtwo) if new_total > 0 else None,
        "weights": reports, "parameters": parameters,
        "tests": {**tests, "local_route_normalizers": normalizers,
            "cyclic_dc_overlaps": cyclic_dc, "truncated_dc_inequality_checks": truncated_dc,
            "own_weight_profile_normalizers": extractions, "source_factor_511_512_boundary": True,
            "frozen512_control_replayed": 512 in weights, "terminal_updates": 32767,
            "every_weight_all_run_and_source_count_sums_exact": True},
        "precision_bits": args.precision, "seconds": time.monotonic()-started,
        "budget_seconds": args.seconds, "peak_process_bytes": peak_memory(),
        "nonnegative_intersections": v.intersections, "generated_bytes": 0,
        "witness_sha256": witness_hash, "hybrid_pin": HYBRID_PIN,
        "directed_pin": hybrid.DIRECTED_PIN, "fine_pin": hybrid.FINE_PIN,
        "base_pin": old.BASE_PIN, "screen_pin": old.SCREEN_PIN,
        "anchor_pin": base.ANCHOR_PIN, "counts_pin": base.COUNT_PIN,
        "adapter_pin": anchor.ADAPTER_PIN, "core_pin": anchor.CORE_PIN,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2))


if __name__ == "__main__":
    main()
