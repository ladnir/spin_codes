"""Directed contiguous source-weight cover from shared cone coefficients.

Compute local caps once and extract every marginal (f,r-2f) for rlo..rhi
from downward-closed rectangles at total degree rhi. No interpolation or
monotonicity in r is assumed. Frozen geometry, caps, alpha-before-routing,
coefficient extraction and the output-only terminal remain unchanged.

Each cell contains integer (r,c), where c is cyclic run count. Scale a
regional matrix/terminal column at (r,f) by
 d=Omega_r^(alpha/256)*rho^(r/256)*sigma^(f/128).
For common-r pairs f0+f1=r-c, take separate entrywise maxima Bbar of
Rhat0*Rhat1 and Cbar of Rhat0*Vhat1. The cell bound is
 z^(-alpha L)*(Bbar^127*Cbar)[0,0]
   *sum_(r,c in cell) run_total(r,c)*rho^-r*sigma^(-(r-c)).
The envelopes may mix maximizing pairs, which only increases the positive
chain. They do not allow actual source weight to change. Degree-feasible
pairs of zero actual multiplicity may safely enlarge the envelopes.

Default: r64..80,z19/20,alpha1/8, singleton weight/run cells, unit marks.
Above256 use the audited JOINT missing-set factor, not ambient padding.
Optional whole-chain-max control isolates loss from the paired envelope.
180s cooperative budget; explicit240s is allowed only when rhi>256: the
frozen256 coefficient cover already took156s, and the larger downset and
all band marginals require additional work. At most33 weights, estimated
512MiB,20GiB free, bounded stdout, no generated files or numerical search.
"""

import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
from itertools import combinations
import json
from math import comb
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True
WEIGHT_PIN = "e32069c849efa8e3abe1d3eee0a252a386e8acb3c161725d9ec83314a7b1678b"
WEIGHT_BYTES = 12607


def load_weight():
    path = Path(__file__).with_name("cheap_lane_cone_weight_witness.py")
    data = path.read_bytes()
    if len(data) != WEIGHT_BYTES or hashlib.sha256(data).hexdigest() != WEIGHT_PIN:
        raise RuntimeError("frozen single-weight source hash/size mismatch")
    spec = importlib.util.spec_from_file_location("cell_weight_witness", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def support_factor(r):
    """Omega_r=beta^r*support_factor(r), before its fractional power."""
    if r <= 256:
        return F(256**2, 255**2*(257-r))
    h = r//256
    c_h = 1-F(r, 256*(h+1))
    assert 0 < h < r and h <= 256 and c_h > 0
    return (1+F(sum(comb(r, j) for j in range(1, h+1)), 255))/(255*c_h)


def pairs(r, c):
    total, maximum = r-c, r//2
    for f0 in range(max(0, total-maximum), min(maximum, total)+1):
        yield f0, total-f0


def rational_tests():
    """Small indexed-cycle union and exact scale cancellation, no floats."""
    def mm(a, b):
        return [[sum((a[i][k]*b[k][j] for k in range(len(b))), F(0))
                 for j in range(len(b[0]))] for i in range(len(a))]

    def mp(a, n):
        value = [[F(i == j) for j in range(2)] for i in range(2)]
        for _ in range(n):
            value = mm(value, a)
        return value

    def maximum(a, b):
        return [[max(x, y) for x, y in zip(row, other)]
                for row, other in zip(a, b)]

    def regular(r, f):
        return [[F(1+f, 20+r), F(2+r-f, 30+r)],
                [F(1+r-f, 40+r), F(2+f, 25+r)]]

    def terminal(r, f):
        return [[F(1+f, 11+r)], [F(3+r-f, 19+r)]]

    cases, supports = 0, 0
    for r in (2, 3, 4):
        hist = {}
        for selected in combinations(range(8), r):
            word = set(selected)
            f0 = sum(2*i in word and 2*i+1 in word for i in range(4))
            f1 = sum(2*i+1 in word and (2*i+2) % 8 in word for i in range(4))
            hist[f0, f1] = hist.get((f0, f1), 0)+1
            supports += 1
        for c in range(1, r+1):
            bb, cc = [[F(0)]*2 for _ in range(2)], [[F(0)], [F(0)]]
            weighted, number = F(0), 0
            for f0, f1 in pairs(r, c):
                a, b, t = regular(r, f0), regular(r, f1), terminal(r, f1)
                product, end = mm(a, b), mm(a, t)
                bb, cc = maximum(bb, product), maximum(cc, end)
                count = hist.get((f0, f1), 0)
                weighted += count*mm(mp(product, 3), end)[0][0]
                number += count
            envelope = mm(mp(bb, 3), cc)[0][0]
            assert weighted <= number*envelope
            # For8 toy regions choose Omega=2^8,rho=3^8,sigma=5^4.
            # Then d(r,f)=2*3^r*5^f and each pair scale is common in(r,c).
            scale = F(4*3**(2*r)*5**(r-c))
            scaled_b = [[scale*x for x in row] for row in bb]
            scaled_c = [[scale*x for x in row] for row in cc]
            got = mm(mp(scaled_b, 3), scaled_c)[0][0]
            assert got == 2**8*(3**8)**r*(5**4)**(r-c)*envelope
            # An additional candidate with zero source multiplicity is safe.
            larger_b = maximum(bb, [[2*x for x in row] for row in bb])
            assert mm(mp(larger_b, 3), cc)[0][0] >= envelope
            cases += 1
        assert sum(hist.values()) == comb(8, r)
    assert support_factor(256) == F(256**2, 255**2)
    assert support_factor(257) == (1+F(257, 255))/(255*F(255, 512))
    return {"rational_cell_unions": cases, "indexed_toy_supports": supports,
            "perfect_power_scaling_and_zero_count_enlargement": "passed"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rlo", type=int, default=64)
    parser.add_argument("--rhi", type=int, default=80)
    parser.add_argument("--z", default="19/20")
    parser.add_argument("--alpha", default="1/8")
    parser.add_argument("--rho", default="1")
    parser.add_argument("--sigma", default="1")
    parser.add_argument("--weight-cell", type=int, choices=(1, 2, 4), default=1)
    parser.add_argument("--run-cell", type=int, choices=(1, 2, 4), default=1)
    parser.add_argument("--control-weight", type=int)
    parser.add_argument("--precision", type=int, choices=(256, 384, 512), default=256)
    parser.add_argument("--seconds", type=int, choices=(180, 240), default=180)
    args = parser.parse_args()
    if not 7 <= args.rlo <= args.rhi <= 384 or args.rhi-args.rlo > 32:
        parser.error("require7<=rlo<=rhi<=384 and at most33 weights")
    if args.seconds == 240 and args.rhi <= 256:
        parser.error("240seconds is reserved for a downset extending above256")
    control = args.control_weight
    if control is None and args.rlo <= 64 <= args.rhi:
        control = 64
    if control is not None and not args.rlo <= control <= args.rhi:
        parser.error("control weight must be in the requested band")
    z, alpha, rho, sigma = map(F, (args.z, args.alpha, args.rho, args.sigma))
    if not 0 < z < 1 or not 0 < alpha <= 1 or min(rho, sigma) <= 0:
        parser.error("require0<z<1,0<alpha<=1,and positive rho,sigma")
    if max(n.bit_length() for x in (z, alpha, rho, sigma)
           for n in (x.numerator, x.denominator)) > 64:
        parser.error("rational numerators/denominators must fit64 bits")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    started = time.monotonic(); deadline = started+args.seconds

    def progress(phase, **values):
        print(json.dumps({"phase": phase, "seconds": time.monotonic()-started, **values}), flush=True)

    weight = load_weight()
    cover = weight.load_cover()
    adapter = cover.load("cheap_lane_cone_interval_tests.py", cover.TEST_PIN, "cell_interval_tests")
    core = adapter.load_verifier()
    assert adapter.PIN == cover.CORE_PIN
    core.Z, core.ALPHA, core.MAX_FULL = z, alpha, args.rhi//2
    packed = cover.load("cheap_lane_cone_packed.py", cover.PACK_PIN, "cell_packed")
    counts = cover.load("cheap_lane_fast_profile_counts.py", cover.COUNT_PIN, "cell_counts")
    counts.DEADLINE = deadline
    from flint import arb, ctx
    ctx.prec = args.precision

    class CellVerifier(weight.make_verifier(core, args.rhi)):
        def guard(self):
            if time.monotonic() > self.deadline:
                raise TimeoutError("bounded shared-coefficient cell budget exceeded")

    tests = rational_tests()
    toy = packed.tiny_tests(core, CellVerifier, adapter.exact_point, deadline)
    count_checks = counts.checks()
    v = CellVerifier(min(32, args.rhi), deadline)
    records, route_checks = core.route_records(min(32, args.rhi))
    records = {key: val for key, val in records.items() if 2*key[0]+key[1] <= args.rhi}
    caps = v.local_caps(sorted({key[2] for key in records}))
    powered = {occ: tuple(v.positive(x)**v.exact(alpha) for x in cap)
               for occ, cap in caps.items()}
    del caps
    progress("local_caps_ready", records=len(records), shapes=len(powered))
    rr, tt, bands = {}, {}, []
    intersections = v.intersections
    for low, high, uscale, vscale in weight.rectangles(args.rhi):
        core.USCALE, core.VSCALE = uscale, vscale
        block = CellVerifier(args.rhi-2*low, deadline); block.packed_full = high
        regular, terminal = block.regions(*cover.local_polynomials(block, core, records, powered))
        for r in range(args.rlo, args.rhi+1):
            for f in range(low, min(high, r//2)+1):
                rr[r, f] = block.extract(regular, f, r-2*f)
                tt[r, f] = block.extract(terminal, f, r-2*f)
        intersections += block.intersections
        bands.append({"full_min": low, "full_max": high, "singleton_max": block.single,
                      "full_scale": str(uscale), "singleton_scale": str(vscale)})
        progress("band_ready", full_min=low, full_max=high, marginal_profiles=len(rr))
        del regular, terminal
    required = {(r, f) for r in range(args.rlo, args.rhi+1) for f in range(r//2+1)}
    assert set(rr) == set(tt) == required
    beta, arho, asigma = v.exact(core.BETA), v.exact(rho), v.exact(sigma)
    omega, rhat, that, run_numbers = {}, {}, {}, {}
    for r in range(args.rlo, args.rhi+1):
        omega[r] = v.exact(support_factor(r))*beta**r
        common = omega[r]**v.exact(alpha/256)*arho**v.exact(F(r, 256))
        for f in range(r//2+1):
            scale = common*asigma**v.exact(F(f, 128))
            rhat[r, f] = [[scale*x for x in row] for row in rr[r, f]]
            that[r, f] = [[scale*x for x in row] for row in tt[r, f]]
        run_numbers[r] = {c: counts.run_total(4096, r, c) for c in range(1, r+1)}
        assert sum(run_numbers[r].values()) == comb(8192, r)
        v.guard()
    zfactor = v.exact(z)**v.exact(-alpha*core.LIMIT)
    logtwo = arb(2).log()

    def envelope_update(left, right):
        return [[x.max(y) for x, y in zip(row, other)] for row, other in zip(left, right)]

    total, cell_count, pair_count = arb(0), 0, 0
    covered, weight_reports, worst = set(), [], None
    control_envelope = arb(0)
    for rstart in range(args.rlo, args.rhi+1, args.weight_cell):
        rend = min(args.rhi, rstart+args.weight_cell-1)
        subtotal = arb(0)
        for cstart in range(1, rend+1, args.run_cell):
            cend = min(rend, cstart+args.run_cell-1)
            bb, cc = [[arb(0)]*2 for _ in range(2)], [[arb(0)], [arb(0)]]
            count_factor, points = arb(0), []
            for r in range(rstart, rend+1):
                for c in range(cstart, min(cend, r)+1):
                    points.append((r, c)); assert (r, c) not in covered
                    covered.add((r, c))
                    count_factor += v.exact(F(run_numbers[r][c])*rho**(-r)*sigma**(-(r-c)))
                    for f0, f1 in pairs(r, c):
                        bb = envelope_update(bb, v.matmul(rhat[r, f0], rhat[r, f1]))
                        cc = envelope_update(cc, v.matmul(rhat[r, f0], that[r, f1]))
                        pair_count += 1
                    v.guard()
            assert points
            value = zfactor*count_factor*v.positive(v.matmul(v.matpow(bb, 127), cc)[0][0])
            total += value; subtotal += value; cell_count += 1
            if args.weight_cell == 1 and rstart == control:
                control_envelope += value
            if worst is None or value.upper() > worst[0].upper():
                worst = value, (rstart, rend), (cstart, cend)
        weight_reports.append({"weights": [rstart, rend], "log2_bound_interval": str(subtotal.log()/logtwo)})
        progress("weight_cells_ready", through_weight=rend, cells=cell_count)
    assert covered == {(r, c) for r in range(args.rlo, args.rhi+1) for c in range(1, r+1)}
    controls = None
    if control is not None:
        r, raw_total = control, arb(0)
        for c in range(1, r+1):
            largest = arb(0)
            for f0, f1 in pairs(r, c):
                largest = largest.max(packed.ordered_chain(v, rr[r, f0], rr[r, f1], tt[r, f1]))
            raw_total += run_numbers[r][c]*largest
            v.guard()
        bound = zfactor*omega[r]**v.exact(alpha)*raw_total
        controls = {"source_weight": r, "kind": "run_total times maximum COMPLETE ordered chain",
                    "log2_bound_interval": str(bound.log()/logtwo),
                    "exact_profile_count_weighting": False}
        if args.weight_cell == 1:
            controls["paired_envelope_log2_interval"] = str(control_envelope.log()/logtwo)
            controls["envelope_minus_control_bits"] = str((control_envelope.log()-bound.log())/logtwo)
    v.guard()
    print(json.dumps({"scope": "all state paths and indexed supports in the entire stated source-weight band",
        "directed": True, "verifier_completed": True, "whole_code_certificate": False,
        "parameters": {"rlo": args.rlo, "rhi": args.rhi, "z": str(z), "alpha": str(alpha),
                       "rho": str(rho), "sigma": str(sigma), "weight_cell": args.weight_cell,
                       "run_cell": args.run_cell, "distance_threshold": core.LIMIT, "bands": bands},
        "source_comparison": "eta_r*beta^r through256; above256 joint missing-set Omega with h=floor(r/256)",
        "missing_factor": "[1+sum(j=1..h)C(r,j)/255]/[255*(1-r/(256*(h+1)))]",
        "log2_bound_interval": str(total.log()/logtwo),
        "proved_joint_below_2_pow_minus_40": bool(total < arb(2)**(-40)),
        "proved_joint_below_2_pow_minus_100": bool(total < arb(2)**(-100)),
        "weight_cell_bounds": weight_reports,
        "worst_cell": {"weights": worst[1], "runs": worst[2], "selection": "largest upper endpoint",
                       "log2_bound_interval": str(worst[0].log()/logtwo)},
        "whole_chain_control": controls,
        "tests": {**tests, **count_checks, "packed_exact_tests": toy,
                  "route_normalizers_before_downset": route_checks, "marginal_profiles": len(required),
                  "run_total_sum_every_weight": "C(8192,r)", "disjoint_cell_points": len(covered),
                  "degree_feasible_ordered_pairs": pair_count, "terminal_updates": 32767},
        "precision_bits": args.precision, "seconds": time.monotonic()-started,
        "runtime_budget_seconds": args.seconds, "estimated_peak_bytes": 512*(1 << 20),
        "nonnegative_intersections": intersections+v.intersections, "generated_bytes": 0,
        "weight_witness_pin": WEIGHT_PIN, "core_pin": cover.CORE_PIN, "counts_pin": cover.COUNT_PIN,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2), flush=True)


if __name__ == "__main__":
    main()
