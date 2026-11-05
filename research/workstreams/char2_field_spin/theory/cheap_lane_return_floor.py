"""A positive block-floor refinement of the fixed-cone return bound.

Retain the pinned disjoint-support argument and the same global cone.
For one RS emission block, let ell be its allowed input count, c=min(4,ell),
d=max(0,13-ell), and Q_ell its EXACT quiet moment. On a fixed nonzero
syndrome, the inside-input output is an affine coset of ker(C_I). For
0<=z<=1 the zero coset maximizes its weight moment. The nonzero RS word
has at least d nonzero coordinates outside I. Thus its nonquiet moment
is bounded by
  D_ell=min(g^ell-Q_ell, (q^c-1) Q_ell z^d), g=z+(1-z)/q.

Combine four blocks positively: (Q,R)->(Q*Q_ell,
R*(Q_ell+D_ell)+Q*D_ell). R counts at least one active block, with no
subtraction of four-block products. The return cap is kappa*q^-16 R,
also capped by the previous independently valid cone return bound.
No new auxiliary independence assumption or change to the encoder.

Fixed z=.85,alpha=.05; FLOAT diagnostic for consecutive256 source groups,
not a directed or whole-code certificate. One48-second job,20GiB free,
128MiB estimated peak, stdout only. No search or generated data files.
"""

from fractions import Fraction
import hashlib
import importlib.util
import itertools
import json
import math
import os
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
PIN = "a3eff28456bad943f64932c8f10219ad316d9baede9220c3fcc380d080361abf"
DUAL_PIN = "84b6b757a3fc4e85f28c0c3e9f59f6063ae8daf70cc3704aca6feb9ebe27600b"


def load_disjoint():
    path = Path(__file__).with_name("cheap_lane_disjoint_return.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != PIN:
        raise RuntimeError("frozen disjoint-return hash mismatch")
    spec = importlib.util.spec_from_file_location("floor_disjoint_return", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def gf4_mul(a, b):
    out = 0
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
        if a & 4:
            a ^= 7  # x^2+x+1, not integer arithmetic modulo4.
    return out


def positive_nonquiet(blocks):
    quiet, active = Fraction(1), Fraction(0)
    for q, d in blocks:
        active, quiet = active*(q+d)+quiet*d, quiet*q
    return quiet, active


def exact_tests():
    """RS[4,2] over GF4 has C=A^T and CA=0; enumerate all cosets."""
    q = 4
    coset_checks = moment_checks = 0
    # Points0..3 are the entire GF4 field for this toy only. The actual
    # GF256 encoder uses its own frozen basis; no GF16 substitution occurs.
    for i, j in itertools.product(range(2), repeat=2):
        value = 0
        for point in range(4):
            value ^= gf4_mul(1 if i == 0 else point, 1 if j == 0 else point)
        assert value == 0
    for z in (Fraction(1, 3), Fraction(2, 3)):
        g = z+(1-z)/q
        for mask in range(16):
            allowed = [i for i in range(4) if mask >> i & 1]
            ell = len(allowed); c = min(2, ell); floor = max(0, 3-ell)
            cosets = {}
            direct_quiet = direct_active = Fraction(0)
            for values in itertools.product(range(q), repeat=ell):
                x = [0]*4
                for i, value in zip(allowed, values):
                    x[i] = value
                a = b = 0
                for i in allowed:
                    a ^= x[i]; b ^= gf4_mul(i, x[i])
                syndrome = a, b
                cosets[syndrome] = cosets.get(syndrome, Fraction(0))+z**sum(v != 0 for v in values)/q**ell
                y = [x[i] ^ a ^ gf4_mul(i, b) for i in range(4)]
                value = z**sum(v != 0 for v in y)/q**ell
                if syndrome == (0, 0):
                    direct_quiet += value
                else:
                    direct_active += value
            quiet = cosets[(0, 0)]
            assert len(cosets) == q**c
            assert all(value <= quiet for value in cosets.values())
            coset_checks += len(cosets)
            formula = sum(Fraction(math.comb(ell, k))*z**k*((1-z)/q)**(ell-k)
                          *Fraction(1, q**min(2, k)) for k in range(ell+1))
            assert quiet == direct_quiet == formula
            cap = min(g**ell-quiet, (q**c-1)*quiet*z**floor)
            assert 0 <= direct_active <= cap
            moment_checks += 1
    blocks = [(Fraction(1, 5), Fraction(2, 7)), (Fraction(1, 3), Fraction(1, 4)),
              (Fraction(1, 8), Fraction(3, 11)), (Fraction(2, 9), Fraction(1, 6))]
    quiet, active = positive_nonquiet(blocks)
    explicit = Fraction(0)
    for choice in itertools.product((0, 1), repeat=4):
        if not any(choice):
            continue
        product = Fraction(1)
        for bit, block in zip(choice, blocks):
            product *= block[bit]
        explicit += product
    assert active == explicit
    assert quiet == math.prod(block[0] for block in blocks)
    assert quiet+active == math.prod(sum(block) for block in blocks)
    assert positive_nonquiet([(Fraction(1), Fraction(0))]*4) == (1, 0)
    return {"exact_GF4_affine_coset_comparisons": coset_checks,
            "exact_GF4_nonquiet_floor_bounds": moment_checks,
            "exact_positive_OR_subsets": 15, "empty_block_terminal_case": "passed"}


def floor_caps(records, z, disjoint, np):
    qb = np.asarray(disjoint.quiet_blocks(z))
    g = z+(1-z)/256
    ell = np.arange(17)
    generic = np.power(g, ell)-qb
    assert np.all(generic >= -1e-15)
    generic = np.maximum(0., generic)
    floor = (np.power(256., np.minimum(4, ell))-1)*qb*np.power(z, np.maximum(0, 13-ell))
    db = np.minimum(generic, floor)
    assert qb[0] == 1. and db[0] == 0.
    occ = np.asarray([record[0] for record in records], dtype=int)
    quiet = np.ones(len(records)); active = np.zeros(len(records))
    for b in range(4):
        q, d = qb[occ[:, b]], db[occ[:, b]]
        active, quiet = active*(q+d)+quiet*d, quiet*q
    kappa, _, _ = disjoint.kappa_exact()
    cap = float(kappa)*256.**-16*active
    old_cap, old_quiet, _ = disjoint.return_caps(records, z, np)
    assert np.allclose(quiet, old_quiet, rtol=3e-14, atol=0.)
    assert np.all(cap <= old_cap*(1+3e-14)+1e-300)
    return cap, quiet, kappa, qb, db


def tilted_counts(a, b, result, np, lse, deadline, pilot):
    path = Path(__file__).with_name("cheap_lane_branch_sensitivity.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != DUAL_PIN:
        raise RuntimeError("frozen positive-dual helper hash mismatch")
    spec = importlib.util.spec_from_file_location("floor_return_duals", path)
    dual = importlib.util.module_from_spec(spec); spec.loader.exec_module(dual)
    dual.np, dual.lse, dual.DEADLINE = np, lse, deadline
    checks = dual.tests(pilot)
    # These route coefficients ALREADY contain local cap**alpha. Seed
    # ordinary binary path counts here, with no second alpha application.
    a, b = dual.seed(a), dual.seed(b, True)
    for _ in range(7):
        b, a = dual.convolution(a, b), dual.convolution(a, a)
    regions, tails = [], []
    for f, s in ((128, 0), (127, 2)):
        normalizer = math.log(math.comb(4096, f)*math.comb(4096-f, s))
        regions.append((a[0][s, f]-normalizer, a[1][s, f]))
        tails.append((b[0][s, f]-normalizer, b[1][s, f]))
    value, means, _ = dual.chain(regions, tails)
    score = result["source_bits_per_orientation"]+value/math.log(2)
    assert abs(score-result["log2_bound"]) < 1e-6
    identity = float(means[0]+means[1]+means[2]-means[3])
    assert abs(identity-32767) < 1e-6
    return {"log2_bound": score, "births": float(means[0]),
            "quiet_regular_steps": float(means[1]), "active_incoming_packets": float(means[2]),
            "terminal_live_probability": float(means[3]),
            "quiet_plus_births_plus_active_minus_terminal": identity,
            "interpretation": "tilted coarse-majorant counts, not actual encoder measurements",
            "tests": {**checks, "scalar_dual_agreement": "passed", "count_identity": "passed"}}


def main():
    if len(sys.argv) != 1:
        raise ValueError("fixed return-floor diagnostic takes no scan arguments")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    started = time.monotonic()
    import numpy as np
    from scipy.special import logsumexp
    disjoint = load_disjoint(); cone = disjoint.load_cone()
    tests = exact_tests(); tests.update(disjoint.exact_tests()); tests.update(cone.exact_tests())
    degree, pilot = cone.load_degree(np, logsumexp, started+48)
    tests.update(degree.tests(pilot))
    singleton, filtered = pilot.imports()
    shapes, index, degrees, weights, records, route_checks = pilot.geometry()
    m, t, refinements = pilot.raw_tables(shapes, index, records, cone.Z, singleton, filtered)
    p = np.asarray([math.comb(16, j)*min(1., 256/255*255.**(-(16-j))) for j in range(1, 17)])
    coarse, terminal = cone.cone_caps(m, t, p, np)
    old_cap, old_quiet, _ = disjoint.return_caps(records, cone.Z, np)
    coarse[:, 1, 0] = np.minimum(coarse[:, 1, 0], old_cap)
    prior = cone.score_polys(*cone.route_polys(coarse, terminal, degrees, weights, cone.ALPHA, pilot), degree, pilot)
    assert abs(prior["log2_bound"]-3736.1967) < .002
    cap, quiet, kappa, qb, db = floor_caps(records, cone.Z, disjoint, np)
    assert np.allclose(quiet, m[:, 0, 0], rtol=3e-14, atol=0.)
    old_return = coarse[:, 1, 0].copy()
    coarse[:, 1, 0] = np.minimum(old_return, cap)
    polynomial, final = cone.route_polys(coarse, terminal, degrees, weights, cone.ALPHA, pilot)
    result = cone.score_polys(polynomial, final, degree, pilot)
    assert result["log2_bound"] <= prior["log2_bound"]+1e-7
    counts = tilted_counts(polynomial, final, result, np, logsumexp, started+48, pilot)
    degree.guard()
    print(json.dumps({"scope": "FLOAT positive block-floor return cone, consecutive256-group supports",
        "whole_code_certificate": False, "directed": False,
        "all_state_masks_covered": True, "all_source_supports_covered": False,
        "parameters": {"z": cone.Z, "alpha": cone.ALPHA, "path_marks": 1},
        "prior_disjoint_return": prior, "block_floor_return": result,
        "improvement_bits": prior["log2_bound"]-result["log2_bound"],
        "return_records_improved": int(np.count_nonzero(cap < old_return)),
        "kappa_exact": str(kappa), "quiet_block_moments": qb.tolist(),
        "nonquiet_block_caps": db.tolist(),
        "tilted_binary_path_counts": counts,
        "tests": {**tests, "exact_route_normalizers": route_checks,
            "quiet_vs_raw_table": "passed", "prior_control": "passed",
            "pointwise_return_cap_monotonicity": "passed"},
        "pair_only_refinements": refinements, "seconds": time.monotonic()-started,
        "estimated_peak_bytes": 128*(1 << 20), "generated_bytes": 0,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "disjoint_return_pin": PIN, "positive_dual_pin": DUAL_PIN}, indent=2))


if __name__ == "__main__":
    main()
