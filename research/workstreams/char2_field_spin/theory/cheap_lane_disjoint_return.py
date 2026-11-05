"""FLOAT actual-support return bound for the fixed live-state cone.

At return u=Cx. First enforce supp(u)=H in the scalar comparison,
paying q/(q-|H|), then sum the DISJOINT supports H. For the fixed cone
f_h=p_h/C(16,h), the resulting cap is
  kappa*q^-16 E_x[z^wt((I+AC)x) 1(Cx!=0)],
  kappa=max_h f_h*q^(16-h)*q/(q-h)=q^3/((q-1)^2*(q-15)).
There is no extra nonzero-coordinate conditioning factor. Since CA=0,
I+AC is invertible. Its restriction to ell allowed input coordinates is
injective, so the expectation is at most g^ell-Q_I(z). Q_I is the EXACT
quiet moment, not an upper bound that may be subtracted.

This bounds the actual-mask cone return, not necessarily the old
cardinality-max return array. Both independently bound that cone sum;
their minimum is valid. Keep old activation/live/terminal caps and the
same fixed cone vector throughout. No encoder or setup change.
Fixed z=.85, alpha=.05; all state masks, consecutive256-group source
supports only. FLOAT diagnostic, NOT a directed or whole-code certificate.
One48-second invocation;20GiB free;128MiB estimated peak; stdout only.
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
PIN = "f450d2bef7520363011e68e7600e9324b176a8fb20bc9350d059b00e2cdfc1d2"


def load_cone():
    path = Path(__file__).with_name("cheap_lane_live_cone.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != PIN:
        raise RuntimeError("frozen cone source hash mismatch")
    spec = importlib.util.spec_from_file_location("disjoint_cone", path)
    cone = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cone)
    return cone


def kappa_exact(q=256, n=16):
    per_mask = [min(Fraction(1), Fraction(q, q-1)*Fraction(1, (q-1)**(n-h)))
                for h in range(1, n+1)]
    factors = [f*q**(n-h)*Fraction(q, q-h) for h, f in enumerate(per_mask, 1)]
    return max(factors), per_mask, factors


def exact_tests():
    # A duplicates each of two coordinates; C=A^T, over GF4. Thus CA=0
    # and S=I+AC is an involution. Enumerate every allowed input set,
    # input and incoming scalar state. Coefficients are binary, so GF4
    # addition alone suffices; no prime-field substitution.
    q, n = 4, 2
    kappa, p, _ = kappa_exact(q, n)
    checked = 0
    for z in (Fraction(1, 3), Fraction(2, 3)):
        g = z+(1-z)/q
        for bits in range(16):
            allowed = [i for i in range(4) if bits >> i & 1]
            ell = len(allowed)
            direct = Fraction(0)
            eliminated = Fraction(0)
            quiet = Fraction(0)
            active = Fraction(0)
            total = Fraction(0)
            for values in itertools.product(range(q), repeat=ell):
                x = [0]*4
                for i, v in zip(allowed, values):
                    x[i] = v
                cx = (x[0] ^ x[2], x[1] ^ x[3])
                y = [x[i] ^ cx[i % 2] for i in range(4)]
                cy = (y[0] ^ y[2], y[1] ^ y[3])
                assert [y[i] ^ cy[i % 2] for i in range(4)] == x
                weight = z**sum(v != 0 for v in y)/q**ell
                total += weight
                support = sum((v != 0) << i for i, v in enumerate(cx))
                if not support:
                    quiet += weight
                else:
                    active += weight
                for mask in range(1, 1 << n):
                    indices = [i for i in range(n) if mask >> i & 1]
                    h = len(indices)
                    correction = Fraction(q, q-h)*p[h-1]
                    if support == mask:
                        eliminated += correction*weight/q**h
                    for uvals in itertools.product(range(1, q), repeat=h):
                        u = [0]*n
                        for i, v in zip(indices, uvals):
                            u[i] = v
                        if tuple(u) == cx:
                            # Denominator is q^h, NOT (q-1)^h.
                            direct += correction*weight/q**h
            assert direct == eliminated
            assert total == quiet+active and total <= g**ell
            assert eliminated <= kappa*active/q**n <= kappa*(g**ell-quiet)/q**n
            checked += 1
    k, _, factors = kappa_exact()
    assert k == Fraction(256**3, 255**2*241)
    assert factors.index(k)+1 == 15
    return {"exact_GF4_disjoint_support_and_injective_moment_tests": checked,
            "exact_kappa_maximizer": 15}


def quiet_blocks(z, q=256):
    # Positive requested-zero expansion for an [16,12] MDS kernel.
    return [sum(math.comb(ell, k)*z**k*((1-z)/q)**(ell-k)
                *q**(-min(4, k)) for k in range(ell+1)) for ell in range(17)]


def return_caps(records, z, np):
    kappa, _, _ = kappa_exact()
    qb = np.asarray(quiet_blocks(z))
    occ = np.asarray([record[0] for record in records], dtype=int)
    ell = occ.sum(axis=1)
    quiet = np.prod(qb[occ], axis=1)
    g = z+(1-z)/256
    active = np.power(g, ell)-quiet
    assert np.all(active >= -1e-15)
    cap = float(kappa)*256.**-16*np.maximum(0., active)
    cap[ell == 0] = 0.
    return cap, quiet, kappa


def main():
    if len(sys.argv) != 1:
        raise ValueError("fixed diagnostic accepts no scan arguments")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free")
    started = time.monotonic()
    import numpy as np
    from scipy.special import logsumexp
    cone = load_cone()
    tests = exact_tests()
    tests.update(cone.exact_tests())
    degree, pilot = cone.load_degree(np, logsumexp, started+48)
    tests.update(degree.tests(pilot))
    singleton, filtered = pilot.imports()
    shapes, index, degrees, weights, records, route_checks = pilot.geometry()
    m, t, refinements = pilot.raw_tables(shapes, index, records, cone.Z, singleton, filtered)
    p = np.asarray([math.comb(16, j)*min(1., 256/255*255.**(-(16-j))) for j in range(1, 17)])
    coarse, terminal = cone.cone_caps(m, t, p, np)
    old = cone.score_polys(*cone.route_polys(coarse, terminal, degrees, weights, cone.ALPHA, pilot), degree, pilot)
    assert abs(old["log2_bound"]-4312.541523826458) < 1e-6
    cap, quiet, kappa = return_caps(records, cone.Z, np)
    assert np.allclose(quiet, m[:, 0, 0], rtol=3e-14, atol=0.)
    old_return = coarse[:, 1, 0].copy()
    coarse[:, 1, 0] = np.minimum(old_return, cap)
    result = cone.score_polys(*cone.route_polys(coarse, terminal, degrees, weights, cone.ALPHA, pilot), degree, pilot)
    assert result["log2_bound"] <= old["log2_bound"]+1e-7
    degree.guard()
    print(json.dumps({"scope": "FLOAT disjoint-support return cone, consecutive256-group supports",
        "whole_code_certificate": False, "directed": False,
        "all_state_masks_covered": True, "all_source_supports_covered": False,
        "parameters": {"z": cone.Z, "alpha": cone.ALPHA, "path_marks": 1},
        "kappa_exact": str(kappa), "kappa": float(kappa),
        "old_cone_control": old, "disjoint_support_return": result,
        "improvement_bits": old["log2_bound"]-result["log2_bound"],
        "return_records_improved": int(np.count_nonzero(cap < old_return)),
        "tests": {**tests, "exact_route_normalizers": route_checks,
            "exact_quiet_formula_vs_raw_table": "passed", "old_cone_control": "passed"},
        "pair_only_refinements": refinements, "seconds": time.monotonic()-started,
        "estimated_peak_bytes": 128*(1 << 20), "generated_bytes": 0,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), "cone_pin": PIN}, indent=2))


if __name__ == "__main__":
    main()
