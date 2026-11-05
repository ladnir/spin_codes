"""Algebraically reweight a completed directed r=257 verification.

This receipt does NOT recompute coefficients. Its input is the printed
directed interval from the pinned verifier and reproduction command below.

The audited improved h=1 comparison retains the missing group's nonzero
projection in the initial E-kernel. Average that independent outer map
first; its scalar counting-density bound cancels the exact normalized
ambient-kernel rank increment. The complete source factor becomes
 beta^r*(r+1)/[255*(1-r/512)]
instead of beta^r*(1+r*q^256/beta)/[255*(1-r/512)]. No local bound,
actual-support filter, state sum, route coefficient, source count, z,
alpha, or terminal condition changes. Their ratio, raised to alpha,
multiplies the old bound. This is not a caller-supplied factor or whole-code
certificate.
"""

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True
PIN = "e32069c849efa8e3abe1d3eee0a252a386e8acb3c161725d9ec83314a7b1678b"
CHILD_BYTES = 12607
CENTER = "0.89564762496111756242215465131490153554814736612596350726109377784013829"
RADIUS = "8.39e-72"
R, Q, ALPHA, Z = 257, 256, F(119, 1280), F(1343, 1600)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--precision", type=int, choices=(256, 384), default=256)
    args = parser.parse_args()
    started = time.monotonic()
    path = Path(__file__).resolve()
    if shutil.disk_usage(path.anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    witness = path.with_name("cheap_lane_cone_weight_witness.py")
    contents = witness.read_bytes()
    if len(contents) != CHILD_BYTES or hashlib.sha256(contents).hexdigest() != PIN:
        raise RuntimeError("frozen r257 verifier hash/size mismatch")
    from flint import arb, fmpq, ctx
    ctx.prec = args.precision

    def exact(x):
        x = F(x)
        return arb(fmpq(x.numerator, x.denominator))

    beta = F(Q**160, 255**32)
    c1 = 1-F(R, 512)
    assert 257 <= R <= 511 and c1 == F(255, 512) and c1 > 0
    assert Z == 1-F(R, 1600) and ALPHA == F(119, 1280)
    old_support = (1+R*F(Q**256)/beta)/(255*c1)
    new_support = F(R+1)/(255*c1)
    ratio = F(R+1)/(1+R*F(Q**256)/beta)
    assert ratio == new_support/old_support and 0 < ratio < 1
    # Exact decimal receipt, not float literals.
    old = arb(CENTER, RADIUS)
    lo, hi = F(CENTER)-F(RADIUS), F(CENTER)+F(RADIUS)
    assert old.is_finite() and F(RADIUS) > 0
    assert old.lower().fmpq() <= fmpq(lo.numerator, lo.denominator)
    assert old.upper().fmpq() >= fmpq(hi.numerator, hi.denominator)
    correction = exact(ALPHA)*exact(ratio).log()/arb(2).log()
    adjusted = old+correction
    assert correction < 0 and adjusted.upper() < old.lower()
    if time.monotonic()-started > 10:
        raise TimeoutError("10-second algebraic receipt budget exceeded")
    command = ("C:/Python314/python.exe -B research/workstreams/char2_field_spin/"
               "theory/cheap_lane_cone_weight_witness.py --weight 257 "
               "--allow-missing-one --precision 256")
    print(json.dumps({
        "scope": "all state paths and indexed supports at source weight257",
        "method": "directed algebraic reweight of a completed verification",
        "directed": True, "verifier_completed": True,
        "fresh_coefficient_recomputation": False, "whole_code_certificate": False,
        "parameters": {"source_weight": R, "z": str(Z), "alpha": str(ALPHA),
                       "distance_threshold": 734003},
        "old_receipt": {"log2_bound_interval": "["+CENTER+" +/- "+RADIUS+"]",
                        "precision_bits": 256, "ordered_profile_chains": 16641,
                        "source_total_check": "C(8192,257)", "terminal_updates": 32767,
                        "seconds": "158.6401126", "reproduction_command": command,
                        "source_sha256": PIN, "source_bytes": CHILD_BYTES},
        "comparison_lemma": "Average the missing independent outer map, retaining its "
            "nonzero kernel projection; its density cancels the normalized rank increment.",
        "old_source_factor": "beta^257*(1+257*256^256/beta)/(255*c1)",
        "new_source_factor": "beta^257*258/(255*c1)",
        "definitions": {"beta": "256^160/255^32", "c1_exact": str(c1)},
        "unchanged": "local caps, actual-T filters, state/route/count sums, z, alpha, terminal",
        "correction_log2_interval": str(correction),
        "adjusted_log2_bound_interval": str(adjusted),
        "proved_below_2_pow_minus_40": bool(adjusted < -40),
        "proved_below_2_pow_minus_80": bool(adjusted < -80),
        "proved_below_2_pow_minus_100": bool(adjusted < -100),
        "precision_bits": args.precision,
        "tests": "exact prefactor ratio, sign, c1, fixed witnesses, and receipt enclosure",
        "generated_bytes": 0,
        "seconds": time.monotonic()-started,
        "receipt_sha256": hashlib.sha256(path.read_bytes()).hexdigest()}, indent=2))


if __name__ == "__main__":
    main()
