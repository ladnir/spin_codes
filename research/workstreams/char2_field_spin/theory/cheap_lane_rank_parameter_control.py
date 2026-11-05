"""Selected r2048 parameter control with the actual residual-rank cap.

Recompute source, local moments and pair corrections at explicit rational
z/alpha witnesses. Keep the frozen profile, indexed run mass, marker fit,
terminal, and FFT controls. FLOAT diagnostic only; not a certificate.
Delegated limits:180s/512MiB/20GiB, stdout only, no generated files.
"""

import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
from pathlib import Path
import sys

sys.dont_write_bytecode = True
CONTROL_PIN = "5b5e18d2b62b8dbbf5614b5b07c258b4aa1fe6e53a554513518906c931212518"
RANK_PIN = "3efa28f205a01cafe47a8058ca4d3a69e065d5e22e75ea00a8fe79ed4a5dcf7c"


def load(name, digest, label):
    path = Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
        raise RuntimeError("frozen source mismatch: " + name)
    spec = importlib.util.spec_from_file_location(label, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--z", type=F, required=True)
    parser.add_argument("--alpha", type=F, default=F(119,1280))
    args = parser.parse_args()
    if not F(1,10) <= args.z <= F(9,10) or not F(1,64) <= args.alpha <= F(1,4):
        parser.error("require 1/10<=z<=9/10 and 1/64<=alpha<=1/4")
    if max(n.bit_length() for x in (args.z,args.alpha) for n in (x.numerator,x.denominator)) > 64:
        parser.error("rational witnesses limited to 64 bits")
    control = load("cheap_lane_pair_profile_calibration.py", CONTROL_PIN, "rank_parameter_profile")
    rank = load("cheap_lane_residual_rank_return.py", RANK_PIN, "rank_parameter_return")
    control.Z, control.ALPHA = args.z, args.alpha

    def augment(caps, z, alpha, guard):
        tests = rank.exact_tests(guard)
        result = rank.augment_caps(caps, z, alpha, guard)
        result["tests"]["residual_rank_exact_tests"] = tests
        return result

    control.main(augment=augment, augmentation_label={
        "name": "explicit z/alpha witness with actual residual-rank return bound",
        "z": str(args.z), "alpha": str(args.alpha),
        "control_pin": CONTROL_PIN, "rank_pin": RANK_PIN,
        "wrapper_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})


if __name__ == "__main__":
    main()
