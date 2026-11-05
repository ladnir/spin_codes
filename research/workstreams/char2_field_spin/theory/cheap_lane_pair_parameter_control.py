"""One explicitly selected z/alpha witness for the fixed r2048 profile.

Reuses the frozen12-proposal marker fit and m256/m512 controls. The
inside-coordinate coset refinement is enabled. Actual source multiplicity,
run512 mass, local moments and fractional powers are recomputed at the
supplied witnesses. No construction, source class, or distribution changes.
Every invocation remains selected-profile FLOAT evidence, not a certificate.
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
REFINEMENT_PIN = "8fd5db13b0aca7ebb4149772827985dfba322a5ca02e17ea6fd2931a8199bd06"


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
        parser.error("this bounded diagnostic requires1/10<=z<=9/10 and1/64<=alpha<=1/4")
    if max(n.bit_length() for x in (args.z,args.alpha) for n in (x.numerator,x.denominator)) > 64:
        parser.error("rational witnesses limited to64 bits")
    control = load("cheap_lane_pair_profile_calibration.py", CONTROL_PIN, "parameter_profile")
    refined = load("cheap_lane_coset_return_refinement.py", REFINEMENT_PIN, "parameter_coset")
    control.Z, control.ALPHA = args.z, args.alpha

    def augment(caps, z, alpha, guard):
        tests = refined.exact_tests(guard)
        result = refined.augment_caps(caps, z, alpha, guard)
        result["tests"]["coset_return_exact_tests"] = tests
        return result

    control.main(augment=augment, augmentation_label={
        "name": "explicit z/alpha witness with coset return bound",
        "z": str(args.z), "alpha": str(args.alpha),
        "control_pin": CONTROL_PIN, "refinement_pin": REFINEMENT_PIN,
        "wrapper_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})


if __name__ == "__main__":
    main()
