"""Matched selected-profile marker control with the coset return bound.

The frozen calibration owns its180s/512MiB/20GiB guards, actualr2048
source factor, fixed run512 mass, twelve-proposal budget, and m256/m512
replays. Only raw E and its local power change before rebuilding the
paired coefficients. This is not a run cover or a distance certificate.
"""

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
    control = load("cheap_lane_pair_profile_calibration.py", CONTROL_PIN, "coset_profile_calibration")
    refined = load("cheap_lane_coset_return_refinement.py", REFINEMENT_PIN, "coset_profile_refinement")

    def augment(caps, z, alpha, guard):
        tests = refined.exact_tests(guard)
        result = refined.augment_caps(caps, z, alpha, guard)
        result["tests"]["coset_return_exact_tests"] = tests
        return result

    control.main(augment=augment, augmentation_label={
        "name": "inside-coordinate coset return bound",
        "control_pin": CONTROL_PIN, "refinement_pin": REFINEMENT_PIN,
        "wrapper_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})


if __name__ == "__main__":
    main()
