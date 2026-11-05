"""Matched FLOAT controls for the actual residual-rank return bound.

Default: selectedr2048 profile768,768 with the frozen marker calibration.
--complete delegates all remaining arguments to the frozen complete-size
screen, for example --complete --weight2048. Both modes preserve source,
routing, terminal and guards; only E changes before paired reconstruction.
No frozen edits, generated files, or directed/global certificate claims.
"""

import hashlib
import importlib.util
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
CONTROL_PIN = "5b5e18d2b62b8dbbf5614b5b07c258b4aa1fe6e53a554513518906c931212518"
SIZES_PIN = "7a7024b78912abb78eaaef9e90f766eefdd7e7f807d8dfe6143f38c86c233afd"
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
    rank = load("cheap_lane_residual_rank_return.py", RANK_PIN, "rank_control_refinement")

    def augment(caps, z, alpha, guard):
        tests = rank.exact_tests(guard)
        result = rank.augment_caps(caps, z, alpha, guard)
        result["tests"]["residual_rank_exact_tests"] = tests
        result["block_report"]["residual_rank_return"]["source_sha256"] = RANK_PIN
        return result

    provenance = {"name": "actual tensor-row residual-rank return bound",
        "rank_pin": RANK_PIN,
        "wrapper_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    if "--complete" not in sys.argv:
        if len(sys.argv) != 1:
            raise ValueError("selected mode has a frozen witness; use --complete for size-screen arguments")
        control = load("cheap_lane_pair_profile_calibration.py", CONTROL_PIN, "rank_selected_control")
        control.main(augment=augment, augmentation_label=provenance)
        return

    sys.argv.remove("--complete")
    sizes = load("cheap_lane_pair_sizes_screen.py", SIZES_PIN, "rank_complete_control")
    original_load = sizes.load

    def with_refinement(name, pin, label):
        module = original_load(name, pin, label)
        if name == "cheap_lane_pair_bank_screen.py":
            original_composition = module.load_composition

            def composition():
                result = original_composition()
                original_build = result.build_caps

                def build(z, alpha, triple, anchor, guard, deadline):
                    return augment(original_build(z, alpha, triple, anchor, guard, deadline), z, alpha, guard)

                result.build_caps = build
                return result

            module.load_composition = composition
        return module

    sizes.load = with_refinement
    print(json.dumps({"phase": "residual_rank_driver_provenance", "sizes_pin": SIZES_PIN,
        **provenance}), flush=True)
    try:
        sizes.main()
    finally:
        sizes.load = original_load


if __name__ == "__main__":
    main()
