"""Matched complete-size FLOAT control changing only the return cap.

Delegate CLI, guards, routing, source factors and all-profile joining to
the frozen size screen. Replace its private composition.build_caps by the
same call followed by the audited coset-return augmentation. Rebuild all
pair corrections from the resulting coherent raw/powered maps. No frozen
file is edited. Exact local tests run within the delegated resource guard.
The final report's block_caps.coset_return records the local change; the
initial provenance record pins this wrapper and both dependencies.
This remains FLOAT reconnaissance, not a directed certificate.
"""

import hashlib
import importlib.util
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
SCREEN_PIN = "7a7024b78912abb78eaaef9e90f766eefdd7e7f807d8dfe6143f38c86c233afd"
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
    screen = load("cheap_lane_pair_sizes_screen.py", SCREEN_PIN, "coset_sizes_baseline")
    refined = load("cheap_lane_coset_return_refinement.py", REFINEMENT_PIN, "coset_sizes_refinement")
    original_load = screen.load

    def with_refinement(name, pin, label):
        module = original_load(name, pin, label)
        if name == "cheap_lane_pair_bank_screen.py":
            original_composition = module.load_composition

            def composition():
                result = original_composition()
                original_build = result.build_caps

                def build(z, alpha, triple, anchor, guard, deadline):
                    tests = refined.exact_tests(guard)
                    caps = original_build(z, alpha, triple, anchor, guard, deadline)
                    caps = refined.augment_caps(caps, z, alpha, guard)
                    caps["tests"]["coset_return_exact_tests"] = tests
                    caps["block_report"]["coset_return"]["source_sha256"] = REFINEMENT_PIN
                    return caps

                result.build_caps = build
                return result

            module.load_composition = composition
        return module

    screen.load = with_refinement
    print(json.dumps({"phase": "coset_return_driver_provenance",
        "screen_pin": SCREEN_PIN, "refinement_pin": REFINEMENT_PIN,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "changed": "raw E and its local power; pair caps rebuilt",
        "routing_source_terminal_unchanged": True}), flush=True)
    try:
        screen.main()
    finally:
        screen.load = original_load


if __name__ == "__main__":
    main()
