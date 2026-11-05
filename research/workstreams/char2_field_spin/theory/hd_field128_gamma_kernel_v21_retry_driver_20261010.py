"""Bounded receipt-only retry driver; all numerical arithmetic is frozen elsewhere."""
from hashlib import sha256
from pathlib import Path
import json
import shutil

import hd_field128_gamma_all_pair_certificate as cert
import hd_field128_gamma_midpoint_weight_bands as bands
import hd_field128_gamma_profile_certificate as output


def run(bank, first, last, folder, witness_source):
    """Retry one explicit assigned band at 2048 nodes; never rebuild the bank."""
    if type(bank) is not cert.AliasBank:
        raise TypeError("a fresh actual all-pair AliasBank is required")
    if (type(first) is not int or type(last) is not int
            or not 1697 <= first <= last <= 1920 or last - first >= 64):
        raise ValueError("expected an explicit assigned band of at most 64 weights")
    here = Path(__file__).resolve()
    driver_bytes = here.read_bytes()
    witness_path = Path(witness_source["path"]).resolve()
    witness_bytes = witness_path.read_bytes()
    if sha256(witness_bytes).hexdigest() != witness_source["sha256"]:
        raise ValueError("witness byte pin mismatch")
    destination = Path(folder).resolve()
    if shutil.disk_usage(destination.parent).free < 750_000_000:
        raise RuntimeError("resource stop: fewer than 750 MB free")

    def check():
        if here.read_bytes() != driver_bytes or witness_path.read_bytes() != witness_bytes:
            raise RuntimeError("driver or witness changed during retry")
        bank.check_pins()

    def write_guarded(path, value):
        if shutil.disk_usage(destination.parent).free < 750_000_000:
            raise RuntimeError("resource stop: fewer than 750 MB free; checkpoints preserved")
        return output.write_new(path, value)

    check()
    destination.mkdir(exist_ok=False)
    try:
        write_guarded(destination / "driver.json", dict(
            schema="gamma-v21-bounded-band-retry-driver-v1",
            driver=dict(path=str(here), sha256=sha256(driver_bytes).hexdigest()),
            witness_source=witness_source, requested_weight_range=[first, last],
            max_nodes=2048, max_prepare_seconds=180, max_cover_seconds=60,
            cell_failure_bits=72, aggregate_failure_bits=54,
            numerical_provider="hd_field128_gamma_midpoint_weight_bands.py",
            fresh_bank_rebuilt=False, old_numerical_bounds_reused=False))
        evaluator = bands.Evaluator(
            bank, first, last, max_seconds=180., witness_source=witness_source,
            progress=lambda row: print(json.dumps(row), flush=True))
        result = evaluator.run_cover(
            max_nodes=2048, max_seconds=60., cell_failure_bits=72, aggregate_bits=54,
            checkpoint=lambda row: write_guarded(
                destination / f"node_{row['sequence']:06d}.json", row),
            context_checkpoint=lambda row: write_guarded(destination / "context.json", row))
        check()
        write_guarded(destination / "summary.json", result)
        print(json.dumps({key: result[key] for key in (
            "status", "requested_weight_range", "evaluated_nodes", "tree_nodes",
            "log2_upper_rational", "meets_requested_band_target", "stop_reason")}), flush=True)
        return result
    finally:
        check()
