"""Explicit resource launcher for the frozen directed regional driver.

The first512MiB pilot stopped at its memory guard without a result. This
launcher permits a caller-selected peak-process budget through1536MiB.
It requires free physical RAM of at least that budget plus2GiB, and keeps
the20GiB disk guard. It changes only the imported driver's MEMORY constant.
Numerical functions, witnesses, precision, deadlines and frozen files are
unchanged. Execution remains serial under the caller's control.

No numerical work occurs on import. Stdout only; no generated files.
"""

import argparse
import ctypes
from ctypes import wintypes
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import sys

sys.dont_write_bytecode = True
DRIVER = "cheap_lane_pair_regions_directed.py"
DRIVER_PIN = "d3a6e858437e19bede76e215a9d9b31f26a49b9b601ecd60705f1a2b99f0f213"
CENTERS = (32, 96, 160, 224, 288, 352, 416, 472)
MIB, GIB = 1 << 20, 1 << 30


def physical_memory():
    """Read Windows physical-memory availability without launching a job."""
    if sys.platform != "win32":
        raise RuntimeError("this launcher requires the audited Windows RAM query")
    class MemoryStatus(ctypes.Structure):
        _fields_ = [("length", wintypes.DWORD), ("load_percent", wintypes.DWORD),
                    *[(name, ctypes.c_ulonglong) for name in (
                        "total_physical", "free_physical", "total_pagefile", "free_pagefile",
                        "total_virtual", "free_virtual", "free_extended_virtual")]]
    status = MemoryStatus()
    status.length = ctypes.sizeof(status)
    query = ctypes.windll.kernel32.GlobalMemoryStatusEx
    query.argtypes = [ctypes.POINTER(MemoryStatus)]
    query.restype = wintypes.BOOL
    if not query(ctypes.byref(status)):
        raise OSError("cannot read free physical RAM")
    return {"total_physical_bytes": int(status.total_physical),
            "free_physical_bytes": int(status.free_physical)}


def preflight(path, budget):
    disk = shutil.disk_usage(path.resolve().anchor).free
    if disk < 20*GIB:
        raise RuntimeError("requires20GiB free on the source volume")
    ram = physical_memory()
    required = budget+2*GIB
    if ram["free_physical_bytes"] < required:
        raise RuntimeError("free physical RAM must exceed the selected budget by at least2GiB")
    return {**ram, "free_disk_bytes": disk,
            "required_free_physical_bytes": required, "required_free_disk_bytes": 20*GIB}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--memory-mib", type=int, choices=(512, 1024, 1536), default=512)
    parser.add_argument("--seconds", type=int, choices=(300, 1200), default=300)
    parser.add_argument("--precision", type=int, choices=(256, 384), default=256)
    parser.add_argument("--test-center", type=int, choices=CENTERS)
    args = parser.parse_args()
    path = Path(__file__).with_name(DRIVER)
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual != DRIVER_PIN:
        raise RuntimeError("frozen regional-driver hash mismatch")
    budget = args.memory_mib*MIB
    resources = preflight(path, budget)
    spec = importlib.util.spec_from_file_location("explicit_memory_pair_regions", path)
    driver = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(driver)
    if driver.MEMORY != 512*MIB:
        raise RuntimeError("unexpected original driver memory budget")
    driver.MEMORY = budget
    # Check again immediately before dispatch. The import above is source-only.
    resources = preflight(path, budget)
    forwarded = [str(path), "--seconds", str(args.seconds), "--precision", str(args.precision)]
    if args.test_center is not None:
        forwarded.extend(("--test-center", str(args.test_center)))
    print(json.dumps({"phase": "explicit_resource_launch", "memory_mib": args.memory_mib,
        "memory_budget_bytes": budget, "seconds": args.seconds, "precision_bits": args.precision,
        "test_center": args.test_center, "resource_preflight": resources,
        "override": "driver.MEMORY only", "driver_pin": DRIVER_PIN,
        "numerical_functions_unchanged": True, "disk_guard_unchanged": True,
        "serial_execution_required": True, "generated_bytes": 0,
        "launcher_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}), flush=True)
    previous = sys.argv
    sys.argv = forwarded
    try:
        driver.main()
    finally:
        sys.argv = previous


if __name__ == "__main__":
    main()
