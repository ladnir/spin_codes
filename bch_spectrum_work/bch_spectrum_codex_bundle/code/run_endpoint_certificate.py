"""Freeze, execute, replay, and check one BCH endpoint experiment.

The stored Windows CNG stream is computational randomness. The binomial
soundness bound is conditional on ideal independent bytes. No retry on hits.
"""
from __future__ import annotations

import argparse
import ctypes
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from fractions import Fraction
from math import comb
from pathlib import Path

from analyze_endpoint_sampling import log_integer, neg_log_one_minus, ratio_record

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "code" / "test_bch_endpoints.cpp"
EXE = ROOT / "generated" / "test_bch_endpoints.exe"
THRESHOLD = ROOT / "generated" / "random_inner_threshold_outward.json"


def sha(path: Path) -> str:
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def write_new(path: Path, value: dict) -> None:
    with path.open("x", encoding="utf-8") as f:
        json.dump(value, f, indent=2)
        f.write("\n")


class ExperimentLock:
    """Keep this task's sampler and replay runs sequential across processes."""
    def __enter__(self):
        self.api = ctypes.WinDLL("kernel32", use_last_error=True)
        self.api.CreateMutexW.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_wchar_p]
        self.api.CreateMutexW.restype = ctypes.c_void_p
        self.api.WaitForSingleObject.argtypes = [ctypes.c_void_p, ctypes.c_uint32]
        self.api.ReleaseMutex.argtypes = [ctypes.c_void_p]
        self.api.CloseHandle.argtypes = [ctypes.c_void_p]
        self.handle = self.api.CreateMutexW(None, False, "Local\\CodexBCH256EndpointExperiment")
        if not self.handle:
            raise ctypes.WinError(ctypes.get_last_error())
        status = self.api.WaitForSingleObject(self.handle, 0)
        if status not in (0, 0x80):
            self.api.CloseHandle(self.handle)
            raise RuntimeError("Another BCH endpoint experiment holds the run lock")
        return self

    def __exit__(self, *_):
        self.api.ReleaseMutex(self.handle)
        self.api.CloseHandle(self.handle)


def verify_trace(directory: Path) -> int:
    # Third implementation: Gaussian elimination on the 19-point Vandermonde
    # matrix, with a field table constructed from carryless multiplication.
    def mul(a: int, b: int) -> int:
        x = 0
        for j in range(8):
            if b & (1 << j):
                x ^= a << j
        for j in range(14, 7, -1):
            if x & (1 << j):
                x ^= 0x14d << (j-8)
        return x

    table = [[mul(a, b) for b in range(256)] for a in range(256)]
    inverse = [0] + [table[a].index(1) for a in range(1, 256)]

    def power(x: int, e: int) -> int:
        v = 1
        for _ in range(e):
            v = table[v][x]
        return v

    target = [power(x, 146) for x in range(256)]
    checked = 0
    previous = -1
    with (directory / "random.bin").open("rb") as tape:
        for line in (directory / "sample.json.audit.jsonl").read_text().splitlines():
            row = json.loads(line)
            assert row["sample"] > previous
            previous = row["sample"]
            tape.seek(row["byte_offset"])
            raw = tape.read(row["next_byte_offset"] - row["byte_offset"])
            nodes = [0]
            for x in raw:
                if x not in nodes:
                    nodes.append(x)
            assert len(nodes) == 19 and nodes == row["nodes"]
            matrix = []
            for x in nodes:
                powers = [1]
                for j in range(18):
                    powers.append(table[powers[-1]][x])
                matrix.append(powers + [1 if x == 0 else target[x]])
            for j in range(19):
                pivot = next(i for i in range(j, 19) if matrix[i][j])
                matrix[j], matrix[pivot] = matrix[pivot], matrix[j]
                inv = inverse[matrix[j][j]]
                matrix[j] = [table[inv][x] for x in matrix[j]]
                for i in range(19):
                    if i == j:
                        continue
                    scale = matrix[i][j]
                    if scale:
                        matrix[i] = [a ^ table[scale][b] for a, b in zip(matrix[i], matrix[j])]
            g = [matrix[i][-1] for i in range(19)]
            assert g == row["coefficients"] and g[0] == 1
            count = 0
            for x in range(1, 256):
                value = 0
                for c in reversed(g):
                    value = table[value][x] ^ c
                count += value == target[x]
            assert count == row["agreements"]
            checked += 1
    return checked


def verify(directory: Path) -> dict:
    manifest_path = directory / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    for relative, expected in manifest["input_sha256"].items():
        assert sha(ROOT / relative) == expected, f"Frozen input changed: {relative}"
    primary = json.loads((directory / "sample.json").read_text())
    replay = json.loads((directory / "replay.json").read_text())
    assert primary["mode"] == "BCryptGenRandom"
    assert replay["mode"] == "independent_replay"
    for key in ("samples", "endpoint_hits", "bytes_consumed", "bytes_loaded",
                "zero_bytes_rejected", "duplicate_bytes_rejected", "extra_agreement_histogram"):
        assert primary[key] == replay[key], key
    n, hits = primary["samples"], primary["endpoint_hits"]
    assert n == manifest["samples"]
    histogram = primary["extra_agreement_histogram"]
    assert len(histogram) == 20 and sum(histogram) == n and histogram[19] == hits
    assert all(isinstance(v, int) and v >= 0 for v in histogram)
    assert primary["bytes_consumed"] == 18*n + primary["zero_bytes_rejected"] + primary["duplicate_bytes_rejected"]
    size = (directory / "random.bin").stat().st_size
    assert size == primary["bytes_loaded"]
    assert 0 <= size - primary["bytes_consumed"] < 1 << 20
    assert sha(directory / "sample.json.audit.jsonl") == sha(directory / "replay.json.audit.jsonl")
    checked = verify_trace(directory)
    cap = int(json.loads(THRESHOLD.read_text())["a38_sufficient_cap"]["certified_integer_cap"])
    assert cap == manifest["a38_cap"]
    first_bad = 19*(cap//3968 + 1)
    p_bad = Fraction(first_bad*comb(37, 18), comb(255, 18))
    decay_lo, decay_hi = neg_log_one_minus(p_bad)
    log2_lo, log2_hi = log_integer(2)
    sufficient = n*decay_lo >= manifest["error_bits"]*log2_hi
    if not manifest["diagnostic_only"]:
        assert sufficient
    accepted = sufficient and hits == 0 and not manifest["diagnostic_only"]
    result = {
        "status": "conditional_statistical_acceptance" if accepted else "diagnostic_only" if manifest["diagnostic_only"] else "inconclusive",
        "a38_cap": cap, "samples": n, "endpoint_hits": hits,
        "full_independent_algorithm_replay_matches": True,
        "independent_python_trace_records_checked": checked,
        "sample_seconds": primary["seconds"], "replay_seconds": replay["seconds"],
        "false_accept_probability_under_IID_at_most": f"2^-{manifest['error_bits']}" if accepted else None,
        "no_hit_log2_probability_upper": ratio_record(-n*decay_lo/log2_hi),
        "randomness_qualification": "Execution uses Windows CNG cryptographic pseudorandom bytes. The statistical bound assumes ideal IID bytes; no unconditional information-theoretic RNG guarantee is claimed. Under computational replacement, add the applicable distinguishing advantage.",
        "scope": "Current RandomStepConv occupation-one BCH weight-38 obligation only; neither full SPIN distance nor a deterministic spectrum proof.",
        "sha256": {name: sha(directory / name) for name in (
            "manifest.json", "random.bin", "sample.json", "replay.json", "sample.json.audit.jsonl", "replay.json.audit.jsonl")}
    }
    return result


def execute(directory: Path, samples: int, diagnostic: bool) -> None:
    with ExperimentLock():
        directory.mkdir(parents=True, exist_ok=False)
        cap = int(json.loads(THRESHOLD.read_text())["a38_sufficient_cap"]["certified_integer_cap"])
        p_bad = Fraction(19*(cap//3968+1)*comb(37, 18), comb(255, 18))
        if not diagnostic:
            assert samples*neg_log_one_minus(p_bad)[0] >= 40*log_integer(2)[1]
        assert samples > 0
        inputs = [SOURCE, EXE, THRESHOLD, Path(__file__), ROOT / "code" / "analyze_endpoint_sampling.py"]
        manifest = {
            "protocol": "fixed-bch256-endpoint-v1", "created_utc": datetime.now(timezone.utc).isoformat(),
            "samples": samples, "error_bits": 40, "diagnostic_only": diagnostic, "a38_cap": cap,
            "field_modulus": "0x14d", "domain_size": 255, "base_size": 18, "endpoint_agreements": 37,
            "acceptance_rule": "exactly zero endpoint hits at the fixed sample count; otherwise inconclusive; no automatic retry",
            "randomness": "BCryptGenRandom, null handle, BCRYPT_USE_SYSTEM_PREFERRED_RNG, 1 MiB blocks retained in full",
            "randomness_model": "Ideal IID bytes for the exact statistical bound; actual CNG execution uses computational pseudorandomness",
            "subset_algorithm": "Read bytes, rejecting zero and within-base duplicates, until 18 distinct nonzero bytes have been accepted",
            "primary_algorithm": "Newton divided differences and AVX2 monomial contribution tables",
            "replay_algorithm": "incremental vanishing-polynomial interpolation and eight-lane scalar Horner evaluation",
            "input_sha256": {str(p.relative_to(ROOT)).replace("\\", "/"): sha(p) for p in inputs}
        }
        write_new(directory / "manifest.json", manifest)
        print("FROZEN", directory / "manifest.json", "sha256=" + sha(directory / "manifest.json"), flush=True)
        for mode, result in (("run", "sample.json"), ("replay", "replay.json")):
            subprocess.run([str(EXE), mode, str(samples), str(directory / "random.bin"),
                            str(directory / result), "fixed-bch256-endpoint-v1"], check=True)
        result = verify(directory)
        write_new(directory / "certificate.json", result)
        print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("operation", choices=("execute", "verify"))
    parser.add_argument("directory", type=Path)
    parser.add_argument("--samples", type=int, default=150361035)
    parser.add_argument("--diagnostic", action="store_true")
    args = parser.parse_args()
    if args.operation == "execute":
        execute(args.directory.resolve(), args.samples, args.diagnostic)
    else:
        print(json.dumps(verify(args.directory.resolve()), indent=2))
