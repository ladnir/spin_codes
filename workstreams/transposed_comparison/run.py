#!/usr/bin/env python3
"""Serial, same-host comparison. Generated observations are deliberately small."""
import argparse
import hashlib
import json
import platform
import statistics
import subprocess
from pathlib import Path

LIBOTE_COMMIT = "fa7bfc3c8a0178fc99b0b666d5f31fc2a9277579"
CRYPTOTOOLS_COMMIT = "643b3ca0b57ee4ec5e3388992c5906e26f3ca72d"
MODES = ("spin", "golay", "rm", "exconv", "raa3", "raa4")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def summarize(rows):
    result = []
    for m in sorted({row["m"] for row in rows}):
        for mode in MODES:
            group = [row for row in rows if row["m"] == m and row["mode"] == mode]
            if not group:
                continue
            assert len({(r["n"], r["k"], r["output_hash"]) for r in group}) == 1
            medians = [statistics.median(row["samples_ms"]) for row in group]
            result.append(dict(mode=mode, m=m, n=group[0]["n"], k=group[0]["k"],
                               median_ms=statistics.median(medians), run_medians_ms=medians))
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("binary", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--m", nargs="+", type=int, default=[16, 18, 20], choices=[16, 18, 20])
    p.add_argument("--trials", type=int, default=101)
    p.add_argument("--runs", type=int, default=3)
    args = p.parse_args()
    if args.output.exists():
        p.error("output exists; use a new results filename")
    if args.trials < 3 or not args.trials % 2 or args.runs < 1:
        p.error("require positive runs and odd trials >= 3")
    binary = args.binary.resolve()
    # Correctness is completed before any measurements. No parallel dispatch.
    subprocess.run([str(binary), "check"], check=True)
    here = Path(__file__).resolve().parent
    repo = here.parent.parent
    sources = sorted(here.glob("*.h")) + sorted(here.glob("*.cpp")) + [here / "CMakeLists.txt", here / "run.py"]
    sources += [repo / "workstreams/bare_bch_rm2sub" / name
                for name in ("Spin.h", "Spin.cpp", "Inner.h", "generated/BchCircuit.h", "generated/BchCircuit.cpp")]
    result = dict(schema=2, status="incomplete", input_policy="no_reset", platform=platform.platform(),
                  processor=subprocess.check_output(["lscpu"], text=True),
                  compiler=subprocess.check_output(["g++", "--version"], text=True),
                  libote_commit=LIBOTE_COMMIT, cryptotools_commit=CRYPTOTOOLS_COMMIT,
                  binary_sha256=digest(binary),
                  source_sha256={str(path.relative_to(repo)): digest(path) for path in sources},
                  exconv_parameters=dict(implementation="ExConvCodeOld", w=33, m=25,
                                         random_taps=24, systematic=False, convolution_passes=1),
                  protocol="one CPU 15; three warmups; native buffers reused with no input reset or copy; "
                           "setup, allocation, guard before/after batch, and hashing excluded; serial processes",
                  compile_commands=json.loads((binary.parent / "compile_commands.json").read_text()),
                  runs=[])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as out:
        json.dump(result, out, indent=2)
    for m in args.m:
        for run in range(args.runs):
            order = MODES[2 * (run % 3):] + MODES[:2 * (run % 3)]
            for mode in order:
                completed = subprocess.run([str(binary), mode, str(m), str(args.trials)],
                                           check=True, capture_output=True, text=True)
                row = json.loads(completed.stdout)
                row["run"] = run
                result["runs"].append(row)
                args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
                print(f"m={m} run={run+1} {mode}: {statistics.median(row['samples_ms']):.3f} ms", flush=True)
    result["summary"] = summarize(result["runs"])
    result["status"] = "complete"
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
