"""Recheck 12 semantic sources locally without changing any live build object."""
import ctypes
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "scripts/map_data"
MODULES = ['SpinCodes.Framework', 'SpinCodes.Distance', 'SpinCodes.Selection', 'SpinCodes.Cover.DenseTail', 'SpinCodes.Structured.Certificate', 'SpinCodes.Structured.Regimes', 'SpinCodes.Structured.Instantiation', 'SpinCodes.Structured.Composition', 'SpinCodes.Structured.Interleaver', 'SpinCodes.Structured.Enumerator', 'SpinCodes.Structured.Schedule', 'SpinCodes.Structured.WeightedNorm']
THEOREMS = ['Spin.Setup.prob_bad_le_cond_sum', 'Spin.one_le_Z_iff', 'Spin.prob_not_good_tendsto_zero', 'Spin.Cover.denseTail', 'Spin.Structured.ScalableCertificate.distance_whp', 'Spin.Structured.total_regime_sum', 'Spin.Structured.Family.distance_whp_native', 'Spin.Structured.prob_two_stage', 'Spin.Structured.card_perm_accWt_eq', 'Spin.Structured.card_tuples_weight', 'Spin.Structured.Nsched_ratio_tendsto_one', 'Spin.rowNormLe_pathProduct']
ALLOWED = {"propext", "Classical.choice", "Quot.sound"}


class MemoryStatus(ctypes.Structure):
    _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong)] + [
        (name, ctypes.c_ulonglong) for name in ["ullTotalPhys", "ullAvailPhys", "ullTotalPageFile",
                                               "ullAvailPageFile", "ullTotalVirtual", "ullAvailVirtual",
                                               "ullAvailExtendedVirtual"]]


def free_gib():
    status = MemoryStatus()
    status.dwLength = ctypes.sizeof(status)
    assert ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status))
    return status.ullAvailPhys / 2**30


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main():
    assert os.name == "nt", "This bounded repair uses local Windows only"
    sys.stdout.reconfigure(encoding="utf-8")
    run_dir = DATA / ("encoder_foundations_isolated_" + time.strftime("%Y%m%d_%H%M%S") + f"_{os.getpid()}")
    objects = run_dir / "objects"
    objects.mkdir(parents=True)
    report_path = DATA / "encoder_foundations_isolated_verification.json"
    report = {"status": "RUNNING", "started_utc": datetime.now(timezone.utc).isoformat(),
              "scope": "Fresh local source checks in isolated outputs; existing external imports reused; live objects unchanged",
              "compiler_cap_mb": 4000, "compiler_jobs": 1, "stop_below_free_gib": 3,
              "initial_free_gib": free_gib(), "run_directory": str(run_dir), "modules": []}

    def save():
        text = json.dumps(report, indent=2) + "\n"
        (run_dir / "report.json").write_text(text, encoding="utf-8")
        report_path.write_text(text, encoding="utf-8")

    save()
    try:
        assert report["initial_free_gib"] >= 3, "Local available memory below 3 GiB; not starting"
        spec = importlib.util.spec_from_file_location("replay", ROOT / "scripts/replay-native-closure.py")
        replay = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(replay)
        selected = set(MODULES)
        order = []
        imported = set()
        for module in MODULES:
            _, chain = replay.graph(module)
            imported.update(chain)
            for module in chain:
                if module in selected and module not in order:
                    order.append(module)
        assert len(order) == len(MODULES)
        live = ROOT / ".lake/build/lib/lean"
        snapshot = {m: {"source": ROOT / (m.replace(".", "/") + ".lean"),
                        "original": live / (m.replace(".", "/") + ".olean")} for m in order}
        for v in snapshot.values():
            v["source_hash"] = sha(v["source"])
            v["original_hash"] = sha(v["original"])
        # Lean resolves a top-level module prefix to one search-path root. Give
        # this isolated root its complete checked import tree. These hardlinks
        # are read-only compiler inputs; none of the selected output targets is linked.
        linked = 0
        for module in sorted(imported - selected):
            original = live / (module.replace(".", "/") + ".olean")
            assert original.exists(), f"Missing checked dependency: {module}"
            for artifact in original.parent.glob(original.name + "*"):
                destination = objects / artifact.relative_to(live)
                destination.parent.mkdir(parents=True, exist_ok=True)
                os.link(artifact, destination)
                linked += 1
        report["reused_dependency_modules"] = len(imported - selected)
        report["read_only_import_hardlinks"] = linked
        lake = shutil.which("lake")
        assert lake
        result = subprocess.run([lake, "env", sys.executable, "-c",
            "import os,json; print(json.dumps({k:os.environ.get(k,'') for k in ['LEAN_PATH','PATH']}))"],
            cwd=ROOT, capture_output=True, text=True, check=True)
        environment = os.environ.copy()
        environment.update(json.loads(result.stdout))
        environment["LEAN_PATH"] = str(objects) + os.pathsep + environment.get("LEAN_PATH", "")
        lean = shutil.which("lean", path=environment["PATH"])
        assert lean
        report["lean_version"] = subprocess.check_output([lean, "--version"], cwd=ROOT, env=environment, text=True).strip()
        assert "version 4.34.0" in report["lean_version"]

        def check(src, output, log):
            assert free_gib() >= 3, "Available memory below 3 GiB before check"
            output.parent.mkdir(parents=True, exist_ok=True)
            command = [lean, "-j1", "-M", "4000", "-o", str(output), src.relative_to(ROOT).as_posix()]
            start, minimum = time.monotonic(), free_gib()
            with log.open("wb") as stream:
                proc = subprocess.Popen(command, cwd=ROOT, env=environment, stdout=stream, stderr=subprocess.STDOUT)
                try:
                    while proc.poll() is None:
                        available = free_gib()
                        minimum = min(minimum, available)
                        if available < 3:
                            proc.terminate()
                            proc.wait()
                            raise RuntimeError("Stopped compiler: local available memory fell below 3 GiB")
                        time.sleep(1)
                finally:
                    if proc.poll() is None:
                        proc.terminate()
                        proc.wait()
            return {"command": command, "exit_code": proc.returncode,
                    "seconds": round(time.monotonic() - start, 3), "minimum_free_gib": minimum,
                    "log": log.relative_to(ROOT).as_posix()}

        for module in order:
            before = snapshot[module]
            assert sha(before["source"]) == before["source_hash"]
            output = objects / (module.replace(".", "/") + ".olean")
            log = run_dir / (module.rsplit(".", 1)[-1] + ".log")
            record = {"module": module, "source": before["source"].relative_to(ROOT).as_posix(),
                      "source_sha256": before["source_hash"], "original_object_sha256": before["original_hash"],
                      "fresh_object": output.relative_to(ROOT).as_posix()}
            report["modules"].append(record)
            record.update(check(before["source"], output, log))
            save()
            assert record["exit_code"] == 0, f"Compiler failed: {module}; see {log}"
            assert sha(before["source"]) == before["source_hash"], f"Source changed: {module}"
            fresh = sha(output)
            record.update(fresh_object_sha256=fresh, fresh_object_bytes_match_original=fresh == before["original_hash"])
            if fresh == before["original_hash"]:
                record["object_sha256"] = fresh  # Only expose a live paired hash when bytes really match.
            record["status"] = "PASS"
            save()
            print(f'PASS {module}; matches live bytes={record["fresh_object_bytes_match_original"]}; '
                  f'min free={record["minimum_free_gib"]:.2f} GiB', flush=True)
        pin = run_dir / "EncoderFoundationsEvidencePin.lean"
        pin.write_text("\n".join("import " + m for m in order) + "\n" +
                       "\n".join("#print axioms " + t for t in THEOREMS) +
                       "\n", encoding="utf-8")
        pin_result = check(pin, run_dir / "EncoderFoundationsEvidencePin.olean", run_dir / "EncoderFoundationsEvidencePin.log")
        report["audit_pin"] = pin_result
        assert pin_result["exit_code"] == 0, "Isolated audit pin failed"
        log = (run_dir / "EncoderFoundationsEvidencePin.log").read_text(encoding="utf-8")
        rows = re.findall(r"'([^']+)' depends on axioms: \[([^]]*)\]", log)
        assert len(rows) == len(THEOREMS), rows
        for name, names in rows:
            assert {a.strip() for a in names.split(",")} <= ALLOWED, (name, names)
        for module, before in snapshot.items():
            assert sha(before["source"]) == before["source_hash"], f"Source changed: {module}"
            assert sha(before["original"]) == before["original_hash"], f"Live object changed: {module}"
        report.update(status="PASS", finished_utc=datetime.now(timezone.utc).isoformat(),
                      all_live_objects_unchanged=True, axiom_audit_count=len(rows),
                      axiom_theorems=[name for name, _ in rows], allowed_axioms=sorted(ALLOWED),
                      all_fresh_objects_match_live=all(r["fresh_object_bytes_match_original"] for r in report["modules"]))
        print("PASS: 12 isolated source checks and 12 axiom audits; live objects preserved.", flush=True)
    except BaseException as error:
        report.update(status="FAIL", error=str(error))
        raise
    finally:
        save()


if __name__ == "__main__":
    main()
