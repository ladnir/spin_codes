"""Wait for the audited native theorem, then compile its separate corollaries."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "scripts/map_data"
BASE = DATA / "native_theorem_verification.json"
REPORT = DATA / "encoder_native_distance_consequences_verification.json"
ALLOWED = {"propext", "Classical.choice", "Quot.sound"}
STEMS = ["ConcreteNativeDistanceConsequences", "ConcreteNativeDistanceConsequencesPin",
         "ConcreteNativeDistanceConsequencesFinal", "ConcreteNativeDistanceConsequencesFinalPin"]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(report):
    tmp = REPORT.with_suffix(".tmp")
    tmp.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    tmp.replace(REPORT)


def audits(log, expected):
    rows = re.findall(r"'([^']+)' depends on axioms: \[([^]]*)\]", log)
    assert len(rows) == expected, (len(rows), expected)
    assert "sorryAx" not in log and "Lean.ofReduceBool" not in log
    for name, names in rows:
        assert {a.strip() for a in names.split(",")} <= ALLOWED, (name, names)
    return [name for name, _ in rows]


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wait", action="store_true")
    args = parser.parse_args()
    report = {"status": "WAITING" if args.wait else "RUNNING",
              "scope": "Separate success-probability and existence corollaries; waits for audited base theorem",
              "checks": []}
    save(report)
    try:
        previous = None
        while True:
            base = json.loads(BASE.read_text(encoding="utf-8")) if BASE.exists() else {"status": "ABSENT"}
            status = base.get("status")
            if status == "PASS":
                break
            assert status not in {"FAIL", "FAILED"}, f"Base theorem checker failed: {base.get('error')}"
            assert args.wait, f"Base theorem not yet PASS: {status}"
            if status != previous:
                print(f"WAITING: base theorem status {status}", flush=True)
                previous = status
            time.sleep(20)
        base_digest = sha(BASE)
        # Check the base record's actual source/object snapshot before and after use.
        def verify_base_files():
            for rel, record in base["modules"].items():
                src = ROOT / rel
                obj = (ROOT / ".lake/build/lib/lean" / rel).with_suffix(".olean")
                assert sha(src) == record["source_sha256"], f"Changed base source: {rel}"
                assert sha(obj) == record["object_sha256"], f"Changed base object: {rel}"
        verify_base_files()
        source_hashes = {s: sha(ROOT / f"SpinCodes/Structured/{s}.lean") for s in STEMS}
        report.update(status="RUNNING", base_report_sha256=base_digest)
        save(report)
        names = []
        for stem in STEMS:
            start = time.monotonic()
            run = subprocess.run([sys.executable, "scripts/peach-lean.py", f"SpinCodes/Structured/{stem}.lean"],
                                 cwd=ROOT, capture_output=True, encoding="utf-8", timeout=1200)
            output = run.stdout + run.stderr
            log = DATA / f"encoder_consequence_check_{stem}.log"
            log.write_text(output, encoding="utf-8")
            report["checks"].append({"module": stem, "exit_code": run.returncode,
                                      "seconds": round(time.monotonic() - start, 3), "log": log.name})
            save(report)
            assert run.returncode == 0, output[-8000:]
            if stem.endswith("Pin"):
                names += audits(output, 5 if stem == "ConcreteNativeDistanceConsequencesPin" else 2)
            print("PASS", stem, flush=True)
        records = []
        for stem in STEMS:
            src = ROOT / f"SpinCodes/Structured/{stem}.lean"
            obj = ROOT / f".lake/build/lib/lean/SpinCodes/Structured/{stem}.olean"
            assert sha(src) == source_hashes[stem], f"Changed corollary source: {stem}"
            assert obj.stat().st_mtime >= src.stat().st_mtime, f"Stale corollary object: {stem}"
            records.append({"module": stem, "source_sha256": sha(src), "olean_sha256": sha(obj)})
        assert sha(BASE) == base_digest, "Base theorem record changed during corollary compilation"
        verify_base_files()
        report.update(status="PASS", modules=records, checked_module_count=len(records),
                      axiom_audit_count=len(names), axiom_theorems=names, allowed_axioms=sorted(ALLOWED))
        print("PASS: actual success probability tends to 1; eventual rate-half realizations have distance >11%.", flush=True)
    except BaseException as error:
        report.update(status="FAIL", error=str(error))
        raise
    finally:
        save(report)


if __name__ == "__main__":
    main()
