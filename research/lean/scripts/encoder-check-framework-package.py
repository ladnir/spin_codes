"""Controlled test of Lake package metadata in isolated Framework output."""
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "scripts/map_data"
spec = importlib.util.spec_from_file_location("runner", ROOT / "scripts/encoder-recheck-foundations-isolated.py")
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


def main():
    comparison = json.loads((DATA / "encoder_framework_comparison.json").read_text())
    assert comparison["status"] == "PASS"
    batch = json.loads((DATA / "encoder_foundations_isolated_verification.json").read_text())
    assert batch["status"] == "PASS"
    run_dir = Path(batch["run_directory"])
    source = ROOT / "SpinCodes/Framework.lean"
    original = ROOT / ".lake/build/lib/lean/SpinCodes/Framework.olean"
    before_source, before_object = runner.sha(source), runner.sha(original)
    setup = run_dir / "Framework.package.setup.json"
    setup.write_text(json.dumps({"name": "SpinCodes.Framework", "package": "spincodes", "isModule": False,
                                 "importArts": {}, "plugins": [], "dynlibs": [], "options": {}}), encoding="utf-8")
    output = run_dir / "Framework.with-package.olean"
    log = run_dir / "Framework.with-package.log"
    environment = os.environ.copy()
    result = subprocess.run([shutil.which("lake"), "env", sys.executable, "-c",
        "import os,json; print(json.dumps({k:os.environ.get(k,'') for k in ['LEAN_PATH','PATH']}))"],
        cwd=ROOT, capture_output=True, text=True, check=True)
    environment.update(json.loads(result.stdout))
    lean = shutil.which("lean", path=environment["PATH"])
    command = [lean, "-j1", "-M", "4000", "--setup", str(setup), "-o", str(output), "SpinCodes/Framework.lean"]
    report = {"status": "RUNNING", "scope": "One isolated recheck with explicit Lake package metadata; existing imports reused", "module": "SpinCodes.Framework",
              "source_sha256": before_source, "original_object_sha256": before_object, "setup": str(setup), "setup_sha256": runner.sha(setup), "command": command}
    (DATA / "encoder_framework_package_verification.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    try:
        minimum = runner.free_gib()
        assert minimum >= 3
        with log.open("wb") as stream:
            proc = subprocess.Popen(command, cwd=ROOT, env=environment, stdout=stream, stderr=subprocess.STDOUT)
            while proc.poll() is None:
                minimum = min(minimum, runner.free_gib())
                if minimum < 3:
                    proc.terminate()
                    proc.wait()
                    report.update(minimum_free_gib=minimum, exit_code=proc.returncode, resource_guard_stopped=True)
                    raise RuntimeError("Free RAM fell below 3 GiB")
                time.sleep(1)
        report.update(exit_code=proc.returncode, minimum_free_gib=minimum, log=str(log), log_sha256=runner.sha(log))
        assert proc.returncode == 0
        assert runner.sha(source) == before_source and runner.sha(original) == before_object
        fresh = runner.sha(output)
        report.update(status="PASS", fresh_object=str(output), fresh_object_sha256=fresh,
                      fresh_object_bytes_match_original=fresh == before_object, all_live_objects_unchanged=True)
        if fresh == before_object:
            report["object_sha256"] = fresh
        print(json.dumps(report, indent=2), flush=True)
    except BaseException as error:
        report.update(status="FAIL", error=str(error))
        raise
    finally:
        report["live_source_unchanged"] = runner.sha(source) == before_source
        report["live_object_unchanged"] = runner.sha(original) == before_object
        (DATA / "encoder_framework_package_verification.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
