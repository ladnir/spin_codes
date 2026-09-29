"""Compare printed Framework declarations in live and isolated import trees."""
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
spec = importlib.util.spec_from_file_location("runner", ROOT / "scripts/encoder-recheck-foundations-isolated.py")
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


def main():
    report_path = DATA / "encoder_framework_comparison.json"
    batch_path = DATA / "encoder_foundations_isolated_verification.json"
    while True:
        batch = json.loads(batch_path.read_text())
        if batch["status"] != "RUNNING":
            break
        time.sleep(3)
    assert batch["status"] == "PASS", "Foundation batch did not finish successfully"
    run_dir = Path(batch["run_directory"])
    source = ROOT / "SpinCodes/Framework.lean"
    names = ["Spin.nonzeroMsgs", "Spin.Setup"]
    names += ["Spin.Setup." + x for x in re.findall(r"^(?:def|lemma|theorem) ([A-Za-z_][A-Za-z_0-9]*)", source.read_text(encoding="utf-8"), re.M)
              if x != "nonzeroMsgs"]
    names += ["Spin.Setup.mk", "Spin.Setup.rec", "Spin.Setup.Pout", "Spin.Setup.Pin", "Spin.Setup.outer", "Spin.Setup.innerWt"]
    pin = run_dir / "EncoderFrameworkCompare.lean"
    pin.write_text("import SpinCodes.Framework\nset_option pp.all true\nset_option pp.proofs true\nset_option format.width 10000\nset_option maxRecDepth 100000\n" +
                   "\n".join("#print " + n + "\n#print axioms " + n for n in names) + "\n", encoding="utf-8")
    environment = os.environ.copy()
    result = subprocess.run([shutil.which("lake"), "env", sys.executable, "-c",
        "import os,json; print(json.dumps({k:os.environ.get(k,'') for k in ['LEAN_PATH','PATH']}))"],
        cwd=ROOT, capture_output=True, text=True, check=True)
    environment.update(json.loads(result.stdout))
    lean = shutil.which("lean", path=environment["PATH"])
    report = {"status": "RUNNING", "scope": "Printed declaration/type/body/axiom comparison, not byte identity or source-output association", "declarations": names, "checks": []}
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    try:
        for mode in ["live", "fresh"]:
            assert runner.free_gib() >= 3
            env = environment.copy()
            if mode == "fresh":
                env["LEAN_PATH"] = str(run_dir / "objects") + os.pathsep + env["LEAN_PATH"]
            log = run_dir / ("FrameworkDeclarations_" + mode + ".log")
            command = [lean, "-j1", "-M", "4000", pin.relative_to(ROOT).as_posix()]
            with log.open("wb") as stream:
                proc = subprocess.Popen(command, cwd=ROOT, env=env, stdout=stream, stderr=subprocess.STDOUT)
                while proc.poll() is None:
                    if runner.free_gib() < 3:
                        proc.terminate()
                        proc.wait()
                        raise RuntimeError("Free RAM fell below 3 GiB")
                    time.sleep(1)
            report["checks"].append({"mode": mode, "exit_code": proc.returncode, "log": log.relative_to(ROOT).as_posix(), "log_sha256": runner.sha(log), "command": command})
            assert proc.returncode == 0
        live, fresh = [ROOT / x["log"] for x in report["checks"]]
        report.update(status="PASS", printed_declarations_identical=live.read_bytes() == fresh.read_bytes())
        print(json.dumps(report, indent=2), flush=True)
    except BaseException as error:
        report.update(status="FAIL", error=str(error))
        raise
    finally:
        report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
