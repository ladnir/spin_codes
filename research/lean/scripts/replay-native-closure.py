"""Optional full project-source replay of the native theorem's import closure.

Dry-run by default. Mathlib and Lean remain pinned external dependencies.
Python schedules Lean processes and records evidence; it proves no theorem.
"""
import argparse
import concurrent.futures
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
TARGET = "SpinCodes.Structured.ConcreteNativeTheoremPin"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def uncomment(text):
    """Strip nested comments without scanning large numeric literals byte by byte."""
    tokens = re.compile(r'/-|-/|--|"(?:\\.|[^"\\])*"')
    comment_tokens = re.compile(r'/-|-/')
    out, position, depth = [], 0, 0
    while match := (comment_tokens if depth else tokens).search(text, position):
        prefix = text[position:match.start()]
        out.append("\n" * prefix.count("\n") if depth else prefix)
        token = match.group()
        if token == "/-":
            depth += 1
            out.append(" ")
        elif token == "-/" and depth:
            depth -= 1
            out.append(" ")
        elif token == "--" and not depth:
            end = text.find("\n", match.end())
            position = len(text) if end < 0 else end
            continue
        elif not depth:
            out.append(token)
        else:
            out.append("\n" * token.count("\n"))
        position = match.end()
    if depth:
        raise ValueError("Unclosed Lean comment")
    out.append(text[position:])
    return "".join(out)


def source(module):
    return ROOT / (module.replace(".", "/") + ".lean")


def graph(target):
    deps, order, active = {}, [], set()

    def visit(module):
        if module in active:
            raise ValueError(f"Import cycle at {module}")
        if module in deps:
            return
        active.add(module)
        text = uncomment(source(module).read_text(encoding="utf-8-sig"))
        imports = []
        for line in text.splitlines():
            line = line.strip()
            if not line or line == "prelude":
                continue
            match = re.fullmatch(r"(?:(?:public|meta)\s+)*import\s+(.+)", line)
            if not match:
                break  # Lean imports precede all ordinary commands.
            imports.append(match.group(1))
        children = set()
        for line in imports:
            for name in line.split():
                if not re.fullmatch(r"[A-Za-z_][A-Za-z_0-9.]*", name):
                    raise ValueError(f"Unsupported import syntax in {module}: {line}")
                if name == "SpinCodes" or name.startswith("SpinCodes."):
                    children.add(name)
        for child in sorted(children):
            visit(child)
        active.remove(module)
        deps[module] = children
        order.append(module)

    visit(target)
    return deps, order


def ensure_idle(backend):
    if backend == "peach":
        spec = importlib.util.spec_from_file_location("peach_helper", ROOT / "scripts/peach-lean.py")
        helper = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(helper)
        command = "ps -eo comm= | awk '$1 == \"lean\" {n++} END {print n+0}'"
        output = subprocess.check_output(helper.SSH + [command], timeout=30, text=True)
        count = int(output.strip())
    elif os.name == "nt":
        output = subprocess.check_output(
            ["tasklist", "/FI", "IMAGENAME eq lean.exe", "/FO", "CSV", "/NH"], text=True)
        count = sum(row and row[0].lower() == "lean.exe" for row in csv.reader(output.splitlines()))
    else:
        output = subprocess.check_output(["ps", "-eo", "comm="], text=True)
        count = sum(line.strip() == "lean" for line in output.splitlines())
    if count:
        raise RuntimeError(f"Refusing replay: {count} Lean process(es) already active on {backend}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", default=TARGET, help="Root project module (default: final theorem pin)")
    parser.add_argument("--execute", action="store_true", help="Actually compile; otherwise write only a plan")
    parser.add_argument("--backend", choices=["peach", "local"], default="peach")
    parser.add_argument("--jobs", type=int, choices=range(1, 5), default=1)
    parser.add_argument("--max-modules", type=int, help="Replay at most this many modules; result stays PARTIAL")
    args = parser.parse_args()
    if args.max_modules is not None and args.max_modules < 1:
        parser.error("--max-modules must be positive")
    if not re.fullmatch(r"SpinCodes(?:\.[A-Za-z_][A-Za-z_0-9]*)*", args.target):
        parser.error("--target must be a SpinCodes module name")
    deps, order = graph(args.target)
    digests = {m: sha(source(m)) for m in order}
    result_dir = ROOT / "scripts/map_data" / ("encoder_native_replay_" + time.strftime("%Y%m%d_%H%M%S") + f"_{os.getpid()}")
    result_dir.mkdir(parents=True)
    state = {"status": "PLAN", "target": args.target, "backend": args.backend, "jobs": args.jobs,
             "scope": "All transitive SpinCodes sources; external Lean/mathlib dependencies reused",
             "pins": {p: sha(ROOT / p) for p in ["lean-toolchain", "lakefile.toml", "lake-manifest.json"]},
             "modules": [{"module": m, "source_sha256": digests[m], "imports": sorted(deps[m])} for m in order],
             "completed": []}
    report = result_dir / "report.json"

    def save():
        report.write_text(json.dumps(state, indent=2) + "\n", encoding="utf-8")

    save()
    print(f"{len(order)} project modules; topological plan: {report}", flush=True)
    if not args.execute:
        print("Dry run only. No Lean processes started.")
        return
    ensure_idle(args.backend)
    lock = ROOT / "scripts/map_data/encoder_native_replay.lock"
    lock_fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    os.write(lock_fd, str(os.getpid()).encode())
    os.close(lock_fd)

    def compile_one(module):
        path = source(module)
        if sha(path) != digests[module]:
            raise RuntimeError(f"Source changed before replay: {module}")
        relative = path.relative_to(ROOT).as_posix()
        obj = ROOT / ".lake/build/lib/lean" / (module.replace(".", "/") + ".olean")
        if args.backend == "peach":
            command = [sys.executable, str(ROOT / "scripts/peach-lean.py"), relative]
        else:
            obj.parent.mkdir(parents=True, exist_ok=True)
            command = ["lake", "env", "lean", "-j1", "-M", "20000", "-o", str(obj), relative]
        log = result_dir / (module.replace(".", "_") + ".log")
        start = time.monotonic()
        with log.open("wb") as stream:
            result = subprocess.run(command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT)
        if result.returncode:
            raise RuntimeError(f"Lean/helper failed for {module}; see {log}")
        if sha(path) != digests[module]:
            raise RuntimeError(f"Source changed during replay: {module}")
        return {"module": module, "source_sha256": digests[module], "olean_sha256": sha(obj),
                "seconds": round(time.monotonic() - start, 3), "log": str(log)}

    try:
        selected = order[:args.max_modules] if args.max_modules else order
        pending, done, running = set(selected), set(), {}
        state["status"] = "RUNNING"
        save()
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
            try:
                while pending or running:
                    ready = [m for m in order if m in pending and deps[m] <= done]
                    for module in ready[:args.jobs - len(running)]:
                        pending.remove(module)
                        running[pool.submit(compile_one, module)] = module
                    if not running:
                        raise RuntimeError("No ready module; incomplete dependency plan")
                    finished, _ = concurrent.futures.wait(running, return_when=concurrent.futures.FIRST_COMPLETED)
                    for future in finished:
                        module = running.pop(future)
                        state["completed"].append(future.result())
                        done.add(module)
                        print(f"PASS {len(done)}/{len(selected)} {module}", flush=True)
                        save()
            except BaseException:
                # Stop submitting work; drain already running compilers before releasing the lock.
                for future in running:
                    future.cancel()
                raise
        for module in order:
            if sha(source(module)) != digests[module]:
                raise RuntimeError(f"Source snapshot changed: {module}")
        for p, digest in state["pins"].items():
            if sha(ROOT / p) != digest:
                raise RuntimeError(f"Dependency pin changed: {p}")
        state["status"] = "PASS" if len(done) == len(order) else "PARTIAL"
        save()
        print(f'{state["status"]}: {report}')
    except BaseException as error:
        state["status"], state["error"] = "FAILED", str(error)
        save()
        raise
    finally:
        lock.unlink()


if __name__ == "__main__":
    main()
