"""Check and replay the paper's finite certificates.

The manifest pins the trusted verifier source tree and every frozen input.
Commands run serially and never invoke a marker-generation subcommand.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import subprocess
import sys
import time
from typing import Any


ROOT = Path(__file__).resolve().parent.parent
DEFAULT_MANIFEST = ROOT / "certificate_manifest.json"
SCHEMA = "expander-codes-certificate-manifest-v1"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def trusted_source_paths() -> list[Path]:
    scripts = ROOT / "scripts"
    return sorted(
        path
        for path in scripts.glob("*.py")
        if not path.name.startswith("test_") and path.name != Path(__file__).name
    )


def source_tree_sha256() -> str:
    digest = hashlib.sha256()
    for path in trusted_source_paths():
        relative = path.relative_to(ROOT).as_posix().encode("utf-8")
        contents = path.read_bytes()
        digest.update(len(relative).to_bytes(8, "big"))
        digest.update(relative)
        digest.update(len(contents).to_bytes(8, "big"))
        digest.update(contents)
    return digest.hexdigest()


def load_manifest(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("schema") != SCHEMA:
        raise ValueError(f"unsupported manifest schema: {data.get('schema')!r}")
    return data


def resolve_artifact(relative: str) -> Path:
    path = (ROOT / relative).resolve()
    path.relative_to(ROOT)
    return path


def check_manifest(data: dict[str, Any], *, strict_versions: bool) -> list[str]:
    errors: list[str] = []
    expected_tree = data.get("trusted_source_tree", {}).get("sha256")
    actual_tree = source_tree_sha256()
    if expected_tree != actual_tree:
        errors.append(
            f"trusted source digest mismatch: expected {expected_tree}, got {actual_tree}"
        )

    artifacts = data.get("artifacts")
    if not isinstance(artifacts, dict):
        errors.append("artifacts must be an object")
        artifacts = {}
    for relative, expected in sorted(artifacts.items()):
        try:
            path = resolve_artifact(relative)
        except ValueError:
            errors.append(f"artifact escapes repository root: {relative}")
            continue
        if not path.is_file():
            errors.append(f"missing artifact: {relative}")
            continue
        actual = sha256_file(path)
        if actual != expected:
            errors.append(
                f"artifact digest mismatch for {relative}: expected {expected}, got {actual}"
            )

    seen: set[str] = set()
    pinned_artifacts = set(artifacts)
    for entry in data.get("entries", []):
        identifier = entry.get("id")
        if not isinstance(identifier, str) or not identifier:
            errors.append("every entry needs a nonempty id")
            continue
        if identifier in seen:
            errors.append(f"duplicate entry id: {identifier}")
        seen.add(identifier)
        if not entry.get("claims"):
            errors.append(f"entry {identifier} has no claim mapping")
        for command in entry.get("commands", []):
            argv = command.get("argv")
            if not isinstance(argv, list) or not argv:
                errors.append(f"entry {identifier} contains an invalid command")
                continue
            if argv[0] != "python":
                errors.append(f"entry {identifier} must invoke Python without a shell")
            for token in argv[1:]:
                if token.endswith((".py", ".json")):
                    try:
                        candidate = resolve_artifact(token)
                    except ValueError:
                        errors.append(f"entry {identifier} path escapes root: {token}")
                        continue
                    if not candidate.exists():
                        errors.append(f"entry {identifier} references missing path: {token}")
                    if token.endswith(".json") and token not in pinned_artifacts:
                        errors.append(
                            f"entry {identifier} references unpinned artifact: {token}"
                        )

    tested = data.get("tested_environment", {})
    for distribution, expected in tested.get("packages", {}).items():
        try:
            actual = importlib.metadata.version(distribution)
        except importlib.metadata.PackageNotFoundError:
            errors.append(f"missing required package: {distribution}")
            continue
        if strict_versions and actual != expected:
            errors.append(
                f"package version mismatch for {distribution}: "
                f"expected {expected}, got {actual}"
            )
    if strict_versions:
        expected_python = tested.get("python")
        actual_python = ".".join(map(str, sys.version_info[:3]))
        if expected_python and actual_python != expected_python:
            errors.append(
                f"Python version mismatch: expected {expected_python}, got {actual_python}"
            )
    return errors


def run_command(identifier: str, command: dict[str, Any]) -> None:
    argv = list(command["argv"])
    argv[0] = sys.executable
    print(f"[{identifier}] {' '.join(command['argv'])}", flush=True)
    started = time.monotonic()
    process = subprocess.Popen(
        argv,
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    assert process.stdout is not None
    output: list[str] = []
    for line in process.stdout:
        print(line, end="", flush=True)
        output.append(line)
    code = process.wait()
    elapsed = time.monotonic() - started
    if code != 0:
        raise RuntimeError(f"{identifier} failed with exit code {code}")
    joined = "".join(output)
    for expected in command.get("expect", []):
        if expected not in joined:
            raise RuntimeError(
                f"{identifier} succeeded but omitted expected output {expected!r}"
            )
    print(f"[{identifier}] verified in {elapsed:.1f} s", flush=True)


def selected_entries(data: dict[str, Any], selector: str) -> list[dict[str, Any]]:
    entries = [entry for entry in data["entries"] if entry.get("kind") == "rigorous"]
    if selector == "all":
        return entries
    selected = [
        entry
        for entry in entries
        if selector == entry["id"] or selector in entry.get("groups", [])
    ]
    if not selected:
        raise ValueError(f"no rigorous entry or group named {selector!r}")
    return selected


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument(
        "--run",
        metavar="GROUP_OR_ID",
        help="run matching rigorous entries serially; use 'all' for every entry",
    )
    parser.add_argument("--strict-versions", action="store_true")
    parser.add_argument("--list", action="store_true")
    args = parser.parse_args()

    data = load_manifest(args.manifest.resolve())
    errors = check_manifest(data, strict_versions=args.strict_versions)
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)

    print(
        f"manifest integrity verified: {len(data['entries'])} claim groups, "
        f"{len(data['artifacts'])} frozen artifacts"
    )
    if args.list:
        for entry in data["entries"]:
            groups = ",".join(entry.get("groups", [])) or "-"
            print(f"{entry['id']}: {entry['kind']} [{groups}]")
    if args.run:
        entries = selected_entries(data, args.run)
        completed_commands: set[tuple[str, ...]] = set()
        for entry in entries:
            for command in entry["commands"]:
                key = tuple(command["argv"])
                if key in completed_commands:
                    print(
                        f"[{entry['id']}] skipped duplicate command already verified",
                        flush=True,
                    )
                    continue
                run_command(entry["id"], command)
                completed_commands.add(key)
        print(f"verified {len(entries)} manifest entries serially")


if __name__ == "__main__":
    main()
