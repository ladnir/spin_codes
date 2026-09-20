#!/usr/bin/env python3
"""Reject bulk experiment artifacts and oversized tracked files."""

from __future__ import annotations

import pathlib
import subprocess
import sys


MAX_TRACKED_BYTES = 5 * 1024 * 1024
FORBIDDEN_DIRECTORIES = {
    "__pycache__",
    ".pytest_cache",
    "old",
    "out",
    "output",
    "previous_draft",
    "receipts",
    "related_work",
    "target",
}
FORBIDDEN_SUFFIXES = {
    ".bin",
    ".dll",
    ".exe",
    ".npy",
    ".npz",
    ".obj",
    ".o",
    ".parquet",
    ".pdb",
    ".prof",
    ".trace",
}
FORBIDDEN_NAME_SUFFIXES = (".stderr.txt", ".stdout.txt")


def tracked_paths() -> list[pathlib.Path]:
    output = subprocess.check_output(["git", "ls-files", "-z"])
    return [pathlib.Path(raw.decode("utf-8")) for raw in output.split(b"\0") if raw]


def main() -> int:
    failures: list[str] = []
    for path in tracked_paths():
        lowered_parts = {part.lower() for part in path.parts}
        lowered_name = path.name.lower()
        if lowered_parts & FORBIDDEN_DIRECTORIES:
            failures.append(f"forbidden bulk-output directory: {path}")
        if path.suffix.lower() in FORBIDDEN_SUFFIXES:
            failures.append(f"forbidden generated/binary suffix: {path}")
        if lowered_name.endswith(FORBIDDEN_NAME_SUFFIXES):
            failures.append(f"forbidden command-output log: {path}")
        if path.is_file() and path.stat().st_size > MAX_TRACKED_BYTES:
            failures.append(
                f"tracked file exceeds 5 MiB ({path.stat().st_size} bytes): {path}"
            )

    if failures:
        print("Repository hygiene check failed:", file=sys.stderr)
        for failure in failures:
            print(f"- {failure}", file=sys.stderr)
        return 1

    print("Repository hygiene check passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
