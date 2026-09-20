#!/usr/bin/env python3
"""Run the four endpoint-aware dimension-33 frontier certificates sequentially."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import run_riffle_dp_2lap_g4_component_dimension32 as base


ROOT = Path(__file__).resolve().parents[1]
EXPLORATIONS = ROOT / "explorations"
FRONTIER = (
    (0x2E, (2, 4, 9, 18)),
    (0x35, (1, 4, 10, 18)),
    (0x4C, (4, 9, 20)),
    (0x53, (1, 2, 10, 20)),
)
DIMENSION = 33


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def receipt_path(
    mode: str, window: int, mask: int, value: int, probe_only: bool
) -> Path:
    stage = "probe_" if probe_only else ""
    return EXPLORATIONS / (
        f"riffle_dp_2lap_g4_component33_w{window}_"
        f"{stage}{mode}_s{mask:02x}_v{value:02d}.json"
    )


def output_argument() -> Path:
    try:
        index = sys.argv.index("--output")
        return Path(sys.argv[index + 1])
    except (ValueError, IndexError) as error:
        raise RuntimeError("dimension-33 runner: missing --output") from error


def main() -> None:
    base.FRONTIER = FRONTIER
    base.DIMENSION = DIMENSION
    base.receipt_path = receipt_path
    output = output_argument()
    exit_code = 0
    try:
        base.main()
    except SystemExit as result:
        exit_code = int(result.code or 0)

    payload = json.loads(output.read_text())
    payload["schema"] = "riffle-dp-2lap-g4-component-dimension33-run-v1"
    payload["runner_source_sha256"] = digest(Path(__file__).resolve())
    payload["base_runner_source_sha256"] = digest(Path(base.__file__).resolve())
    payload["scope_limitation"] = (
        "This run covers only the four dimension-33 frontier supports. "
        "The conclusion through dimension 33 also uses the audited "
        "dimension-at-most-32 certificate."
    )
    output.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"dimension33_output={output}")
    raise SystemExit(exit_code)


if __name__ == "__main__":
    main()
