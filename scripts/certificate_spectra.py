"""Load, structurally validate, and fingerprint finite-certificate spectra."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _validate(
    rows: tuple[tuple[int, int], ...],
    *,
    name: str,
    length: int,
    dimension: int,
    minimum_distance: int,
) -> tuple[tuple[int, int], ...]:
    if not rows:
        raise SystemExit(f"{name}: empty spectrum")
    weights = [weight for weight, _count in rows]
    if weights != sorted(set(weights)):
        raise SystemExit(f"{name}: weights are not unique and increasing")
    if any(weight < 0 or weight > length or count <= 0 for weight, count in rows):
        raise SystemExit(f"{name}: invalid weight or count")
    table = dict(rows)
    if table.get(0) != 1:
        raise SystemExit(f"{name}: zero-word multiplicity is not one")
    if sum(table.values()) != 1 << dimension:
        raise SystemExit(f"{name}: multiplicities do not sum to 2^{dimension}")
    observed_minimum = min(weight for weight, count in rows if weight > 0 and count)
    if observed_minimum != minimum_distance:
        raise SystemExit(
            f"{name}: minimum positive weight {observed_minimum} != {minimum_distance}"
        )
    if any(table.get(length - weight, 0) != count for weight, count in rows):
        raise SystemExit(f"{name}: spectrum is not complement symmetric")
    return rows


def load_rm512_spectrum(path: Path) -> tuple[tuple[int, int], ...]:
    with path.open(newline="") as handle:
        rows = tuple(
            (int(row["weight"]), int(row["count"])) for row in csv.DictReader(handle)
        )
    return _validate(
        rows,
        name="RM(4,9) spectrum",
        length=512,
        dimension=256,
        minimum_distance=32,
    )


def load_ebch128_spectrum(path: Path) -> tuple[tuple[int, int], ...]:
    rows: list[tuple[int, int]] = []
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        weight, count = stripped.split()
        rows.append((int(weight), int(count)))
    return _validate(
        tuple(rows),
        name="extended BCH [128,64,22] spectrum",
        length=128,
        dimension=64,
        minimum_distance=22,
    )
