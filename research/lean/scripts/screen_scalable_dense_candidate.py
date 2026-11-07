"""Screen rational Scalable dense proposals before Lean checks them.

This catches bad generated candidates early; it is not part of the proof.
"""

from __future__ import annotations

import argparse
from fractions import Fraction as Q
import json
from pathlib import Path
import sys

from flint import ctx

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from research.workstreams.turbo_spin_migration.shared_inner_dense_certificate import (
    kl, num, upper,
)


SPECTRUM = ((0, 1), (16, 2604), (24, 291648), (28, 888832),
            (32, 1828134), (36, 888832), (40, 291648),
            (48, 2604), (64, 1))
LIVE = tuple(w for w, _ in SPECTRUM if w)


def exact_transfer_check(data: dict) -> bool:
    p, y, z, rho = map(Q, (data["p"], data["y"], data["z"], data["rho"]))
    v = tuple(map(Q, data["vector"]))
    beta = p * y
    g0 = 1 - beta + beta * z
    g1 = beta + (1 - beta) * z
    h1 = 1 - beta - beta * z
    weight = lambda w: g0 ** (64 - w) * g1**w
    H = sum(count * g0 ** (64 - w) * h1**w
            for w, count in SPECTRUM) / 4194304
    G = g0**64 - H
    F = max(weight(w) for w in LIVE)
    U = (sum(count * weight(w) for w, count in SPECTRUM)
         - weight(0)) / 4194303
    u = 1 - sum(count * (1 - 2 * beta)**w
                for w, count in SPECTRUM) / 4194304
    rows = (H * v[0] + G * v[1],
            min(F, u) / 4194303 * v[0] + F * v[2],
            min(U, u) / 4194303 * v[0] + U * v[2])
    return all(row <= rho * vi for row, vi in zip(rows, v))


def exponent_bound(box: dict, witness: dict) -> Q:
    p, y, z, rho = map(Q, (witness["p"], witness["y"],
                            witness["z"], witness["rho"]))
    slope, intercept = map(Q, box["outer_line"])
    return max(upper(num(alpha * (slope * x + intercept))
                     + kl(alpha, p) + num(alpha) * kl(x, y)
                     + num(rho).log() / 64
                     - num(Q(11, 100)) * num(z).log())
               for alpha in map(Q, box["alpha"])
               for x in map(Q, box["row_density"]))


def screen(path: Path, precision: int) -> dict:
    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["status"] == "candidate_complete"
    ctx.prec = precision
    witnesses = payload["witnesses"]
    bad_transfer = [key for key, data in witnesses.items()
                    if not exact_transfer_check(data)]
    bad_prefactor = [key for key, data in witnesses.items()
                     if any(10_000_000 * Q(v) < 1
                            for v in data["vector"])]
    target = -Q(1, 10_000_000)
    worst = None
    worst_index = None
    bad_vertices = []
    for index, box in enumerate(payload["boxes"]):
        bound = exponent_bound(box, witnesses[str(box["witness"])])
        if worst is None or bound > worst:
            worst, worst_index = bound, index
        if bound > target:
            bad_vertices.append(index)
    return {"boxes": len(payload["boxes"]),
            "witnesses": len(witnesses),
            "bad_transfer": bad_transfer,
            "bad_prefactor": bad_prefactor,
            "bad_vertices_count": len(bad_vertices),
            "bad_vertices_first": bad_vertices[:10],
            "worst_box": worst_index,
            "worst_exponent_upper": str(worst),
            "target": str(target)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--precision", type=int, default=256)
    args = parser.parse_args()
    print(json.dumps(screen(args.candidate, args.precision)))


if __name__ == "__main__":
    main()
