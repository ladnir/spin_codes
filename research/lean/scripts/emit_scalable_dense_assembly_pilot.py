"""Emit two independently checked chunks and their Lean split claim."""

from __future__ import annotations

import argparse
from fractions import Fraction as Q
import json
from pathlib import Path

from emit_scalable_dense_geometry_pilot import emit, rat


def make_pilot(candidate: Path, lean_root: Path,
               max_child_leaves: int) -> dict:
    payload = json.loads(candidate.read_text(encoding="utf-8"))
    assert payload["status"] == "candidate_complete"
    nodes = payload["nodes"]
    counts = [None] * len(nodes)

    def count(index: int) -> int:
        if counts[index] is None:
            node = nodes[index]
            counts[index] = (1 if node["kind"] == "leaf" else
                             count(node["left"]) + count(node["right"]))
        return counts[index]

    best = None

    def search(index: int, rect: tuple[Q, ...]) -> None:
        nonlocal best
        node = nodes[index]
        if node["kind"] == "leaf":
            return
        left_id, right_id = node["left"], node["right"]
        left_count, right_count = count(left_id), count(right_id)
        if (left_count <= max_child_leaves and
                right_count <= max_child_leaves):
            score = left_count + right_count
            if best is None or score > best[0]:
                best = (score, index, rect)
            return
        alo, ahi, xlo, xhi = rect
        cut = Q(node["cut"])
        if node["kind"] == "alpha":
            search(left_id, (alo, cut, xlo, xhi))
            search(right_id, (cut, ahi, xlo, xhi))
        else:
            search(left_id, (alo, ahi, xlo, cut))
            search(right_id, (alo, ahi, cut, xhi))

    alo, ahi = map(Q, payload["domain"]["alpha"])
    for root in payload["roots"]:
        search(root["node"], (alo, ahi, Q(root["xlo"]), Q(root["xhi"])))
    assert best is not None
    size, parent_id, rect = best
    parent = nodes[parent_id]
    cut = Q(parent["cut"])
    generated = lean_root / "SpinCodes" / "UniformMixer" / "Generated"
    scratch = lean_root / ".lake"
    for suffix, node_id in (("Left", parent["left"]),
                            ("Right", parent["right"])):
        emit(candidate, scratch / f"Geometry{suffix}.lean", max_child_leaves,
             generated / f"ScalableDensePilot{suffix}.lean", True,
             node_id)
    rect_text = ", ".join(rat(str(x)) for x in rect)
    left = "Spin.UniformMixer.ScalableDensePilotLeft"
    right = "Spin.UniformMixer.ScalableDensePilotRight"
    combine = ("denseClaim_alpha" if parent["kind"] == "alpha"
               else "denseClaim_density")
    source = f"""import SpinCodes.UniformMixer.ScalableDenseCoverCore
import SpinCodes.UniformMixer.Generated.ScalableDensePilotLeft
import SpinCodes.UniformMixer.Generated.ScalableDensePilotRight

/-! Two independently checked numerical chunks compose across one exact cut. -/

namespace Spin.UniformMixer.ScalableDenseAssemblyPilot

open Spin.Structured.DenseGeometry

def root : Rect := ⟨{rect_text}⟩

theorem root_claim :
    DenseClaim (DenseTarget (1 / 10000000) 10000000) root := by
  apply {combine} (DenseTarget (1 / 10000000) 10000000) root {rat(str(cut))}
  · simpa [root, {left}.root] using {left}.root_claim_global
  · simpa [root, {right}.root] using {right}.root_claim_global

#print axioms root_claim

end Spin.UniformMixer.ScalableDenseAssemblyPilot
"""
    destination = lean_root / "SpinCodes" / "UniformMixer" / \
        "ScalableDenseAssemblyPilot.lean"
    destination.write_bytes(source.encode("utf-8"))
    return {"parent_node": parent_id, "leaves": size,
            "left_leaves": count(parent["left"]),
            "right_leaves": count(parent["right"]),
            "kind": parent["kind"], "destination": str(destination)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("lean_root", type=Path)
    parser.add_argument("--max-child-leaves", type=int, default=8)
    args = parser.parse_args()
    print(json.dumps(make_pilot(args.candidate, args.lean_root,
                                args.max_child_leaves)))


if __name__ == "__main__":
    main()
