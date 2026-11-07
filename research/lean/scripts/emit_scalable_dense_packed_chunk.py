"""Emit a bounded compact Scalable dense-cover chunk for Lean checking.

The JSON and this generator are untrusted. Lean evaluates `denseChunkCheck`,
which verifies the split geometry and every referenced numerical leaf.
"""

from __future__ import annotations

import argparse
from fractions import Fraction as Q
import json
from pathlib import Path

from emit_scalable_dense_geometry_pilot import emit as choose_subtree, qinput, rat


def emit(candidate: Path, destination: Path, max_leaves: int,
         selected_node: int | None = None, shared: bool = False,
         payload: dict | None = None,
         selected_bounds: tuple[Q, ...] | None = None) -> dict:
    if payload is None:
        payload = json.loads(candidate.read_text(encoding="utf-8"))
    assert payload["status"] == "candidate_complete"
    nodes = payload["nodes"]
    boxes = payload["boxes"]
    if selected_bounds is None:
        scratch = candidate.parent / "ScalableDensePackedGeometryPilot.lean"
        selection = choose_subtree(candidate, scratch, max_leaves,
                                   selected_node=selected_node)
        root_id = selection["root_node"]
    else:
        assert selected_node is not None
        root_id = selected_node

        def leaf_count(index: int) -> int:
            node = nodes[index]
            return (1 if node["kind"] == "leaf" else
                    leaf_count(node["left"]) + leaf_count(node["right"]))

        selection = {"leaves": leaf_count(root_id)}
    domain_a0, domain_a1 = map(Q, payload["domain"]["alpha"])
    bounds = selected_bounds

    def locate(index: int, rect: tuple[Q, ...]) -> None:
        nonlocal bounds
        if index == root_id:
            bounds = rect
            return
        node = nodes[index]
        if node["kind"] == "leaf":
            return
        a0, a1, x0, x1 = rect
        cut = Q(node["cut"])
        if node["kind"] == "alpha":
            locate(node["left"], (a0, cut, x0, x1))
            locate(node["right"], (cut, a1, x0, x1))
        else:
            locate(node["left"], (a0, a1, x0, cut))
            locate(node["right"], (a0, a1, cut, x1))

    if bounds is None:
        for root in payload["roots"]:
            locate(root["node"], (domain_a0, domain_a1,
                                  Q(root["xlo"]), Q(root["xhi"])))
    assert bounds is not None
    local_boxes = []

    def tree(index: int) -> str:
        node = nodes[index]
        if node["kind"] == "leaf":
            local_boxes.append(boxes[node["box"]])
            return f"(.leaf {len(local_boxes) - 1})"
        constructor = ".alpha" if node["kind"] == "alpha" else ".density"
        return (f"({constructor} {rat(node['cut'])} "
                f"{tree(node['left'])} {tree(node['right'])})")

    tree_term = tree(root_id)
    assert len(local_boxes) == selection["leaves"]
    line = tuple(local_boxes[0]["outer_line"])
    assert all(tuple(box["outer_line"]) == line for box in local_boxes)
    witness_ids = sorted({box["witness"] for box in local_boxes})
    local_witness = {wid: i for i, wid in enumerate(witness_ids)}

    def witness_term(wid: int) -> str:
        w = payload["witnesses"][str(wid)]
        vector = ", ".join(qinput(v) for v in w["vector"])
        return (f"⟨{qinput(w['p'])}, {qinput(w['y'])}, "
                f"{qinput(w['rho'])}, {qinput(w['z'])}, ![{vector}]⟩")

    def ref_term(box: dict) -> str:
        a0, a1 = box["alpha"]
        x0, x1 = box["row_density"]
        return (f"⟨{local_witness[box['witness']]}, "
                f"⟨{qinput(a0)}, {qinput(a1)}, "
                f"{qinput(x0)}, {qinput(x1)}⟩⟩")

    name = destination.stem
    witness_text = ",\n".join("  " + witness_term(wid) for wid in witness_ids)
    ref_text = ",\n".join("  " + ref_term(box) for box in local_boxes)
    leaf_check = "denseRefBoxCheck" if shared else "denseLeafRefCheck"
    leaf_proofs = "\n".join(
        f"theorem leaf{i}_checked :\n"
        f"    {leaf_check} line witnesses {ref_term(box)} = true := by\n"
        "  decide +kernel\n"
        for i, box in enumerate(local_boxes))
    leaf_names = ", ".join(f"leaf{i}_checked" for i in range(len(local_boxes)))
    if shared:
        witness_proofs = "\n".join(
            f"theorem witness{i}_checked :\n"
            f"    denseWitnessCheck {witness_term(wid)} = true := by\n"
            "  decide +kernel\n"
            for i, wid in enumerate(witness_ids))
        witness_names = ", ".join(
            f"witness{i}_checked" for i in range(len(witness_ids)))
        preliminary = f"""theorem line_checked : decide line.Valid = true := by
  decide +kernel

{witness_proofs}

theorem witnesses_checked : witnesses.all denseWitnessCheck = true := by
  simp [witnesses, List.all, {witness_names}]
"""
        checked_body = """theorem checked :
    denseChunkCheckShared line witnesses refs tree root = true := by
  simpa only [denseChunkCheckShared, Bool.and_eq_true] using
    (show ((check (refs.map (fun r => r.box.rect)) tree root = true ∧
      decide line.Valid = true) ∧ witnesses.all denseWitnessCheck = true) ∧
      refs.all (denseRefBoxCheck line witnesses) = true from
      ⟨⟨⟨geometry_checked, line_checked⟩, witnesses_checked⟩,
        numerical_checked⟩)
"""
        chunk_check = "denseChunkCheckShared"
        chunk_sound = "denseChunkCheckShared_sound"
    else:
        preliminary = ""
        checked_body = """theorem checked :
    denseChunkCheck line witnesses refs tree root = true := by
  simpa only [denseChunkCheck, Bool.and_eq_true] using
    (show check (refs.map (fun r => r.box.rect)) tree root = true ∧
      refs.all (denseLeafRefCheck line witnesses) = true from
      ⟨geometry_checked, numerical_checked⟩)
"""
        chunk_check = "denseChunkCheck"
        chunk_sound = "denseChunkCheck_sound"
    root_text = ", ".join(rat(str(v)) for v in bounds)
    source = f"""import SpinCodes.UniformMixer.ScalableDensePackedLeaf

namespace Spin.UniformMixer.{name}

open Spin.Structured.DenseGeometry

set_option maxHeartbeats 1000000
set_option maxRecDepth 100000

def line : DenseLine := ⟨{qinput(line[0])}, {qinput(line[1])}⟩
def witnesses : List DenseWitness := [
{witness_text}
]
def refs : List DenseLeafRef := [
{ref_text}
]
def tree : Tree := {tree_term}
def root : Rect := ⟨{root_text}⟩

theorem geometry_checked :
    check (refs.map (fun r => r.box.rect)) tree root = true := by
  decide +kernel

{preliminary}

{leaf_proofs}

theorem numerical_checked :
    refs.all ({leaf_check} line witnesses) = true := by
  simp [refs, List.all, {leaf_names}]

{checked_body}

theorem root_claim :
    DenseClaim (DenseTarget (1 / 10000000) 10000000) root :=
  {chunk_sound} line witnesses refs tree root checked

#print axioms root_claim

end Spin.UniformMixer.{name}
"""
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(source.encode("utf-8"))
    return {"leaves": len(local_boxes), "witnesses": len(witness_ids),
            "root_node": root_id, "source_bytes": len(source.encode("utf-8")),
            "shared": shared, "check": chunk_check,
            "destination": str(destination)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--max-leaves", type=int, default=8)
    parser.add_argument("--node-id", type=int)
    parser.add_argument("--shared", action="store_true")
    args = parser.parse_args()
    print(json.dumps(emit(args.candidate, args.destination,
                          args.max_leaves, args.node_id, args.shared)))


if __name__ == "__main__":
    main()
