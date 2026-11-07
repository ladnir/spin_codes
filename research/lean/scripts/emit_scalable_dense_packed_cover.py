"""Partition the full proposed dense cover and emit Lean claim assembly.

The JSON and this generator are untrusted. Every generated part is checked by
`denseChunkCheckShared`; the assembly uses exact split lemmas. The manifest
records source paths and a digest of the candidate for reproducibility.
"""

from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

from emit_scalable_dense_geometry_pilot import rat
from emit_scalable_dense_packed_chunk import emit as emit_chunk


TARGET = "(DenseTarget (1 / 10000000) 10000000)"


def rect_term(bounds: tuple[Q, ...]) -> str:
    return "⟨" + ", ".join(rat(str(v)) for v in bounds) + "⟩"


def generate(candidate: Path, output: Path, leaf_limit: int,
             vertex_limit: int, emit_parts: bool,
             selected_parts: set[int] | None = None) -> dict:
    raw = candidate.read_bytes()
    payload = json.loads(raw)
    assert payload["status"] == "candidate_complete"
    nodes = payload["nodes"]
    boxes = payload["boxes"]
    roots = payload["roots"]
    a0, a1 = map(Q, payload["domain"]["alpha"])
    counts = [None] * len(nodes)

    def count(index: int) -> int:
        if counts[index] is None:
            node = nodes[index]
            counts[index] = (1 if node["kind"] == "leaf" else
                             count(node["left"]) + count(node["right"]))
        return counts[index]

    parts: list[dict] = []

    def partition(index: int, bounds: tuple[Q, ...]):
        node = nodes[index]
        start = len(parts)
        if node["kind"] == "leaf":
            box = boxes[node["box"]]
            vertices = {(box["witness"], a, x)
                        for a in box["alpha"] for x in box["row_density"]}
            children = None
        else:
            lo_a, hi_a, lo_x, hi_x = bounds
            cut = Q(node["cut"])
            if node["kind"] == "alpha":
                left_bounds = (lo_a, cut, lo_x, hi_x)
                right_bounds = (cut, hi_a, lo_x, hi_x)
            else:
                left_bounds = (lo_a, hi_a, lo_x, cut)
                right_bounds = (lo_a, hi_a, cut, hi_x)
            left_set, left_tree = partition(node["left"], left_bounds)
            right_set, right_tree = partition(node["right"], right_bounds)
            vertices = left_set | right_set
            children = (left_tree, right_tree)
        if count(index) <= leaf_limit and len(vertices) <= vertex_limit:
            del parts[start:]
            part_id = len(parts)
            parts.append({"id": part_id, "node": index,
                          "bounds": [str(v) for v in bounds],
                          "leaves": count(index),
                          "vertices": len(vertices)})
            return vertices, ("part", part_id, bounds)
        assert children is not None
        return vertices, (node["kind"], Q(node["cut"]), bounds,
                          children[0], children[1])

    segments = []
    for number, root in enumerate(roots):
        bounds = (a0, a1, Q(root["xlo"]), Q(root["xhi"]))
        _, tree = partition(root["node"], bounds)
        segments.append((number, bounds, tree))

    def ids(tree) -> list[int]:
        if tree[0] == "part":
            return [tree[1]]
        return ids(tree[3]) + ids(tree[4])

    def claim_expr(tree, prefix: str) -> str:
        if tree[0] == "part":
            return f"h{prefix}{tree[1]}"
        kind, cut, bounds, left, right = tree
        combine = "denseClaim_alpha" if kind == "alpha" else "denseClaim_density"
        return (f"({combine} {TARGET} {rect_term(bounds)} {rat(str(cut))} "
                f"({claim_expr(left, prefix)}) ({claim_expr(right, prefix)}))")

    source_dir = output / "sources"
    source_dir.mkdir(parents=True, exist_ok=True)
    module_prefix = "SpinCodes.UniformMixer.Generated.FullCover"
    namespace_prefix = "Spin.UniformMixer"
    module_paths = {}
    for part in parts:
        part_id = part["id"]
        module_paths[part_id] = f"{module_prefix}.Part{part_id}"
        part["source"] = str(source_dir / f"Part{part_id}.lean")
        if emit_parts or (selected_parts is not None and part_id in selected_parts):
            result = emit_chunk(candidate, Path(part["source"]), leaf_limit,
                                selected_node=part["node"], shared=True,
                                payload=payload,
                                selected_bounds=tuple(map(Q, part["bounds"])))
            assert result["leaves"] == part["leaves"]

    segment_paths = []
    for number, bounds, tree in segments:
        path = source_dir / f"Segment{number}.lean"
        segment_paths.append(path)
        part_ids = ids(tree)
        imports = "\n".join(f"import {module_paths[i]}" for i in part_ids)
        claims = []
        for i in part_ids:
            namespace = f"{namespace_prefix}.Part{i}"
            part_bounds = tuple(map(Q, parts[i]["bounds"]))
            claims.append(
                f"  have hPart{i} : DenseClaim {TARGET} "
                f"{rect_term(part_bounds)} := by\n"
                f"    simpa [{namespace}.root] using {namespace}.root_claim")
        claims_text = "\n".join(claims)
        source = f"""import SpinCodes.UniformMixer.ScalableDenseCoverCore
{imports}

namespace Spin.UniformMixer.Segment{number}

open Spin.Structured.DenseGeometry
set_option maxRecDepth 100000
set_option maxHeartbeats 1000000

def root : Rect := {rect_term(bounds)}

theorem root_claim : DenseClaim {TARGET} root := by
{claims_text}
  exact {claim_expr(tree, 'Part')}

#print axioms root_claim

end Spin.UniformMixer.Segment{number}
"""
        path.write_bytes(source.encode("utf-8"))

    def global_expr(lo: int, hi: int) -> str:
        if hi - lo == 1:
            return f"hSegment{lo}"
        mid = (lo + hi) // 2
        bounds = (a0, a1, Q(roots[lo]["xlo"]), Q(roots[hi - 1]["xhi"]))
        cut = Q(roots[mid]["xlo"])
        return (f"(denseClaim_density {TARGET} {rect_term(bounds)} "
                f"{rat(str(cut))} ({global_expr(lo, mid)}) "
                f"({global_expr(mid, hi)}))")

    global_imports = "\n".join(
        f"import {module_prefix}.Segment{i}" for i in range(len(segments)))
    global_claims = "\n".join(
        f"  have hSegment{i} : DenseClaim {TARGET} "
        f"{rect_term(bounds)} := by\n"
        f"    simpa [Spin.UniformMixer.Segment{i}.root] using "
        f"Spin.UniformMixer.Segment{i}.root_claim"
        for i, bounds, _ in segments)
    global_bounds = (a0, a1, Q(roots[0]["xlo"]), Q(roots[-1]["xhi"]))
    global_source = f"""import SpinCodes.UniformMixer.ScalableDenseCoverCore
{global_imports}

namespace Spin.UniformMixer.ScalableDenseFullCover

open Spin.Structured.DenseGeometry
set_option maxRecDepth 100000
set_option maxHeartbeats 1000000

def root : Rect := {rect_term(global_bounds)}

theorem root_claim : DenseClaim {TARGET} root := by
{global_claims}
  exact {global_expr(0, len(segments))}

theorem rates : ScalableDenseRates (1 / 10000000) 10000000 := by
  apply scalableDenseRates_of_claim
  simpa [root] using root_claim

#print axioms rates

end Spin.UniformMixer.ScalableDenseFullCover
"""
    (source_dir / "ScalableDenseFullCover.lean").write_bytes(
        global_source.encode("utf-8"))
    manifest = {"candidate_sha256": hashlib.sha256(raw).hexdigest(),
                "status": "source_emitted" if emit_parts else
                          "selected_sources_emitted" if selected_parts else "plan_only",
                "leaf_limit": leaf_limit, "vertex_limit": vertex_limit,
                "segments": len(segments), "parts": parts,
                "total_leaves": sum(p["leaves"] for p in parts),
                "assembly_sources": [str(p) for p in segment_paths],
                "global_source": str(source_dir / "ScalableDenseFullCover.lean")}
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8")
    return {"parts": len(parts), "segments": len(segments),
            "total_leaves": manifest["total_leaves"],
            "candidate_sha256": manifest["candidate_sha256"],
            "status": manifest["status"], "output": str(output)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--leaf-limit", type=int, default=256)
    parser.add_argument("--vertex-limit", type=int, default=400)
    parser.add_argument("--emit-parts", action="store_true")
    parser.add_argument("--emit-part", type=int, action="append", default=[])
    args = parser.parse_args()
    print(json.dumps(generate(args.candidate, args.output, args.leaf_limit,
                              args.vertex_limit, args.emit_parts,
                              set(args.emit_part))))


if __name__ == "__main__":
    main()
