"""Diagnose a proposed Scalable dense cover before Lean checks it.

This script checks tree bookkeeping and reports sharing opportunities. Its
result is not a certificate: Lean must check the generated tree and leaves.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from statistics import median


def q(value: str) -> Fraction:
    return Fraction(value)


def inspect(path: Path, max_chunk_leaves: int = 32,
            max_chunk_vertices: int = 80) -> dict:
    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["status"] == "candidate_complete"
    nodes = payload["nodes"]
    boxes = payload["boxes"]
    roots = payload["roots"]
    witnesses = payload["witnesses"]
    assert all(node is not None for node in nodes)

    alpha_lo, alpha_hi = map(q, payload["domain"]["alpha"])
    x_lo, x_hi = map(q, payload["domain"]["row_density"])
    assert alpha_lo < alpha_hi and x_lo < x_hi
    assert q(roots[0]["xlo"]) == x_lo
    assert q(roots[-1]["xhi"]) == x_hi
    assert all(q(a["xhi"]) == q(b["xlo"])
               for a, b in zip(roots, roots[1:]))

    seen_nodes: set[int] = set()
    seen_boxes: set[int] = set()
    unique_vertices: set[tuple] = set()
    alpha_terms: set[tuple] = set()
    density_terms: set[tuple] = set()
    depth = 0

    def walk(index: int, rect: tuple[Fraction, ...], level: int) -> None:
        nonlocal depth
        assert 0 <= index < len(nodes) and index not in seen_nodes
        seen_nodes.add(index)
        depth = max(depth, level)
        node = nodes[index]
        kind = node["kind"]
        alo, ahi, xlo, xhi = rect
        if kind == "leaf":
            box_index = node["box"]
            assert 0 <= box_index < len(boxes) and box_index not in seen_boxes
            seen_boxes.add(box_index)
            box = boxes[box_index]
            actual = (*map(q, box["alpha"]), *map(q, box["row_density"]))
            assert actual == rect
            witness = box["witness"]
            assert str(witness) in witnesses
            for a in (alo, ahi):
                alpha_terms.add((witness, a))
                for x in (xlo, xhi):
                    unique_vertices.add((witness, tuple(box["outer_line"]), a, x))
            for x in (xlo, xhi):
                density_terms.add((witness, x))
            return
        assert kind in ("alpha", "density")
        cut = q(node["cut"])
        if kind == "alpha":
            assert alo < cut < ahi
            left = (alo, cut, xlo, xhi)
            right = (cut, ahi, xlo, xhi)
        else:
            assert xlo < cut < xhi
            left = (alo, ahi, xlo, cut)
            right = (alo, ahi, cut, xhi)
        walk(node["left"], left, level + 1)
        walk(node["right"], right, level + 1)

    for root in roots:
        walk(root["node"],
             (alpha_lo, alpha_hi, q(root["xlo"]), q(root["xhi"])), 0)

    assert len(seen_nodes) == len(nodes)
    assert len(seen_boxes) == len(boxes)
    assert len(witnesses) == len({box["witness"] for box in boxes})
    for item in witnesses.values():
        assert 0 < q(item["p"]) < 1
        assert 0 < q(item["y"]) < 1
        assert 0 < q(item["z"]) <= 1
        assert q(item["rho"]) > 0
        assert len(item["vector"]) == 3
        assert all(q(v) > 0 for v in item["vector"])
    counts: list[int | None] = [None] * len(nodes)

    def count(index: int) -> int:
        if counts[index] is None:
            node = nodes[index]
            counts[index] = (1 if node["kind"] == "leaf" else
                             count(node["left"]) + count(node["right"]))
        return counts[index]

    chunks: list[int] = []

    def cut(index: int) -> None:
        if count(index) <= max_chunk_leaves:
            chunks.append(index)
            return
        node = nodes[index]
        cut(node["left"])
        cut(node["right"])

    for root in roots:
        cut(root["node"])

    def chunk_vertices(index: int) -> set[tuple]:
        node = nodes[index]
        if node["kind"] != "leaf":
            return chunk_vertices(node["left"]) | chunk_vertices(node["right"])
        box = boxes[node["box"]]
        return {(box["witness"], tuple(box["outer_line"]), a, x)
                for a in box["alpha"] for x in box["row_density"]}

    chunk_vertex_counts = [len(chunk_vertices(index)) for index in chunks]

    def vertex_cut(index: int) -> tuple[set[tuple], list[int]]:
        node = nodes[index]
        if node["kind"] == "leaf":
            return chunk_vertices(index), [index]
        left_set, left_chunks = vertex_cut(node["left"])
        right_set, right_chunks = vertex_cut(node["right"])
        combined = left_set | right_set
        if len(combined) <= max_chunk_vertices:
            return combined, [index]
        return combined, left_chunks + right_chunks

    vertex_chunks = []
    for root in roots:
        vertex_chunks.extend(vertex_cut(root["node"])[1])
    vertex_chunk_sizes = [len(chunk_vertices(index)) for index in vertex_chunks]

    def dual_cut(index: int) -> tuple[set[tuple], list[int]]:
        node = nodes[index]
        if node["kind"] == "leaf":
            return chunk_vertices(index), [index]
        left_set, left_chunks = dual_cut(node["left"])
        right_set, right_chunks = dual_cut(node["right"])
        combined = left_set | right_set
        if count(index) <= max_chunk_leaves and len(combined) <= max_chunk_vertices:
            return combined, [index]
        return combined, left_chunks + right_chunks

    dual_chunks = []
    for root in roots:
        dual_chunks.extend(dual_cut(root["node"])[1])
    return {"roots": len(roots), "nodes": len(nodes),
            "boxes": len(boxes), "witnesses": len(witnesses),
            "unique_vertices": len(unique_vertices), "max_tree_depth": depth,
            "unique_alpha_kl_terms": len(alpha_terms),
            "unique_density_kl_terms": len(density_terms),
            "bytes": path.stat().st_size,
            "chunk_leaves_limit": max_chunk_leaves,
            "chunks": len(chunks),
            "median_vertices_per_chunk": median(chunk_vertex_counts),
            "max_vertices_per_chunk": max(chunk_vertex_counts),
            "chunk_vertices_limit": max_chunk_vertices,
            "vertex_chunks": len(vertex_chunks),
            "median_vertices_in_vertex_chunk": median(vertex_chunk_sizes),
            "max_leaves_in_vertex_chunk": max(count(i) for i in vertex_chunks),
            "dual_chunks": len(dual_chunks),
            "max_vertices_in_dual_chunk": max(len(chunk_vertices(i))
                                               for i in dual_chunks),
            "heaviest_dual_node": max(dual_chunks,
                                       key=lambda i: len(chunk_vertices(i)))}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--max-chunk-leaves", type=int, default=32)
    parser.add_argument("--max-chunk-vertices", type=int, default=80)
    args = parser.parse_args()
    assert args.max_chunk_leaves > 0 and args.max_chunk_vertices >= 4
    print(json.dumps(inspect(args.candidate, args.max_chunk_leaves,
                             args.max_chunk_vertices)))


if __name__ == "__main__":
    main()
