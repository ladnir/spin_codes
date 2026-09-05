#!/usr/bin/env python3
"""Utilities for published BCH weight-distribution tables.

The Okayama/Kusaka tables store distributions as whitespace-separated
weight/count pairs.  This helper keeps the first affine-coset validation rung
reproducible: read an extended primitive BCH table, derive the punctured
primitive spectrum using the extension relation, and compare with a published
primitive table.

For an extended binary primitive BCH code of length N=n+1 whose automorphism
group is transitive on coordinates, the primitive length-n distribution A is
recovered from the extended distribution B by

    A_{w-1} = B_w * w / N,
    A_w     = B_w * (N - w) / N

for even extended weights w, with the all-zero term handled separately.  This
is the relation used by the length-128 Desaki--Fujiwara--Kasami tables.
"""

from __future__ import annotations

import argparse
import csv
import html
import math
import re
import ssl
import urllib.request
from pathlib import Path


def read_text(path_or_url: str, insecure_url: bool = False) -> str:
    if path_or_url.startswith(("http://", "https://")):
        context = ssl._create_unverified_context() if insecure_url else None
        with urllib.request.urlopen(path_or_url, timeout=30, context=context) as response:
            return response.read().decode("utf-8")
    return Path(path_or_url).read_text()


def parse_wd_table(text: str) -> tuple[str, str, dict[int, int]]:
    words = text.split()
    if not words or words[0] != "#":
        raise ValueError("expected table to start with '#'")
    if len(words) < 4:
        raise ValueError("table header is too short")
    name = words[1]
    source = words[2]
    pairs = words[3:]
    if len(pairs) % 2 != 0:
        raise ValueError("weight/count payload has odd length")
    spectrum: dict[int, int] = {}
    for i in range(0, len(pairs), 2):
        weight = int(pairs[i])
        count = int(pairs[i + 1])
        if count:
            spectrum[weight] = spectrum.get(weight, 0) + count
    return name, source, spectrum


def extended_to_primitive(extended: dict[int, int], length: int) -> dict[int, int]:
    out: dict[int, int] = {}
    for weight, count in sorted(extended.items()):
        if weight == 0:
            out[0] = out.get(0, 0) + count
            continue
        left_num = count * weight
        right_num = count * (length - weight)
        if left_num % length != 0 or right_num % length != 0:
            raise ValueError(f"nonintegral split at extended weight {weight}")
        left = left_num // length
        right = right_num // length
        if left:
            out[weight - 1] = out.get(weight - 1, 0) + left
        if right:
            out[weight] = out.get(weight, 0) + right
    return out


def summarize_spectrum(spectrum: dict[int, int]) -> dict[str, str]:
    total = sum(spectrum.values())
    nonzero = [w for w, c in spectrum.items() if w > 0 and c]
    max_weight = max(spectrum) if spectrum else 0
    dmin = min(nonzero) if nonzero else None
    return {
        "length": str(max_weight),
        "total_log2": f"{math.log2(total):.6f}" if total > 0 else "-inf",
        "minimum_weight": "" if dmin is None else str(dmin),
        "nonzero_weights": str(sum(1 for c in spectrum.values() if c)),
    }


def write_csv(path: Path, spectrum: dict[int, int]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["weight", "count"])
        for weight in sorted(spectrum):
            count = spectrum[weight]
            if count:
                writer.writerow([weight, count])


def compare_spectra(left: dict[int, int], right: dict[int, int]) -> list[tuple[int, int, int]]:
    weights = sorted(set(left) | set(right))
    return [(w, left.get(w, 0), right.get(w, 0)) for w in weights if left.get(w, 0) != right.get(w, 0)]


def section_text(index_html: str, section_title: str) -> str:
    match = re.search(rf"<h1>\s*{re.escape(section_title)}\s*</h1>", index_html, re.IGNORECASE)
    if match is None:
        raise ValueError(f"section not found: {section_title}")
    next_section = re.search(r"<h1>", index_html[match.end() :], re.IGNORECASE)
    if next_section is None:
        return index_html[match.end() :]
    return index_html[match.end() : match.end() + next_section.start()]


def audit_index(index_html: str, section_title: str, length: int) -> list[tuple[int, bool, str]]:
    section = section_text(index_html, section_title)
    linked: dict[int, str] = {}
    for href, n_s, k_s in re.findall(r'<a\s+href="([^"]+)">\s*\((\d+),(\d+)\)', section, re.IGNORECASE):
        if int(n_s) == length:
            linked[int(k_s)] = html.unescape(href)

    all_dims = set(linked)
    for n_s, k_s in re.findall(r"\((\d+),(\d+)\)", section):
        if int(n_s) == length:
            all_dims.add(int(k_s))
    return [(k, k in linked, linked.get(k, "")) for k in sorted(all_dims)]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extended", default=None, help="Extended BCH .wd path or URL")
    parser.add_argument("--primitive", default=None, help="Primitive BCH .wd path or URL to compare")
    parser.add_argument("--length", type=int, default=None, help="Extended length N; default max table weight")
    parser.add_argument("--out-csv", type=Path, default=None, help="Optional derived primitive CSV")
    parser.add_argument("--insecure-url", action="store_true", help="Disable TLS verification for public table URLs")
    parser.add_argument("--index-url", default=None, help="Audit linked dimensions on a weight-table index page")
    parser.add_argument("--section", default="Extended Binary Primitive BCH Codes")
    args = parser.parse_args()

    if args.index_url is not None:
        if args.length is None:
            raise SystemExit("--length is required with --index-url")
        rows = audit_index(read_text(args.index_url, args.insecure_url), args.section, args.length)
        print("BCH weight-distribution index audit")
        print(f"index_url={args.index_url}")
        print(f"section={args.section}")
        print(f"length={args.length}")
        print("dimension,linked,href")
        for dim, linked, href in rows:
            print(f"{dim},{int(linked)},{href}")
        return 0

    if args.extended is None:
        raise SystemExit("--extended is required unless --index-url is used")

    ext_name, ext_source, ext_spectrum = parse_wd_table(read_text(args.extended, args.insecure_url))
    length = args.length if args.length is not None else max(ext_spectrum)
    primitive = extended_to_primitive(ext_spectrum, length)

    print("BCH weight-distribution table tool")
    print(f"extended_table={ext_name}, source={ext_source}, length={length}")
    for key, value in summarize_spectrum(ext_spectrum).items():
        print(f"extended_{key}={value}")
    for key, value in summarize_spectrum(primitive).items():
        print(f"derived_primitive_{key}={value}")

    if args.primitive is not None:
        prim_name, prim_source, prim_spectrum = parse_wd_table(read_text(args.primitive, args.insecure_url))
        mismatches = compare_spectra(primitive, prim_spectrum)
        print(f"primitive_table={prim_name}, source={prim_source}")
        print(f"mismatches={len(mismatches)}")
        if mismatches:
            print("weight,derived,published")
            for weight, derived, published in mismatches[:20]:
                print(f"{weight},{derived},{published}")
        print(f"matches_published={len(mismatches) == 0}")

    if args.out_csv is not None:
        write_csv(args.out_csv, primitive)
        print(f"wrote_csv={args.out_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
