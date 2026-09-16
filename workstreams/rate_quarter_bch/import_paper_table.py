"""Extract Table 7's k=71 column from the user-provided 2021 paper.

Optional dependency: pypdf. This reads but never copies or modifies the PDF.
The published k=63 column is cross-checked against the retained author table.
Only numerical data and a provenance receipt are saved in this workstream.
"""
import argparse
import hashlib
import json
from pathlib import Path
from pypdf import PdfReader
import reconstruct as r

SOURCE = "Fujiwara-Kusaka 2021, DOI 10.1587/transfun.2020EAP1119, Table 7, p.1326"


def extract(reader):
    r.require(len(reader.pages) == 8, "expected eight-page final paper")
    first = reader.pages[0].extract_text()
    r.require("10.1587/transfun.2020EAP1119" in first, "wrong paper DOI")
    text = reader.pages[5].extract_text()
    r.require("Table 7" in text and "Table 8" in text, "table anchors missing")
    table = text.split("Table 7", 1)[1].split("Table 8", 1)[0]
    r.require("k = 71" in table and "t = 29" in table, "wrong BCH table")
    half, small_half = {}, {}
    for line in table.splitlines():
        columns = line.split()
        if not columns or not columns[0].isdigit():
            continue
        r.require(all(v.isdigit() for v in columns), "noninteger table row")
        weight = int(columns[0])
        r.require(weight not in half, "duplicate weight")
        half[weight] = int(columns[-1])
        small_half[weight] = int(columns[-2]) if len(columns) >= 3 else 0
    r.require(set(half) == {0, *range(62, 129, 2)}, "unexpected table row support")
    parent = [0] * 257
    small = r.read_spectrum(r.HERE / "sources/EBCH256_63.wd")
    for w, v in half.items():
        parent[w] = parent[256-w] = v  # central weight 128 is assigned once
        r.require(small[w] == small_half[w], f"published k=63 mismatch at weight {w}")
    r.audit_spectrum(parent, 71, 60)
    selected = r.reconstruct(small, parent)
    r.audit_spectrum(selected, 64, 60)
    return parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdf", type=Path)
    parser.add_argument("--check", action="store_true", help="compare, do not write")
    args = parser.parse_args()
    parent = extract(PdfReader(args.pdf))
    text = ("# " + SOURCE + "; k=71, t=29\n"
            "# Header states A_w = A_(256-w); blanks are zero; w=128 is not doubled.\n"
            + "".join(f"{w} {v}\n" for w, v in enumerate(parent) if v))
    receipt = dict(schema=1, source=SOURCE, pdf_sha256=hashlib.sha256(args.pdf.read_bytes()).hexdigest(),
                   pdf_bytes=args.pdf.stat().st_size, pdf_page=6, printed_page=1326,
                   table=7, parent_dimension=71, parent_t=29,
                   parent_designed_distance=59, smaller_dimension=63, smaller_t=30,
                   smaller_designed_distance=61, lower_half_rows=35,
                   smaller_column_matches_public_table=True,
                   coefficients_sha256=hashlib.sha256(text.encode()).hexdigest(),
                   pdf_redistributed=False)
    outputs = {r.HERE / "sources/EBCH256_71.wd": text,
               r.HERE / "sources/PAPER_TABLE_PROVENANCE.json": json.dumps(receipt, indent=2)+"\n"}
    for path, content in outputs.items():
        if args.check:
            r.require(path.read_text(encoding="utf-8") == content, f"stale extraction: {path.name}")
        else:
            path.write_text(content, encoding="utf-8", newline="\n")
    print("TABLE 7 CHECK PASSED: 35 rows; both BCH columns cross-checked; exact audits passed")


if __name__ == "__main__":
    main()
