"""Reproduce source fingerprints and paper/Lean IMT basis comparison."""
import hashlib
import json
import pathlib
import re

root = pathlib.Path(__file__).resolve().parents[1]
paper = root.parent / "paper"
files = [
    "Schedule", "ConcreteOuter", "ConcreteOuterNative", "ConcreteRoutePermutation",
    "ConcreteSerialization", "Transvection", "ConcreteMapData", "ConcreteMaps",
    "ConcreteEncoder", "ConcreteEncoderOutput", "ConcreteRoutedMoment",
    "ConcreteNativeFamily", "ConcreteNativeTotal", "ConcreteNativeCodeword",
    "ConcreteNativeLinearCode", "ConcreteNativeLinearCodeDistance",
    "ConcreteOuterTailClosure", "ConcreteNativeDenseRates",
    "ConcreteNativeFixedLargeClosure", "ConcreteFixedInsertionNorm",
]
paths = [root / "SpinCodes" / "Structured" / (n + ".lean") for n in files]
paths += [paper / (n + ".tex") for n in
          ["structured_spin", "structured_proof", "structured_imt_appendix"]]
data = (root / "SpinCodes/Structured/ConcreteMapData.lean").read_text(encoding="utf-8")
appendix = (paper / "structured_imt_appendix.tex").read_text(encoding="utf-8")
table = re.findall(r"^\s*(\d+) & \\texttt\{([0-9a-f]+)\} & \\texttt\{([0-9a-f]+)\}", appendix, re.M)
assert [int(t[0]) for t in table] == list(range(19))
results = {}
for name, col in [("aRows", 1), ("cTransposeRows", 2)]:
    raw = re.search(r"def " + name + r" : List Nat :=\s*\[(.*?)\]", data, re.S).group(1)
    lean = [int(v, 16) for v in re.findall(r"0x([0-9a-f]+)", raw)]
    expected = [int(t[col], 16) for t in table]
    assert lean == expected, name
    results[name] = {"entries": len(lean), "paper_equal": True}
snapshot = {
    "scope": "Source semantic review snapshot; not a replacement for Lean kernel checks",
    "map_basis_comparison": results,
    "files": [{"path": str(p.relative_to(root.parent)).replace("\\", "/"),
               "sha256": hashlib.sha256(p.read_bytes()).hexdigest()} for p in paths],
}
out = root / "scripts/map_data/encoder_native_semantic_snapshot.json"
out.write_text(json.dumps(snapshot, indent=2) + "\n", encoding="utf-8")
print("PASS: all 38 paper/Lean basis words equal; source snapshot:", out)
