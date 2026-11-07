"""Generate untrusted packed-map spectrum witnesses for Lean kernel replay.

The generated Lean modules check the production row masks and every claimed
weight. This script's arithmetic is not part of the trusted certificate.
Use --pilot for the first and high-offset blocks; --all emits all 4096 blocks.
"""

import argparse
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
MAP_SOURCE = ROOT / "SpinCodes/UniformMixer/Map22Packed.lean"
OUT = ROOT / "SpinCodes/UniformMixer/Map22SpectrumData"
BLOCK_SIZE = 1024
BLOCK_COUNT = (1 << 22) // BLOCK_SIZE


def rows_from_source():
    source = MAP_SOURCE.read_text(encoding="utf-8")
    match = re.search(r"def productionPackedRows : List Nat := \[(.*?)\]", source, re.S)
    if match is None:
        raise ValueError("productionPackedRows was not found")
    rows = [int(word, 16) for word in re.findall(r"0x[0-9a-f]+", match.group(1))]
    if len(rows) != 22 or any(row >= 1 << 64 for row in rows):
        raise ValueError("expected 22 rows of at most 64 bits")
    return rows


def encode(rows, bits):
    result = 0
    while bits:
        low = bits & -bits
        result ^= rows[low.bit_length() - 1]
        bits ^= low
    return result


def lean_list(values, formatter=str):
    return "[" + ", ".join(map(formatter, values)) + "]"


def basis_module(rows):
    low = [encode(rows[:10], q) for q in range(BLOCK_SIZE)]
    return "\n".join([
        "import SpinCodes.UniformMixer.Map22PackedBridge",
        "namespace Spin.UniformMixer.Map22SpectrumData",
        "set_option maxRecDepth 1000000",
        "set_option maxHeartbeats 0",
        f"def lowImages : List Nat := {lean_list(low, hex)}",
        "theorem lowImages_checked :",
        "    (List.range 1024).map (Spin.Structured.PackedMap.eval",
        "      (productionPackedRows.take 10)) = lowImages := by",
        "  decide +kernel",
        "end Spin.UniformMixer.Map22SpectrumData",
        "",
    ])


def block_module(rows, block):
    high = encode(rows[10:], block)
    weights = [(low ^ high).bit_count()
               for low in (encode(rows[:10], q) for q in range(BLOCK_SIZE))]
    hist = [0] * 129
    for weight in weights:
        hist[weight] += 1
    stem = f"block{block}"
    return "\n".join([
        "import SpinCodes.UniformMixer.Map22SpectrumData.Basis",
        "namespace Spin.UniformMixer.Map22SpectrumData",
        "set_option maxRecDepth 1000000",
        "set_option maxHeartbeats 0",
        f"def {stem}High : Nat := {hex(high)}",
        f"theorem {stem}_high_checked :",
        f"    Spin.Structured.PackedMap.eval (productionPackedRows.drop 10) {block} =",
        f"      {stem}High := by decide +kernel",
        f"def {stem}Weights : List Nat := {lean_list(weights)}",
        f"theorem {stem}_fast_checked :",
        f"    lowImages.map (fun q => Spin.Structured.MapSpectrum.fastWeight 128",
        f"      (q ^^^ {stem}High)) = {stem}Weights := by decide +kernel",
        f"theorem {stem}_weights_checked :",
        f"    lowImages.map (fun q => Spin.Structured.PackedMap.weight 128",
        f"      (q ^^^ {stem}High)) = {stem}Weights := by",
        f"  simpa only [Spin.Structured.MapSpectrum.fastWeight_eq] using {stem}_fast_checked",
        f"def {stem}Histogram : List Nat := {lean_list(hist)}",
        f"theorem {stem}_hist_checked :",
        f"    {stem}Weights.foldl (fun h w => Spin.Structured.MapSpectrum.bump w h)",
        f"      (List.replicate 129 0) = {stem}Histogram := by decide +kernel",
        f"theorem {stem}_checked :",
        f"    Spin.Structured.MapSpectrum.block productionPackedRows ({block} * 1024) 1024 =",
        f"      {stem}Histogram := by",
        f"  calc",
        f"    _ = {stem}Weights.foldl",
        f"        (fun h w => Spin.Structured.MapSpectrum.bump w h)",
        f"        (List.replicate 129 0) := by",
        f"          apply Spin.Structured.MapSpectrum.block_split_weights",
        f"            productionPackedRows 10 {block} (by decide) lowImages {stem}High",
        f"            {stem}Weights lowImages_checked {stem}_high_checked",
        f"            {stem}_weights_checked",
        f"    _ = {stem}Histogram := {stem}_hist_checked",
        "end Spin.UniformMixer.Map22SpectrumData",
        "",
    ])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument("--pilot", action="store_true")
    choice.add_argument("--all", action="store_true")
    args = parser.parse_args()
    rows = rows_from_source()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "Basis.lean").write_text(basis_module(rows), encoding="utf-8")
    blocks = [0, BLOCK_COUNT // 2] if args.pilot else range(BLOCK_COUNT)
    for block in blocks:
        (OUT / f"Block{block}.lean").write_text(block_module(rows, block), encoding="utf-8")
    print(f"Generated {len(blocks)} untrusted block witness modules")


if __name__ == "__main__":
    main()
