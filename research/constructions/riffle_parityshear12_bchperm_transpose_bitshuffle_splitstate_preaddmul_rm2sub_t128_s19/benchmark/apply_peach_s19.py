#!/usr/bin/env python3
"""Apply the isolated Peach s=19 benchmark substitutions with count checks."""

from pathlib import Path


ROOT = Path("/tmp/riffle_parityshear_s19_v1")


def replace_exact(text: str, old: str, new: str, count: int) -> str:
    actual = text.count(old)
    if actual != count:
        raise RuntimeError(f"expected {count} copies of {old!r}, found {actual}")
    return text.replace(old, new)


header_path = ROOT / "libOTe/Tools/RiffleCode/RiffleExactPermFieldCheckpoint.h"
header = header_path.read_text()
header = replace_exact(
    header, '#include "RiffleRm2SubS16.h"', '#include "RiffleRm2SubS19.h"', 1
)
header_path.write_text(header)

bench_path = ROOT / "RiffleFieldCheckpoint_Bench.cpp"
bench = bench_path.read_text()
bench = replace_exact(bench, "verifyRm2SubS16Transpose", "verifyRm2SubS19Transpose", 2)
bench = replace_exact(bench, "RiffleRm2SubS16Transpose", "RiffleRm2SubS19Transpose", 3)
bench = replace_exact(bench, "rm2SubS16Only", "rm2SubS19Only", 2)
bench = replace_exact(bench, "rm2sub-s16", "rm2sub-s19", 2)
bench = replace_exact(bench, "RM2SubS16", "RM2SubS19", 1)
bench = replace_exact(bench, "optimized RM2Sub s16", "optimized RM2Sub s19", 2)
bench = replace_exact(bench, '<< " t=128 s=16 epochs="', '<< " t=128 s=19 epochs="', 1)

start = bench.index("\tvoid verifyRm2SubS19Transpose()")
end = bench.index("\n#ifndef _WIN32", start)
prefix, verify, suffix = bench[:start], bench[start:end], bench[end:]
verify = replace_exact(verify, "std::array<block, 16>", "std::array<block, 19>", 3)
verify = replace_exact(
    verify, "u16 column = detail::Rm2SubColumns[point];",
    "u32 column = detail::Rm2Sub19Columns[point];", 2
)
verify = replace_exact(
    verify, "column &= static_cast<u16>(column - 1);", "column &= column - 1;", 2
)
verify = replace_exact(verify, "row < 16", "row < 19", 1)
verify = replace_exact(verify, "u16 mask = masks[row];", "u32 mask = masks[row];", 1)
verify = replace_exact(
    verify, "mask &= static_cast<u16>(mask - 1);", "mask &= mask - 1;", 1
)
bench_path.write_text(prefix + verify + suffix)
print("applied s19 benchmark substitutions")
