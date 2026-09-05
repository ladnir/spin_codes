#!/usr/bin/env python3
"""Move the isolated fanout benchmark from 29x33 to 31x33."""

from pathlib import Path


ROOT = Path("/tmp/riffle_parityfanout31x33_s19_v1")


def replace_count(text: str, old: str, new: str, count: int) -> str:
    actual = text.count(old)
    if actual != count:
        raise RuntimeError(f"expected {count} copies of {old!r}, found {actual}")
    return text.replace(old, new)


header_path = ROOT / "libOTe/Tools/RiffleCode/RiffleExactPermFieldCheckpoint.h"
header = header_path.read_text()
header = replace_count(header, "ParityFanout29x33Schedule", "ParityFanout31x33Schedule", 3)
header = replace_count(header, "std::array<u8, 29> sources", "std::array<u8, 31> sources", 1)
header = replace_count(header, "std::array<u8, 62> sampled", "std::array<u8, 64> sampled", 1)
header_path.write_text(header)

bench_path = ROOT / "RiffleFieldCheckpoint_Bench.cpp"
bench = bench_path.read_text()
bench = replace_count(bench, "fanout29x33", "fanout31x33", 2)
bench = replace_count(bench, "ParityFanout29x33", "ParityFanout31x33", 1)
bench_path.write_text(bench)
print("applied ParityFanout-31x33 benchmark substitutions")
