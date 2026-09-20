#!/usr/bin/env python3
"""Finish the shear insertion after the guarded first pass modified the header."""

from pathlib import Path


ROOT = Path("/tmp/riffle_parityshear_s19_v1")
header_path = ROOT / "libOTe/Tools/RiffleCode/RiffleExactPermFieldCheckpoint.h"
header = header_path.read_text()

for receipt in (
    "struct ParityShear12Schedule",
    "std::vector<ParityShear12Schedule> mParityShears",
    "mParityShears.resize(outerBlocks)",
):
    if header.count(receipt) != 1:
        raise RuntimeError(f"partial header does not contain exactly one {receipt!r}")

custom_marker = "const auto& shear = mParityShears[outerBase + offset + lane];"
custom_start = header.index(custom_marker)
pair_anchor = '''\t\t\t\tfor (u64 offset = 0; offset < tileOuterBlocks; offset += 2)
\t\t\t\t{
\t\t\t\t\tExtendedBch256x128Eq3::transposeBlock2(
\t\t\t\t\t\ttile + offset * outerLength,
'''
pair_replacement = '''\t\t\t\tfor (u64 offset = 0; offset < tileOuterBlocks; offset += 2)
\t\t\t\t{
\t\t\t\t\tfor (u64 lane = 0; lane < 2; ++lane)
\t\t\t\t\t{
\t\t\t\t\t\tblock* word = tile + (offset + lane) * outerLength;
\t\t\t\t\t\tconst auto& shear = mParityShears[outerBase + offset + lane];
\t\t\t\t\t\tconst block pivot = word[shear.target];
\t\t\t\t\t\tfor (const u8 source : shear.sources) word[source] ^= pivot;
\t\t\t\t\t}
\t\t\t\t\tExtendedBch256x128Eq3::transposeBlock2(
\t\t\t\t\t\ttile + offset * outerLength,
'''
prefix, custom = header[:custom_start], header[custom_start:]
if pair_anchor not in custom:
    raise RuntimeError("custom-inner pair loop not found")
header_path.write_text(prefix + custom.replace(pair_anchor, pair_replacement, 1))

bench_path = ROOT / "RiffleFieldCheckpoint_Bench.cpp"
bench = bench_path.read_text()
staged_anchor = '''\t\t\tinner.emitReverse(source.data(), N,
\t\t\t\t[&](u64 index, block value) { innerWord[index] = value; });
\t\t\tstd::cerr << "rm2sub: full inner materialization PASS\\n";
\t\t\tschedule.optimizedGatherAndOuterTranspose(
\t\t\t\tinnerWord.data(), stagedMessage.data(), tiledWorkspace.data());
\t\t\tstd::cerr << "rm2sub: staged complete PASS\\n";
'''
staged_replacement = '''\t\t\tinner.emitReverse(source.data(), N,
\t\t\t\t[&](u64 index, block value) { innerWord[index] = value; });
\t\t\tstd::cerr << "rm2sub: full inner materialization PASS\\n";
\t\t\tfor (u64 innerIndex = 0; innerIndex < N; ++innerIndex)
\t\t\t\touterWord[promoted.innerToOuter()[innerIndex]] = innerWord[innerIndex];
\t\t\tfor (u64 outer = 0; outer < OuterBlocks; outer += 2)
\t\t\t{
\t\t\t\tfor (u64 lane = 0; lane < 2; ++lane)
\t\t\t\t{
\t\t\t\t\tblock* word = outerWord.data() + (outer + lane) * OuterLength;
\t\t\t\t\tconst auto& shear = promoted.parityShears()[outer + lane];
\t\t\t\t\tconst block pivot = word[shear.target];
\t\t\t\t\tfor (const u8 sourceIndex : shear.sources) word[sourceIndex] ^= pivot;
\t\t\t\t}
\t\t\t\tExtendedBch256x128Eq3::transposeBlock2(
\t\t\t\t\touterWord.data() + outer * OuterLength,
\t\t\t\t\touterWord.data() + (outer + 1) * OuterLength,
\t\t\t\t\tstagedMessage.data() + outer * OuterDimension,
\t\t\t\t\tstagedMessage.data() + (outer + 1) * OuterDimension);
\t\t\t}
\t\t\tstd::cerr << "rm2sub: staged shear12 complete PASS\\n";
'''
if bench.count(staged_anchor) != 1:
    raise RuntimeError("staged benchmark anchor changed")
bench = bench.replace(staged_anchor, staged_replacement, 1)
for old, new in (
    (
        "construction=Riffle-BCHPerm-BitShuffle-RM2SubS19",
        "construction=Riffle-ParityShear12-BCHPerm-BitShuffle-RM2SubS19",
    ),
    (
        "correctness=dense-inner+staged-complete-PASS",
        "correctness=dense-inner+staged-shear12-complete-PASS",
    ),
):
    if bench.count(old) != 1:
        raise RuntimeError(f"benchmark label anchor {old!r} changed")
    bench = bench.replace(old, new, 1)
bench_path.write_text(bench)
print("finished ParityShear-12 insertion")
