#!/usr/bin/env python3
"""Fuse ParityShear-12 into the isolated Peach s=19 benchmark."""

from pathlib import Path


ROOT = Path("/tmp/riffle_parityshear_s19_v1")


def replace_once(text: str, old: str, new: str) -> str:
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"expected one copy of insertion anchor, found {count}")
    return text.replace(old, new, 1)


header_path = ROOT / "libOTe/Tools/RiffleCode/RiffleExactPermFieldCheckpoint.h"
header = header_path.read_text()

workspace_anchor = '''\t\tstruct ExperimentalWorkspace : Workspace
\t\t{
\t\t\tstd::vector<u32> bucketOffsets;

\t\t\tExperimentalWorkspace()
\t\t\t\t: bucketOffsets(codeBlocks)
\t\t\t{
\t\t\t}
\t\t};
'''
workspace_replacement = workspace_anchor + '''
\t\tstruct ParityShear12Schedule
\t\t{
\t\t\tu8 target;
\t\t\tstd::array<u8, 12> sources;
\t\t};
'''
header = replace_once(header, workspace_anchor, workspace_replacement)

accessor_anchor = '''\t\tconst std::vector<u32>& innerToOuter() const noexcept { return mInnerToOuter; }
\t\tconst std::vector<u64>& coefficients() const noexcept { return mCoefficients; }
'''
accessor_replacement = accessor_anchor + '''\t\tconst std::vector<ParityShear12Schedule>& parityShears() const noexcept
\t\t{
\t\t\treturn mParityShears;
\t\t}
'''
header = replace_once(header, accessor_anchor, accessor_replacement)

member_anchor = '''\t\tstd::vector<u8> mBucketOffsetsPair5;
\t\tstd::vector<u64> mCoefficients;
'''
member_replacement = member_anchor + '''\t\tstd::vector<ParityShear12Schedule> mParityShears;
'''
header = replace_once(header, member_anchor, member_replacement)

schedule_anchor = '''\t\t\tfor (u64 pair = 0; pair < codeBlocks / 2; ++pair)
\t\t\t{
\t\t\t\tconst u64 slots = static_cast<u64>(mInnerToSlot[2 * pair]) |
\t\t\t\t\t(static_cast<u64>(mInnerToSlot[2 * pair + 1]) << 21);
\t\t\t\tauto* slotBytes = mInnerToSlotPair6.data() + 6 * pair;
\t\t\t\tfor (u64 byte = 0; byte < 6; ++byte)
\t\t\t\t\tslotBytes[byte] = static_cast<u8>(slots >> (8 * byte));

\t\t\t\tconst u64 offsets = static_cast<u64>(mBucketOffsets[2 * pair]) |
\t\t\t\t\t(static_cast<u64>(mBucketOffsets[2 * pair + 1]) << 19);
\t\t\t\tauto* offsetBytes = mBucketOffsetsPair5.data() + 5 * pair;
\t\t\t\tfor (u64 byte = 0; byte < 5; ++byte)
\t\t\t\t\toffsetBytes[byte] = static_cast<u8>(offsets >> (8 * byte));
\t\t\t}
'''
schedule_replacement = schedule_anchor + '''
\t\t\tmParityShears.resize(outerBlocks);
\t\t\tstd::uniform_int_distribution<unsigned> coordinateDistribution(0, outerLength - 1);
\t\t\tfor (auto& shear : mParityShears)
\t\t\t{
\t\t\t\tstd::array<u8, 13> sampled;
\t\t\t\tfor (unsigned pick = 0; pick < sampled.size(); ++pick)
\t\t\t\t{
\t\t\t\t\tu8 candidate;
\t\t\t\t\tdo candidate = static_cast<u8>(coordinateDistribution(random));
\t\t\t\t\twhile (std::find(sampled.begin(), sampled.begin() + pick, candidate) != sampled.begin() + pick);
\t\t\t\t\tsampled[pick] = candidate;
\t\t\t\t}
\t\t\t\tshear.target = sampled[0];
\t\t\t\tstd::copy(sampled.begin() + 1, sampled.end(), shear.sources.begin());
\t\t\t}
'''
header = replace_once(header, schedule_anchor, schedule_replacement)

avx4_anchor = '''\t\t\t\tfor (u64 offset = 0; offset < tileOuterBlocks; offset += 4)
\t\t\t\t{
\t\t\t\t\tExtendedBch256x128Eq3::transposeBlock4(
\t\t\t\t\t\ttile + offset * outerLength,
'''
avx4_replacement = '''\t\t\t\tfor (u64 offset = 0; offset < tileOuterBlocks; offset += 4)
\t\t\t\t{
\t\t\t\t\tfor (u64 lane = 0; lane < 4; ++lane)
\t\t\t\t\t{
\t\t\t\t\t\tblock* word = tile + (offset + lane) * outerLength;
\t\t\t\t\t\tconst auto& shear = mParityShears[outerBase + offset + lane];
\t\t\t\t\t\tconst block pivot = word[shear.target];
\t\t\t\t\t\tfor (const u8 source : shear.sources) word[source] ^= pivot;
\t\t\t\t\t}
\t\t\t\t\tExtendedBch256x128Eq3::transposeBlock4(
\t\t\t\t\t\ttile + offset * outerLength,
'''
# Only modify the custom-inner method, which is the final occurrence.
avx4_count = header.count(avx4_anchor)
if avx4_count < 1:
    raise RuntimeError("AVX4 loop anchor missing")
head, sep, tail = header.rpartition(avx4_anchor)
header = head + avx4_replacement + tail

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
# The custom-inner pair loop is immediately after the modified AVX4 loop.
custom_start = header.index(avx4_replacement)
before, custom = header[:custom_start], header[custom_start:]
if pair_anchor not in custom:
    raise RuntimeError("custom-inner pair loop anchor missing")
custom = custom.replace(pair_anchor, pair_replacement, 1)
header = before + custom
header_path.write_text(header)

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
bench = replace_once(bench, staged_anchor, staged_replacement)
bench = replace_once(
    bench,
    'construction=Riffle-BCHPerm-BitShuffle-RM2SubS19',
    'construction=Riffle-ParityShear12-BCHPerm-BitShuffle-RM2SubS19',
)
bench = replace_once(
    bench,
    'correctness=dense-inner+staged-complete-PASS',
    'correctness=dense-inner+staged-shear12-complete-PASS',
)
bench_path.write_text(bench)
print("fused ParityShear-12 into s19 benchmark")
