#!/usr/bin/env python3
"""Replace ParityShear-12 by ParityFanout-29x33 in the isolated benchmark."""

from pathlib import Path


ROOT = Path("/tmp/riffle_parityfanout29x33_s19_v1")


def replace_exact(text: str, old: str, new: str, count: int = 1) -> str:
    actual = text.count(old)
    if actual != count:
        raise RuntimeError(f"expected {count} copies of anchor, found {actual}")
    return text.replace(old, new)


header_path = ROOT / "libOTe/Tools/RiffleCode/RiffleExactPermFieldCheckpoint.h"
header = header_path.read_text()
header = replace_exact(
    header,
    '''\t\tstruct ParityShear12Schedule
\t\t{
\t\t\tu8 target;
\t\t\tstd::array<u8, 12> sources;
\t\t};''',
    '''\t\tstruct ParityFanout29x33Schedule
\t\t{
\t\t\tstd::array<u8, 33> targets;
\t\t\tstd::array<u8, 29> sources;
\t\t};''',
)
header = replace_exact(
    header,
    "const std::vector<ParityShear12Schedule>& parityShears() const noexcept",
    "const std::vector<ParityFanout29x33Schedule>& parityFanouts() const noexcept",
)
header = replace_exact(header, "return mParityShears;", "return mParityFanouts;")
header = replace_exact(
    header,
    "std::vector<ParityShear12Schedule> mParityShears;",
    "std::vector<ParityFanout29x33Schedule> mParityFanouts;",
)

setup_old = '''\t\t\tmParityShears.resize(outerBlocks);
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
\t\t\t}'''
setup_new = '''\t\t\tmParityFanouts.resize(outerBlocks);
\t\t\tstd::uniform_int_distribution<unsigned> coordinateDistribution(0, outerLength - 1);
\t\t\tfor (auto& fanout : mParityFanouts)
\t\t\t{
\t\t\t\tstd::array<u8, 62> sampled;
\t\t\t\tfor (unsigned pick = 0; pick < sampled.size(); ++pick)
\t\t\t\t{
\t\t\t\t\tu8 candidate;
\t\t\t\t\tdo candidate = static_cast<u8>(coordinateDistribution(random));
\t\t\t\t\twhile (std::find(sampled.begin(), sampled.begin() + pick, candidate) != sampled.begin() + pick);
\t\t\t\t\tsampled[pick] = candidate;
\t\t\t\t}
\t\t\t\tstd::copy(sampled.begin(), sampled.begin() + 33, fanout.targets.begin());
\t\t\t\tstd::copy(sampled.begin() + 33, sampled.end(), fanout.sources.begin());
\t\t\t}'''
header = replace_exact(header, setup_old, setup_new)

hot_old = '''\t\t\t\t\t\tconst auto& shear = mParityShears[outerBase + offset + lane];
\t\t\t\t\t\tconst block pivot = word[shear.target];
\t\t\t\t\t\tfor (const u8 source : shear.sources) word[source] ^= pivot;'''
hot_new = '''\t\t\t\t\t\tconst auto& fanout = mParityFanouts[outerBase + offset + lane];
\t\t\t\t\t\tblock pivot = word[fanout.targets[0]];
\t\t\t\t\t\tfor (unsigned index = 1; index < fanout.targets.size(); ++index)
\t\t\t\t\t\t\tpivot ^= word[fanout.targets[index]];
\t\t\t\t\t\tfor (const u8 source : fanout.sources) word[source] ^= pivot;'''
header = replace_exact(header, hot_old, hot_new, 2)
header_path.write_text(header)

bench_path = ROOT / "RiffleFieldCheckpoint_Bench.cpp"
bench = bench_path.read_text()
staged_old = '''\t\t\t\t\tconst auto& shear = promoted.parityShears()[outer + lane];
\t\t\t\t\tconst block pivot = word[shear.target];
\t\t\t\t\tfor (const u8 sourceIndex : shear.sources) word[sourceIndex] ^= pivot;'''
staged_new = '''\t\t\t\t\tconst auto& fanout = promoted.parityFanouts()[outer + lane];
\t\t\t\t\tblock pivot = word[fanout.targets[0]];
\t\t\t\t\tfor (unsigned index = 1; index < fanout.targets.size(); ++index)
\t\t\t\t\t\tpivot ^= word[fanout.targets[index]];
\t\t\t\t\tfor (const u8 sourceIndex : fanout.sources) word[sourceIndex] ^= pivot;'''
bench = replace_exact(bench, staged_old, staged_new)
bench = replace_exact(
    bench, "staged shear12 complete PASS", "staged fanout29x33 complete PASS"
)
bench = replace_exact(
    bench,
    "Riffle-ParityShear12-BCHPerm-BitShuffle-RM2SubS19",
    "Riffle-ParityFanout29x33-BCHPerm-BitShuffle-RM2SubS19",
)
bench = replace_exact(
    bench,
    "dense-inner+staged-shear12-complete-PASS",
    "dense-inner+staged-fanout29x33-complete-PASS",
)
bench_path.write_text(bench)
print("applied ParityFanout-29x33 benchmark substitutions")
