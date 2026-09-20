#!/usr/bin/env python3
"""Interleave the two fanout maps consumed by one paired BCH transpose."""

from pathlib import Path


path = Path(
    "/tmp/riffle_parityfanout31x33_s19_v1/"
    "libOTe/Tools/RiffleCode/RiffleExactPermFieldCheckpoint.h"
)
text = path.read_text()
old = '''\t\t\t\tfor (u64 offset = 0; offset < tileOuterBlocks; offset += 2)
\t\t\t\t{
\t\t\t\t\tfor (u64 lane = 0; lane < 2; ++lane)
\t\t\t\t\t{
\t\t\t\t\t\tblock* word = tile + (offset + lane) * outerLength;
\t\t\t\t\t\tconst auto& fanout = mParityFanouts[outerBase + offset + lane];
\t\t\t\t\t\tblock pivot = word[fanout.targets[0]];
\t\t\t\t\t\tfor (unsigned index = 1; index < fanout.targets.size(); ++index)
\t\t\t\t\t\t\tpivot ^= word[fanout.targets[index]];
\t\t\t\t\t\tfor (const u8 source : fanout.sources) word[source] ^= pivot;
\t\t\t\t\t}
\t\t\t\t\tExtendedBch256x128Eq3::transposeBlock2(
\t\t\t\t\t\ttile + offset * outerLength,
\t\t\t\t\t\ttile + (offset + 1) * outerLength,
\t\t\t\t\t\toutput + (outerBase + offset) * outerDimension,
\t\t\t\t\t\toutput + (outerBase + offset + 1) * outerDimension);
\t\t\t\t}'''
new = '''\t\t\t\tfor (u64 offset = 0; offset < tileOuterBlocks; offset += 2)
\t\t\t\t{
\t\t\t\t\tblock* word0 = tile + offset * outerLength;
\t\t\t\t\tblock* word1 = word0 + outerLength;
\t\t\t\t\tconst auto& fanout0 = mParityFanouts[outerBase + offset];
\t\t\t\t\tconst auto& fanout1 = mParityFanouts[outerBase + offset + 1];
\t\t\t\t\tblock pivot0 = word0[fanout0.targets[0]];
\t\t\t\t\tblock pivot1 = word1[fanout1.targets[0]];
\t\t\t\t\tfor (unsigned index = 1; index < fanout0.targets.size(); ++index)
\t\t\t\t\t{
\t\t\t\t\t\tpivot0 ^= word0[fanout0.targets[index]];
\t\t\t\t\t\tpivot1 ^= word1[fanout1.targets[index]];
\t\t\t\t\t}
\t\t\t\t\tfor (unsigned index = 0; index < fanout0.sources.size(); ++index)
\t\t\t\t\t{
\t\t\t\t\t\tword0[fanout0.sources[index]] ^= pivot0;
\t\t\t\t\t\tword1[fanout1.sources[index]] ^= pivot1;
\t\t\t\t\t}
\t\t\t\t\tExtendedBch256x128Eq3::transposeBlock2(
\t\t\t\t\t\tword0,
\t\t\t\t\t\tword1,
\t\t\t\t\t\toutput + (outerBase + offset) * outerDimension,
\t\t\t\t\t\toutput + (outerBase + offset + 1) * outerDimension);
\t\t\t\t}'''
if text.count(old) != 1:
    raise RuntimeError(f"expected one pair loop, found {text.count(old)}")
path.write_text(text.replace(old, new, 1))
print("interleaved paired fanout maps")
