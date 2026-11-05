# Compact7: exact Toeplitz4 key compression

Status: default-off implementation; scalar/source and native correctness checks
pass. Matched full-encoder timing is effectively tied. This is
an implementation change to the existing byte-field Toeplitz4 map, not a new
outer distribution or a new probability certificate.

`HDGROUP_TOEPLITZ_COMPACT7=1` requires `HDGROUP_TOEPLITZ=1`,
`HDGROUP_BYTE_KEYS=1`, `HDGROUP_SUBFIELD=1`, and group width four. The driver
restricts dispatch and correctness checks to that width. The flag defaults to zero.
No other native source or proof helper was changed.

## Exact map

Fix seven bytes `t[0],...,t[6]`. The sampled matrix has entry
`T[i,j] = t[3+i-j]`. The established nine-product circuit prepares

```
A0 = t3             A1 = t3 ^ t2             A2 = t3 ^ t4
V0 = t5 ^ t3        V1 = V0 ^ t4 ^ t2        V2 = V0 ^ t6 ^ t4
W0 = t1 ^ t3        W1 = W0 ^ t0 ^ t2        W2 = W0 ^ t2 ^ t4
```

Compact7 stores `(A0,A1,A2,V0,V2,W0,W1)`. The hot kernel reconstructs
`K = A1 ^ A2`, `V1 = V0 ^ K`, and `W2 = W0 ^ K` after byte-key expansion.
Expansion replicates each of four independent bytes within its own 128-bit
record. It is binary-linear, so reconstruction commutes with expansion.

The kernel then applies the same nine products and combines them in the same
order. For every matrix, including singular matrices, and every input column,
its result is exactly `T*x`. In particular, the original rejection sampler is
unchanged: it makes the same draws and accepts the same invertible matrices.
RS transforms, output layout, routing, swap bits, and inner setup are unchanged.
The existing probability argument therefore needs no new distributional lemma.
The native checks below validate the compiled SIMD implementation separately.

## Static cost and layout

One SIMD call processes one column of four independent groups. Counts below
are source-level sites, not measured instructions or latency.

| Per column call | Existing | Compact7 |
|---|---:|---:|
| Packed four-byte key loads/expansions | 9 | 7 |
| Key byte shuffles | 9 | 7 |
| Variable GFNI products | 9 | 9 |
| Data XORs | 15 | 15 |
| Extra coefficient XORs | 0 | 3 |

There is no expanded-key array, runtime dispatch, callback, or additional data
pass. The explicit straight-line circuit retains the existing four-group layout.
The stored-coefficient count is separate from the nine-product count.

Each group has 64 columns. Hot keys decrease from 576 to 448 bytes per group.
At K20, 8192 groups use **3.5 MiB instead of 4.5 MiB**, saving exactly 1 MiB.
The 131072 column calls remove 262144 key expansions and add 393216 coefficient
XOR sites. Setup retains the original scalar matrices for its independent oracle;
the stated saving concerns the hot prepared-key table, not all setup storage.

Register allocation, shuffle-port pressure, and key-cache behavior determine the
actual result. The compiler may combine coefficient XORs or alter their schedule.
No speedup, and specifically no ten-percent speedup, is claimed.

## Checks completed

`theory/hd_toeplitz_compact7_tests.py` ran seven tests in 0.079 seconds:

- all 56 coefficient-bit directions reconstruct the old nine keys;
- all 56 by 32 coefficient/input basis pairs agree with a dense AES-field oracle;
- random dense columns and the zero singular matrix agree;
- byte replication preserves lane order and commutes with XOR;
- the packed seven-key layout covers its table exactly, without overlaps;
- the native source has seven expansion sites, nine GFNI sites, and 18 XOR sites;
- removing only Compact7 preprocessor branches reproduces both baseline source
  files, apart from blank separators.

The bilinear basis check proves the scalar identity for all matrix/input bytes.
The source-preservation check is not a binary-code comparison. It confirms that
flag-off source branches retain their previous content.

The flag-on native `--check` now additionally includes 1792 coefficient/input
basis cases, 4096 individual SIMD-bit cases, and 16 random columns against a
dense scalar oracle. The full-path checks reconstruct all nine prepared keys
from every stored seven-key record before comparing with the unchanged dense
encoder oracle. These native checks have now passed, including the K20 check.

Reproduce the completed scalar/source tests from the theory directory:

```powershell
C:/Python314/python.exe -m unittest hd_toeplitz_compact7_tests -v
```

## Review and validation plan

1. Review the two native diffs against `compact7_baseline/`. Its snapshots are
   byte-identical to the pre-edit header and driver.
2. Build the certified Toeplitz4/paired-route profile with the flag off and on,
   changing no other options. Confirm the off binary's hot-path code remains
   unchanged; avoid attributing unrelated compiler changes to Compact7.
3. Run `--check` and `--large-check` for the compact binary. Compare full-output
   checksums and setup-rejection counts with the baseline at the same seed.
4. Inspect the outer kernel for seven key expansions, unchanged nine products,
   no new spills, and no outlined helper in the column loop. Confirm the 448-byte
   group stride and the flag-on metadata.
5. Only after correctness passes, run isolated, serial outer/full measurements
   against the frozen baseline. Inner-only timing is an unchanged-path control.

Keep the flag off unless a later experiment establishes a reproducible full-path
gain. The current test establishes the key-memory saving, not a useful speedup.

## Native result, 2026-10-09

Root reviewed the source changes, built both variants with GCC15 and the same
AVX512/GFNI flags, and ran the candidate checks before timing. The flag-off
executable is byte-identical to the earlier `inline-toeplitz` executable.
The candidate passes 1792 coefficient/input basis cases, 4096 SIMD-bit cases,
30720 inner input/state directions, all pair masks and cyclic wrap cases,
small/natural-size full encodes, and the independent K20 oracle.
K20 retains checksum `7701295210980077712` and 2060 setup rejections.

Disassembly of each two-column coset loop has 14 key-shuffle instructions
instead of 18, the same 22 GFNI instructions including its final RS stage,
and no vector stack references. It adds no outlined mixer helper. The existing
RS-prefix call remains. These are static loop-body counts, not total timings.

Three serial, alternating-order comparisons each use 501 complete encodes:

| Variant | Three run medians, ms | Median of run medians, ms |
| --- | --- | ---: |
| Existing nine-key layout | 3.685541, 3.694597, 3.673448 | 3.685541 |
| Compact7 | 3.679881, 3.678037, 3.674790 | 3.678037 |

The difference is about 0.20%, and one paired round reverses its sign.
This does not establish a useful performance win or meet the 3.4-ms goal.
The absolute times are not a speedup against a previous day's 3.798-ms run.
Both variants ran on Peach CPU2 under both benchmark locks; no benchmarks
were run concurrently by this task. Other users' shared-host load is not
controlled.

Scratch directory: `/tmp/hd-compact7-20261009-uot619`.
The local `compact7_build.sh`, `compact7_bench.sh`, and `compact7_asm.py`
record the build, serial timing sequence, and read-only disassembly check.
The bundled WinSCP could not read the current key format, so source transfer
used the existing pinned-key plink connection instead. Frozen reference
directories were not changed.

Native hashes:

- control: `0b31e23d66948665135485da064868ad769d007bb5cc8a2c4c4570698f1638bc`;
- Compact7: `dd8067632562c60b639b279360e885fef08c70efb81d74fe471a7b1a42143588`;
- build/check log: `1b06d94104b36a0b958a1af617efef72a943542cf478c7f7bfddcfbdd49c7055`;
- matched timing log: `8eff72f04a95ad9628df72109242004950bf61b0f7eddab90737e9f0741ff8bc`.

## Source pins

Pre-edit snapshots:

- `compact7_baseline/HdGroupedOuter.h`:
  `f9d4aa8e934533ef2dd0c0c0ab39d5d1f3b64a3d14697e0dabcd1ab53856f962`
- `compact7_baseline/hd_grouped_forward.cpp`:
  `2223df3595805e54fe548ce9beac96915f758b639794387736467facb1fb837e`

Candidate sources:

- `spin/experiments/char2_field_spin/HdGroupedOuter.h`:
  `3a107366e33d5acf2f9b86933344e439c44fb8361358098fc0588b47cbb43e8f`
- `spin/experiments/char2_field_spin/hd_grouped_forward.cpp`:
  `628211f9f002d8f0cfb3d3226c4bae87b96411e1517bf1898c8ce39f05ccd9db`
- `theory/hd_toeplitz_compact7_tests.py`:
  `3a56f15411f5a8f6faa6072572d40a9a7d7fd82953d470f645fb02d31bb823f5`
