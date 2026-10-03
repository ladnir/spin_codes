# K16 packed-randomizer performance record

Measured on Peach, Ryzen 7950X, GCC 15.2.0, CPU15, normal pages, AVX-512/VBMI/GFNI.
The build uses the same Release flags as the frozen wider24 experiment,
including `-mtune=znver4`, `-fno-strict-aliasing`, and no LTO. K=65536,
N=131072, 128-bit XOR elements. Every number below includes inner, route,
randomizer, RS outer, and external packing/unpacking. Setup, table creation,
buffer allocation, and correctness checks occur before timing.

One benchmark process runs at a time, under the shared benchmark locks.
Each process performs five warmup calls, then 2001 timed calls. Stage clocks
are disabled. Four fresh seeds give two run medians per candidate per seed.
Mode order alternates between `0 4 7 20 20 7 4 0` and
`20 7 4 0 0 4 7 20`. No run or sample was removed from the reported medians.
Raw clocks have occasional slow runs; claims smaller than this variation
are not used to select a new construction.

## Holdout: every run median in microseconds

| Seed | Frozen, mode0 | Flat9, mode4 | Flat9 + identity symbol, mode7 | MDS sandwich, mode20 |
|---|---|---|---|---|
| 271 | 97.682, 97.692 | 93.855, 93.133 | 95.058, 98.093 | 96.971, 96.691 |
| 419 | 98.354, 97.943 | 93.525, 93.665 | 93.805, 93.836 | 96.430, 97.422 |
| 557 | 97.592, 97.282 | 93.985, 94.146 | 93.946, 93.204 | 96.590, 95.559 |
| 863 | 98.013, 101.119 | 94.005, 93.935 | 93.414, 93.274 | 97.653, 98.324 |
| **Median of all eight** | **97.8175** | **93.8950** | **93.8205** | **96.8310** |

The exact-map winner reduces time by 4.0100%, a 1.0418x throughput ratio.
Checksums for modes0 and4 match at every seed. Changed-family checksums need
not match the old family; each is checked against its own literal scalar
forward and transpose maps before timing.

## Exploratory screens

Screens used seeds41/113 and then53, with 1001 or1501 calls. They decide
which variants entered the holdout; their minima are not headline results.

| Mode | Change | Screen result / decision |
|---|---|---|
| 1 | One payload half at a time, 4 KiB workspace | about106 us; reject |
| 2 | Parity-first, full payload, 4 KiB workspace | about97–99 us; no clear gain |
| 3 | Parity-first + half payload, 2 KiB workspace | about106–107 us; reject |
| 4 | Nine byte products with flattened XOR circuit | about94 us; holdout winner |
| 5 | Same flat9, precomputed GFNI affine matrices | about95.7 us; larger metadata |
| 6 | Direct sixteen GFNI affine blocks | about100–102 us; reject |
| 7 | Flat9 with symbol0 map omitted | near mode4; no stable extra win |
| 8 | Parity-first + flat9 | 98.834 and107.100 us; reject |
| 9 | Parity-first + flat9 + identity symbol | 96.149 and96.380 us; reject |
| 10–15 | Write-prefetch 1,2,4,8,16,32 reverse epochs ahead | about101–109 us; reject |
| 20 | MDS sandwich, fixed GFNI middle | four-xtime version about98 us; three-xtime version in holdout |
| 21 | MDS sandwich, nibble shifts/masks | three-xtime version about108–111 us; reject |

Reduced workspace variants preserve the arithmetic and do not inherently
remove aggregate temporary traffic. No intermediate representation conversion
was removed: the frozen implementation already remains byte-packed between
the inner and outer.

## Measured identities

Final holdout binary SHA256:

```
568d8c237042e9722b6b0ce57734bf395d9d77178077b752a0c8259bbd423ff1
```

The local source archive `tmp/packet8-opt-combined-source.tar.gz` contains the
experiment sources used for that binary, before these reporting documents:

```
ad835ee6313570ae1867959c14ee1d57b4494c30c1fcc1fd1360d119e4c9d688
```

The baseline sources and their transitive research dependencies are the frozen
iteration5 files authenticated by `whole_wider24_v1.json`. The same prebuilt
SPIN package as that campaign supplies buffers and capability detection.
This is a kernel benchmark, not a measurement of a newly promoted public API.

Local ignored archive `tmp/packet8-opt-native-results.tar.gz` retains every
screen, all four holdout logs, and the final binary:

```
d1804df2c95eaa8179b21ea64e679fd321d3e98759afc81fd85ed18df16e1b41
```

The serial runner prints the binary digest into every log. Earlier binaries
and their source archives correspond to the exploratory stages only. The
remote campaign directory was `/tmp/spin-packet8opt-lH41eb`; reproduction
does not require that directory to survive.
