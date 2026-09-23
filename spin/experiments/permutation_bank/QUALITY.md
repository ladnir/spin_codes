# Screening the fast composed family

The fast cyclic-row-shift candidate meets the median timing target, but its
bank-adaptive overlap tails differ sharply from uniform permutations.
The stronger row-map prototype improves the tested tails. Its final optimized
version meets the median timing target; see [TARGET_RESULT.md](TARGET_RESULT.md).
Neither family is a certified replacement. Timing sections below retain the
earlier measurements, including failed optimizations and short-run variability.

## What remains exact for one active row

Use the family defined in `COMPOSED.md`. Fix any valid bank, any outer row,
and any support of size w in that row. The support may depend on the bank,
but not on the fresh setup. Assume the specified independent uniform setup
sampling: a shared region permutation sigma, independent cyclic position
shifts b_g, and row maps independent of sigma and the position shifts.

The resulting inner support has exactly this distribution:

1. Select a uniform w-subset of the 256 actual regions.
2. Independently select a uniform position among the 2,048 positions of
   each selected region.

To see this, condition on the row maps. They identify a fixed set of w stored
regions. The inverse of uniform sigma sends that set to a uniform w-subset.
Condition also on sigma and all selected region tables. Each active region
contains a single fixed table position for the active row. Subtracting its
independent uniform b_g makes the actual position uniform. Those shifts are
independent across regions, proving the stated distribution.

This matches the one-active-row routing distribution of the original
independently shuffled structured route. It does not require the bank itself
to remain random after it is published. This is a statement about a fixed
word's marginal distribution, not independence between different words.
It also does not transfer any multi-row bound or certify the whole code.
A one-active-row first-moment calculation that uses only this marginal and
the unchanged independent inner can therefore be reused.

## Why the cyclic-row-shift mean was misleading

Fix two rows and 32 stored regions. After reading the bank, choose the 32
columns in each row that reach those regions before fresh row shifts.
Those supports are then fixed for the experiment.

`composed_screen.cpp` enumerates all 256^2 pairs of row shifts. A shared
sigma does not change the size of the intersection, so this enumeration
also accounts for every region relabeling. Across bank seeds 913 through
920, the exact distribution for these supports gave:

| Statistic | Cyclic row shifts | Independent uniform 32-subsets |
|---|---:|---:|
| Mean intersection | 4 | 4 |
| Probability of intersection at least 16 | 2^-16 | approximately 1.233e-8 |
| Probability of intersection 32 | 2^-16 | approximately 1.717e-41 |

The mean is exactly 4 for any fixed pair of 32-bit row supports: each
region is reached with probability 32/256 in each independently shifted row.
The mean therefore cannot establish that the tails are suitable.

The experiment also checks interval and regularly spaced column supports.
Those did not exhibit intersection at least 16 in these enumerations.
That does not repair the bank-adaptive counterexample. These arbitrary
supports are not shown to be BCH codewords, so this is not a SPIN attack.

## Stronger row-map prototype

`row-affine1` replaces each cyclic column shift with an independently sampled
map on eight-bit column indices:

```
F_r(c) = ((a_r * c + b_r) mod 256) XOR e_r,
```

where a_r is uniform over the 128 odd residues, and b_r,e_r are independent
uniform bytes. Odd multiplication makes every map bijective. The region
bank, shared region permutation, region position shifts, and IMT masks are
otherwise unchanged. The prototype keeps a 32-bit parameter per row.

`affine_tail.cpp` computes exact complete-overlap probabilities for selected
bank-adaptive support pairs. Two algebraic identities reduce the enumeration:

```
(a,b,e) ~ (a,b+128,e XOR 128)
(a,b,e) ~ (256-a,255-b,e XOR 255).
```

All arithmetic in b is modulo 256. Restricting odd a and b to values below
128 leaves 2^21 equally weighted representatives. Each represents four
original keys. The program checks these identities exhaustively on all
column indices, enumerates support images, sorts them, and counts equal pairs.
It uses about 128 MiB for the two support-image arrays and writes only the
small result record.

For the original supports in rows 0 and 1, both bank seeds 913 and 914 gave
complete-overlap probability exactly 2^-42. This is not a universal bound.
A more adaptive construction finds a 32-element common orbit under the two
rows' conjugates of addition by 128. Its supports have extra symmetries.
The exact results are:

| Bank seed | Second row | Support selection | Equal representative pairs | Probability |
|---|---:|---|---:|---:|
| 913 | 1 | First 32 stored regions | 1 | 2^-42 |
| 914 | 1 | First 32 stored regions | 1 | 2^-42 |
| 913 | 35 | Common 32-element orbit | 64 | 2^-36 |
| 914 | 10 | Common 32-element orbit | 64 | 2^-36 |

The denominator is 2^42 in every row. These are exact results for the
specified support pair and bank, not security margins or worst-case bounds.
Odd multiplication, addition, and XOR all preserve the position of the
lowest differing bit. This remaining structure motivates testing a byte
rotation or another inexpensive non-preserving operation next.
In particular, every map in this prototype satisfies
`F(c+128) = F(c)+128` modulo 256. Adding more choices of odd multiplier
does not remove this identity.

## Adding a fresh byte rotation

`row-rotate1` samples an independent rotation d_r from 0 through 7 for each
row and uses

```
F_r(c) = ((a_r * rotr8(c,d_r) + b_r) mod 256) XOR e_r.
```

The remaining parameters and routing stages are unchanged. Rotation breaks
the common addition-by-128 identity across the family. It does not make
the family uniform over all byte permutations.

For each of the two common-orbit support pairs above, `affine_tail.cpp`
enumerates all 64 pairs of rotations as well as the affine representatives.
Both calculations give 64 matching representative pairs out of 2^48:
complete overlap has probability exactly 2^-42, rather than 2^-36.
Only one rotation pair contributes in each calculation. The enumeration
uses about 576 MiB of temporary memory and retains a small result record.
These are the same two selected supports, not a worst-case bound for
supports adapted to the new family. In particular, 42 here is not a code
security margin.

Nine matched serial runs on Ryzen gave the following process medians:

| Flow | Fresh setup | Encode | Total | Ratio to fixed |
|---|---:|---:|---:|---:|
| Original-family fixed code | amortized | 1.492354 ms | 1.492354 ms | 1.00000 |
| Odd multiply, add, XOR | 0.006462 ms | 1.630292 ms | 1.636764 ms | 1.09677 |
| Rotation, odd multiply, add, XOR | 0.006482 ms | 1.657353 ms | 1.663233 ms | 1.11450 |

The rotation candidate remains above the 10% target. Ratios use total time,
including fresh masks and route parameters, but exclude reusable bank and
workspace preparation. Every timed fresh mode first passes same-map checks
against the materialized-route encoder for two seeds.

Preparing only 128 addresses per inner step, instead of 2,048 per region,
raised overhead to 16.5% in nine matched runs. Forcing that preparation
inline still gave 15.6% in five runs. Prefetching output stores 32 positions
ahead gave 19.9% in five runs. These variants remain opt-in experiments;
none changes the permutation family or the installed encoder defaults.

Removing table lookahead gave 11.6% overhead; disabling loop unrolling gave
13.3%. Unrolling two or eight vectors gave 11.1% and 11.3%, respectively.
These five-run screens do not establish an improvement over the four-vector
baseline.

`SPIN_COMPOSED_ROW_INDEX=1` stores the row index in unused high bits of a
redundant bank array. This removes index reconstruction before each SIMD
lookup. Nine matched runs gave 1.659020 ms fresh versus 1.494263 ms fixed,
or 11.03% overhead. The prototype retains its original array as a scalar
reference, so C=1 uses 4 MiB of reusable bank storage instead of 2 MiB.
Only the indexed array feeds the shifted SIMD route; per-code storage does
not grow.

`SPIN_COMPOSED_ROTATE_KEYS=1` additionally repacks each row key so byte
rotation and multiplication need one fewer mask operation. Nine runs gave
1.655262 ms fresh versus 1.487449 ms fixed, or 11.28% overhead. The paired
median ratio was 1.11089. This does not establish a further improvement.
The indexed layout and repacked keys preserve the sampled maps and seeds.
Both pass focused AVX2 ASan/UBSan routing tests. The combined native
AVX-512 build also passes, including scalar/SIMD equivalence and 128-address
versus full-region preparation. Full-encoder checksums match across the
original rotated implementation, indexed layout, and repacked keys.
All 24 registered portable
bank tests pass, including the rotated full encoder. Those measurements
preceded the successful byte-offset implementation in
[TARGET_RESULT.md](TARGET_RESULT.md).

## Region-local screens

For one actual region, uniform sigma selects one of 256 stored regions.
The region table index and cyclic shift then determine the row-label map.
The row-column maps do not change those row labels.

`composed_screen.cpp` enumerates this conditional marginal exactly, using
eight banks, 12 position gaps, and four four-point shapes. For C=1, the mean
XOR collision excess over gaps/banks was 0.00424, with maximum 0.00834.
For C=4, the corresponding values were 0.00106 and 0.00202. These are genuine
finite-family biases, not Monte Carlo estimation noise.

Four-point XOR-zero and common-quarter probabilities were close to their
uniform reference values in this screen. That does not imply joint uniformity.
When two fresh codes use the same stored region/table, cyclic position shifts
preserve the circular ordering of all row labels. For C=1, every stored
region's row ordering is reused up to a shift somewhere in every fresh code.
Public region relabeling moves that ordering; it does not remove it.

## Performance and decision

The richer map uses packed-coordinate arithmetic, 16-bit SIMD multiplication,
unrolled lookups, and the exact vectorized IMT mask generator. Nine serial
matched runs on the same Ryzen setup as `COMPOSED.md` gave:

| Flow | Fresh setup | Encode | Total | Ratio to fixed |
|---|---:|---:|---:|---:|
| Original-family fixed code | amortized | 1.489719 ms | 1.489719 ms | 1.00000 |
| Cyclic row shifts | 0.006452 ms | 1.617878 ms | 1.624400 ms | 1.09041 |
| Odd multiply, add, XOR | 0.006482 ms | 1.634650 ms | 1.641192 ms | 1.10168 |

The stronger map narrowly misses the 10% target. The faster candidate meets
it but has the larger adaptive tail above. These were interim measurements;
the final byte-offset implementation and independent confirmations are in
[TARGET_RESULT.md](TARGET_RESULT.md). No public default has changed.

All 23 portable bank tests pass, including same-map encoding checks for the
new prototype. Focused AVX2 and native AVX-512 ASan/UBSan routing checks also
pass. The tail calculations above are separate finite enumerations, not
claims implied by those correctness tests.

## Reproduce

```sh
g++ -O3 -std=c++20 spin/experiments/permutation_bank/composed_screen.cpp \
  -o out/composed-screen
out/composed-screen > out/composed-screen.csv
g++ -O3 -std=c++20 spin/experiments/permutation_bank/affine_tail.cpp \
  -o out/affine-tail
out/affine-tail 913
out/affine-tail 913 cycle
out/affine-tail 914
out/affine-tail 914 cycle
out/affine-tail 913 cycle rotate
out/affine-tail 914 cycle rotate
```

The performance build uses the flags in `COMPOSED.md`, plus `-mavx512bw`
and `-DSPIN_BANK_AFFINE16=1`. Use comparison modes
`fixed,row-shift1-epoch,row-affine1` and nine passes. The raw result directory
is `results-bank-optimize-affine16`; data remain outside source control.

For the rotated family use the same native flags and comparison modes
`fixed,row-affine1,row-rotate1`. Its nine-run result directory is
`results-bank-optimize-rotate`. The scheduling experiments use, respectively,
`SPIN_COMPOSED_EPOCH_ROUTING=1`, that flag plus
`SPIN_COMPOSED_INLINE_ROUTE=1`, or `SPIN_COMPOSED_OUTPUT_PREFETCH=32`.
Their result directories end in `epochroute`, `inline`, and `prefetch32`.
The indexed layout uses the same native flags plus
`-DSPIN_COMPOSED_ROW_INDEX=1`; its results end in `index`. Adding
`-DSPIN_COMPOSED_ROTATE_KEYS=1` gives the `rotkeys` experiment. The loop
unroll screens use `SPIN_COMPOSED_UNROLL_VECTORS=2` or `8` with
`SPIN_COMPOSED_SHIFT_UNROLL=1` and produce `unroll2` and `unroll8`.
