# Fixed permutation banks with lightweight randomization

The statistical experiments evaluate permutation families, without searching
SPIN codewords. A separate [full-encoder comparison](FLOWS.md) implements
fully precomputed and fresh-offset flows. Neither changes the public library,
libOTe, or any certificate.

The [optimization report](OPTIMIZATION.md) tracks the full fresh-code flow
against the target of 10% overhead over fully precomputed encoding.
The [final performance and tradeoff report](TARGET_RESULT.md) records the
screened heuristic candidate, confirmations, and reproduction command.

The main question is whether inexpensive transformations give useful diversity
when a small bank of permutations is sampled once and then reused. Regenerating
the bank for every observation would obscure precisely the correlations we
want to measure.

## Families

The domain has n=2^d positions. A bank contains B independent permutations
P_j, generated with Fisher--Yates and rejection sampling. The experiment uses
the existing seeded word generator for reproducible pseudorandom draws.
The bank remains fixed throughout each experiment.

Each instance selects j uniformly from the bank and independent d-bit words
a, b, c, e. Addition below is modulo n. The three families are:

```
XOR:         P_j[x XOR a] XOR c
XOR+add:     (P_j[(x XOR a) + b] XOR c) + e
Rotate+add:  rot_d(P_j[rot_d(x XOR a, ri) + b] XOR c, ro) + e
```

The rotation variant samples ri and ro independently from 0,...,d-1 for every
instance. These are rotations of d-bit coordinates, not 32-bit rotations.
The rotation counts are not fixed constants. Every family is bijective and
has an explicit inverse. Unused parameters in the simpler families have no
effect on the permutation.

The statistical screen uses ordinary unsigned tables for simplicity. The
separate address benchmark uses packed byte row tables and 16-bit region
tables. Neither changes the production encoder. Full-encoder measurements
are reported separately in [FLOWS.md](FLOWS.md).

## An exact limitation of XOR and XOR+add

Let v(x,y) be the position of the lowest bit at which distinct coordinates
x and y differ. Equivalently, v(x,y) is the number of trailing zeros of
x XOR y. Both XOR by a constant and addition by a constant modulo 2^d
preserve v(x,y). For addition, this follows because x-y is divisible by
2^v but not 2^(v+1), and the difference is unchanged modulo 2^d.

In particular, an input pair separated by the top bit remains such a pair
after the input wrapper. The output wrapper also preserves v. For a fixed
bank, the distribution of v(pi(0),pi(n/2)) under XOR+add is therefore exactly
the distribution obtained by enumerating

```
v(P_j[x], P_j[x XOR (n/2)])
```

over all bank entries j and all x. Offsets do not erase this bank-specific
distribution. A uniform permutation instead gives probability
2^(d-v-1)/(n-1) for valuation v. We calculate both distributions directly
and compare the prediction against sampled wrappers.

XOR-only has a stronger limitation: its entire XOR-difference histogram,
for each input difference, is unchanged by the wrappers. With one bank entry,
each nonzero input difference can produce at most n/2 distinct output XOR
differences, because inputs come in pairs. At least n/2-1 nonzero output
differences are impossible. This is a property of the conditional family,
not evidence that the bank permutation itself was sampled incorrectly.

Rotations break the lowest-differing-bit invariant. They do not imply full
independence or uniformly distributed permutations.

## Experiments

The broad screen uses domains of 256 and 2,048 positions, bank sizes
1, 4, 16, and 64, and eight independently seeded banks for each size/domain.
Bank sizes use nested prefixes for each seed, so comparisons across bank
sizes are paired, not independent. It samples 32,768 instances per family
and fixed bank. Input pairs are (0,2^i), (0,n-1), and (0,n/2-1).

The deeper 256-position screen uses the same eight banks, 262,144 instances
per family/bank, and all pairs (0,dx) for dx=1,...,255. This covers all input
XOR differences. The random input-XOR mask makes each fixed input pair's
distribution depend only on its XOR difference, although the measured output
statistics do not describe the full joint output distribution.

For each pair we record output XOR and modular-difference histograms. We
estimate collision excess over the uniform distribution on n-1 differences:

```
C = (n-1) * sum_y p_y^2 - 1.
```

The estimator uses counts c_y as
(n-1)*sum_y c_y(c_y-1)/(S(S-1))-1. It removes the ordinary histogram collision
bias for each fixed pair. Taking the maximum over input differences still
introduces selection bias. Matched uniform controls undergo the same maximum
selection, with the same sample count. C is not statistical distance or a
failure probability.

Three fixed supports each contain 32 positions: a contiguous interval, the
subspace of evenly spaced multiples of n/32, and the preimage of positions
0,...,31 under the first bank entry. The last is selected after seeing the
bank; the uniform control uses that same support. We measure concentration
in the four high-bit quarters and four low-bit residue classes.

For cross-instance tests, two independent parameter sets are forced to share
the same bank entry. We measure the intersection of their images of each
support. Under independent entry selection this event occurs with probability
1/B; these conditional measurements must not be reported as unconditional
rates. Uniform controls use two independent uniform permutations, sampled
only on the queried coordinates via rejection sampling without replacement.

## Results

The rotation variant is the most promising of the three. XOR-only retains
large conditional pair biases. XOR+add improves them but retains the exact
valuation obstruction. Random input and output rotations substantially reduce
the measured biases. A one-entry bank still has detectable residual bias;
16 and 64 entries deserve further study.

In the deeper 256-position screen, the mean over eight banks of the maximum
XOR-difference collision excess was:

| Bank entries | Uniform control | XOR | XOR+add | Rotate+add |
|---|---:|---:|---:|---:|
| 1 | 0.000246 | 2.617213 | 0.169916 | 0.002163 |
| 4 | 0.000257 | 0.627834 | 0.047766 | 0.000566 |
| 16 | 0.000261 | 0.157457 | 0.009730 | 0.000309 |
| 64 | 0.000270 | 0.040584 | 0.002521 | 0.000269 |

The corresponding maximum modular-difference excesses for 16 entries were
0.000282 (uniform), 0.025504 (XOR), 0.016015 (XOR+add), and 0.000374
(Rotate+add). For 64 entries they were 0.000265, 0.006790, 0.003924, and
0.000262, respectively. Thus rotation is not merely hiding an XOR bias in
the modular-difference statistic. The 64-entry results are close to the
control's measured sampling/maximum-selection noise level, not a proof of
uniform pair distributions.

The exact bank-specific valuation excess for the 16-entry banks averaged
0.002550. The sampled XOR+add excess averaged 0.002559. This agreement supports
the predicted mechanism rather than just an unexplained empirical bias.

The broad 2,048-position screen also strongly separated XOR-only from the
other families. At 16 entries, its maximum XOR-difference excess averaged
0.136974, compared with 0.004408 for XOR+add, 0.003637 for Rotate+add, and
0.003187 for uniform. That screen has fewer samples and only 13 input
differences; its resolution is insufficient to declare the latter families
equivalent.

Cross-instance correlations matter even when average overlaps look normal.
For two 32-position subspace supports in the 256-position domain, forced to
use the same bank entry, the probability of an overlap of at least 12 was:

| Bank entries | Uniform control | XOR | XOR+add | Rotate+add |
|---|---:|---:|---:|---:|
| 16 | 0.0000882 | 0.0011392 | 0.0003848 | 0.0000744 |
| 64 | 0.0000729 | 0.0011821 | 0.0003653 | 0.0000777 |

Larger banks do not remove the conditional same-entry issue: XOR-only remains
roughly 13--16 times above the control here, and XOR+add roughly 4--5 times.
The rotation variant shows no comparable inflation in this screen. These
are empirical tail frequencies, not security failure probabilities.

For support concentration, Rotate+add at 16 entries gave probabilities
0.003779, 0.003950, and 0.003854 that an interval, subspace, or bank-selected
preimage placed at least 16 of its 32 outputs in one quarter. The matched
controls were 0.003801, 0.003797, and 0.003813. The simpler families can show
both excess and deficient tails: a low tail count alone is not evidence of
a uniform family.

Release tests and ASan/UBSan tests passed for bijectivity and inverses in
256-, 2,048-, and 8,192-position domains. The exact identities above were
also checked exhaustively on 256-position domains. These tests validate the
implementation and identities, not the distance of a resulting SPIN code.

## Four-point and adaptive follow-up

The next screen tests 128 four-point shapes per bank. Shapes include binary
planes spanned by pairs of coordinate bits, arithmetic progressions, and
preimages of two fixed four-point shapes under up to 16 bank entries. Random
distinct four-point shapes fill the remaining slots. All shapes are fixed
before sampling the lightweight instance parameters.

For each shape we measure whether the four outputs XOR to zero, whether
y0+y3=y1+y2 modulo n, and whether all four occupy the same high-bit quarter.
We use eight banks of each size, both domain sizes, and 262,144 instances
per bank/family. Each uniform control uses the same shapes and maximum
selection. A fixed four-point shape has XOR-zero probability 1/(n-3) under
a uniform permutation: the first three distinct images determine one valid
fourth image among the n-3 remaining positions.

Mean per-bank maxima over the 128 shapes were:

| Domain; bank entries | Family | XOR zero | Modular relation | One quarter |
|---|---|---:|---:|---:|
| 256; 16 | Uniform | 0.004307 | 0.004275 | 0.015130 |
| 256; 16 | XOR+add | 0.005135 | 0.005465 | 0.016863 |
| 256; 16 | Rotate+add | 0.004303 | 0.004309 | 0.015176 |
| 256; 64 | Uniform | 0.004306 | 0.004259 | 0.015115 |
| 256; 64 | Rotate+add | 0.004313 | 0.004315 | 0.015178 |
| 2,048; 16 | Uniform | 0.000601 | 0.000596 | 0.016080 |
| 2,048; 16 | Rotate+add | 0.000603 | 0.000598 | 0.016087 |
| 2,048; 64 | Uniform | 0.000596 | 0.000608 | 0.016147 |
| 2,048; 64 | Rotate+add | 0.000607 | 0.000608 | 0.016143 |

The rotated variant remains close to the matched controls on these statistics.
These maxima should not be compared directly to the fixed-shape theoretical
mean without accounting for selection across shapes.

The adaptive search operates on the 256-position domain. For each fixed bank
and family, it samples 4,096 full permutations and counts, for every pair of
inputs, how often their outputs share a high-bit quarter. It selects a support
of size 32 that maximizes the sum of these pair counts. Sixteen initial
supports undergo greedy swaps, capped at 256 improving swaps per restart.
Initial supports include bank preimages and random supports. This is a local
search for average concentration, not exhaustive optimization of rare tails.

The selected support is then fixed. An independent stream supplies holdout
instances. Uniform controls undergo the same training and selection process,
so their training scores also benefit from overfitting. The initial holdout
used 262,144 instances per bank. A small 64-entry Rotate+add tail difference
prompted a separate confirmation run with 1,048,576 fresh holdout instances
per bank, without changing the selected supports.

Confirmation results, averaged over eight banks:

| Entries | Family | Holdout same-quarter pair rate | Holdout peak >=16 |
|---|---|---:|---:|
| 16 | Uniform | 0.2470607 | 0.0037795 |
| 16 | XOR+add | 0.2476040 | 0.0045064 |
| 16 | Rotate+add | 0.2470658 | 0.0038010 |
| 64 | Uniform | 0.2470527 | 0.0038180 |
| 64 | XOR+add | 0.2471378 | 0.0039061 |
| 64 | Rotate+add | 0.2470616 | 0.0038041 |

For comparison, a uniform permutation gives same-quarter pair probability
63/255 = 0.2470588. The rotated variant's apparent initial tail inflation
did not persist in the larger independent holdout. XOR+add retains a visible
effect, especially at 16 entries. The search supplies evidence about these
families, not a guarantee against all bank-adaptive supports.

## Packed routing cost

`Routing.h` implements the same Rotate+add family with packed inverse tables.
It samples parameters independently for 2,048 rows and 256 regions. Given an
inner position, it inverts the region permutation to obtain the row, then
inverts that row's permutation to obtain the original coordinate. Thus one
complete routing address needs two bank lookups, not one.

The row domain uses native 8-bit rotations. The 11-bit region domain uses
masked shifts. A 16-address batch separates lookup stages to expose independent
work. No dynamic dispatch or allocations occur inside the address loop.

Address-only results at K=2^18, Ryzen 7950X CPU 15, GCC 15.2, `-O3 -mavx2
-mtune=znver4`:

| Routing | Persistent bank | Fresh instance parameters | Fresh setup | Scalar addresses | Batched addresses |
|---|---:|---:|---:|---:|---:|
| Bank, 16 entries | 68 KiB | 27 KiB | 0.01638 ms | 2.01160 ms | 1.31696 ms |
| Bank, 64 entries | 272 KiB | 27 KiB | 0.01640 ms | 2.07796 ms | 1.36649 ms |
| Six-round Feistel | None | 264 KiB | 0.02032 ms | 3.12746 ms | 1.99321 ms |

Each address pass computes and writes all 524,288 inner-to-outer indices to
a preallocated 2 MiB array, visiting batches in reverse order. The table shows
warm passes. First-pass medians for the batched methods were 1.30962, 1.36271,
and 1.94127 ms, respectively. These timings exclude IMT masks, input elements,
the inner and BCH kernels, and full encoder memory traffic. They are not
on-demand encoding timings or predicted total encoding latencies.

The 16-entry bank reduces this address-generation work by about 34% relative
to batched six-round Feistel; the 64-entry bank reduces it by about 31%.
The 64-entry version is about 4% slower than the 16-entry version in this test.
One-time bank construction, including the generic temporary tables used by
this prototype, took approximately 0.23 ms and 0.80 ms. It is excluded from
fresh-instance setup because the bank is reused. Only inverse tables are
retained by this benchmark.

Measurements are medians of three process medians, reversing family order
on the middle pass. Each process uses 11 fresh instance seeds, discards two,
and measures one first pass followed by four warm passes. Compiler fences
prevent elimination of repeated passes. All benchmarks run serially under
the shared encoder locks. Untimed checks verify bijectivity and exact equality
between scalar and batched addresses. The test suite also compares every
address against the original generic bank implementation. Release and
ASan/UBSan routing checks passed.

## Recommended next experiments

Keep 16- and 64-entry Rotate+add as candidates. Both survived the additional
family screens substantially better than XOR+add. The fused on-demand encoder
is now measured in [FLOWS.md](FLOWS.md); its total speedup is not the same as
the address-only speedup. Further work should target that fused path. No
full-SPIN minimum-distance search is needed for this phase.

For reference, tightly packed forward row tables would occupy 4 KiB for
16 entries or 16 KiB for 64 entries. Forward tables on 2,048 positions, with
16-bit entries, would occupy 64 KiB or 256 KiB. Inverse tables would double
these figures if both directions are retained. The current screen uses
32-bit entries and stores both directions; these packed sizes are estimates,
not memory measurements of an optimized encoder.

## Reproduce

```sh
cmake -S spin -B out/bank -DCMAKE_BUILD_TYPE=Release -DSPIN_BUILD_EXPERIMENTS=ON
cmake --build out/bank --target spin_bank_test spin_bank_screen spin_bank_higher spin_bank_addresses -j2
ctest --test-dir out/bank -R '^spin_bank$' --output-on-failure
out/bank/spin_bank_screen 32768 8 > out/bank-screen.csv
out/bank/spin_bank_screen 262144 8 8 0 all-pairs > out/bank-screen-deep.csv
python3 spin/experiments/permutation_bank/summarize.py out/bank-screen-deep.csv
out/bank/spin_bank_higher 262144 8 > out/bank-higher.csv
out/bank/spin_bank_higher 1048576 8 adaptive 1 > out/bank-adaptive-confirm.csv
python3 spin/experiments/permutation_bank/summarize_higher.py out/bank-higher.csv
python3 spin/experiments/permutation_bank/summarize_higher.py out/bank-adaptive-confirm.csv
```

Arguments are sample count, number of bank seeds, optional domain bits,
optional bank size, and optional `all-pairs`. A zero filter selects all
configured domains/sizes. The tests check inverses, bijectivity, and the
exact XOR-histogram and XOR+add valuation identities. Raw experiment data
belongs under ignored `out/`, not in source control.

For Ryzen address timings, configure with `-DSPIN_TUNE=znver4`. The serial
script expects the executables under `COMPARISON_ROOT/build`:

```sh
bash spin/experiments/permutation_bank/compare_addresses.sh COMPARISON_ROOT 15
python3 spin/experiments/permutation_bank/summarize_addresses.py COMPARISON_ROOT/results-bank-addresses
```
