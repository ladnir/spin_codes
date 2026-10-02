# Two Routes Toward 6 ms and 10 Percent Distance

2026-09-30. This bounded exploration tested a change to BCH evaluation and
a new local packet mixer. Neither currently meets both targets. The fused
BCH evaluation gives no speedup. The packet mixer improves the support
bounds, but its optimized prototype takes 7.315 ms and lacks a complete
distance certificate. Existing code distributions and certificates remain
unchanged; these are research options, not production defaults.

## Measured Results

All headline times concern the precomputed in-place transposed encoder at
K=2^20, rate 1/2, with 128-bit elements on Peach's Ryzen 7950X. Setup and
allocation are excluded. Each headline is the median of four process
medians, using two seeds and balanced serial order. Each process uses
three warmups and 101 timed calls. No benchmarks ran concurrently.

| Implementation | Time | Proof status |
|---|---:|---|
| Existing shared GF16, two updates | 5.823 ms | Existing complete >7% / >44.61-bit certificate |
| Existing independent-row GF16, four updates | 7.491 ms | Existing complete >10% / >48.65-bit certificate |
| Same independent-row code, fused BCH preparation | 7.548 ms | Same code and proof; no speed improvement |
| Same independent-row code, direct cold gathers | 8.994 ms | Same code and proof; slower |
| New local pair mixer, two updates | 7.315 ms | Two dense points pass at 9.5%; no complete certificate |
| New local pair mixer, four updates | 8.041 ms | No completed proof diagnostic for this update count |

The shared four-update control in the mixer batch took 6.491 ms. A bounded
one-update mixer screen took 6.993 ms; this is one 31-call process, not a
repeated headline result. Neither dropping an update nor removing the
canonical BCH temporary reaches the target in these prototypes.

The [engineering report](fused_bch_REPORT.md) records ranges, phase timings,
correctness checks, sanitizer scope, hashes, and reproduction commands.
Its matched fresh times should not be conflated with the older retained
5.863 ms and 7.506 ms measurements.

## Why the Fused Loads Did Not Help

The new BCH modes reuse an existing generated helper that gathers the
independently permuted columns while preparing the GFNI inputs. With the
sequential tile copy retained, the repacking-plus-BCH phase stays near
3.7 ms. Reading scattered coordinates directly from the cold routing
buffer raises that phase to about 5.2 ms. These measurements do not isolate
individual cache or instruction effects, but they reject this proposed
fusion as an end-to-end speedup on the tested machine.

## The Corrected Pair Mixer

For each four-row BCH group, uniformly match its 256 packet positions.
For a matched pair of GF16 values x,y, independently sample nonzero a,b,e
and set U=a*x+b*y and V=e*(a*x+2*b*y). Arithmetic uses the polynomial basis
modulo X^4+X+1. An independent uniform column shuffle follows the mixer;
the regional routing and IMT maps are unchanged. Setup is fixed for all
messages.

The third scalar repairs a dependency in the initial proposal. Without e,
input (x,0) gives equal output labels. With e, common-scaling symmetry
makes the first active output uniform, and e independently randomizes the
second. Conditional on the entire support, all active packet labels are
independent and uniform nonzero. The existing conditional inner interface
therefore applies to this new outer distribution.

The local support transition is exact: no active input gives no active
output; one gives two; two give one with probability 2/15 and two otherwise.
If u positions are active and J matched pairs contain two active positions,
the new support has size 2u-2J-Z, where Z is binomial with parameters J and
2/15. The matching distribution and this transition are evaluated rationally.

For example, expected support grows from 114 to 174.115. More importantly,
fresh coupled BCH bounds and the support transition reduce the log2 upper
expected count at support threshold 114 from 309.83 to 242.90. At threshold
128 it falls from 325.32 to 270.15. These comparisons use ordinary coupled
counts, not the six stronger retained counting witnesses from the 7% proof.

The support kernel is not stochastically monotone near full support. The
implementation uses a suffix envelope when transporting upper CDFs; it
does not assume their increments are actual shell counts. Eight exact
tests cover the field distribution, conditional independence, small matching
enumerations, safe CDF transport, and cutoff retargeting.

## Distance Diagnostics

The one-layer, two-update screen regenerates coupled BCH counts and checks
the positive comparison measure, then evaluates two dense points with
256-bit outward arithmetic. The mean column refers to the comparison
distribution, not a message weight or relative distance.

| Comparison mean | log2 upper at 10% | At 9.5% | At 9% |
|---|---:|---:|---:|
| .032 | +1154.495 | -219.708 | -1593.910 |
| .096 | +2944.022 | -746.845 | -4437.712 |

The lower-distance entries reuse the same witnesses and rigorously adjust
the cutoff factor. They do not independently replay the original bounds.
Positive entries fail even a subunit probability target; they do not prove
that low-weight codewords exist. Negative entries cover only their two
points. No sparse cover or continuous dense cover was attempted, so there
is no whole-code 9.5% claim.

Two mixer layers improve the CDF caps further at thresholds 114 and 128,
but no inner diagnostic or performance implementation was completed for
two layers. A GF256 pair-randomizer control has a simpler monotone support
kernel, but weaker CDF caps at every checked threshold. It was not benchmarked.

The initial in-place mixer spent about 11.36 ms processing the cold routed
buffer. Copying each 16 KiB tile into aligned local storage, mixing it there,
and feeding BCH without a copy back reduces the time to 7.315 ms. The
mixer-plus-BCH phase still takes about 4.27 ms versus 2.67 ms for the
unmodified shared-route BCH phase. The support gain therefore has a
substantial implementation cost, even without another global routing pass.

## Reproduction and Next Step

Proof sources are [pair_mixer.py](gf16_packets/pair_mixer.py),
[pair_mixer_retarget.py](gf16_packets/pair_mixer_retarget.py), and
[test_pair_mixer.py](gf16_packets/test_pair_mixer.py). From the repository root:

```text
python -B research/workstreams/permutation_locality/gf16_packets/pair_mixer.py --layers 1 --ceil-cdf --means .032 .096 --distance .1 --output tmp/pair-mixer/fresh-layer1-screen.json
python -B research/workstreams/permutation_locality/gf16_packets/pair_mixer_retarget.py tmp/pair-mixer/fresh-layer1-screen.json --output tmp/pair-mixer/fresh-retarget.json
python -B -m unittest discover -s research/workstreams/permutation_locality/gf16_packets -p test_pair_mixer.py -v
```

Use new output paths. The retained screen conservatively rounds expected
CDF bounds upward to integers; --ceil-cdf reproduces that choice. The
current default retains exact fractions. Receipts remain ignored under
tmp/pair-mixer; raw engineering logs remain under tmp/fused-bch-results.
The stopped four-update ablation produced no result and supplies no evidence.

Preserve these experiments, but do not start a full certificate search for
the mixer solely to pursue 6 ms. At similar runtime, the independent-row
reference already has a 10% certificate. The next recommended diagnostic
is to keep the fast shared route and its outer comparison fixed, then
compare its failing classes under actual and idealized inner mixing. That
can identify whether a cheaper inner redesign is worth pursuing; an
idealized result would be a diagnostic, not a certificate for the code.
