# Larger outers with the cheap GF256 inner

The question is whether changing the outer alone can raise the rate to
three quarters while retaining 10% distance and a 40-bit failure bound.
The current answer is negative for the **existing two-state proof bound**.
It is not a counterexample to the large-outer constructions.

## Construction and scope

Each outer row uses a fixed GF256 MDS code of length `b` and dimension `k`.
There are `L = 64 ceil(2^20 / (64 k))` rows and `N = b L` output symbols.
One regional route independently permutes each of the `b` columns.
The inner processes eight symbols with a four-symbol state. It uses the
repaired cubic maps, `eta=2`, and independent nonzero input/state diagonals.
Its proved return minima are `(8,6,5,4,3,2,2,1)`.

The inner, routing distribution, distance target, and setup distribution
remain unchanged. Only the outer dimensions vary. These are mathematical
screens, not performance measurements or completed certificates.

The new screen reuses the audited positive transfer. It generalizes the
source parity count, region count, and geometry to `b<=256`.
Length 256 can use evaluation at every GF256 point; no larger alphabet is
needed. An efficient implementation of that outer is a separate question.

## Larger blocks do not close the present proof

At `K=2^20`, the following table reports optimized single-weight margins
in bits. A negative entry means that this upper bound is uninformative.
No entry is a bound on the construction's actual failure from below.

| Outer | Nominal rate | Requested rate | Density .03 | .06 | .10 | .14 | .20 | .24 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| [256,176] | 68.750% | 68.085% | -4,083 | -329 | 3,213 | 3,203 | 169,365 | 362,732 |
| [256,184] | 71.875% | 71.111% | -13,116 | -18,905 | -28,472 | -41,963 | 81,448 | 237,712 |
| [256,192] | 75.000% | 74.419% | -22,067 | -37,252 | -59,628 | -86,163 | -1,326 | 121,147 |
| [256,200] | 78.125% | 78.049% | -31,035 | -55,572 | -90,602 | -129,911 | -80,510 | 11,980 |

The .06 optimization for [256,176] stopped on a line-search condition.
The other displayed optimizations reported success. These numerical
optimizers provide candidate witnesses, not proofs of global optimality.

At nominal rate 3/4, increasing the outer length gives only a modest gain:

| Outer | Density .10 | Density .14 |
|---|---:|---:|
| [128,96] | -62,763 | -89,761 |
| [192,144] | -60,326 | -86,913 |
| [256,192] | -59,628 | -86,163 |

The exact source calculation excludes impossible nonempty MDS supports
when the source tilt permits the established homogenization argument.
The screen also uses the central-mode conditioning factor for each region.
The large deficits therefore persist after those existing refinements.

Sparse classes fail as well. The entries below cover one active-row count,
not a sum over active-row counts.

| Outer | A=1 | A=8 | A=32 | A=64 | A=128 |
|---|---:|---:|---:|---:|---:|
| [256,176] | 418 | 663 | -72 | -690 | -1,513 |
| [256,184] | 354 | 401 | -838 | -2,160 | -4,430 |
| [256,192] | 290 | 152 | -1,566 | -3,560 | -7,208 |
| [256,200] | 226 | -83 | -2,258 | -4,891 | -9,849 |

## What is being lost?

The relaxed row-source factor has two positive terms: an empty row and an
occupied row. Replacing their log-sum by their maximum removes at most one
bit per outer row. This replacement is an **optimistic diagnostic**, not
a valid finite upper bound. It approximates the large-block source envelope.

For [256,192], optimizing this diagnostic still gives roughly -55,365 bits
at density .10 and -81,275 bits at density .14. Thus the finite entropy of
choosing occupied outer rows does not explain most of the current gap.
Simply making blocks larger cannot be expected to supply the missing gain
through that term alone.

At the [256,192], density-.14 witness, the rank tilt is about 11.72, rather
than the available maximum 255. Raising it would reward outer parity more
strongly, but makes the inner's rank-potential bound more expensive.
This identifies a coupled source/inner proof loss. It does not establish
that the actual four-symbol-state code needs the same loss.

An MDS outer already maximizes the rank of every coordinate restriction
allowed by the dimensions. Replacing it by another deterministic linear
outer cannot improve that worst-support rank interface. A genuinely
different outer would need a different proof interface, such as an
ensemble-average argument that retains more information than support rank.

For comparison, the GF256 Hamming-ball exponent at relative radius .10 is
`H2(.10) + .10 log2(255) = 1.2684309373` bits per symbol. The corresponding
random-linear-code asymptotic rate is about 84.145%. Three-quarter rate is
therefore not ruled out by this reference bound. The present short-state
transfer is much more restrictive than that benchmark.

## Distinguishing a real window obstruction

There is a separate, deterministic failure test. Fix one packet-aligned
window in each region, with total length at most `ceil(N/10)-1`. Let `h_i`
be the number of outer-row coordinates landing in these windows and put
`p=b-k`. The shortened MDS source space has dimension

`D = sum_i max(0, h_i-p)`.

Requiring the inner state to be zero after each window imposes at most
`4b` linear equations. Consequently, `D>4b` gives a nonzero codeword
confined to the windows, for every fixed choice of the diagonal values.
This is an actual distance obstruction for the realized route.

For independent uniform regional permutations and fixed equal-width
windows of occupancy `a`, each `h_i` is Binomial(`b,a`). Hence

`E[D] = L sum_{h=p+1}^b (h-p) binom(b,h) a^h (1-a)^(b-h)`.

This expectation concerns fixed windows, not windows chosen after seeing
the route. It neither certifies the code nor excludes an adaptive attack.

At nominal rate 3/4, packet-aligned windows just below 10% give:

| Outer | Fixed-window E[D] | Exit-state constraints | D for route seeds 1,2,3 |
|---|---:|---:|---:|
| [32,24] | 187.78 | 128 | 176, 163, 176 |
| [64,48] | 3.9922 | 256 | 8, 5, 3 |
| [128,96] | 0.00380 | 512 | 0, 0, 0 |
| [256,192] | 7.15e-9 | 1024 | 0, 0, 0 |

The [32,24] routes are concrete obstructions. Their guaranteed kernel
dimensions are 48, 35, and 48. The seeds specify realized routes; three
samples are not a quantified failure-probability claim. The larger-outers'
failed first-moment screens have no such counterexample here.

## Reproduction and next step

Run from the repository root:

```powershell
python research/workstreams/char2_field_spin/theory/gf8_high_rate_outer_tests.py
python research/workstreams/char2_field_spin/theory/gf8_high_rate_outer_screen.py
python research/workstreams/char2_field_spin/theory/gf8_high_rate_outer_screen.py --mode sparse
python research/workstreams/char2_field_spin/theory/gf8_high_rate_outer_screen.py --outer 256 192 --source bulk --density .10 .14 --iterations 100
python research/workstreams/char2_field_spin/theory/gf8_high_rate_outer_screen.py --mode window --outer 32 24 --window-seed 1
```

Seven tests pass. They cover length-128 agreement with the existing middle
and sparse screens, padded geometry, the log-sum diagnostic, exact window
expectations, and the realized short-outer obstruction.

The next useful step is the rank-history refinement of the same cheap
inner, or a local map/state change that allows a substantially larger rank
tilt. Another one-dimension outer sweep is unlikely to produce the desired
rate jump. No frozen certificate or production implementation was changed.
