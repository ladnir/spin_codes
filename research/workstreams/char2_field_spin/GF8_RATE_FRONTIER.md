# GF256 rate frontier at 10% distance

The K20 rate climb now has a complete certificate at nominal rate `37/44`
(84.0909%), 10% symbol distance, and 183 bits of setup-failure margin.
The GF256 random-linear-code benchmark at this distance is 84.1446%.
This is a proof checkpoint, not a production profile. The subsequent
[native comparison](GF8_RATE_NATIVE.md) implements three higher-rate profiles
and recommends 84% for its better rate--time tradeoff.

## Complete checkpoints

All rows below target requested dimension `K=2^20`. The certificate covers
the entire padded code and therefore its shortened requested-dimension subcode.
Rate means input symbols divided by retained output symbols. Weight counts
nonzero GF256 symbols, not binary bits.

| Outer GRS | Retained outputs per eight-input packet | Nominal rate | Padded dimension | Output length | Requested rate | Whole-code margin |
|---|---|---:|---:|---:|---:|---:|
| [128,84] | 7 | 75.0000% | 1,053,696 | 1,404,928 | 74.6356% | 116 bits |
| [128,80] | 6 | 83.3333% | 1,049,600 | 1,259,520 | 83.2520% | 264 bits |
| [128,67] | 5 | 83.7500% | 1,050,560 | 1,254,400 | 83.5918% | 368 bits |
| [128,84] | repeat (6,6,6,7) | 84.0000% | 1,053,696 | 1,254,400 | 83.5918% | 115 bits |
| [128,74] | repeat (5,6) | 84.0909% | 1,051,392 | 1,250,304 | 83.8657% | 183 bits |

The first row is the preserved [three-quarter checkpoint](GF8_PUNCTURED_WHOLE_CERTIFICATE.md).
The other four rows were replayed through the new whole-code aggregator.
These margins bound the probability, over the once-sampled setup, that any
nonzero message produces an output below the stated distance. They are not
per-message probabilities or measured empirical failure rates.

Each 128-bit application record carries sixteen lanes of the same sampled
GF256 code. A nonzero lane's support is contained in the record support,
so the certificate also gives 10% nonzero-record distance. This extension
requires no sixteen-fold setup-failure union.

## Construction changes

The outer is a fixed GRS[128,d] constituent, repeated over padded source rows.
One independent uniform permutation of row labels is sampled per outer
coordinate. Together these form the single regional routing stage.
The inner keeps the frozen full-Cauchy eight-input/eight-state maps.
Input and state diagonals have independent uniform nonzero GF256 entries.

Only the outer dimension and retained-output schedule change across the table.
An entry `h` means retaining the first `h` of the packet's eight output symbols.
A periodic schedule is public and deterministic; it adds no randomness or
permutation. The recursive state continues across packets and period boundaries.
The proof multiplies packet transfers in chronological order, not their average.

The shared map digest is
`6e3cf3e5d6b5876838afc5b2d7230f931974f363f9f67b402843f677521450a1`.
The aggregator checks the full Cauchy factorization and the causal pivot for
every retained-output count. Puncturing is not assumed to preserve dimension:
the whole-code bad event includes nonzero messages that would map to zero.

## What limited the climb

Keeping seven outputs and weakening the outer stalls in the current sparse
bound around 78--79% rate. The failed bounds are not low-distance witnesses.
Keeping fewer outputs permits a stronger outer at the same final rate;
this moves the complete frontier past 83%. Periodic retention supplies finer
rate choices without another routing stage.

Near 84%, the dense bound becomes important. Its old MDS dimension bound
overcounts full-support words by a factor involving `(256/255)^W`.
That loss prevents the alternating five/six candidate from closing.
The [exact-spectrum refinement](GF8_EXACT_SPECTRUM_DENSE.md) removes this
overcount. It separates a constant membership baseline from nonnegative
excess terms and bounds unequal regional densities with AM-GM and Maclaurin.
This is a proof improvement for the same sampled construction.

For the highest-rate checkpoint, the three complete components are:

| Component | Covered class | Margin |
|---|---|---:|
| Sparse | 1 through 991 active outer rows | 312 bits |
| Middle | raw weights 54,559 through 1,636,761 | 32,550 bits |
| Dense | raw weights 1,636,762 through 1,818,624 | 184 bits |

An active outer row has weight at least 55. Thus every word outside the sparse
range has weight at least `55*992=54,560`, within the middle/dense cover.
The dense proof includes all omitted-output completions, including those
with retained output zero. Adding the three outward-rounded failure bounds
as exact rationals gives the conservative whole-code margin of 183 bits.

The sparse, middle, and dense partition digests are, respectively:

```text
bff62cfe77bc98c86f283dd6d5e3406e5afb0d12c023c659bc508e1d3b01b9ee
aa679e45f42e36af5a41ba4d386dd207c1a7737e7005da429c08bdc7db8cbd44
61c3a29a382a27cb266b31f1ba857cd50ad45d99e65eab93226a505f82d83136
```

See the [sparse record](GF8_RATE_CLIMB_SPARSE.md) and
[middle/dense record](GF8_RATE_CLIMB_COMPONENTS.md) for the other profiles.

## How close to the random-code benchmark?

For GF256, the random-linear-code/Gilbert--Varshamov reference at distance
`delta=0.1` is `1-H_256(delta)=0.8414461328406041`.
Here `H_q(delta)=(delta*log2(q-1)-delta*log2(delta)-(1-delta)*log2(1-delta))/log2(q)`.
Random linear codes achieve this asymptotic reference; it is not a universal
upper bound on all structured codes. See
[Guruswami's coding-theory notes](https://www.cs.cmu.edu/~venkatg/teaching/codingtheory/notes/notes2.pdf).

Our nominal rate is lower by **0.053704 percentage points**. For exactly K20,
the current 2,816 padding symbols lower the requested rate to `1024/1221`
(83.8657%). A separate directed random-generator first-moment calculation
gives a sufficient 40-bit reference at output length 1,246,165 for exactly K20.
Our 1,250,304 outputs are about 0.332% more. That comparison does not transfer
a random-generator proof to SPIN.

The remaining improvements are correspondingly small. The current dense
output-ball bound uses 1.27 bits per output instead of the entropy value
1.2684309373. Regional conditioning and interval endpoints also cost bits.
The constant 1.27 alone limits this dense argument to an asymptotic rate
below 84.125%. A sharper ball bound is therefore the first proof refinement
before pursuing the final few hundredths of a percentage point.
None of these observations proves that the construction has reached its
actual maximum rate.

## Reproduction and implementation follow-up

Run from the adjacent `theory` directory, with the existing Python proof
dependencies (`python-flint`, NumPy, SciPy):

```text
python -B gf8_rate_climb_certificate.py --profile keep6
python -B gf8_rate_climb_certificate.py --profile keep5
python -B gf8_rate_climb_certificate.py --profile periodic84
python -B gf8_rate_climb_certificate.py --profile periodic74
python -B gf8_random_reference.py
```

The full aggregator recomputes every component. Synthetic receipt tests only
test aggregation logic; they are not substitutes for these replays.
Separate exact tests cover prescribed-support MDS counts, unequal regional
profiles, chronological packet products, interval domination, and rejected
coverage/metadata/completion-factor mutations.

No native timings accompanied this proof checkpoint. The subsequent
[matched comparison](GF8_RATE_NATIVE.md) measures 4.807--4.849 ms for 84%
and 5.318--5.350 ms for 84.09%, including compact retained-input loads, routing,
outer encoding, and stores. The 84% design processes 1,605,632 raw symbols,
versus 1,818,624 for 84.09%: almost equal retained lengths do not imply
equal encoding work. Frozen certificates and production code are preserved.
