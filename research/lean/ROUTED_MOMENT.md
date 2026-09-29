# Concrete routed IMT moment

> **Historical snapshot — archived 2026-09-28.** Status statements and task recommendations below describe an earlier stage.
> See [STATUS.md](STATUS.md) and [NATIVE_RESULT.md](NATIVE_RESULT.md) for the closed native distance/rate result and verification scope.

The finite routed encoder now has a proved moment bound. Fix an `L` by `b`
matrix of input bits, represented by `rows : Fin L → Finset (Fin b)`. Let
`Q` count its nonempty rows, `w` be its total weight, and `q = w/(Lb)`.
Independently shuffle each row, transpose to `b` regions of length `L`, and
independently shuffle each region. Serialize the regions in increasing order
and split the stream into `R = Lb/128` blocks.

For each block `X`, the encoder emits `X ∆ Aset state` and updates the state
to `act u v state ∆ Cset X`. It starts at the empty state and carries that
state across every block and region boundary. It emits no final tail.
Each round samples an independent uniform valid pair `(u,v)` with `u ≠ ∅`
and even intersection cardinality. The maps are the paper's checked concrete
19-to-128 and 128-to-19 maps.

## Proved bound

`ConcreteRoutedEncoder.stream_moment_bound` proves

    E[z^W] ≤ (b+1)^Q (L+1)^b · (T_occ(q,z)^R e_Z).total

for positive `L,b`, `128 ∣ b*L`, `0 < w < L*b`, and `0 ≤ z ≤ 1`.
Here `W` is the actual emitted weight. The expectation is over the product
of all row permutations, region permutations, and transvection pairs.
The numerical occupation matrix is already linked to the actual maps and
all 129 fixed input weights. This theorem has no route-identification,
one-step domination, cancellation-table, or finite-certificate premise.

`moment_bound` permits any fixed bijection from region positions to block
positions. `regionMajor_position` checks the concrete serialization by
the identity

    regionOffset + L*region = blockOffset + 128*round.

`sparse_moment_bound` additionally proves

    E[(1-(8/5)α)^W] ≤ (b+1)^Q (L+1)^b · 2048 (1-96α)^R

when `0 < α ≤ 1/10000` and **the supplied matrix has density `(4/5)α`**.
This is a specialization in the input density. It does not identify `α`
with the fraction `Q/L` of active rows in an arbitrary outer profile.

## Proof interfaces

| Modules | Checked connection |
| --- | --- |
| `FiniteLaw`, `FiniteProductLaw` | Pushforwards, conditional sampling, product laws, and preservation of domination |
| `ConcreteShuffle`, `ConcreteRowLaw` | Exact uniform weight-layer law and row cost `1` for an empty row, `b+1` otherwise |
| `ConcreteShuffleMixture`, `ConcreteShufflePermutation` | Independent uniform permutations produce the layer laws; shuffling preserves weight probabilities |
| `ConcreteRouteTranspose`, `ConcreteRoute`, `ConcreteRouteDomination` | The transpose of the independent comparison cells and the complete `(b+1)^Q (L+1)^b` domination bound |
| `ConcreteRoutePermutation` | The algebraic route law equals the explicit independent-permutation experiment |
| `ConcreteEncoder`, `ConcreteEncoderMoment` | The explicit recurrence's generating function equals the iterated Bernoulli kernel and obeys the numerical matrix bound |
| `ConcreteReshape`, `ConcreteSerialization` | Fixed rewiring preserves independent Bernoulli bits; region-major order is explicit |
| `ConcreteRoutedMoment` | The combined finite routed-encoder theorem |

`shuffleSupport` uses the pullback convention for a permutation. The paper's
coordinate-forward convention is obtained by inverting each seed; uniform
permutations have the same law under inversion. The implemented route and
its uniform law are explicit, rather than supplied as identification
hypotheses. The encoder uses uniform valid transvection pairs; the existing
equal-cardinality perpendicular fibers give the paper's equivalent sequential
sampling interpretation.

## Verification

Run `C:/Python314/python.exe scripts/check-routed-bridge.py` from this directory.
The script sequentially replays all 19 new modules on Peach, checks the 26
consolidated theorem axiom closures, checks `ConcreteRoutedPin.lean` on
Windows, and runs the default build and `scripts/check.sh`. It also checks
the previous 76-module low-weight checkpoint hashes, preserves the original
`SpinCodes/Pin.lean`, and verifies that the paper files remain unchanged.

The report is `scripts/map_data/routed_bridge_verification.json`; only
`status: PASS` records completed verification. Earlier numerical certificates
are reused. This replay does not rebuild every spectrum or polynomial block.
The allowed axiom closure is `propext`, `Classical.choice`, and `Quot.sound`.
These modules are checked by the supplemental replay; they are not added to
the default `SpinCodes` imports.

The current report is **PASS**: all 19 modules, all 26 axiom audits, the
Windows statement check, and the default project checks succeeded. No
verification job remains running.

## Distance-proof obligations at this checkpoint

The theorem fixes the supplied input matrix. It does not construct or
average over the paper's shared outer-code realization. At this checkpoint,
the full distance theorem remained conditional on the following connections:

1. Construct the native-schedule `Family`, its shared Golay-BAA realization,
   occupation statistic, and exact `floor(0.11N)` threshold.
2. Connect the finite outer spectra to the uniform asymptotic bounds that
   make selection succeed with probability tending to one.
3. Prove the sparse fair-row message-counting comparison and marked-position
   conditioning. The generic route cost `(L+1)^b` is too costly for this
   regime; the density specialization above does not discharge it.
4. Connect fixed occupation to isolated-impulse kernels and their uniform
   finite-to-continuum remainder.
5. Connect positive occupation to density changes, KL costs, profile counts,
   certified transfer witnesses, and uniform remainders.

These connections were intended to supply the conditional first moments
required by `Family.distance_whp_native`. The planned parallel milestone
separated the family/selection assembly, sparse counting-and-conditioning
bridge, and fixed-occupation kernel bridge.
