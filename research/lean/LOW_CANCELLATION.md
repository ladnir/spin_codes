# Low-input cancellation and the Bernoulli occupation step

> **Historical snapshot — archived 2026-09-28.** Status statements and task recommendations below describe an earlier stage.
> See [STATUS.md](STATUS.md) and [NATIVE_RESULT.md](NATIVE_RESULT.md) for the closed native distance/rate result and verification scope.

The numerical transfer uses additional cancellation envelopes for input
weights 1 and 2. Their exact identification is implemented in
`ConcreteLowCancellation.lean`; the replay report records verification status.

## Certificate and coverage

`scripts/emit-low-cancellation.py` reads the hexadecimal map tables and
enumerates all 128 singleton inputs and 8,128 unordered pairs. It groups
inputs by their feedback syndrome. The singleton layer has 128 groups;
the pair layer has 7,909 groups.

Each group records its syndrome, the weight of that syndrome's expansion
under A, and each input's packed word and emitted weight. Lean checks every
record against the packed A/C evaluators. It also checks that each group's
emitted exponents match a pattern in the frozen numerical array.

Coverage does not rely on the generator's enumeration. Lean checks that
the flattened input list has length `choose(128,j)` and has no duplicates.
Each listed input has weight j. The generic `inputs_layer` theorem therefore
identifies this list's supports with the entire weight-j layer.

The duplicate check sorts the packed words and checks strict adjacency.
`LowCancellationSort` proves that the executable sort preserves the list's
multiset. Its recursion is structural so the kernel can evaluate it.
The certificate separately checks distinct syndrome keys and all five
shell histograms against the frozen `lowShells` arrays.

`LowCancellationMoments` identifies the actual fiber sum with the group sum.
Distinct syndrome keys imply that at most one group contributes to any
state. Pattern membership then gives the required maximum bound.
`LowCancellationShells` proves exact equality of the shell-average moments
with the five stored polynomials, including their normalization factors.

The generator and Python cross-checks are untrusted. The Lean statements
require kernel evaluation of all finite certificates.

## Resulting transfer statements

`fixedStep_liveDominates` in `ConcreteTransferAll.lean` supplies the actual fixed-weight
one-step domination for all 129 weights and `0 ≤ z ≤ 1`. It assumes only a
nonnegative coordinate budget and the established `LiveDominates` invariant.
There is no exceptional-weight or cancellation-table premise.

`bernoulliRow` in `ConcreteBernoulli.lean` sums directly over all 128-bit inputs, with
weight `β^card(x) (1-β)^(128-card(x))`, emitted factor
`z^card(x ∆ Aset q)`, and the exact transvection law for the next state.
`bernoulliStep_eq_mixture` proves that this kernel equals the binomial mixture
of the fixed-weight kernels.

`bernoulliStep_liveDominates` in `ConcreteOccupation.lean` applies the numerical
occupation matrix to this actual transition for `0 ≤ β,z ≤ 1`.
`bernoulli_moment_bound` bounds every finite iterate from the zero state.
At the paper's sparse parameters, `sparse_bernoulli_moment` gives

    sum_q (bernoulliStep((4/5)α, 1-(8/5)α)^R δ_empty)(q)
      ≤ 2048 (1-96α)^R,          0 < α ≤ 1/10000.

The routing and encoder experiment is now linked to this finite-step
theorem in `ConcreteRoutedMoment`; see `ROUTED_MOMENT.md`. At this checkpoint,
the concrete outer family's conditional first moments and the final occupation
and selection estimates were still required. The full distance theorem was open.

## Replay

From `research/lean`, with the prepared Peach workspace:

```sh
python scripts/emit-low-cancellation.py
python scripts/check-low-cancellation.py --blocks --assembly
python scripts/check-low-final.py
```

The first checker runs 63 certificate blocks and 13 assembly/pin modules,
then audits 16 theorem closures. It records source and object hashes in
`scripts/map_data/low_cancellation_verification.json`.
The final checker verifies those hashes, checks the fetched objects locally,
and runs the default build and `scripts/check.sh`. Its report is
`scripts/map_data/low_final_verification.json`.
Only reports marked `PASS` denote completed checks. Earlier spectrum,
Fourier, sparse-polynomial, and mathlib dependencies are reused.

Both current reports are `PASS`: all 76 modules and 16 axiom audits succeeded.
The Windows import/axiom check and the default build/invariant check passed.
Only the three standard axioms occur in the audited theorem closures.

The original distance statement pins and paper are unchanged.
