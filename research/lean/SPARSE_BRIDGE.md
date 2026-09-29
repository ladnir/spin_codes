# Sparse bridge: complete numerical-matrix contraction

> **Historical snapshot — archived 2026-09-28.** Status statements and task recommendations below describe an earlier stage.
> See [STATUS.md](STATUS.md) and [NATIVE_RESULT.md](NATIVE_RESULT.md) for the closed native distance/rate result and verification scope.

The target remains the paper's concrete distance theorem. The complete sparse
polynomial program and its comparison with all seven coordinates of the
occupation matrix are now proved. `Sparse.sparse_collatz` proves T3a for the
frozen numerical matrix; `Sparse.sparse_iterate_bound` gives its repeated-step
bound with prefactor 2048.

## Definitions and real-valued proofs

`Structured/Occupation.lean` defines the paper's fixed-weight matrix and its
binomial mixture. Its shell row preserves the entering shell when the live
syndrome mass is zero. The three column-action identities, nonnegativity of
the fixed matrix, and preservation of column bounds under the binomial
mixture are proved.

`SparseWitness.lean` fixes the exact shell counts and affine witness. It proves
the zero shell-weighted correction, the constant live average `1/1024`, and
the lower bound `1/2048` on every coordinate throughout the closed interval.

`SparseModel.lean` defines the numerical matrix using the actual maxima and
minima from the appendix. Its kernel counts and fiber caps are still inputs;
the frozen instances are in `SparseModelData.lean`. The theorem
`SparsePolynomial.eval_zero_contribution` identifies each zero-coordinate
polynomial contribution with its weighted matrix action.

`SparsePolynomial.lean` identifies the executable hypergeometric moments,
binomial probabilities, and affine coordinates with their real-valued
formulas. The executable factorial quotient is proved equal to `Nat.choose`.

## Exact arithmetic

`PolyIdentityDefs.lean` and `PolyIdentity.lean` provide integer coefficient
arithmetic with common positive denominators and prove evaluation soundness.
Addition uses an LCM; equality ignores trailing zeros. The LCM quotients are
computed in `Nat` before conversion to `Int`.

`SparsePowers.lean` supplies the coefficient lists for the powers of
`1-x/6250`, `x/12500`, and `1-x/12500` through degree 128. All 384 successive
multiplications are kernel-checked. `SparsePowersSound.lean` identifies the
tables with the corresponding real powers for every argument.

`SparsePolynomialDefs.lean` is the executable sparse verifier, with explicit
branch choices. Each generated `SparseBridge/Weight*.lean` checks its seven
action polynomials and its binomial probability polynomial. There are 129
weight modules and 1,032 identities. These are coefficient identities in
Lean, not hash comparisons.

`PolySignIdentity.lean` proves the nonpositive coefficient certificate after
removing any initial zeros. The generated `Dominance*.lean` modules check
the 693 maximum comparisons and the moment/pattern identities they use.
`SparseDataChecks.lean` checks all 129 kernel-count bounds and choice-index
ranges. The kernel-count bound is needed when relating integer subtraction
in the program to natural subtraction in the fiber-cap formula.

`SparseMaxima.lean` and `SparseMaximumAll.lean` turn all maximum comparisons
into real inequalities on the whole closed interval, for every input weight.
The additional pattern maxima at weights 1 and 2 are included. The latter
file also defines `Sparse.numericalMatrix` from the frozen finite data.

`PolyPacked.lean` proves a faster exact check for sums of polynomial products.
The radix exceeds the sum of absolute coefficients of the difference. A zero
integer encoding then implies the zero polynomial. Every denominator and
scaling factor is checked. This replaces the prohibitively expensive direct
kernel convolutions for the degree-257 sums.

All 119 product-sum chunks and all seven final residual identities pass kernel
checking. `SparseContributionSound.lean` connects the chunks to the original
program contributions. `SparseResidualsSound.lean` restores each initial factor
of `x`, divides by its checked positive denominator, and applies the existing
strict sign certificate. `SparseProgram.lean` proves that the grouping covers
`List.range 129` exactly and obtains

    programResidual_neg (i : Fin 7) (0 < x) (x ≤ 1).

Here `programResidual` is the sum of the 129 original executable contributions
minus `(1 - 96*(x/10000)) * witness_i(x)`. No coefficient identity or sign
condition is assumed. The generator also cross-checks this transcription
against `SPARSE_EXACT.json`, but those Python checks are not used as Lean proofs.

Reports are in `scripts/sparse_data/`. A report with status `RUNNING` is not a
completed replay. `SparsePin.lean` and `SparseBridge/Pin.lean` supplement the
original statement pins; the original pins have not been changed. The initial
weight replay had 27 missing-import failures caused by a concurrent dependency
build. The successful retries are retained separately, and
`kernel_weights_complete.json` combines the 129 successful checks.

## Completed matrix connection

`SparseCancellation.lean` proves positive denominators and converts the
checked finite choices into bounds for the cancellation and fresh-state
minima. It uses the kernel bound to identify natural subtraction in the real
shell-cap formula with integer subtraction in the polynomial program.

`SparseColumnBounds.lean` evaluates all six live polynomial actions and bounds
the corresponding fixed-weight matrix columns. It verifies the equivalence
of the zero-live-mass branch and `kernel = choose(128,j)`; both endpoint input
weights therefore retain the correct lazy shell behavior.

`SparseContraction.lean` combines those bounds with the zero-row identity,
the nonnegative binomial mixture, exact range-sum conversion, and all seven
strict program residual inequalities. Its theorem `Sparse.sparse_collatz`
states the paper's matrix inequality directly in `α`, for `0 < α ≤ 1/10000`.

`SparseNonneg.lean` proves nonnegativity of the fixed and mixed matrices.
`SparseIteration.lean` then proves `e_Z T_occ^R 1 ≤ 2048 (1-96α)^R` for every
natural `R`. The supplemental `SparseBridge/Pin.lean` pins both conclusions.

The actual maps, spectra, and fiber counts are now connected to the frozen
data. The actual fixed-weight step is proved for all 129 input weights,
including both special cancellation tables. `ConcreteOccupation` applies
the sparse iterate bound to the direct Bernoulli-input transition, with no
one-step premise. At this checkpoint, the concrete family, route linkage,
and final conditional first-moment estimates were still required for the
full distance theorem.
See `LOW_CANCELLATION.md` and `CLOSURE_AUDIT.md`.

## Replay

From `research/lean`:

```sh
bash scripts/check-sparse.sh
bash scripts/check.sh
```

For an unchanged numerical certificate, replay just the new semantic bridge
with `python -X utf8 -B scripts/check-sparse-semantic.py`. Its report records
the source hash, exit code, and log of each module in dependency order.

Directly evaluating a recursively expanded polynomial power in a kernel
decision repeats substantial work. Keep the explicit power tables and their
checked recurrences. Likewise, `polyEq` must remain structurally recursive;
the first formulation introduced an opaque accessibility proof that stopped
kernel reduction even on small inputs. Use the packed checker for the large
product sums. Finish all dependency builds before launching numerical workers;
the runner now records and checks dependency hashes.
