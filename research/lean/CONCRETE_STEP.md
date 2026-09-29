# Concrete fixed-weight transfer

> **Historical snapshot — archived 2026-09-28.** Status statements and task recommendations below describe an earlier stage.
> See [STATUS.md](STATUS.md) and [NATIVE_RESULT.md](NATIVE_RESULT.md) for the closed native distance/rate result and verification scope.

`fixedStep_liveDominates` in `ConcreteTransferAll.lean` proves the paper's
one-step domination for the actual maps, for every input weight `j` from
0 through 128, and every real `z` with `0 ≤ z ≤ 1`.
It applies to any nonnegative coordinate budget `c` and any state measure
dominated by `c` under `LiveDominates`. Its conclusion uses the existing
numerical fixed-weight matrix, with no row-domination premise.

The special tables at weights 1 and 2 are now identified with the actual
maps by certificates covering all 8,256 low-weight inputs. See
`LOW_CANCELLATION.md` for coverage, pattern maxima, and exact shell sums.
The earlier `fixedStep_liveDominates_generic` theorem remains available.

## The actual transition

For an entering state `q`, the input `x` is uniform among the 128-bit words
of weight `j`. The emitted factor is `z ^ card(x ∆ Aset q)`. The next state
is `act p.1 p.2 q ∆ Cset x`, where `p` is uniform over the previously defined
transvection pairs. `ConcreteStep.stepRow` uses the exact finite transvection
law to assign weighted mass to each target. `fixedStep` applies this kernel
linearly to an entering measure.

`ConcreteMoments` proves that the emitted moment is exactly the
hypergeometric formula in the numerical matrix. It also bounds each emitted
factor by the power determined by the difference between the two weights.
`ConcreteCancellation` derives the general pointwise and shell-average
cancellation bounds from the actual kernel counts, fiber caps, and shells.

`ConcreteStep` proves the exact lazy/refresh decomposition for every nonzero
state. It proves the exact zero row and the refresh zero-column bound by
the minimum of the emitted moment and the nonzero-syndrome probability.
`ConcreteRows` and `ConcreteShellRows` turn these facts into coordinate
budgets. `ConcreteZeroRow` and `ConcreteClosedShells` handle the case where
every input in a layer has zero syndrome: the lazy branch then preserves
the entering shell, as required by the matrix.

`LiveKernel.kernelApply_liveDominates` combines the zero, diffuse, and shell
row bounds for an arbitrary dominated measure. It uses the diffuse witness's
zero support condition and nonnegative coordinate budgets explicitly.
`ConcreteCancellationRows` checks that only weights 1 and 2 use special
tables and proves the remaining cancellation inequalities. `ConcreteTransfer`
then supplies the actual row bounds to the kernel theorem.

`fixedStep_liveDominates_of_cancellation` also handles any input weight once
its pointwise and shell-average cancellation inequalities are supplied.
`ConcreteLowCancellation` supplies those inequalities at all weights, and
`ConcreteTransferAll` discharges the cancellation premises.

## Verification and work remaining at this checkpoint

Run `python scripts/check-concrete-step.py` from the local Lean directory
after preparing the Peach workspace described in `CONCRETE_MAPS.md`.
The checker recompiles eleven semantic/pin modules sequentially on Peach,
downloads their objects, and audits seventeen theorem closures. Existing
numerical dependencies are reused. The source/object hashes and per-module
logs are recorded in `scripts/map_data/concrete_step_verification.json`.
Only a `PASS` report denotes a completed replay.
The current report is `PASS`: all eleven modules and all seventeen axiom
audits succeeded, with only standard axioms in their closures.
The local import/pin check and the default project invariant check also
passed; `step_final_verification.json` records those checks and confirms
that the local source/object hashes still match the semantic replay.

The additional low-weight and Bernoulli connection uses
`scripts/check-low-cancellation.py --blocks --assembly`, followed by
`scripts/check-low-final.py`. Its status and hashes are recorded in
`low_cancellation_verification.json` and `low_final_verification.json`.

`ConcreteBernoulli` identifies the binomial mixture with the direct transition
for independent Bernoulli input bits. `ConcreteOccupation` supplies the
one-step law and the explicit sparse finite-iterate bound. `ConcreteEncoder`
and `ConcreteRoutedMoment` now connect them to the explicit routed encoder;
see `ROUTED_MOMENT.md`. At this checkpoint, the concrete outer family and
final selection/occupation first moments and asymptotics were still required
for the unconditional distance theorem.
