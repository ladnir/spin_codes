# The invariant needed by the concrete transfer law

The seven coordinates bound a weighted state measure `μ`: `Z` bounds its
mass at zero, `D` bounds a nonnegative measure supported on nonzero states,
and each `Sᵢ` bounds a multiple of the uniform measure on one nonzero shell.
The concrete matrix propagates nonnegative coordinate budgets.

The earlier `Dominates` definition and `imt_moment` theorem are valid as
written, but their interface is insufficient for this application in two
ways. First, the one-step premise of `imt_moment` quantifies over negative
budgets too. `Spin.Imt.unrestricted_step_forces_rowZ_D_nonpos`
proves that this premise forces `T.rowZ.D ≤ 0`:
the permitted budget `(-1,0,0)` would otherwise propagate to a negative
diffuse budget. This conflicts with a transfer that moves positive mass
from zero into the diffuse coordinate.

Second, `Dominates` does not require its diffuse witness to vanish at zero.
It therefore allows `D` to cover zero mass, contrary to the representation
used to derive the matrix's rows. Merely requiring nonnegative budgets
does not repair this support condition.

`LiveInduction.lean` adds `LiveDominates sys c μ`, defined as the conjunction
of `Dominates sys c μ` and `μ ∅ ≤ c.Z`. The theorem
`liveDominates_iff_witness` proves that this is exactly domination with a
nonnegative diffuse witness that vanishes at zero. The predicate inherits
the terminal mass bound and is preserved by nonnegative scaling, finite
sums, and coordinatewise increases.

`live_moment_eZ` propagates this invariant together with coordinate
nonnegativity. Its one-step premise applies only to nonnegative budgets
satisfying `LiveDominates`. `Spin.Imt.Occupation.Sparse.sparse_live_moment_of_step`
then gives the paper's
`2048 (1-96α)^R` bound from that premise and the completed numerical matrix
certificate. `sparse_bernoulli_moment` in `ConcreteOccupation.lean` now
discharges that premise for the actual Bernoulli-input kernel.

`LiveKernel.kernelApply_liveDominates` now lifts zero, diffuse, and shell
row bounds through a nonnegative transition kernel. `ConcreteTransfer`
supplies these bounds for the actual fixed-weight step. The extension in
`ConcreteTransferAll` includes weights 1 and 2 using the checked low-input
tables. `ConcreteOccupation` combines all 129 weights and iterates the
resulting Bernoulli-input law. See `CONCRETE_STEP.md` and `LOW_CANCELLATION.md`.

The original definitions, induction theorem, and statement pins remain
unchanged. No alteration of the paper's distance claim is involved.
