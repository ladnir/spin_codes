# Native construction semantic review

> **Historical snapshot — archived 2026-09-28.** Status statements and task recommendations below describe an earlier stage.
> See [STATUS.md](STATUS.md) and [NATIVE_RESULT.md](NATIVE_RESULT.md) for the closed native distance/rate result and verification scope.

Reviewed 2026-09-28 against `../paper/structured_proof.tex`,
`../paper/structured_spin.tex`, and `../paper/structured_imt_appendix.tex`.
This is a source-level semantic review, separate from the kernel checks.
No construction mismatch was found for the native minimum-distance and
rate claims. The paper's additional encoder-work claim remains outside the
reviewed formal theorem.

The probability to be closed is over one shared Golay-BAA outer seed and
independent route/IMT seeds. It concerns the minimum Hamming distance of the
actual emitted binary linear code, at the exact threshold
`floor (0.11 * Nsched m)`. It is not a theorem about an abstract surrogate
matrix or about a code conditioned on a successful setup.

| Paper object | Lean definition or theorem | Finding |
| --- | --- | --- |
| Native schedule and dimensions | `Schedule.Lsched`, `SchedOK`, `bsched`, `Nsched`; `ConcreteOuterNative.native_dimension`, `native_width`, `native_round_divisibility` | `Lsched m = 128*(m+1)` is the paper's positive index shifted by one. `bsched` uses `Nat.find` for the least admissible positive multiple of 24, with the exact coefficient `39/4`. Input dimension is `N/2`, output length is `N`, and all regions contain an integral number of 128-bit rounds. |
| Extended Golay direct sum followed by two permute–accumulate stages | `ConcreteOuter.encode`, `seedLaw`; `ConcreteOuterGolay.golay_injective` | `encode seed x = accF ((accF (outerWord Golay.enc x ∘ seed.1)) ∘ seed.2)`. The two permutations are independent and uniform. |
| Same sampled outer in every diagonal position | `ConcreteOuter.rowSupports`; `ConcreteOuterNative.nativeRows`, `nativeSeedLaw` | One `seed` is passed to every row. There is no independently resampled outer per row. Counting uses the realized shared spectrum. |
| Independent row permutations, transpose, independent region permutations | `ConcreteRoutePermutation.Seeds`, `seedLaw`, `routeEval`, `permutation_law_eq` | The exact permutation experiment is defined first; its equality to the support-subset law is proved. `ConcreteShufflePermutation.shuffleSupport` uses the pullback convention `S.map σ.symm`, consistently with coordinate precomposition. |
| Region-major serialization | `ConcreteSerialization.regionMajor`, `regionMajor_position`, `streamWiring` | The proved index identity is `outputBit + 128*round = withinRegionBit + L*region`. It matches the paper's consecutive regions. |
| IMT randomness | `ConcreteEncoder.Transvection`, `transvectionLaw`; `Transvection.transPairs`, `card_perp_eq`; `ConcreteRoutedMoment.experimentLaw` | A round samples uniformly from pairs with nonzero `u` and even overlap of `u,v`; `v=0` is allowed. Equal orthogonal-fiber cardinalities make this the same distribution as uniform nonzero `u`, then uniform `v` in its orthogonal subspace. Products give independence across rounds and from the route; `nativeSetup.Pout` and `Pin` are independent in the joint law. |
| Fixed maps | `ConcreteMapData.aRows`, `cTransposeRows`, `cRows`; `ConcreteMaps.Aset`, `Cset` | All 19 expansion basis words and all 19 feedback-transpose basis words agree numerically with the paper table, including its leading-zero conventions. Feedback is obtained by transposing the specified feedback-transpose basis. |
| Inner recurrence | `ConcreteEncoder.nextState`, `outputWeight`; `ConcreteEncoderOutput.outputBlocks` | The emitted block is `X Δ A(Q)`; the next state is `act u v Q Δ C(X)`. Feedback consumes raw input. The initial state is empty, state persists across region boundaries, and the final state is discarded without a flush. |
| Actual emitted word and linear code | `ConcreteNativeCodeword.codeword`, `codeword_weight`, `codeword_injective`; `ConcreteNativeLinearCode.encoderLinearMap`, `realizedCode` | The finite experiment is connected to an explicit length-`N` binary word. Injectivity and linearity hold for every seed realization. |
| Rate and distance event | `ConcreteNativeLinearCode.realizedCode_finrank`, `realizedCode_rate`; `ConcreteNativeLinearCodeDistance.minimumDistance_le_iff_pair`, `concrete_probBad_minimumDistance`; `ConcreteNativeFamily.threshold_floor` | Every realized code has exact rate `1/2`. The first-moment bad event is exactly minimum pairwise Hamming distance at most the natural floor of `0.11*N`; it is not merely a message-weight proxy. |

`selectedGood` in `ConcreteNativeTotal` uses the whole space only when the
raw selection event has zero probability, making conditioning total at every
finite index. This does not alter the encoder or its unconditional bad-event
probability: `concrete_probBad_eq` proves independence from that auxiliary
choice. `ConcreteOuterTailClosure.native_selection_failure_tendsto` and
`native_good_positive_eventually` remove this small-index device
asymptotically. No rejection-sampled setup is substituted for the paper's
setup.

The remaining assumptions should be read from the final theorem, rather than
from earlier closure interfaces. In the reviewed source snapshot,
`ConcreteNativeDenseRates.minimum_distance_of_fixed_and_dense_rates` still
accepts fixed-occupation limits and a `DenseRates` certificate with positive
constants. The older `ConcreteNativeFixedLargeClosure.native_large_EZ_tendsto`
also exposes a `FixedLargeNorm` premise. Concurrent final assembly now supplies
that norm through `ConcreteFixedLargeNorm.fugacityContinuum_paper_rowNorm`,
and `ConcreteNativeFixedAll.native_fixed_EZ_tendsto` combines the three fixed
cases. Their kernel status must be taken from the owning workers' audit logs;
this review does not promote newly written source to a checked theorem.
The scalar-to-matrix interface itself is checked: the six-module
`encoder_fixed_insertion_norm_verification.json` records seven axiom audits,
all limited to `propext`, `Classical.choice`, and `Quot.sound`.

The principal scope limitation is computational work, which `STATUS.md`
already identifies as outside the previously agreed formalization scope. The reviewed native
code formalizes algebraic encoding, injectivity, linearity, and distance, but
does not specify a bit-operation cost semantics or prove a linear bound for
an executable encoder. It also does not formalize the paper's reverse-time
transposed-encoder recurrence. Consequently, completing the distance theorem
and combining it with `realizedCode_rate` certifies the distance and rate
clauses of the paper theorem; it does not yet certify its separate ordinary
and transposed `O(N)` bit-work clause. The all-length wrapper after the native
theorem is likewise a separate scope from this native review.

To reproduce the concrete map comparison and the reviewed source identities,
run `C:/Python314/python.exe scripts/encoder_native_semantic_review.py` from
this Lean directory. It checks all 38 basis words as integers and writes
`scripts/map_data/encoder_native_semantic_snapshot.json`, containing SHA-256
hashes of the reviewed paper and Lean sources. It does not modify them.

The recommendation at this checkpoint was to finish the fixed/dense
substitutions, audit an unconditional theorem naming `minimumDistance` and
the full product seed law, and present its rate corollary. Those steps are
now complete; the current entry points are in [NATIVE_RESULT.md](NATIVE_RESULT.md).
A formal encoder-work proof remains a separate follow-up objective.
