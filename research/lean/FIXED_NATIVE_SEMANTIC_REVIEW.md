# Fixed occupations in the native family

> **Historical snapshot — archived 2026-09-28.** Status statements and task recommendations below describe an earlier stage.
> See [STATUS.md](STATUS.md) and [NATIVE_RESULT.md](NATIVE_RESULT.md) for the closed native distance/rate result and verification scope.

Read-only semantic review, 2026-09-28. No blocking mismatch was found between the fixed-occupation argument and the actual native first moment. This review concerns statement correspondence and quantifiers; it neither recompiles the dependencies nor certifies completion of the remaining dense replay and final assembly.

The endpoint is `ConcreteNativeFamily.native_fixed_EZ_tendsto` in `ConcreteNativeFixedAll.lean`: for each fixed natural number `Q > 0`, `concreteFamily.EZ m Q` tends to zero as `m` tends to infinity. Its specialization `native_fixed_range_tendsto` covers exactly `Q ∈ Ico 1 4096`. Both statements have no spectrum, continuum-norm, placement-law, or route-identification premise.

## The event, probability space, and native parameters

`Family.EZ` in `Instantiation.lean` is the expected number of nonzero messages of occupation `Q` whose encoded weight is at most the family's threshold. The outer seed is drawn once from its law conditioned on the selection event; every row uses that same seed. The inner law supplies the actual routing permutations and encoder randomness. This is the first moment used by the final union bound, not a surrogate with independently resampled outer constituents.

`ConcreteNativeTotal.selectedGood` uses `nativeGood` whenever that event has positive probability, with a whole-space fallback otherwise. `native_good_positive_eventually` proves eventual positivity for the actual outer law. Consequently the fallback disappears before the asymptotic bounds are applied. `ConcreteNativeDenseIntegration.EZ_le_tuple_layer_bound` averages a pointwise tuple bound valid for every selected shared seed. Conditional expectation of this constant bound is the same bound; an additional inverse selection probability is neither needed nor silently discarded. The probability of outer selection failure is handled separately in the final bad-event bound.

The route transfer uses `R = m + 1`, region length `L = 128 R = Lsched m`, row length `b = bsched m`, and `b L / 128` encoder rounds. `tupleWiring_eq_regionMajor` and `native_fixed_profile_routed` supply the exact serialization and dimension equalities. The threshold remains `threshold m = ⌊0.11 Nsched m⌋`; replacing it by `0.11 L b` only enlarges an upper bound. For `z = exp(-θ/L)`, the low-weight bound divides the moment by `z^d`. Since `0 < z ≤ 1`, the direction is correct: weight at most `d` implies `z^weight ≥ z^d`.

## One and two active rows

`ConcreteFixedFairProduct.fairRegionKernel_eq_fugacity` identifies the actual fair marked-region kernel with the fugacity kernel at all fugacities equal to one, divided by `2^Q`. Its product and serialization identities connect this to the actual routed input stream. The continuum norm transfers with a fixed positive margin, uniformly over every number `b` of regions once `R` is large enough.

`actual_one_mark_moment` and `actual_two_mark_moment` give respectively

\[
1000(1011/2000)^b,\qquad 1000(1011/4000)^b,
\]

at tilts `θ = 2` and `θ = 4`. The routed one- and two-row modules apply these bounds to the same fair experiment used by `occupation_failure_le_fair`. That counting comparison retains the active-position factor `choose(L,Q)` and the realized spectrum factor `(2^b B)^Q` for a fixed shared seed. The native bounds use `choose(L,1)=L` and `choose(L,2)≤L²`; these costs are included in the subsequent exponents.

`native_one_EZ_tendsto` and `native_two_EZ_tendsto` follow from eventual bounds `1000 exp(-b/50)` and `1000 exp(-b/20)`. These are deliberately weaker than the paper's intermediate numerical exponents: the appendix uses `θ = Q log(2)/p₀` for these two cases. The checked Lean argument proves the required limits with different tilts and conservative constants; it does not establish those sharper intermediate estimates.

## At least three active rows

The placement distribution is derived from the actual shuffle. `shuffledRegionKernel` is the shuffle-law expectation of the finite encoder endpoint kernel. `shuffledRegionKernel_uniform_limit` compares it with the lifted continuum kernel uniformly over supports of fixed cardinality and endpoint states. The continuum kernel includes the factorial normalization of the ordered simplex.

Two normalization checks are material here. `ConcreteFixedSubsetMarginal.picked_subset_average` starts with a uniform fixed-size subset and uniform label permutation and proves the exact smaller-subset marginal, including its binomial and factorial denominators. `ConcreteFixedInsertionLimit.fugacityContinuum_eq_insertion` then identifies the fugacity continuum with the insertion integral. The ordered-site factorial cancels the permutation-average factorial; no unaccounted factor `Q!` remains.

`ConcreteFixedOrderStatistic.orderedSite_orderStatistic_integral` proves the ordered-coordinate density directly. Gap symmetry and the endpoint cases feed `ConcreteFixedLargeLaplaceActual.liveDuration_laplace_le`, which supplies the actual spacing bound without a marginal-law premise. The enlarged positive insertion matrix is used through an entrywise upper bound, not treated as a probability transition matrix.

`ConcreteFixedLargeNorm.fugacityContinuum_paper_rowNorm` yields the actual weighted norm with

\[
\sigma=127/250,\quad v=3/1600,\quad r=1/524287,
\quad \tau=133/125,\quad p_0=262144/524287,
\quad \theta=Q\tau/p_0.
\]

Its cost is `exp(Q * paperLargeExponent) * ∏ᵢ mv(σ,v,r,uᵢ)`. The scalar exponent uses the proved `Q≥3` comparison with `G₃`. `fixedLargeNorm_proved` discharges the formerly explicit norm premise before `native_fixed_large_tendsto` is exposed.

`ConcreteFixedProfileRouted.actual_profile_routed_probability` bounds the actual experiment by the transferred kernel cost divided by

\[
\min(1,v)\;\prod_i {b\choose w_i}\;\prod_i u_i^{w_i}\;\exp(-\theta/L)^d.
\]

The row-law/profile identities justify the binomial denominators; nonnegative coefficient extraction justifies the fugacity denominators. Row embedding inserts zeros outside the selected active positions and is identified with the actual route. No independence of the shared outer seed is introduced.

The native spectrum bound and the finite rational fugacity family control every selected row density in the closed interval `[0.104,0.896]`. For fixed `Q`, the type of fugacity tuples is finite. `fixed_profile_uniform` takes an eventual bound for all such tuples before choosing a tuple from the row profile. Thus the eventual threshold does not depend on the later message, profile, or active-position set.

The 62 endpoint checks cover the actual refined majorant and give a gap `0.0008679`. The implementation reserves `0.0001` for the continuum transfer, leaving `0.0007679`. The profile sum retains `(b+1)^Q`, the selection cost retains `b^(2Q)`, and the position sum retains `choose(L,Q)≤L^Q`. `ConcreteNativeFixedLargeLimit.fixedLargeBound` explicitly includes the resulting logarithmic errors and prefactor `1600/3`. The native schedule and vanishing outer-envelope remainder give the eventual bound

\[
\mathrm{EZ}_m(Q)\le(1600/3)\exp(-Qb/2000).
\]

This proves the required fixed-`Q` limit. It does not reproduce the paper's sharper `-0.0008679 Qb + o_Q(b)` intermediate exponent. Using a fixed norm margin also avoids needing an error that remains `exp(o(1))` after multiplying an increasing number of region kernels.

## Summation and scope

`Regimes.fixed_regime` applies `tendsto_finsetSum` to the fixed set `Ico 1 4096`, containing 4095 occupations. Each occupation may have its own eventual threshold. This is sufficient for their finite sum to tend to zero; no uniform statement for a growing `Q(m)` is inferred. The later regime partition starts the sparse range at 4096 and uses consecutive half-open intervals, so the fixed cutoff is neither omitted nor counted twice.

The reviewed chain matches the fixed-occupation clauses in `paper/structured_imt_appendix.tex` (the fixed-occupation continuum argument) and `paper/structured_proof.tex` (the final regime sum), with the conservative constants described above. The then-pending global dense certificates, the full minimum-distance assembly, and construction complexity were outside this review. The recommended completion check was to combine the unconditional fixed theorem with the independently verified sparse and dense bounds, without strengthening its fixed-occupation quantifiers. The native distance/rate assembly is now complete; see [NATIVE_RESULT.md](NATIVE_RESULT.md).
