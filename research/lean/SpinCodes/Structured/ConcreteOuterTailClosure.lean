import SpinCodes.Structured.ConcreteOuterTailBound
import SpinCodes.Structured.ConcreteOuterTailDenseLimit
import SpinCodes.Structured.ConcreteOuterDefectEnvelope

noncomputable section
attribute [local instance] Classical.propDecidable
namespace Spin.Structured.ConcreteOuter
open Filter Finset

lemma outerTail_nonneg (k : ℕ) : 0 ≤ outerTail k := by
  apply sum_nonneg
  intro w hw
  exact FinPMF.expect_nonneg _ (fun seed => by positivity)

/-- The actual finite Golay-BAA outer has vanishing expected count outside the
paper's retained weight window. Sparse and dense messages are both included. -/
theorem outerTail_tendsto : Tendsto outerTail atTop (nhds 0) := by
  have h := sparse_outer_tails_tendsto.add (dense_nativeBlocks_exp_tendsto 5)
  simp only [add_zero] at h
  apply squeeze_zero' (Eventually.of_forall outerTail_nonneg) _ h
  filter_upwards [eventually_ge_atTop 1] with k hk
  exact outerTail_le hk

end Spin.Structured.ConcreteOuter

namespace Spin.Structured.ConcreteNativeFamily
open Filter Finset ConcreteOuter

lemma nativeBlocks_tendsto : Tendsto nativeBlocks atTop atTop := by
  apply tendsto_atTop.2
  intro n
  filter_upwards [bsched_tendsto.eventually (eventually_ge_atTop (n*24))] with m hm
  have hw := native_width m
  omega

/-- Unconditional outer-tail closure for the concrete native family. -/
theorem nativeTail_tendsto : Tendsto nativeTail atTop (nhds 0) := by
  have h := outerTail_tendsto.comp nativeBlocks_tendsto
  apply h.congr
  intro m
  simp only [Function.comp_def,outerTail,nativeTail,Abar,nativeSeedLaw,native_width]

/-- The actual one-sample selection event has positive probability eventually. -/
theorem native_good_positive_eventually :
    ∀ᶠ m in atTop, 0 < (nativeSeedLaw m).prob (nativeGood m) :=
  raw_selection_positive_eventually nativeTail_tendsto

/-- The outer good event fails with vanishing probability, without analytic premises. -/
theorem native_selection_failure_tendsto :
    Tendsto concreteFamily.probNotGood atTop (nhds 0) :=
  concrete_selection_failure_tendsto nativeTail_tendsto

/-- The actual native family's entire growing sparse range is now unconditional. -/
theorem native_sparse_eventually :
    ∀ᶠ m in atTop, ∀ Q∈Ico 4096 (nativeCut m+1),
      concreteFamily.EZ m Q ≤ Real.exp (-(0.006:ℝ)*(Q:ℝ)*bsched m) :=
  concrete_sparse_eventually_of_tail nativeTail_tendsto

end Spin.Structured.ConcreteNativeFamily

