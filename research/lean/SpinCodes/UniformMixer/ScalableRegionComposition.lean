import SpinCodes.UniformMixer.ScalableFullMoment

/-!
# Serial region boundaries in the actual encoder

The flattened step stream is the concatenation of the regions in their
encoder order. Hence the exact full-state kernel composes across region
boundaries, with no state reset.
-/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteOuter

theorem scalableFlatInput_stream (m : ℕ)
    (regions : Fin (bsched m) → Finset (Fin (Lsched m))) :
    List.ofFn (scalableFlatInput m regions) =
      (List.ofFn (fun j : Fin (bsched m) =>
        List.ofFn (scalableRegionInput m (regions j)))).flatten := by
  rw [List.ofFn_mul]
  apply congrArg List.flatten
  apply congrArg List.ofFn
  funext j
  apply congrArg List.ofFn
  funext k
  have heq : (⟨j.val * scalableStepsPerRegion m + k.val,
      by
        have hfirst : j.val * scalableStepsPerRegion m + k.val <
            (j.val + 1) * scalableStepsPerRegion m := by
          rw [add_mul, one_mul]
          omega
        exact hfirst.trans_le
          (Nat.mul_le_mul_right _ (Nat.succ_le_of_lt j.isLt))⟩ :
      Fin (bsched m * scalableStepsPerRegion m)) =
      finProdFinEquiv (j, k) := by
    ext
    ac_rfl
  rw [heq]
  simp [scalableFlatInput]

def kernelListProduct : List (State22 → State22 → ℝ) →
    State22 → State22 → ℝ
  | [] => identity
  | K :: Ks => compose K (kernelListProduct Ks)

theorem inputListKernel_flatten {b : ℕ} (z : ℝ)
    (blocks : Fin b → List Input64) :
    inputListKernel z (List.ofFn blocks).flatten =
      kernelListProduct (List.ofFn (fun j => inputListKernel z (blocks j))) := by
  induction b with
  | zero => rfl
  | succ b ih =>
      simp only [List.ofFn_succ, List.flatten_cons, kernelListProduct]
      rw [inputListKernel_append]
      congr 1
      exact ih (Fin.tail blocks)

theorem scalable_kernel_region_product (m : ℕ) (z : ℝ)
    (regions : Fin (bsched m) → Finset (Fin (Lsched m))) :
    inputListKernel z (List.ofFn (scalableFlatInput m regions)) =
      kernelListProduct (List.ofFn (fun j : Fin (bsched m) =>
        chronologicalRegionKernel z
          (scalableRegionSupport m (regions j)))) := by
  rw [scalableFlatInput_stream, inputListKernel_flatten]
  rfl

end Spin.UniformMixer
