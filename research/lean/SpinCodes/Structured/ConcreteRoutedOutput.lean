import SpinCodes.Structured.ConcreteEncoderOutput
import SpinCodes.Structured.ConcreteRoutedMoment

noncomputable section
namespace Spin.Structured.ConcreteRoute
open Finset Routing
open scoped symmDiff

theorem shuffleSupport_injective {n : ℕ} (σ : Equiv.Perm (Fin n)) :
    Function.Injective (shuffleSupport σ) := Finset.map_injective _

theorem routeEval_injective {L b : ℕ} (seed : Seeds L b) : Function.Injective (fun rows => routeEval rows seed) := by
  intro rows rows' h
  have ht : transpose (fun i => shuffleSupport (seed.1 i) (rows i)) =
      transpose (fun i => shuffleSupport (seed.1 i) (rows' i)) := by
    funext j
    exact shuffleSupport_injective (seed.2 j) (congrFun h j)
  have hr := (transposeEquiv L b).injective ht
  funext i
  exact shuffleSupport_injective (seed.1 i) (congrFun hr i)

theorem routeEval_xor {L b : ℕ} (rows rows' : Fin L → Finset (Fin b)) (seed : Seeds L b) :
    routeEval (fun i => rows i ∆ rows' i) seed = fun j => routeEval rows seed j ∆ routeEval rows' seed j := by
  funext j
  ext i
  simp only [routeEval, mem_shuffleSupport, mem_transpose, mem_symmDiff]

theorem reshape_xor {a n c m : ℕ} (e : (Fin a × Fin n) ≃ (Fin c × Fin m))
    (x y : Fin a → Finset (Fin n)) :
    reshape e (fun i => x i ∆ y i) = fun j => reshape e x j ∆ reshape e y j := by
  funext j
  ext i
  simp only [mem_reshape, mem_symmDiff]

end Spin.Structured.ConcreteRoute

namespace Spin.Structured.ConcreteBinaryEncoding
open Finset
open scoped symmDiff

theorem bitEquiv_xor (a b : Bool) : bitEquiv (xor a b) = bitEquiv a + bitEquiv b := by
  cases a <;> cases b <;> decide

theorem rowsToWord_xor {L b : ℕ} (x y : Fin L → Finset (Fin b)) :
    rowsToWord (fun i => x i ∆ y i) = rowsToWord x + rowsToWord y := by
  funext k
  simp only [rowsToWord, Pi.add_apply]
  have h : decide ((finProdFinEquiv.symm k).2 ∈
      x (finProdFinEquiv.symm k).1 ∆ y (finProdFinEquiv.symm k).1) =
      xor (decide ((finProdFinEquiv.symm k).2 ∈ x (finProdFinEquiv.symm k).1))
          (decide ((finProdFinEquiv.symm k).2 ∈ y (finProdFinEquiv.symm k).1)) := by
    by_cases hx : (finProdFinEquiv.symm k).2 ∈ x (finProdFinEquiv.symm k).1 <;>
      by_cases hy : (finProdFinEquiv.symm k).2 ∈ y (finProdFinEquiv.symm k).1 <;>
      simp only [mem_symmDiff, hx, hy, not_true_eq_false, not_false_eq_true,
        true_and, false_and, or_self, or_false, false_or, decide_true, decide_false]
    all_goals rfl
  rw [h, bitEquiv_xor]

end Spin.Structured.ConcreteBinaryEncoding

namespace Spin.Structured.ConcreteRoutedEncoder
open ConcreteRoute ConcreteEncoder ConcreteBinaryEncoding
open scoped symmDiff

/-- The binary word emitted by the sampled route and IMT, in block-major order. -/
def word {L b R : ℕ} (e : (Fin b × Fin L) ≃ (Fin R × Fin 128))
    (rows : Fin L → Finset (Fin b)) (ω : Seeds L b × (Fin R → Transvection)) : Fin (R*128) → ZMod 2 :=
  outputWord (reshape e (routeEval rows ω.1)) ω.2 ∅

theorem word_weight {L b R : ℕ} (e : (Fin b × Fin L) ≃ (Fin R × Fin 128))
    (rows : Fin L → Finset (Fin b)) (ω : Seeds L b × (Fin R → Transvection)) :
    binaryWeight (word e rows ω) = weight e rows ω := outputWord_weight _ _ _

theorem word_injective {L b R : ℕ} (e : (Fin b × Fin L) ≃ (Fin R × Fin 128))
    (ω : Seeds L b × (Fin R → Transvection)) : Function.Injective (fun rows => word e rows ω) :=
  (outputWord_injective ω.2 ∅).comp ((reshapeEquiv e).injective.comp (routeEval_injective ω.1))

theorem word_xor {L b R : ℕ} (e : (Fin b × Fin L) ≃ (Fin R × Fin 128))
    (rows rows' : Fin L → Finset (Fin b)) (ω : Seeds L b × (Fin R → Transvection)) :
    word e (fun i => rows i ∆ rows' i) ω = word e rows ω + word e rows' ω := by
  unfold word outputWord
  rw [routeEval_xor, reshape_xor]
  have h := outputBlocks_xor (reshape e (routeEval rows ω.1))
    (reshape e (routeEval rows' ω.1)) ω.2 ∅ ∅
  rw [symmDiff_self] at h
  exact (congrArg rowsToWord h).trans (rowsToWord_xor _ _)

end Spin.Structured.ConcreteRoutedEncoder
