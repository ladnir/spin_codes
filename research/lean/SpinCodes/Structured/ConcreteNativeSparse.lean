import SpinCodes.Structured.ConcreteNativeSelection
import SpinCodes.Structured.ConcreteOuterCountingSelected

noncomputable section
attribute [local instance] Classical.propDecidable
namespace Spin.Structured.ConcreteNativeFamily
open Finset ConcreteOuter ConcreteRoute ConcreteBinaryEncoding

/-- The message bijection used by the actual binary native setup. -/
def tupleEquiv (m : ℕ) : (Fin (Nsched m / 2) → ZMod 2) ≃ NativeMessage m :=
  (wordEquiv (Nsched m / 2)).symm.trans (nativeMessageEquiv m)

@[simp] theorem tupleEquiv_zero (m : ℕ) : tupleEquiv m 0 = 0 := rfl

theorem tupleEquiv_ne_zero (m : ℕ) (x : Fin (Nsched m / 2) → ZMod 2) :
    tupleEquiv m x ≠ 0 ↔ x ≠ 0 := by
  rw [← tupleEquiv_zero m, (tupleEquiv m).injective.ne_iff]

theorem occ_eq_tuple (m : ℕ) (seed : NativeSeed m) (x : Fin (Nsched m / 2) → ZMod 2) :
    occ m seed x = occupation (tupleEquiv m x) := native_activeRows m seed _

/-- Reindexing an equal finite width does not change the routed experiment. -/
theorem failureProbability_width_cast {L b b' R : ℕ} (h : b = b')
    (e : (Fin b' × Fin L) ≃ (Fin R × Fin 128))
    (x : Fin L → Finset (Fin b)) (d : ℕ) :
    failureProbability e (fun i => (x i).map (finCongr h).toEmbedding) d =
      failureProbability ((Equiv.prodCongr (finCongr h) (Equiv.refl (Fin L))).trans e) x d := by
  subst b'
  simp
  rfl

def tupleWiring (m : ℕ) :
    (Fin (nativeBlocks m * 24) × Fin (Lsched m)) ≃ (Fin (rounds m) × Fin 128) :=
  (Equiv.prodCongr (finCongr (native_width m)) (Equiv.refl (Fin (Lsched m)))).trans
    (streamWiring (native_round_divisibility m))

theorem qd_eq_tuple_failure (m : ℕ) (seed : NativeSeed m)
    (x : Fin (Nsched m / 2) → ZMod 2) (d : ℕ) :
    (nativeSetup m).qd d ((nativeSetup m).outer seed x) =
      failureProbability (tupleWiring m) (rowSupports seed (tupleEquiv m x)) d := by
  have hcast := failureProbability_width_cast (native_width m)
    (streamWiring (native_round_divisibility m)) (rowSupports seed (tupleEquiv m x)) d
  trans failureProbability (streamWiring (native_round_divisibility m)) (rows m seed x) d
  · simp only [Setup.qd, nativeSetup, wordToRows_rowsToWord, failureProbability, rounds]
    exact prob_decidable _ _ _ _
  · exact hcast

theorem qd_sum_eq_tuple (m : ℕ) (seed : NativeSeed m) (Q d : ℕ) :
    (∑ x ∈ (nonzeroMsgs (Fin (Nsched m / 2) → ZMod 2)).filter (fun x => occ m seed x = Q),
      (nativeSetup m).qd d ((nativeSetup m).outer seed x)) =
    ∑ x ∈ (nonzeroMsgs (NativeMessage m)).filter (fun x => occupation x = Q),
      failureProbability (tupleWiring m) (rowSupports seed x) d := by
  simp only [nonzeroMsgs, Finset.sum_filter, qd_eq_tuple_failure, occ_eq_tuple]
  have h := (tupleEquiv m).sum_comp (fun x => if x ≠ 0 then
    (if occupation x = Q then failureProbability (tupleWiring m) (rowSupports seed x) d else 0) else 0)
  simpa only [tupleEquiv_ne_zero] using h

end Spin.Structured.ConcreteNativeFamily

