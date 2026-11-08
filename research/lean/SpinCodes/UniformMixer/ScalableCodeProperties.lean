import SpinCodes.UniformMixer.ScalableRealization

/-!
# Vector realization of the sampled Scalable inner

The state and output recurrences are exactly those used by
`sampledOutputWeight`. The output retains chronological 64-bit blocks.
Injectivity is causal and requires no probabilistic assumption or selected
outer. The old IMT encoder is not imported as an encoder conclusion.
-/
noncomputable section
namespace Spin.UniformMixer

def sampledOutput : {R : ℕ} → (Fin R → Input64) →
    (Fin R → StateMap) → State22 → (Fin R → Input64)
  | 0, _, _, _ => Fin.elim0
  | _R + 1, xs, gs, q =>
      Fin.cons (fun j => coordMap j q + xs 0 j)
        (sampledOutput (Fin.tail xs) (Fin.tail gs)
          ((gs 0) q + stepFeedback (xs 0)))

def sampledWordWeight {R : ℕ} (ys : Fin R → Input64) : ℕ :=
  ∑ i, ∑ j, if ys i j = 1 then 1 else 0

theorem sampledOutput_weight {R : ℕ} (xs : Fin R → Input64)
    (gs : Fin R → StateMap) (q : State22) :
    sampledWordWeight (sampledOutput xs gs q) = sampledOutputWeight xs gs q := by
  induction R generalizing q with
  | zero => simp [sampledWordWeight, sampledOutputWeight]
  | succ R ih =>
    rw [sampledWordWeight, Fin.sum_univ_succ]
    change stepWeight (xs 0) q + sampledWordWeight
      (sampledOutput (Fin.tail xs) (Fin.tail gs)
        ((gs 0) q + stepFeedback (xs 0))) = _
    rw [ih]
    rfl

theorem stepFeedback_add (x y : Input64) :
    stepFeedback (x + y) = stepFeedback x + stepFeedback y := by
  funext i
  simp [stepFeedback, Pi.add_apply, add_mul, Finset.sum_add_distrib]

theorem sampledOutput_add {R : ℕ} (xs ys : Fin R → Input64)
    (gs : Fin R → StateMap) (q r : State22) :
    sampledOutput (xs + ys) gs (q + r) =
      sampledOutput xs gs q + sampledOutput ys gs r := by
  induction R generalizing q r with
  | zero => funext i; exact Fin.elim0 i
  | succ R ih =>
    funext i
    refine Fin.cases ?_ (fun k => ?_) i
    · funext j
      simp [sampledOutput, map_add, Pi.add_apply]
      abel
    · simp only [sampledOutput, Fin.cons_succ, Pi.add_apply]
      have hs : (gs 0) (q + r) + stepFeedback (xs 0 + ys 0) =
          ((gs 0) q + stepFeedback (xs 0)) +
            ((gs 0) r + stepFeedback (ys 0)) := by
        simp only [map_add, Pi.add_apply, stepFeedback_add]
        abel
      rw [hs]
      exact congrFun (ih (Fin.tail xs) (Fin.tail ys) (Fin.tail gs) _ _) k

theorem sampledOutput_injective {R : ℕ} (gs : Fin R → StateMap) (q : State22) :
    Function.Injective (fun xs : Fin R → Input64 => sampledOutput xs gs q) := by
  induction R generalizing q with
  | zero => intro xs ys _; funext i; exact Fin.elim0 i
  | succ R ih =>
    intro xs ys h
    have hhead : xs 0 = ys 0 := by
      funext j
      have hh := congrFun (congrFun h 0) j
      simp only [sampledOutput, Fin.cons_zero] at hh
      exact add_left_cancel hh
    have htail : Fin.tail xs = Fin.tail ys := by
      apply ih (Fin.tail gs) ((gs 0) q + stepFeedback (xs 0))
      have hh := congrArg Fin.tail h
      simpa only [sampledOutput, Fin.tail_cons, hhead] using hh
    funext i
    exact Fin.cases hhead (fun k => congrFun htail k) i

theorem sampledOutput_zero {R : ℕ} (gs : Fin R → StateMap) :
    sampledOutput 0 gs 0 = 0 := by
  have h := sampledOutput_add (0 : Fin R → Input64) 0 gs 0 0
  simp only [add_zero] at h
  exact (add_left_cancel (h.symm.trans (add_zero _).symm) : _)

def sampledEncoderLinearMap {R : ℕ} (gs : Fin R → StateMap) :
    (Fin R → Input64) →ₗ[ZMod 2] (Fin R → Input64) where
  toFun := fun xs => sampledOutput xs gs 0
  map_add' xs ys := by simpa only [add_zero] using sampledOutput_add xs ys gs 0 0
  map_smul' c xs := by
    have hc : c = 0 ∨ c = 1 := (by decide : ∀ c : ZMod 2, c = 0 ∨ c = 1) c
    rcases hc with rfl | rfl <;> simp [sampledOutput_zero]

#print axioms sampledOutput_weight
#print axioms sampledOutput_add
#print axioms sampledOutput_injective
#print axioms sampledEncoderLinearMap
end Spin.UniformMixer
