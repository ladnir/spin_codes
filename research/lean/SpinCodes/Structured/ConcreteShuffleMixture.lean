import SpinCodes.Structured.ConcreteShuffle
import SpinCodes.FiniteLaw

namespace Spin.Structured.Routing

noncomputable def shuffledLaw {n : ℕ} (P : FinPMF (Finset (Fin n))) :=
  P.bind shuffleLaw

lemma shuffledLaw_fiber_uniform {n : ℕ} (P : FinPMF (Finset (Fin n)))
    (T U : Finset (Fin n)) (h : T.card = U.card) :
    (shuffledLaw P).p T = (shuffledLaw P).p U := by
  simp only [shuffledLaw, FinPMF.bind_p, shuffleLaw_apply, h]

lemma shuffledLaw_prob_card {n : ℕ} (P : FinPMF (Finset (Fin n))) (k : ℕ) :
    (shuffledLaw P).prob (fun T => T.card = k) = P.prob (fun T => T.card = k) := by
  classical
  rw [shuffledLaw, FinPMF.bind_prob]
  simp only [shuffleLaw_prob_card]
  exact (FinPMF.prob_eq_expect_indicator P _).symm

end Spin.Structured.Routing
