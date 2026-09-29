import SpinCodes.Structured.ConcreteRouteTranspose
import SpinCodes.FiniteProductLaw

noncomputable section
attribute [local instance] Classical.propDecidable
namespace Spin.Structured.ConcreteMarked
open Finset

def coin (p : ℝ) (hp0 : 0 ≤ p) (hp1 : p ≤ 1) : FinPMF Bool where
  p := fun b => if b then p else 1 - p
  nonneg := fun b => by cases b <;> simp <;> linarith
  total := by simp

def support {n : ℕ} (x : Fin n → Bool) : Finset (Fin n) :=
  univ.filter (fun i => x i = true)

@[simp] theorem mem_support {n : ℕ} (x : Fin n → Bool) (i : Fin n) :
    i ∈ support x ↔ x i = true := by simp [support]

def supportEquiv (n : ℕ) : (Fin n → Bool) ≃ Finset (Fin n) where
  toFun := support
  invFun := fun S i => decide (i ∈ S)
  left_inv := fun x => by funext i; simp [support]
  right_inv := fun S => by ext i; simp [support]

def iidBits (n : ℕ) (p : ℝ) (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :=
  poissonBinom (fun _ : Fin n => p) (fun _ => hp0) (fun _ => hp1)

theorem support_coin {n : ℕ} (p : ℝ) (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    (piPMF (fun _ : Fin n => coin p hp0 hp1)).map support = iidBits n p hp0 hp1 := by
  ext S
  change ((piPMF (fun _ : Fin n => coin p hp0 hp1)).map (supportEquiv n)).p S = _
  rw [FinPMF.map_equiv_apply]
  simp only [piPMF_apply, coin, supportEquiv, iidBits, ConcreteRoute.poissonBinom_eq_prod,
    decide_eq_true_eq]
  change (∏ i : Fin n, if decide (i ∈ S) = true then p else 1-p) = _
  simp

/-- Independent fair bits are supplied at all positions, then masked by the marks. -/
def fairOnMarks {n : ℕ} (S : Finset (Fin n)) : FinPMF (Finset (Fin n)) :=
  (piPMF (fun _ : Fin n => coin (1/2) (by norm_num) (by norm_num))).map
    (fun bits => support (fun i => decide (i ∈ S) && bits i))

theorem fairOnMarks_support {n : ℕ} (marks : Fin n → Bool) :
    fairOnMarks (support marks) =
      (piPMF (fun i => (coin (1/2) (by norm_num) (by norm_num)).map
        (fun bit => marks i && bit))).map support := by
  rw [← piPMF_map, FinPMF.map_comp]
  unfold fairOnMarks
  congr 1
  funext bits
  congr 1
  funext i
  simp

theorem coin_thinning (p : ℝ) (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    (coin p hp0 hp1).bind (fun mark =>
      (coin (1/2) (by norm_num) (by norm_num)).map (fun bit => mark && bit)) =
        coin (p/2) (by positivity) (by linarith) := by
  ext b
  cases b <;> simp [FinPMF.bind_p, FinPMF.map_p, coin] <;> ring

/-- Independently retaining a fair bit at each Bernoulli-p mark gives iid Bernoulli-p/2. -/
theorem iid_thinning (n : ℕ) (p : ℝ) (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    (iidBits n p hp0 hp1).bind fairOnMarks =
      iidBits n (p/2) (by positivity) (by linarith) := by
  rw [← support_coin p hp0 hp1, FinPMF.bind_map]
  simp only [fairOnMarks_support]
  rw [← FinPMF.map_bind]
  rw [piPMF_bind (fun _ : Fin n => coin p hp0 hp1)
    (fun (_ : Fin n) mark => (coin (1/2) (by norm_num) (by norm_num)).map
      (fun bit => mark && bit))]
  simp only [coin_thinning]
  exact support_coin _ _ _

end Spin.Structured.ConcreteMarked


