import SpinCodes.Structured.ConcreteRouteDomination
import SpinCodes.FiniteProductLaw

noncomputable section
namespace Spin.Structured.ConcreteRoute
open Finset

/-- Rewire a bit matrix along a bijection of its positions. -/
def reshape {a n c m : ℕ} (e : (Fin a × Fin n) ≃ (Fin c × Fin m))
    (x : Fin a → Finset (Fin n)) : Fin c → Finset (Fin m) :=
  fun j => univ.filter (fun k => (e.symm (j, k)).2 ∈ x (e.symm (j, k)).1)

@[simp] theorem mem_reshape {a n c m : ℕ}
    (e : (Fin a × Fin n) ≃ (Fin c × Fin m))
    (x : Fin a → Finset (Fin n)) (j : Fin c) (k : Fin m) :
    k ∈ reshape e x j ↔ (e.symm (j, k)).2 ∈ x (e.symm (j, k)).1 := by
  simp [reshape]

@[simp] theorem reshape_inverse {a n c m : ℕ}
    (e : (Fin a × Fin n) ≃ (Fin c × Fin m)) (x : Fin a → Finset (Fin n)) :
    reshape e.symm (reshape e x) = x := by
  funext j
  ext k
  simp only [mem_reshape, Equiv.symm_symm]
  simpa using Iff.rfl

def reshapeEquiv {a n c m : ℕ} (e : (Fin a × Fin n) ≃ (Fin c × Fin m)) :
    (Fin a → Finset (Fin n)) ≃ (Fin c → Finset (Fin m)) where
  toFun := reshape e
  invFun := reshape e.symm
  left_inv := reshape_inverse e
  right_inv := reshape_inverse e.symm

/-- Independent identically distributed bits are invariant under any fixed rewiring. -/
theorem reshape_iid {a n c m : ℕ} (e : (Fin a × Fin n) ≃ (Fin c × Fin m))
    (q : ℝ) (hq0 : 0 ≤ q) (hq1 : q ≤ 1) :
    (iidLaw n a q hq0 hq1).map (reshape e) = iidLaw m c q hq0 hq1 := by
  ext y
  change ((iidLaw n a q hq0 hq1).map (reshapeEquiv e)).p y = _
  rw [FinPMF.map_equiv_apply]
  change (iidLaw n a q hq0 hq1).p (reshape e.symm y) = _
  simp only [iidLaw, piPMF_apply, poissonBinom_eq_prod, mem_reshape,
    Equiv.symm_symm]
  have h := Fintype.prod_equiv e
    (fun x => if (e x).2 ∈ y (e x).1 then q else 1 - q)
    (fun x => if x.2 ∈ y x.1 then q else 1 - q) (fun _ => rfl)
  simpa only [Fintype.prod_prod_type] using h

end Spin.Structured.ConcreteRoute
