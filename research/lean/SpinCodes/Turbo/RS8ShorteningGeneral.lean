import SpinCodes.Turbo.RS8Code

/-!
# General RS8 shortening

The four-column information-set equivalence reduces zero constraints on
the first `m ≤ 4` selected columns to `4-m` free columns. The selected
positions may be any distinct positions; no adjacency is assumed.
-/

noncomputable section

namespace Spin.Turbo

variable {F : Type*} [Field F] [Fintype F]

def rs8ZeroOnPrefix (m : ℕ) (hm : m ≤ 4)
    (points : Fin 8 → F) (positions : Fin 4 → Fin 8)
    (rows : RS8Group F) : Prop :=
  ∀ i : Fin m,
    rs8Column points rows (positions ⟨i.val, by omega⟩) = 0

noncomputable instance (m : ℕ) (hm : m ≤ 4)
    (points : Fin 8 → F) (positions : Fin 4 → Fin 8) :
    Fintype {rows : RS8Group F // rs8ZeroOnPrefix m hm points positions rows} :=
  Fintype.ofFinite _

private def extendFree (m : ℕ) (hm : m ≤ 4)
    (v : Fin (4-m) → Fin 8 → F) : Fin 4 → Fin 8 → F :=
  fun i => if h : i.val < m then 0 else v ⟨i.val-m, by omega⟩

/-- The `4-m` unforced columns parameterize the shortened RS8 group. -/
noncomputable def rs8ShorteningEquiv
    (m : ℕ) (hm : m ≤ 4)
    (points : Fin 8 → F) (hpoints : Function.Injective points)
    (positions : Fin 4 → Fin 8) (hpositions : Function.Injective positions) :
    {rows : RS8Group F // rs8ZeroOnPrefix m hm points positions rows} ≃
      (Fin (4-m) → Fin 8 → F) := by
  let e := rs8_four_eval_equiv points hpoints positions hpositions
  refine ⟨(fun x j => e x.val ⟨m+j.val, by omega⟩),
    (fun v => ⟨e.symm (extendFree m hm v), ?_⟩), ?_, ?_⟩
  · intro i
    have hz : e (e.symm (extendFree m hm v)) ⟨i.val, by omega⟩ = 0 := by
      rw [e.apply_symm_apply]
      have hi : i.val < m := i.isLt
      simp [extendFree, hi]
    funext r
    exact (rs8_four_eval_equiv_apply points hpoints positions
      hpositions _ _ r).symm.trans (congrFun hz r)
  · intro x
    apply Subtype.ext
    apply e.injective
    funext i
    by_cases hi : i.val < m
    · have hz : e x.val i = 0 := by
        let j : Fin m := ⟨i.val, hi⟩
        have hzero := x.property j
        have hieq : i = ⟨j.val, by omega⟩ := Fin.ext rfl
        rw [hieq]
        funext r
        exact (rs8_four_eval_equiv_apply points hpoints positions
          hpositions x.val _ r).trans (congrFun hzero r)
      change e (e.symm (extendFree m hm (fun j => e x.val ⟨m+j.val, by omega⟩))) i =
        e x.val i
      rw [e.apply_symm_apply]
      simp [extendFree, hi, hz]
    · let j : Fin (4-m) := ⟨i.val-m, by omega⟩
      have hij : (⟨m+j.val, by omega⟩ : Fin 4) = i := Fin.ext (by simp [j]; omega)
      change e (e.symm (extendFree m hm (fun j => e x.val ⟨m+j.val, by omega⟩))) i =
        e x.val i
      rw [e.apply_symm_apply]
      simp [extendFree, ← hij]
  · intro v
    funext j
    change e (e.symm (extendFree m hm v)) ⟨m+j.val, by omega⟩ = v j
    rw [e.apply_symm_apply]
    have hj : ¬ m + j.val < m := by omega
    simp [extendFree, hj]

/-- Zeroing `m ≤ 4` selected columns leaves `|F|^(8*(4-m))` groups. -/
theorem rs8Shortening_card [DecidableEq F]
    (m : ℕ) (hm : m ≤ 4)
    (points : Fin 8 → F) (hpoints : Function.Injective points)
    (positions : Fin 4 → Fin 8) (hpositions : Function.Injective positions) :
    Fintype.card {rows : RS8Group F //
      rs8ZeroOnPrefix m hm points positions rows} =
      (Fintype.card F) ^ (8 * (4-m)) := by
  classical
  rw [Fintype.card_congr (rs8ShorteningEquiv m hm points hpoints positions hpositions)]
  simp [← pow_mul, mul_comm]

#print axioms rs8ShorteningEquiv
#print axioms rs8Shortening_card

end Spin.Turbo
