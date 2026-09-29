import SpinCodes.Structured.ConcreteRouteDomination
import SpinCodes.Structured.ConcreteShufflePermutation

noncomputable section
namespace Spin.Structured.ConcreteRoute

open Routing

abbrev Seeds (L b : ℕ) :=
  (Fin L → Equiv.Perm (Fin b)) × (Fin b → Equiv.Perm (Fin L))

/-- The route uses independent uniform row and region permutations. -/
def seedLaw (L b : ℕ) : FinPMF (Seeds L b) :=
  (piPMF (fun _ : Fin L => FinPMF.uniform (Equiv.Perm (Fin b)))).prod
    (piPMF (fun _ : Fin b => FinPMF.uniform (Equiv.Perm (Fin L))))

/-- The deterministic action of the sampled route permutations. -/
def routeEval {L b : ℕ} (rows : Fin L → Finset (Fin b)) (seed : Seeds L b) :
    Fin b → Finset (Fin L) :=
  fun j => shuffleSupport (seed.2 j)
    (transpose (fun i => shuffleSupport (seed.1 i) (rows i)) j)

theorem uniform_rows {L b : ℕ} (rows : Fin L → Finset (Fin b)) :
    (piPMF (fun _ : Fin L => FinPMF.uniform (Equiv.Perm (Fin b)))).map
      (fun σ i => shuffleSupport (σ i) (rows i)) = rowLaw rows := by
  rw [piPMF_map _ (fun i σ => shuffleSupport σ (rows i))]
  simp only [uniform_permutation_shuffleLaw, rowLaw]

theorem uniform_regions {L b : ℕ} (regions : Fin b → Finset (Fin L)) :
    (piPMF (fun _ : Fin b => FinPMF.uniform (Equiv.Perm (Fin L)))).map
      (fun σ j => shuffleSupport (σ j) (regions j)) = regionKernel regions := by
  rw [piPMF_map _ (fun j σ => shuffleSupport σ (regions j))]
  simp only [uniform_permutation_shuffleLaw, regionKernel]

/-- The algebraic route law is exactly the stated uniform-permutation experiment. -/
theorem permutation_law_eq {L b : ℕ} (rows : Fin L → Finset (Fin b)) :
    (seedLaw L b).map (routeEval rows) = law rows := by
  apply FinPMF.eq_of_expect_eq
  intro F
  rw [FinPMF.expect_map, seedLaw, FinPMF.expect_prod]
  rw [law, FinPMF.expect_bind, FinPMF.expect_map, ← uniform_rows rows,
    FinPMF.expect_map]
  apply congrArg
  funext σ
  rw [← uniform_regions, FinPMF.expect_map]
  rfl

/-- Pointwise domination for the explicit permutation experiment. -/
theorem permutation_dominates {L b : ℕ} (hL : 0 < L) (hb : 0 < b)
    (rows : Fin L → Finset (Fin b)) (hq0 : 0 < density rows) (hq1 : density rows < 1) :
    Dominates ((seedLaw L b).map (routeEval rows))
      (iidLaw L b (density rows) hq0.le hq1.le)
      (((b : ℝ) + 1) ^ activeRows rows * ((L : ℝ) + 1) ^ b) := by
  rw [permutation_law_eq]
  exact law_dominates hL hb rows hq0 hq1

/-- The paper's route inequality, with every random permutation explicit. -/
theorem permutation_expect_le {L b : ℕ} (hL : 0 < L) (hb : 0 < b)
    (rows : Fin L → Finset (Fin b)) (hq0 : 0 < density rows) (hq1 : density rows < 1)
    (F : (Fin b → Finset (Fin L)) → ℝ) (hF : ∀ x, 0 ≤ F x) :
    (seedLaw L b).expect (fun seed => F (routeEval rows seed)) ≤
      (((b : ℝ) + 1) ^ activeRows rows * ((L : ℝ) + 1) ^ b) *
        (iidLaw L b (density rows) hq0.le hq1.le).expect F := by
  rw [← FinPMF.expect_map, permutation_law_eq]
  exact expect_le hL hb rows hq0 hq1 F hF





/-- Every choice of permutations preserves the supplied total Hamming weight. -/
theorem routeEval_totalWeight {L b : ℕ} (rows : Fin L → Finset (Fin b))
    (seed : Seeds L b) : totalWeight (routeEval rows seed) = totalWeight rows := by
  simp only [totalWeight, routeEval, shuffleSupport_card]
  change totalWeight (transpose (fun i => shuffleSupport (seed.1 i) (rows i))) = _
  rw [totalWeight_transpose]
  simp only [totalWeight, shuffleSupport_card]

/-- The route inequality under natural-number weight conditions alone. -/
theorem permutation_expect_le_of_weight {L b : ℕ} (hL : 0 < L) (hb : 0 < b)
    (rows : Fin L → Finset (Fin b)) (hw0 : 0 < totalWeight rows)
    (hw1 : totalWeight rows < L * b)
    (F : (Fin b → Finset (Fin L)) → ℝ) (hF : ∀ x, 0 ≤ F x) :
    (seedLaw L b).expect (fun seed => F (routeEval rows seed)) ≤
      (((b : ℝ) + 1) ^ activeRows rows * ((L : ℝ) + 1) ^ b) *
        (iidLaw L b (density rows) (density_nonneg rows)
          (density_le_one hL hb rows)).expect F := by
  have hq0 : 0 < density rows := by
    unfold density
    exact div_pos (by exact_mod_cast hw0) (by positivity)
  have hq1 : density rows < 1 := by
    unfold density
    apply (div_lt_one (by positivity)).mpr
    exact_mod_cast hw1
  exact permutation_expect_le hL hb rows hq0 hq1 F hF

end Spin.Structured.ConcreteRoute


