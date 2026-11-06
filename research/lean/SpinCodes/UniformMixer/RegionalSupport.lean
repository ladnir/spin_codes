import SpinCodes.UniformMixer.RegionalPlacement
import Mathlib.Data.Finset.Sort

/-!
# From a shuffled support to separated marked steps

We sort the active bit positions in one region. Each position determines
its 64-bit step and its coordinate within that step. The preceding empty
gap and final empty gap are then determined by subtraction. The good
placement event makes every such gap nonempty.
-/

noncomputable section

namespace Spin.UniformMixer

/-- Empty gap and local coordinate for each selected bit, in bit order. -/
def supportMarks (R : ℕ) :
    List (Fin (64 * R)) → ℕ → List (ℕ × Fin 64)
  | [], _ => []
  | x :: xs, prev =>
      (x.val / 64 - prev, ⟨x.val % 64, Nat.mod_lt _ (by omega)⟩) ::
        supportMarks R xs (x.val / 64 + 1)

/-- Number of empty steps after the last selected bit. -/
def supportLastGap (R : ℕ) : List (Fin (64 * R)) → ℕ → ℕ
  | [], prev => R - prev
  | x :: xs, _ => supportLastGap R xs (x.val / 64 + 1)

theorem supportMarks_length (R : ℕ) (xs : List (Fin (64 * R)))
    (prev : ℕ) : (supportMarks R xs prev).length = xs.length := by
  induction xs generalizing prev with
  | nil => rfl
  | cons x tail ih => simp [supportMarks, ih]

/-- The ordered support has a nonempty empty-step gap before each bit and
after its last bit. -/
inductive ValidSupport (R : ℕ) : ℕ → List (Fin (64 * R)) → Prop
  | nil {prev} (hprev : prev < R) : ValidSupport R prev []
  | cons {prev x xs}
      (hgap : prev < x.val / 64)
      (hinterior : x.val / 64 + 1 < R)
      (htail : ValidSupport R (x.val / 64 + 1) xs) :
      ValidSupport R prev (x :: xs)

theorem validSupport_mark_gaps {R prev : ℕ}
    {xs : List (Fin (64 * R))} (h : ValidSupport R prev xs) :
    (∀ p ∈ supportMarks R xs prev, 0 < p.1 ∧ p.1 ≤ R) ∧
      0 < supportLastGap R xs prev ∧ supportLastGap R xs prev ≤ R := by
  induction h with
  | nil hprev =>
      simp [supportMarks, supportLastGap]
      omega
  | cons hgap hinterior htail ih =>
      rcases ih with ⟨hmarks, hlast⟩
      constructor
      · intro p hp
        simp only [supportMarks, List.mem_cons] at hp
        rcases hp with rfl | hp
        · constructor <;> omega
        · exact hmarks p hp
      · exact hlast

theorem validSupport_of_pairwise {R prev : ℕ}
    (xs : List (Fin (64 * R)))
    (hprevR : prev < R)
    (hpair : xs.Pairwise (· < ·))
    (hint : ∀ x ∈ xs,
      1 ≤ x.val / 64 ∧ x.val / 64 + 1 < R)
    (hsep : ∀ x ∈ xs, ∀ y ∈ xs, x ≠ y →
      2 ≤ Nat.dist (x.val / 64) (y.val / 64))
    (hprev : ∀ x ∈ xs, prev < x.val / 64) :
    ValidSupport R prev xs := by
  induction xs generalizing prev with
  | nil => exact ValidSupport.nil hprevR
  | cons x tail ih =>
      obtain ⟨hhead, htailPair⟩ := List.pairwise_cons.mp hpair
      have hxmem : x ∈ x :: tail := by simp
      have hxgap := hprev x hxmem
      have hxint := (hint x hxmem).2
      apply ValidSupport.cons hxgap hxint
      apply ih hxint htailPair
      · intro y hy
        exact hint y (by simp [hy])
      · intro y hy z hz hyz
        exact hsep y (by simp [hy]) z (by simp [hz]) hyz
      · intro y hy
        have hxy : x < y := hhead y hy
        have hstep : x.val / 64 ≤ y.val / 64 := by
          have hval : x.val < y.val := hxy
          omega
        have hd := hsep x hxmem y (by simp [hy]) (ne_of_lt hxy)
        unfold Nat.dist at hd
        omega

/-- Every support outside the discarded event has the gap structure used
by `markedRegionKernel`. -/
theorem goodRegionalSites_validSupport {R : ℕ} (hR : 0 < R)
    (T : Finset (Fin (64 * R))) (hgood : ¬ badRegionalSites T) :
    ValidSupport R 0 (T.sort (· ≤ ·)) := by
  have hpair : (T.sort (· ≤ ·)).Pairwise (· < ·) :=
    T.sortedLT_sort.pairwise
  have hint (x : Fin (64 * R)) (hx : x ∈ T.sort (· ≤ ·)) :
      1 ≤ x.val / 64 ∧ x.val / 64 + 1 < R :=
    goodRegionalSites_interior T hgood x ((Finset.mem_sort (· ≤ ·)).mp hx)
  have hsep (x : Fin (64 * R)) (hx : x ∈ T.sort (· ≤ ·))
      (y : Fin (64 * R)) (hy : y ∈ T.sort (· ≤ ·)) (hxy : x ≠ y) :
      2 ≤ Nat.dist (x.val / 64) (y.val / 64) :=
    goodRegionalSites_step_separation T hgood x y
      ((Finset.mem_sort (· ≤ ·)).mp hx)
      ((Finset.mem_sort (· ≤ ·)).mp hy) hxy
  exact validSupport_of_pairwise (T.sort (· ≤ ·)) hR hpair hint hsep
    (fun x hx => Nat.zero_lt_of_lt (hint x hx).1)

theorem goodRegionalSites_mark_gaps {R : ℕ} (hR : 0 < R)
    (T : Finset (Fin (64 * R))) (hgood : ¬ badRegionalSites T) :
    (∀ p ∈ supportMarks R (T.sort (· ≤ ·)) 0,
      0 < p.1 ∧ p.1 ≤ R) ∧
    0 < supportLastGap R (T.sort (· ≤ ·)) 0 ∧
      supportLastGap R (T.sort (· ≤ ·)) 0 ≤ R :=
  validSupport_mark_gaps (goodRegionalSites_validSupport hR T hgood)

/-- The kernel indexed by a support after sorting its active positions.
The equality to the stepwise encoder kernel is proved separately. -/
def supportRegionKernel {R : ℕ} (z : ℝ)
    (T : Finset (Fin (64 * R))) : State22 → State22 → ℝ :=
  markedRegionKernel z (supportMarks R (T.sort (· ≤ ·)) 0)
    (supportLastGap R (T.sort (· ≤ ·)) 0)

theorem good_supportRegion_coarse_upper {R : ℕ} (hR : 0 < R)
    (T : Finset (Fin (64 * R))) (hgood : ¬ badRegionalSites T)
    (θ : ℝ) (hθ : 0 ≤ θ) (q r : State22) :
    let L : ℝ := 64 * R
    let C := ((64 * θ + 64 * θ ^ 2) / L) / Real.exp (-θ)
    supportRegionKernel (Real.exp (-(θ / L))) T q r ≤
      (1 + C) ^ (T.card + 1) *
        (liftM * coarseRegionProduct
          (fun g => Real.exp (-((θ / L) * liveWeightAverage coordMap)) ^ g)
          (supportMarks R (T.sort (· ≤ ·)) 0)
          (supportLastGap R (T.sort (· ≤ ·)) 0) * projectM) q r := by
  dsimp [supportRegionKernel]
  have hgap := goodRegionalSites_mark_gaps hR T hgood
  have hL : (0 : ℝ) < (64 * R : ℕ) := by
    exact_mod_cast Nat.mul_pos (by omega : 0 < 64) hR
  have hcap (g : ℕ) (hg : g ≤ R) :
      (g : ℝ) ≤ ((64 * R : ℕ) : ℝ) / 64 := by
    norm_num [Nat.cast_mul]
    exact hg
  have hmarks : ∀ p ∈ supportMarks R (T.sort (· ≤ ·)) 0,
      0 < p.1 ∧ (p.1 : ℝ) ≤ ((64 * R : ℕ) : ℝ) / 64 := by
    intro p hp
    exact ⟨(hgap.1 p hp).1, hcap p.1 (hgap.1 p hp).2⟩
  have hlast : 0 < supportLastGap R (T.sort (· ≤ ·)) 0 ∧
      (supportLastGap R (T.sort (· ≤ ·)) 0 : ℝ) ≤
        ((64 * R : ℕ) : ℝ) / 64 :=
    ⟨hgap.2.1, hcap _ hgap.2.2⟩
  have hbound := production_markedRegion_coarse_upper θ
    ((64 * R : ℕ) : ℝ) hθ hL
    (supportMarks R (T.sort (· ≤ ·)) 0)
    (supportLastGap R (T.sort (· ≤ ·)) 0)
    hmarks hlast q r
  rw [supportMarks_length, Finset.length_sort] at hbound
  simpa only [Nat.cast_mul, Nat.cast_ofNat] using hbound

end Spin.UniformMixer
