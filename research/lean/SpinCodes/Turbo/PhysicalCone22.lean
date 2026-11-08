import SpinCodes.Turbo.PhysicalThinning22
import SpinCodes.UniformMixer.QuadraticSpectrumSupport

/-!
# Semantic cone for the Turbo 11-coordinate physical envelope

The retained numerical producer tracks mass at state zero, arbitrary
nonzero mass, uniformly dominated nonzero mass, and eight nonzero
expansion-weight classes. A cone element is an explicit decomposition
of a nonnegative state measure into those eleven parts. This file
defines the exact semantic interface before checking any numerical row.
-/

noncomputable section

namespace Spin.Turbo

open Finset Spin.UniformMixer

def physicalClassWeight22 (i : Fin 11) : ℕ :=
  if i.val = 3 then 16 else
  if i.val = 4 then 24 else
  if i.val = 5 then 28 else
  if i.val = 6 then 32 else
  if i.val = 7 then 36 else
  if i.val = 8 then 40 else
  if i.val = 9 then 48 else
  if i.val = 10 then 64 else 0

def physicalClassAllows22 (i : Fin 11) (q : State22) : Prop :=
  if i.val = 0 then q = 0
  else if i.val = 1 ∨ i.val = 2 then q ≠ 0
  else coordinateWeight coordMap q = physicalClassWeight22 i

/-- The eight retained positive weight classes cover every nonzero
production state. This connects the cone's class indices to the
checked production-map spectrum, rather than to a data-file label. -/
theorem physicalClass22_covers_nonzero (q : State22)
    (hq : q ≠ 0) :
    ∃ i : Fin 11, 3 ≤ i.val ∧ physicalClassAllows22 i q := by
  have hw := productionWeight_mem_spectrumWeights q
  simp only [spectrumWeights, Finset.mem_insert,
    Finset.mem_singleton] at hw
  rcases hw with h | h | h | h | h | h | h | h | h
  · exact (hq (productionWeight_zero_only q h)).elim
  · exact ⟨3, by norm_num,
      by simpa [physicalClassAllows22, physicalClassWeight22] using h⟩
  · exact ⟨4, by norm_num,
      by simpa [physicalClassAllows22, physicalClassWeight22] using h⟩
  · exact ⟨5, by norm_num,
      by simpa [physicalClassAllows22, physicalClassWeight22] using h⟩
  · exact ⟨6, by norm_num,
      by simpa [physicalClassAllows22, physicalClassWeight22] using h⟩
  · exact ⟨7, by norm_num,
      by simpa [physicalClassAllows22, physicalClassWeight22] using h⟩
  · exact ⟨8, by norm_num,
      by simpa [physicalClassAllows22, physicalClassWeight22] using h⟩
  · exact ⟨9, by norm_num,
      by simpa [physicalClassAllows22, physicalClassWeight22] using h⟩
  · exact ⟨10, by norm_num,
      by simpa [physicalClassAllows22, physicalClassWeight22] using h⟩

/-- A piece of state mass represented by one envelope coordinate.
Coordinate 2 additionally has pointwise density at most `b/(2^22-1)`
on each nonzero state. Keeping the mass cap as well makes the abstract
cone sound even before using the exact state count. -/
def physicalPiece22 (i : Fin 11) (b : ℝ)
    (μ : State22 → ℝ) : Prop :=
  (∀ q, 0 ≤ μ q) ∧
  (∀ q, μ q ≠ 0 → physicalClassAllows22 i q) ∧
  (∑ q : State22, μ q) ≤ b ∧
  (i.val = 2 → ∀ q : State22,
    μ q ≤ b / ((2 : ℝ) ^ 22 - 1))

/-- The eleven coordinates represent a complete nonnegative measure,
with a separately bounded piece for each coordinate. -/
def physicalCone22 (v : Fin 11 → ℝ)
    (μ : State22 → ℝ) : Prop :=
  ∃ parts : Fin 11 → State22 → ℝ,
    (∀ q, μ q = ∑ i : Fin 11, parts i q) ∧
      ∀ i, physicalPiece22 i (v i) (parts i)

theorem physicalPiece22_nonneg {i : Fin 11} {b : ℝ}
    {μ : State22 → ℝ} (h : physicalPiece22 i b μ) :
    ∀ q, 0 ≤ μ q := h.1

theorem physicalPiece22_mass_le {i : Fin 11} {b : ℝ}
    {μ : State22 → ℝ} (h : physicalPiece22 i b μ) :
    (∑ q : State22, μ q) ≤ b := h.2.2.1

theorem physicalCone22_nonneg {v : Fin 11 → ℝ}
    {μ : State22 → ℝ} (h : physicalCone22 v μ) :
    ∀ q, 0 ≤ μ q := by
  obtain ⟨parts, hsum, hparts⟩ := h
  intro q
  rw [hsum q]
  exact Finset.sum_nonneg (fun i _ => (hparts i).1 q)

theorem physicalCone22_mass_le {v : Fin 11 → ℝ}
    {μ : State22 → ℝ} (h : physicalCone22 v μ) :
    (∑ q : State22, μ q) ≤ ∑ i : Fin 11, v i := by
  obtain ⟨parts, hsum, hparts⟩ := h
  simp_rw [hsum]
  rw [Finset.sum_comm]
  exact Finset.sum_le_sum (fun i _ => (hparts i).2.2.1)

theorem physicalPiece22_sum {ι : Type*} [Fintype ι]
    (j : Fin 11) (bounds : ι → ℝ)
    (parts : ι → State22 → ℝ)
    (hparts : ∀ i, physicalPiece22 j (bounds i) (parts i)) :
    physicalPiece22 j (∑ i, bounds i)
      (fun q => ∑ i, parts i q) := by
  constructor
  · intro q
    exact Finset.sum_nonneg (fun i _ => (hparts i).1 q)
  constructor
  · intro q hq
    by_contra hallow
    have hz (i : ι) : parts i q = 0 := by
      by_contra hne
      exact hallow ((hparts i).2.1 q hne)
    apply hq
    simp only [hz, Finset.sum_const_zero]
  constructor
  · calc
      (∑ q : State22, ∑ i, parts i q) =
          ∑ i, ∑ q : State22, parts i q := Finset.sum_comm
      _ ≤ ∑ i, bounds i :=
        Finset.sum_le_sum (fun i _ => (hparts i).2.2.1)
  · intro hj q
    calc
      (∑ i, parts i q) ≤
          ∑ i, bounds i / ((2 : ℝ) ^ 22 - 1) :=
        Finset.sum_le_sum (fun i _ => (hparts i).2.2.2 hj q)
      _ = (∑ i, bounds i) / ((2 : ℝ) ^ 22 - 1) := by
        rw [Finset.sum_div]

def physicalKernelAction22 (K : Matrix State22 State22 ℝ)
    (μ : State22 → ℝ) : State22 → ℝ :=
  fun r => ∑ q : State22, μ q * K q r

def physicalEnvelopeApply22 (v : Fin 11 → ℝ)
    (B : Matrix (Fin 11) (Fin 11) ℝ) : Fin 11 → ℝ :=
  fun j => ∑ i : Fin 11, v i * B i j

theorem physicalKernelAction22_mul
    (K₁ K₂ : Matrix State22 State22 ℝ)
    (μ : State22 → ℝ) :
    physicalKernelAction22 (K₁ * K₂) μ =
      physicalKernelAction22 K₂
        (physicalKernelAction22 K₁ μ) := by
  funext r
  simp only [physicalKernelAction22, Matrix.mul_apply,
    Finset.mul_sum, Finset.sum_mul]
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro q _
  apply Finset.sum_congr rfl
  intro s _
  ring

theorem physicalEnvelopeApply22_mul
    (v : Fin 11 → ℝ)
    (B₁ B₂ : Matrix (Fin 11) (Fin 11) ℝ) :
    physicalEnvelopeApply22 v (B₁ * B₂) =
      physicalEnvelopeApply22
        (physicalEnvelopeApply22 v B₁) B₂ := by
  funext j
  simp only [physicalEnvelopeApply22, Matrix.mul_apply,
    Finset.mul_sum, Finset.sum_mul]
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro i _
  apply Finset.sum_congr rfl
  intro k _
  ring

/-- A local 11-coordinate matrix is sound if each complete source
piece can be split into sound output pieces with the stated whole row.
This condition is deliberately on complete rows; entrywise minima of
different candidate decompositions do not satisfy it automatically. -/
def physicalRowSound22 (K : Matrix State22 State22 ℝ)
    (B : Matrix (Fin 11) (Fin 11) ℝ) : Prop :=
  ∀ (i : Fin 11) (b : ℝ) (μ : State22 → ℝ),
    physicalPiece22 i b μ →
      ∃ parts : Fin 11 → State22 → ℝ,
        (∀ q, physicalKernelAction22 K μ q =
          ∑ j : Fin 11, parts j q) ∧
        ∀ j, physicalPiece22 j (b * B i j) (parts j)

/-- A sound complete-row matrix transfers any represented state
measure through one exact full-state kernel. -/
theorem physicalRowSound22_transfer
    {K : Matrix State22 State22 ℝ}
    {B : Matrix (Fin 11) (Fin 11) ℝ}
    (hB : physicalRowSound22 K B)
    {v : Fin 11 → ℝ} {μ : State22 → ℝ}
    (hμ : physicalCone22 v μ) :
    physicalCone22 (physicalEnvelopeApply22 v B)
      (physicalKernelAction22 K μ) := by
  classical
  obtain ⟨parts, hsum, hparts⟩ := hμ
  choose out hout hshape using
    (fun i : Fin 11 => hB i (v i) (parts i) (hparts i))
  refine ⟨(fun j q => ∑ i : Fin 11, out i j q), ?_, ?_⟩
  · intro r
    unfold physicalKernelAction22
    simp_rw [hsum]
    simp only [Finset.sum_mul]
    calc
      (∑ q : State22, ∑ i : Fin 11, parts i q * K q r) =
          ∑ i : Fin 11, physicalKernelAction22 K (parts i) r := by
        rw [Finset.sum_comm]
        rfl
      _ = ∑ i : Fin 11, ∑ j : Fin 11, out i j r := by
        simp only [hout]
      _ = ∑ j : Fin 11, ∑ i : Fin 11, out i j r :=
        Finset.sum_comm
  · intro j
    exact physicalPiece22_sum j (fun i => v i * B i j)
      (fun i q => out i j q) (fun i => hshape i j)

def physicalInitialEnvelope22 (i : Fin 11) : ℝ :=
  if i.val = 0 then 1 else 0

def physicalInitialMass22 (q : State22) : ℝ :=
  if q = 0 then 1 else 0

theorem physicalPiece22_zero (i : Fin 11) (b : ℝ)
    (hb : 0 ≤ b) :
    physicalPiece22 i b (fun _ => 0) := by
  constructor
  · intro q
    exact le_refl _
  constructor
  · intro q h
    exact (h rfl).elim
  constructor
  · simp only [Finset.sum_const_zero]
    exact hb
  · intro _ q
    have hden : (0 : ℝ) < (2 : ℝ) ^ 22 - 1 := by norm_num
    exact div_nonneg hb hden.le

theorem physicalInitialMass22_piece :
    physicalPiece22 (0 : Fin 11) 1 physicalInitialMass22 := by
  constructor
  · intro q
    simp only [physicalInitialMass22]
    split_ifs <;> norm_num
  constructor
  · intro q h
    change q = 0
    by_contra hn
    exact h (by simp [physicalInitialMass22, hn])
  constructor
  · simp only [physicalInitialMass22,
      Finset.sum_ite_eq', Finset.mem_univ, ite_true]
    exact le_refl 1
  · intro h
    omega

theorem physicalInitialCone22 :
    physicalCone22 physicalInitialEnvelope22
      physicalInitialMass22 := by
  refine ⟨(fun i q => if i = 0 then physicalInitialMass22 q else 0), ?_, ?_⟩
  · intro q
    simp only [Finset.sum_ite_eq', Finset.mem_univ, ite_true]
  · intro i
    by_cases hi : i = 0
    · subst i
      simpa only [physicalInitialEnvelope22, Fin.val_zero,
        ↓reduceIte] using physicalInitialMass22_piece
    · have hval : i.val ≠ 0 := by
        intro h
        exact hi (Fin.ext h)
      simpa only [hi, ite_false, physicalInitialEnvelope22,
        hval] using physicalPiece22_zero i 0 (le_refl 0)

#print axioms physicalCone22_nonneg
#print axioms physicalClass22_covers_nonzero
#print axioms physicalCone22_mass_le
#print axioms physicalRowSound22_transfer
#print axioms physicalKernelAction22_mul
#print axioms physicalEnvelopeApply22_mul
#print axioms physicalInitialCone22

end Spin.Turbo
