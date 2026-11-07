import SpinCodes.UniformMixer.QuadraticProductionSpectrum
import SpinCodes.UniformMixer.Map22Rank

/-!
# Support and mass of the production-map weight spectrum

Every production word lies in one of the nine weights admitted by the
quadratic radical classes. Summing the corresponding multiplicities counts
each production state once.
-/

namespace Spin.UniformMixer

def spectrumWeights : Finset Nat := {0, 16, 24, 28, 32, 36, 40, 48, 64}

theorem quadraticWordWeight_mem_spectrumWeights
    (b : ZMod 2) (a : State6) (c : Fin 15 → ZMod 2) :
    quadraticWordWeight b a c ∈ spectrumWeights := by
  rcases quadraticRadical_four_classes c with h64 | h16 | h4 | h1
  · have h := quadraticWordWeight_class b a c 64 0 64
      (by omega) (by omega) (fun ha => walsh_value_radical_64 c h64 a ha)
    rcases h with h | h | h <;> simp [spectrumWeights, h]
  · have h := quadraticWordWeight_class b a c 32 16 48
      (by omega) (by omega) (fun ha => walsh_value_radical_16 c h16 a ha)
    rcases h with h | h | h <;> simp [spectrumWeights, h]
  · have h := quadraticWordWeight_class b a c 16 24 40
      (by omega) (by omega) (fun ha => walsh_value_radical_4 c h4 a ha)
    rcases h with h | h | h <;> simp [spectrumWeights, h]
  · have h := quadraticWordWeight_class b a c 8 28 36
      (by omega) (by omega) (fun ha => walsh_value_radical_1 c h1 a ha)
    rcases h with h | h | h <;> simp [spectrumWeights, h]

theorem productionWeight_mem_spectrumWeights (q : State22) :
    coordinateWeight coordMap q ∈ spectrumWeights := by
  rw [productionWeight_as_triple]
  exact quadraticWordWeight_mem_spectrumWeights
    (productionTriple q).1 (productionTriple q).2.1
    (productionTriple q).2.2

theorem productionWeightCount_total_mass :
    (∑ w ∈ spectrumWeights, productionWeightCount w) = 2 ^ 22 := by
  simp only [productionWeightCount]
  rw [Finset.sum_comm]
  have hinner (q : State22) :
      (∑ w ∈ spectrumWeights,
        if coordinateWeight coordMap q = w then 1 else 0) = 1 := by
    have hmem := productionWeight_mem_spectrumWeights q
    simp only [eq_comm, Finset.sum_ite_eq', hmem, ite_true]
  simp_rw [hinner]
  norm_num [Fintype.card_fun]

theorem productionWeightCount_32 : productionWeightCount 32 = 1828134 := by
  have h := productionWeightCount_total_mass
  change spectrumWeights.sum productionWeightCount = 2 ^ 22 at h
  simp [spectrumWeights] at h
  rw [productionWeightCount_0, productionWeightCount_16,
    productionWeightCount_24, productionWeightCount_28,
    productionWeightCount_36, productionWeightCount_40,
    productionWeightCount_48, productionWeightCount_64] at h
  generalize hk : productionWeightCount 32 = k at *
  omega

theorem productionWeightCount_other_zero (w : Nat)
    (hw : w ∉ spectrumWeights) : productionWeightCount w = 0 := by
  unfold productionWeightCount
  apply Finset.sum_eq_zero
  intro q _
  have hne : coordinateWeight coordMap q ≠ w := by
    intro h
    exact hw (h ▸ productionWeight_mem_spectrumWeights q)
  simp [hne]

def productionSpectrumMultiplicity (w : Nat) : Nat :=
  if w = 0 then 1 else
  if w = 16 then 2604 else
  if w = 24 then 291648 else
  if w = 28 then 888832 else
  if w = 32 then 1828134 else
  if w = 36 then 888832 else
  if w = 40 then 291648 else
  if w = 48 then 2604 else
  if w = 64 then 1 else 0

theorem productionWeightCount_exact (w : Nat) :
    productionWeightCount w = productionSpectrumMultiplicity w := by
  by_cases hmem : w ∈ spectrumWeights
  · simp only [spectrumWeights, Finset.mem_insert,
      Finset.mem_singleton] at hmem
    rcases hmem with h | h | h | h | h | h | h | h | h <;>
      subst w <;>
      simp [productionSpectrumMultiplicity, productionWeightCount_0,
        productionWeightCount_16, productionWeightCount_24,
        productionWeightCount_28, productionWeightCount_32,
        productionWeightCount_36, productionWeightCount_40,
        productionWeightCount_48, productionWeightCount_64]
  · rw [productionWeightCount_other_zero w hmem]
    unfold productionSpectrumMultiplicity
    have h0 : w ≠ 0 := by intro h; exact hmem (h.symm ▸ by decide)
    have h16 : w ≠ 16 := by intro h; exact hmem (h.symm ▸ by decide)
    have h24 : w ≠ 24 := by intro h; exact hmem (h.symm ▸ by decide)
    have h28 : w ≠ 28 := by intro h; exact hmem (h.symm ▸ by decide)
    have h32 : w ≠ 32 := by intro h; exact hmem (h.symm ▸ by decide)
    have h36 : w ≠ 36 := by intro h; exact hmem (h.symm ▸ by decide)
    have h40 : w ≠ 40 := by intro h; exact hmem (h.symm ▸ by decide)
    have h48 : w ≠ 48 := by intro h; exact hmem (h.symm ▸ by decide)
    have h64 : w ≠ 64 := by intro h; exact hmem (h.symm ▸ by decide)
    simp [h0, h16, h24, h28, h32, h36, h40, h48, h64]

theorem productionWeight_zero_only (q : State22)
    (hw : coordinateWeight coordMap q = 0) : q = 0 := by
  have hpoint (j : Fin 64) :
      (if coordMap j q = 1 then 1 else 0 : Nat) = 0 := by
    have hsum :
        (∑ j : Fin 64, if coordMap j q = 1 then 1 else 0 : Nat) = 0 := hw
    exact (Finset.sum_eq_zero_iff_of_nonneg
      (fun _ _ => Nat.zero_le _)).mp hsum j (Finset.mem_univ j)
  have hmap : (fun j : Fin 64 => coordMap j q) =
      (fun j : Fin 64 => coordMap j (0 : State22)) := by
    funext j
    have hne : coordMap j q ≠ 1 := by
      intro he
      have hp := hpoint j
      simp [he] at hp
    have hz := (zmod2_zero_or_one (coordMap j q)).resolve_right hne
    simp [hz]
  exact expansion_injective hmap

theorem productionWeight_nonzero_lower_bound (q : State22) (hq : q ≠ 0) :
    16 ≤ coordinateWeight coordMap q := by
  have hw := productionWeight_mem_spectrumWeights q
  have hzero : coordinateWeight coordMap q ≠ 0 := by
    intro h
    exact hq (productionWeight_zero_only q h)
  simp only [spectrumWeights, Finset.mem_insert,
    Finset.mem_singleton] at hw
  omega

theorem productionWeight_sixteen_exists :
    ∃ q : State22, q ≠ 0 ∧ coordinateWeight coordMap q = 16 := by
  have hcard :
      (Finset.univ.filter (fun q : State22 =>
        coordinateWeight coordMap q = 16)).card = 2604 := by
    rw [← productionWeightCount_eq_card_filter, productionWeightCount_16]
  obtain ⟨q, hq⟩ := Finset.card_pos.mp (by rw [hcard]; omega)
  have hw : coordinateWeight coordMap q = 16 := (Finset.mem_filter.mp hq).2
  refine ⟨q, ?_, hw⟩
  intro hzero
  subst q
  simp [coordinateWeight_zero] at hw

end Spin.UniformMixer
