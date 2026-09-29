import SpinCodes.Structured.LiveInduction
import SpinCodes.Structured.Occupation

/-! Lift zero, diffuse, and shell row bounds through a nonnegative kernel. -/

namespace Spin.Imt
open Finset
variable {s k : ℕ}

def Coords.add (c d : Coords k) : Coords k := ⟨c.Z + d.Z, c.D + d.D, fun i => c.S i + d.S i⟩

theorem LiveDominates.of_le {sys : ShellSystem s k} {c : Coords k}
    {μ ν : Finset (Fin s) → ℝ} (h : LiveDominates sys c ν) (hm : ∀ q, μ q ≤ ν q) :
    LiveDominates sys c μ := by
  obtain ⟨⟨v, hv, ht, hb⟩, hz⟩ := h
  exact ⟨⟨v, hv, ht, fun q => (hm q).trans (hb q)⟩, (hm ∅).trans hz⟩

theorem liveDominates_zero (sys : ShellSystem s k) :
    LiveDominates sys ⟨0, 0, fun _ => 0⟩ (fun _ => 0) := by
  classical
  refine ⟨⟨fun _ => 0, fun _ => le_rfl, by simp, ?_⟩, le_rfl⟩
  intro q
  simp

theorem LiveDominates.add {sys : ShellSystem s k} {c d : Coords k}
    {μ ν : Finset (Fin s) → ℝ} (h : LiveDominates sys c μ) (h' : LiveDominates sys d ν) :
    LiveDominates sys (c.add d) (fun q => μ q + ν q) := by
  classical
  obtain ⟨⟨v, hv, ht, hb⟩, hz⟩ := h
  obtain ⟨⟨v', hv', ht', hb'⟩, hz'⟩ := h'
  refine ⟨⟨fun q => v q + v' q, fun q => add_nonneg (hv q) (hv' q), ?_, ?_⟩,
    add_le_add hz hz'⟩
  · simpa only [Finset.sum_add_distrib, Coords.add] using add_le_add ht ht'
  · intro q
    have hh := add_le_add (hb q) (hb' q)
    have he : (∑ i, if q ∈ sys.shell i then (c.S i + d.S i) / ((sys.shell i).card : ℝ) else 0) =
        (∑ i, if q ∈ sys.shell i then c.S i / ((sys.shell i).card : ℝ) else 0) +
        (∑ i, if q ∈ sys.shell i then d.S i / ((sys.shell i).card : ℝ) else 0) := by
      rw [← Finset.sum_add_distrib]
      apply Finset.sum_congr rfl
      intro i _
      split_ifs <;> ring
    change μ q + ν q ≤ (if q = ∅ then c.Z + d.Z else 0) + (v q + v' q) + _
    dsimp only [Coords.add]
    rw [he]
    split_ifs at * <;> linarith

def kernelApply (K : Finset (Fin s) → Finset (Fin s) → ℝ)
    (μ : Finset (Fin s) → ℝ) (r : Finset (Fin s)) : ℝ := ∑ q, μ q * K q r

theorem kernelApply_liveDominates (sys : ShellSystem s k) (T : Transfer k)
    (K : Finset (Fin s) → Finset (Fin s) → ℝ) (hK : ∀ q r, 0 ≤ K q r)
    (hD : T.rowD.Nonneg)
    (hZrow : LiveDominates sys T.rowZ (K ∅))
    (hDrow : ∀ q ≠ ∅, LiveDominates sys T.rowD (K q))
    (hSrow : ∀ i, LiveDominates sys (T.rowS i)
      (fun r => (∑ q ∈ sys.shell i, K q r) / ((sys.shell i).card : ℝ)))
    (c : Coords k) (μ : Finset (Fin s) → ℝ) (hc : c.Nonneg)
    (hμ : LiveDominates sys c μ) : LiveDominates sys (T.apply c) (kernelApply K μ) := by
  classical
  obtain ⟨v, hv, hv0, hvt, hb⟩ := (liveDominates_iff_witness sys c μ).mp hμ
  have hd_each (q : Finset (Fin s)) :
      LiveDominates sys (Coords.smul (v q) T.rowD) (fun r => v q * K q r) := by
    by_cases hq : q = ∅
    · subst q
      simpa only [hv0, Coords.smul, zero_mul] using liveDominates_zero sys
    · exact liveDominates_smul (hDrow q hq) (hv q)
  have hd := liveDominates_sum hd_each
  have he : Coords.sum (fun q => Coords.smul (v q) T.rowD) = Coords.smul (∑ q, v q) T.rowD := by
    apply Coords.ext
    · simp only [Coords.sum, Coords.smul, Finset.sum_mul]
    · simp only [Coords.sum, Coords.smul, Finset.sum_mul]
    · funext i
      simp only [Coords.sum, Coords.smul, Finset.sum_mul]
  rw [he] at hd
  have hd' : LiveDominates sys (Coords.smul c.D T.rowD) (fun r => ∑ q, v q * K q r) :=
    hd.mono (mul_le_mul_of_nonneg_right hvt hD.1)
      (mul_le_mul_of_nonneg_right hvt hD.2.1) (fun i => mul_le_mul_of_nonneg_right hvt (hD.2.2 i))
  have hz := liveDominates_smul hZrow hc.1
  have hs := liveDominates_sum (fun i => liveDominates_smul (hSrow i) (hc.2.2 i))
  have hall := (hz.add hd').add hs
  have hcoords : ((Coords.smul c.Z T.rowZ).add (Coords.smul c.D T.rowD)).add
      (Coords.sum fun i => Coords.smul (c.S i) (T.rowS i)) = T.apply c := rfl
  rw [hcoords] at hall
  apply hall.of_le
  intro r
  have heZ : (∑ q : Finset (Fin s), (if q = ∅ then c.Z else 0) * K q r) = c.Z * K ∅ r := by
    simp only [ite_mul, zero_mul, Finset.sum_ite_eq', Finset.mem_univ, ite_true]
  have heS : (∑ q : Finset (Fin s),
      (∑ i, if q ∈ sys.shell i then c.S i / ((sys.shell i).card : ℝ) else 0) * K q r) =
      ∑ i, c.S i * ((∑ q ∈ sys.shell i, K q r) / ((sys.shell i).card : ℝ)) := by
    simp only [Finset.sum_mul]
    rw [Finset.sum_comm]
    apply Finset.sum_congr rfl
    intro i _
    simp only [ite_mul, zero_mul, Finset.sum_ite_mem, Finset.univ_inter]
    rw [← Finset.mul_sum]
    ring
  calc
    kernelApply K μ r ≤ ∑ q, ((if q = ∅ then c.Z else 0) + v q +
        ∑ i, if q ∈ sys.shell i then c.S i / ((sys.shell i).card : ℝ) else 0) * K q r :=
      Finset.sum_le_sum fun q _ => mul_le_mul_of_nonneg_right (hb q) (hK q r)
    _ = _ := by
      simp only [add_mul, Finset.sum_add_distrib]
      rw [heZ, heS]

end Spin.Imt
