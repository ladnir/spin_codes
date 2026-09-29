import SpinCodes.Structured.ConcreteNativeLinearCode

noncomputable section
namespace Spin.Structured.ConcreteNativeFamily
open ConcreteOuter ConcreteEncoder Finset
attribute [local instance] Classical.propDecidable

lemma binaryWeight_sub {n : ℕ} (x y : Fin n→ZMod 2) :
    binaryWeight (x-y)=hammingDist x y := by
  simp only [binaryWeight,hammingDist,Pi.sub_apply,sub_ne_zero]

lemma codewords_zero_mem (m : ℕ) (out : NativeSeed m) (inner : InnerSeed m) :
    (0:Fin (Nsched m)→ZMod 2)∈codewords m out inner :=
  (mem_realizedCode m out inner _).mp (Submodule.zero_mem _)

lemma codewords_sub_mem {m : ℕ} {out : NativeSeed m} {inner : InnerSeed m}
    {x y : Fin (Nsched m)→ZMod 2} (hx : x∈codewords m out inner) (hy : y∈codewords m out inner) :
    x-y∈codewords m out inner :=
  (mem_realizedCode m out inner _).mp ((realizedCode m out inner).sub_mem
    ((mem_realizedCode m out inner _).mpr hx) ((mem_realizedCode m out inner _).mpr hy))

/-- In this actual realized linear code, low nonzero weight and low pairwise
Hamming distance are exactly the same event. -/
theorem low_weight_iff_close_pair (m : ℕ) (out : NativeSeed m) (inner : InnerSeed m) (d : ℕ) :
    (∃ c∈codewords m out inner,c≠0 ∧ binaryWeight c≤d) ↔
    ∃ x∈codewords m out inner,∃ y∈codewords m out inner,x≠y ∧ hammingDist x y≤d := by
  constructor
  · rintro ⟨c,hc,hn,hw⟩
    refine ⟨c,hc,0,codewords_zero_mem m out inner,hn,?_⟩
    simpa only [←binaryWeight_sub,sub_zero] using hw
  · rintro ⟨x,hx,y,hy,hne,hd⟩
    exact ⟨x-y,codewords_sub_mem hx hy,sub_ne_zero.mpr hne,(binaryWeight_sub x y).symm ▸ hd⟩

def nonzeroCodewords (m : ℕ) (out : NativeSeed m) (inner : InnerSeed m) :
    Finset (Fin (Nsched m)→ZMod 2) := (codewords m out inner).filter (fun c => c≠0)

lemma nonzeroCodewords_nonempty (m : ℕ) (out : NativeSeed m) (inner : InnerSeed m) :
    (nonzeroCodewords m out inner).Nonempty := by
  have hK : 0<Nsched m/2 := by have := native_dimension_twice m; have := Nsched_pos m; omega
  let x : Fin (Nsched m/2)→ZMod 2 := fun _ => 1
  have hx : x≠0 := by
    intro h
    have hh := congrFun h ⟨0,hK⟩
    norm_num [x] at hh
  refine ⟨codeword m out inner x,mem_filter.mpr ⟨mem_image.mpr ⟨x,mem_univ _,rfl⟩,?_⟩⟩
  intro hz
  exact hx (codeword_injective m out inner (hz.trans (codeword_zero m out inner).symm))

/-- The minimum Hamming distance of the actual realized binary linear code. -/
def minimumDistance (m : ℕ) (out : NativeSeed m) (inner : InnerSeed m) : ℕ :=
  (nonzeroCodewords m out inner).inf' (nonzeroCodewords_nonempty m out inner) binaryWeight

theorem minimumDistance_le_iff (m : ℕ) (out : NativeSeed m) (inner : InnerSeed m) (d : ℕ) :
    minimumDistance m out inner≤d ↔ ∃ c∈codewords m out inner,c≠0 ∧ binaryWeight c≤d := by
  simp only [minimumDistance,inf'_le_iff,nonzeroCodewords,mem_filter]
  aesop

/-- The named minimum distance equals the minimum over distinct emitted pairs. -/
theorem minimumDistance_le_iff_pair (m : ℕ) (out : NativeSeed m) (inner : InnerSeed m) (d : ℕ) :
    minimumDistance m out inner≤d ↔
    ∃ x∈codewords m out inner,∃ y∈codewords m out inner,x≠y ∧ hammingDist x y≤d :=
  (minimumDistance_le_iff m out inner d).trans (low_weight_iff_close_pair m out inner d)

theorem minimumDistance_gt_iff (m : ℕ) (out : NativeSeed m) (inner : InnerSeed m) (d : ℕ) :
    d < minimumDistance m out inner ↔
      ∀ x∈codewords m out inner,∀ y∈codewords m out inner,x≠y → d < hammingDist x y := by
  rw [←not_le,minimumDistance_le_iff_pair]
  push_neg
  rfl

/-- The Family bad probability is exactly failure of the actual minimum-distance threshold. -/
theorem concrete_probBad_minimumDistance (m : ℕ) : concreteFamily.probBad m=
    ((nativeSeedLaw m).prod
      (ConcreteRoutedEncoder.experimentLaw (Lsched m) (bsched m) (rounds m))).prob
      (fun ω => minimumDistance m ω.1 ω.2≤threshold m) := by
  rw [concrete_probBad_codewords]
  apply Spin.FinPMF.prob_congr
  intro ω
  exact (minimumDistance_le_iff m ω.1 ω.2 (threshold m)).symm

end Spin.Structured.ConcreteNativeFamily

