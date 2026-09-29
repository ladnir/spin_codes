import SpinCodes.Structured.ConcreteRoutedOutput
import SpinCodes.Structured.ConcreteNativeTotal

noncomputable section
namespace Spin.Structured.ConcreteEncoder
open Finset

theorem binaryWeight_reindex {n m : ℕ} (e : Fin n ≃ Fin m) (x : Fin m → ZMod 2) :
    binaryWeight (fun i => x (e i)) = binaryWeight x := by
  simp only [binaryWeight, card_eq_sum_ones, sum_filter]
  exact Equiv.sum_comp e (fun j => if x j ≠ 0 then (1 : ℕ) else 0)

end Spin.Structured.ConcreteEncoder

namespace Spin.Structured.ConcreteNativeFamily
open ConcreteOuter ConcreteRoute ConcreteEncoder ConcreteBinaryEncoding Finset

theorem output_length (m : ℕ) : rounds m*128 = Nsched m := by
  rw [rounds, Nat.div_mul_cancel (native_round_divisibility m), Nsched_eq, Nat.mul_comm]

/-- The actual emitted binary word on exactly the paper's native output coordinates. -/
def codeword (m : ℕ) (out : NativeSeed m) (inner : InnerSeed m)
    (x : Fin (Nsched m/2) → ZMod 2) : Fin (Nsched m) → ZMod 2 :=
  fun i => ConcreteRoutedEncoder.word (streamWiring (native_round_divisibility m))
    (rows m out x) inner (finCongr (output_length m).symm i)

theorem codeword_weight (m : ℕ) (out : NativeSeed m) (inner : InnerSeed m)
    (x : Fin (Nsched m/2) → ZMod 2) :
    binaryWeight (codeword m out inner x) =
      ConcreteRoutedEncoder.weight (streamWiring (native_round_divisibility m)) (rows m out x) inner := by
  exact (binaryWeight_reindex (finCongr (output_length m).symm)
    (ConcreteRoutedEncoder.word (R := rounds m) (streamWiring (native_round_divisibility m))
      (rows m out x) inner)).trans (ConcreteRoutedEncoder.word_weight _ _ _)

theorem codeword_injective (m : ℕ) (out : NativeSeed m) (inner : InnerSeed m) :
    Function.Injective (codeword m out inner) := by
  let e : Fin (Nsched m) ≃ Fin (rounds m*128) := finCongr (output_length m).symm
  have he : Function.Injective (fun v : Fin (rounds m*128) → ZMod 2 => fun i => v (e i)) := by
    intro v w h
    funext j
    have hh := congrFun h (e.symm j)
    simpa only [Equiv.apply_symm_apply] using hh
  exact he.comp ((ConcreteRoutedEncoder.word_injective (R := rounds m)
    (streamWiring (native_round_divisibility m)) inner).comp
      ((native_rows_injective m out).comp (wordEquiv (Nsched m/2)).symm.injective))

theorem rows_zero (m : ℕ) (out : NativeSeed m) : rows m out 0 = fun _ => ∅ := by
  funext i
  simp [rows, nativeRows, rowSupports, support]

theorem codeword_zero (m : ℕ) (out : NativeSeed m) (inner : InnerSeed m) : codeword m out inner 0 = 0 := by
  let e := streamWiring (native_round_divisibility m)
  let v := ConcreteRoutedEncoder.word e (fun _ => ∅) inner
  have h := ConcreteRoutedEncoder.word_xor e (fun _ => ∅) (fun _ => ∅) inner
  simp only [symmDiff_self] at h
  have hh : v+v=v+0 := by
    change v = v+v at h
    simpa only [add_zero] using h.symm
  have hv : v=0 := add_left_cancel hh
  unfold codeword
  rw [rows_zero]
  change (fun i => v (finCongr (output_length m).symm i)) = 0
  rw [hv]
  rfl

/-- The realized code as a concrete set of emitted binary words. -/
def codewords (m : ℕ) (out : NativeSeed m) (inner : InnerSeed m) :
    Finset (Fin (Nsched m) → ZMod 2) := univ.image (codeword m out inner)

theorem codewords_card (m : ℕ) (out : NativeSeed m) (inner : InnerSeed m) :
    (codewords m out inner).card = 2^(Nsched m/2) := by
  rw [codewords, card_image_of_injective _ (codeword_injective m out inner), card_univ]
  simp

/-- The Family bad event is exactly existence of a nonzero emitted codeword of low weight. -/
theorem codeword_bad_iff (m : ℕ) (out : NativeSeed m) (inner : InnerSeed m) :
    (∃ x : Fin (Nsched m/2) → ZMod 2, x ≠ 0 ∧
      ConcreteRoutedEncoder.weight (streamWiring (native_round_divisibility m)) (rows m out x) inner ≤ threshold m) ↔
    ∃ c ∈ codewords m out inner, c ≠ 0 ∧ binaryWeight c ≤ threshold m := by
  constructor
  · rintro ⟨x, hx, hw⟩
    refine ⟨codeword m out inner x, mem_image.mpr ⟨x, mem_univ x, rfl⟩, ?_, ?_⟩
    · intro hz
      exact hx (codeword_injective m out inner (hz.trans (codeword_zero m out inner).symm))
    · rwa [codeword_weight]
  · rintro ⟨c, hc, hn, hw⟩
    obtain ⟨x, _, rfl⟩ := mem_image.mp hc
    refine ⟨x, ?_, ?_⟩
    · intro hx
      exact hn (hx ▸ codeword_zero m out inner)
    · rwa [codeword_weight] at hw

theorem concrete_probBad_codewords (m : ℕ) : concreteFamily.probBad m =
    ((nativeSeedLaw m).prod
      (ConcreteRoutedEncoder.experimentLaw (Lsched m) (bsched m) (rounds m))).prob
      (fun ω => ∃ c ∈ codewords m ω.1 ω.2, c ≠ 0 ∧ binaryWeight c ≤ threshold m) := by
  rw [concrete_probBad_explicit]
  apply FinPMF.prob_congr
  intro ω
  exact codeword_bad_iff m ω.1 ω.2

end Spin.Structured.ConcreteNativeFamily
