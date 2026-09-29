import SpinCodes.Structured.MapSpectrumSumDefs
import SpinCodes.Structured.MapSpectrum

namespace Spin.Structured.MapSpectrum

theorem getD_zipWith_add (a b : List ℕ) {i : ℕ}
    (ha : i < a.length) (hb : i < b.length) :
    (List.zipWith (· + ·) a b).getD i 0 = a.getD i 0 + b.getD i 0 := by
  induction a generalizing b i with
  | nil => simp at ha
  | cons x a ih =>
    cases b with
    | nil => simp at hb
    | cons y b =>
      cases i with
      | zero => rfl
      | succ i => exact ih b (by simpa using ha) (by simpa using hb)

theorem sumHistograms_length (hs : List (List ℕ)) (h : ∀ a ∈ hs, a.length = 129) :
    (sumHistograms hs).length = 129 := by
  induction hs with
  | nil => simp [sumHistograms]
  | cons a hs ih =>
    change (List.zipWith (· + ·) a (sumHistograms hs)).length = 129
    rw [List.length_zipWith, h a (by simp), ih (fun a ha => h a (by simp [ha]))]
    rfl

theorem getD_sumHistograms (hs : List (List ℕ)) (h : ∀ a ∈ hs, a.length = 129)
    {i : ℕ} (hi : i < 129) :
    (sumHistograms hs).getD i 0 = (hs.map fun a => a.getD i 0).sum := by
  induction hs with
  | nil => simp only [sumHistograms, List.foldr_nil, List.getD_replicate 0 hi,
      List.map_nil, List.sum_nil]
  | cons a hs ih =>
    have ht : ∀ b ∈ hs, b.length = 129 := fun b hb => h b (by simp [hb])
    change (List.zipWith (· + ·) a (sumHistograms hs)).getD i 0 = _
    rw [getD_zipWith_add _ _ (by simpa [h a (by simp)] using hi)
      (by simpa [sumHistograms_length hs ht] using hi), ih ht]
    rfl

end Spin.Structured.MapSpectrum
