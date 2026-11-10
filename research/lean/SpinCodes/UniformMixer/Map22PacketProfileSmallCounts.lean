import SpinCodes.UniformMixer.Map22PacketProfileSmallDefs

/-! Kernel-checked weight histograms on only 2048, 128, and 16 small quadratics. -/
noncomputable section
namespace Spin.UniformMixer
open Spin.Structured

def profileSmallPackedWeight (width : Nat) (rows : List Nat) (q : Nat) : Nat :=
  MapSpectrum.fastWeight width (PackedMap.eval rows q)

def profileSmallPackedHistogram (width bits : Nat) (rows : List Nat) : List Nat :=
  (List.range (2^bits)).foldl
    (fun h q => MapSpectrum.bump (profileSmallPackedWeight width rows q) h)
    (List.replicate (width+1) 0)

def profileQuad4WeightCounts : List Nat := [1,0,0,0,140,0,448,0,870,0,448,0,140,0,0,0,1]
def profileQuad3WeightCounts : List Nat := [1,0,28,0,70,0,28,0,1]
def profileQuad2WeightCounts : List Nat := [1,4,6,4,1]

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem profileQuad4_histogram_checked :
    profileSmallPackedHistogram 16 11 profileQuad4Rows = profileQuad4WeightCounts := by decide +kernel

set_option maxRecDepth 20000 in
set_option maxHeartbeats 0 in
theorem profileQuad3_histogram_checked :
    profileSmallPackedHistogram 8 7 profileQuad3Rows = profileQuad3WeightCounts := by decide +kernel

set_option maxRecDepth 20000 in
theorem profileQuad2_histogram_checked :
    profileSmallPackedHistogram 4 4 profileQuad2Rows = profileQuad2WeightCounts := by decide +kernel

theorem profileSmallPackedHistogram_count (width bits : Nat) (rows : List Nat)
    (w : Fin (width+1)) :
    (profileSmallPackedHistogram width bits rows).getD w.val 0 =
      (Finset.univ.filter (fun q : Fin (2^bits) => profileSmallPackedWeight width rows q.val = w.val)).card := by
  classical
  unfold profileSmallPackedHistogram
  rw [MapSpectrum.foldl_bump_count _ _ _ _ (by simpa using w.isLt)]
  simp only [List.getD_replicate 0 w.isLt, Nat.zero_add]
  rw [MapSpectrum.countP_range_eq_card]
  simp

theorem profileQuad4PackedWeight_eq (q : Nat) :
    profileSmallPackedWeight 16 profileQuad4Rows q = profileQuad4Weight (profilePackedBits 11 q) := by
  rw [profileSmallPackedWeight, MapSpectrum.fastWeight_eq, PackedMap.weight_eq_sum_bits]
  unfold profileQuad4Weight
  apply Finset.sum_congr rfl
  intro v _
  rw [packedBit_weight, packed_eval_coordinate]
  rfl

theorem profileQuad3PackedWeight_eq (q : Nat) :
    profileSmallPackedWeight 8 profileQuad3Rows q = profileQuad3Weight (profilePackedBits 7 q) := by
  rw [profileSmallPackedWeight, MapSpectrum.fastWeight_eq, PackedMap.weight_eq_sum_bits]
  unfold profileQuad3Weight
  apply Finset.sum_congr rfl
  intro v _
  rw [packedBit_weight, packed_eval_coordinate]
  rfl

theorem profileQuad2PackedWeight_eq (q : Nat) :
    profileSmallPackedWeight 4 profileQuad2Rows q = profileQuad2Weight (profilePackedBits 4 q) := by
  rw [profileSmallPackedWeight, MapSpectrum.fastWeight_eq, PackedMap.weight_eq_sum_bits]
  unfold profileQuad2Weight
  apply Finset.sum_congr rfl
  intro v _
  rw [packedBit_weight, packed_eval_coordinate]
  rfl

theorem profileQuad4_weight_count (w : Fin 17) :
    (∑ q : ProfileQuad4, if profileQuad4Weight q = w.val then 1 else 0) =
      profileQuad4WeightCounts.getD w.val 0 := by
  classical
  calc
    _ = ∑ q : Fin (2^11), if profileQuad4Weight (profilePackedBits 11 q.val) = w.val then 1 else 0 := by
      symm
      apply Fintype.sum_equiv (profilePackedBitsEquiv 11)
      intro q; rfl
    _ = (Finset.univ.filter (fun q : Fin (2^11) => profileSmallPackedWeight 16 profileQuad4Rows q.val = w.val)).card := by
      rw [Finset.card_filter]
      simp only [profileQuad4PackedWeight_eq]
    _ = _ := by rw [← profileSmallPackedHistogram_count 16 11 profileQuad4Rows w,
        profileQuad4_histogram_checked]

theorem profileQuad3_weight_count (w : Fin 9) :
    (∑ q : ProfileQuad3, if profileQuad3Weight q = w.val then 1 else 0) =
      profileQuad3WeightCounts.getD w.val 0 := by
  classical
  calc
    _ = ∑ q : Fin (2^7), if profileQuad3Weight (profilePackedBits 7 q.val) = w.val then 1 else 0 := by
      symm
      apply Fintype.sum_equiv (profilePackedBitsEquiv 7)
      intro q; rfl
    _ = (Finset.univ.filter (fun q : Fin (2^7) => profileSmallPackedWeight 8 profileQuad3Rows q.val = w.val)).card := by
      rw [Finset.card_filter]
      simp only [profileQuad3PackedWeight_eq]
    _ = _ := by rw [← profileSmallPackedHistogram_count 8 7 profileQuad3Rows w,
        profileQuad3_histogram_checked]

theorem profileQuad2_weight_count (w : Fin 5) :
    (∑ q : ProfileQuad2, if profileQuad2Weight q = w.val then 1 else 0) =
      profileQuad2WeightCounts.getD w.val 0 := by
  classical
  calc
    _ = ∑ q : Fin (2^4), if profileQuad2Weight (profilePackedBits 4 q.val) = w.val then 1 else 0 := by
      symm
      apply Fintype.sum_equiv (profilePackedBitsEquiv 4)
      intro q; rfl
    _ = (Finset.univ.filter (fun q : Fin (2^4) => profileSmallPackedWeight 4 profileQuad2Rows q.val = w.val)).card := by
      rw [Finset.card_filter]
      simp only [profileQuad2PackedWeight_eq]
    _ = _ := by rw [← profileSmallPackedHistogram_count 4 4 profileQuad2Rows w,
        profileQuad2_histogram_checked]

#print axioms profileQuad4_weight_count
#print axioms profileQuad3_weight_count
#print axioms profileQuad2_weight_count
end Spin.UniformMixer
