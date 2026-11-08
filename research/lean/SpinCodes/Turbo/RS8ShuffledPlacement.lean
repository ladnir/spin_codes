import SpinCodes.Turbo.PacketShuffleMixture

/-!
# Exact one-group Turbo packet placement

For each fixed nonzero RS8 group message, the sampled column multipliers
and an independent uniform row permutation induce a complete 64-packet
word law. Summing its point mass over group messages connects the
checked expected packet counts `C_v` to the actual value-carrying route.
-/

noncomputable section

namespace Spin.Turbo

theorem sampledRS8RowLaw_weight_prob
    (points : Fin 8 → Field16) (rows : RS8Group Field16)
    (v : ℕ) :
    (sampledRS8RowLaw points rows).prob
      (fun y => (packetWordSupport y).card = v) =
    (piPMF (fun _ : Fin 8 =>
      FinPMF.uniform {c : Field32 // c ≠ 0})).prob
        (fun coeff => rs8SampledPacketWeight points rows coeff = v) := by
  rw [sampledRS8RowLaw, FinPMF.map_prob,
    sampledRS8PacketColumnsLaw, FinPMF.map_prob]
  apply FinPMF.prob_congr
  intro coeff
  exact (by
    change (packetWordSupport (sampledRS8RowWord points rows coeff)).card = v ↔
      rs8SampledPacketWeight points rows coeff = v
    rw [sampledRS8RowWord_weight])

/-- The complete target atom for one fixed message depends only on its
packet weight. The coefficient tuple is sampled once per setup. -/
theorem sampledRS8Row_shuffled_apply
    (points : Fin 8 → Field16) (rows : RS8Group Field16)
    (y : PacketWord64) :
    (shuffledPacketLaw (sampledRS8RowLaw points rows)).p y =
    (piPMF (fun _ : Fin 8 =>
      FinPMF.uniform {c : Field32 // c ≠ 0})).prob
        (fun coeff => rs8SampledPacketWeight points rows coeff =
          (packetWordSupport y).card) /
      ((Nat.choose 64 (packetWordSupport y).card : ℝ) *
        (15 : ℝ) ^ (packetWordSupport y).card) := by
  rw [shuffledPacketLaw_apply _
    (sampledRS8RowLaw_fiber_uniform points rows) y,
    sampledRS8RowLaw_weight_prob]

/-- Expected number of nonzero group messages whose shuffled packet row
equals the specified complete target word. -/
def rs8ExpectedShuffledWordCount
    (points : Fin 8 → Field16) (y : PacketWord64) : ℝ :=
  ∑ rows ∈ nonzeroMsgs (RS8Group Field16),
    (shuffledPacketLaw (sampledRS8RowLaw points rows)).p y

/-- The exact `C_v` placement formula, including the fifteen choices
for each active four-bit packet label. -/
theorem rs8_expected_shuffled_word_count
    (points : Fin 8 → Field16) (y : PacketWord64) :
    rs8ExpectedShuffledWordCount points y =
      rs8ExpectedPacketCount points (packetWordSupport y).card /
        ((Nat.choose 64 (packetWordSupport y).card : ℝ) *
          (15 : ℝ) ^ (packetWordSupport y).card) := by
  unfold rs8ExpectedShuffledWordCount rs8ExpectedPacketCount
  simp_rw [sampledRS8Row_shuffled_apply]
  rw [Finset.sum_div]

#print axioms sampledRS8RowLaw_weight_prob
#print axioms sampledRS8Row_shuffled_apply
#print axioms rs8_expected_shuffled_word_count

end Spin.Turbo
