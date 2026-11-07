import SpinCodes.Turbo.PacketProductMoment

/-!
# The one-column packet generating function

The coefficient of `z^v` counts 32-bit columns with `v` active four-bit
packets. The binomial identity is over a fixed eight-packet column and
is checked in Lean before connecting it to the probability law.
-/

noncomputable section

namespace Spin.Turbo

theorem packet_binomial8 (z : ℝ) :
    (∑ v ∈ Finset.range 9,
      (Nat.choose 8 v : ℝ) * 15 ^ v * z ^ v) =
      (1 + 15 * z) ^ 8 := by
  norm_num [Finset.sum_range_succ, Nat.choose]
  ring

private theorem packet_zero_weight_impossible
    (x : {x : PacketColumn // x ≠ 0}) :
    (packetColumnSupport x.val).card ≠ 0 := by
  intro hz
  have hempty : packetColumnSupport x.val = ∅ := Finset.card_eq_zero.mp hz
  apply x.property
  funext i
  by_contra hn
  have hi : i ∈ packetColumnSupport x.val :=
    Finset.mem_filter.mpr ⟨Finset.mem_univ _, hn⟩
  rw [hempty] at hi
  simp at hi

theorem uniform_nonzero_packet_zero_prob :
    (FinPMF.uniform {x : PacketColumn // x ≠ 0}).prob
      (fun x => (packetColumnSupport x.val).card = 0) = 0 := by
  classical
  have heq : (Finset.univ.filter (fun x :
      {x : PacketColumn // x ≠ 0} =>
      (packetColumnSupport x.val).card = 0)) = ∅ := by
    ext x
    constructor
    · intro hx
      exact (packet_zero_weight_impossible x (Finset.mem_filter.mp hx).2).elim
    · intro hx
      simp at hx
  unfold FinPMF.prob
  rw [heq]
  simp

/-- The one-column packet-weight generating function is exactly the
paper's `P(z)/(2^32-1)`, where `P(z)=(1+15z)^8-1`. -/
theorem uniform_nonzero_packet_gf (z : ℝ) :
    (FinPMF.uniform {x : PacketColumn // x ≠ 0}).expect
      (fun x => z ^ (packetColumnSupport x.val).card) =
      ((1 + 15 * z) ^ 8 - 1) / (2 ^ 32 - 1 : ℝ) := by
  let P := FinPMF.uniform {x : PacketColumn // x ≠ 0}
  have hwt (x : {x : PacketColumn // x ≠ 0}) :
      (packetColumnSupport x.val).card ∈ Finset.range 9 := by
    simp only [Finset.mem_range]
    have hle := Finset.card_le_card
      (Finset.filter_subset (p := fun i => x.val i ≠ 0)
        (Finset.univ : Finset (Fin 8)))
    simp only [Finset.card_univ, Fintype.card_fin] at hle
    change (packetColumnSupport x.val).card ≤ 8 at hle
    omega
  calc
    P.expect (fun x => z ^ (packetColumnSupport x.val).card) =
        ∑ v ∈ Finset.range 9,
          P.prob (fun x => (packetColumnSupport x.val).card = v) * z ^ v :=
      FinPMF.expect_comp P (fun x => (packetColumnSupport x.val).card)
        (Finset.range 9) hwt (fun v => z ^ v)
    _ = _ := by
      have hzero : P.prob
          (fun x => packetColumnSupport x.val = ∅) = 0 := by
        simpa only [Finset.card_eq_zero] using
          uniform_nonzero_packet_zero_prob
      norm_num [Finset.sum_range_succ,
        uniform_nonzero_packet_weight_prob, Nat.choose, P]
      rw [hzero]
      ring

/-- For `h` independent active columns, the packet-weight generating
function is the `h`th power of the one-column generating function. -/
theorem uniform_active_columns_gf (h : ℕ) (z : ℝ) :
    (piPMF (fun _ : Fin h =>
      FinPMF.uniform {x : PacketColumn // x ≠ 0})).expect
      (fun cols => z ^ (∑ i : Fin h,
        (packetColumnSupport (cols i).val).card)) =
      (((1 + 15 * z) ^ 8 - 1) / (2 ^ 32 - 1 : ℝ)) ^ h := by
  have hpow (cols : Fin h → {x : PacketColumn // x ≠ 0}) :
      z ^ (∑ i : Fin h, (packetColumnSupport (cols i).val).card) =
      ∏ i : Fin h, z ^ (packetColumnSupport (cols i).val).card :=
    (Finset.prod_pow_eq_pow_sum Finset.univ
      (fun i : Fin h => (packetColumnSupport (cols i).val).card) z).symm
  simp_rw [hpow]
  calc
    _ = ∏ i : Fin h,
        (FinPMF.uniform {x : PacketColumn // x ≠ 0}).expect
          (fun col => z ^ (packetColumnSupport col.val).card) :=
      piPMF_expect_prod
        (fun _ : Fin h => FinPMF.uniform {x : PacketColumn // x ≠ 0})
        (fun _ col => z ^ (packetColumnSupport col.val).card)
    _ = _ := by simp [uniform_nonzero_packet_gf, div_pow]

#print axioms packet_binomial8
#print axioms uniform_nonzero_packet_gf
#print axioms uniform_active_columns_gf

end Spin.Turbo
