import SpinCodes.Turbo.ColumnAdjointField
import SpinCodes.Turbo.RS8SpectrumNumeric
import SpinCodes.Framework

/-!
# The sampled packet moment of one RS8 group

The eight RS rows determine eight 32-bit columns. Independent nonzero
GF(2^32) coefficients act on those columns through their binary adjoints.
This file connects that sampled law to the actual RS column support; it
does not replace the selected setup by independent codewords.
-/

noncomputable section

namespace Spin.Turbo

abbrev Field16 := GaloisField 2 4

noncomputable instance : Fintype Field16 := Fintype.ofFinite Field16
noncomputable instance : DecidableEq Field16 := Classical.decEq Field16

theorem field16_card : Fintype.card Field16 = 16 := by
  rw [Fintype.card_eq_nat_card]
  simpa [Field16] using (GaloisField.card 2 4 (by decide))

/-- An additive identification of one eight-symbol RS column with 32
binary coordinates. The expected packet moment is independent of the
chosen basis. -/
def rs8ColumnBits : (Fin 8 → Field16) ≃ₗ[ZMod 2] Bits32 :=
  LinearEquiv.ofFinrankEq (Fin 8 → Field16) Bits32 (by
    rw [Module.finrank_pi_fintype, Module.finrank_fin_fun]
    rw [GaloisField.finrank 2 (n := 4) (by decide)]
    norm_num)

theorem rs8ColumnBits_ne_zero {col : Fin 8 → Field16} :
    rs8ColumnBits col ≠ 0 ↔ col ≠ 0 := by
  constructor
  · intro h hc
    apply h
    simp [hc]
  · intro h hc
    apply h
    exact rs8ColumnBits.injective (by simpa using hc)

/-- The 32-bit view of an actual RS output column. -/
def rs8Bits (points : Fin 8 → Field16)
    (rows : RS8Group Field16) (j : Fin 8) : Bits32 :=
  rs8ColumnBits (rs8Column points rows j)

theorem rs8Bits_ne_zero_iff (points : Fin 8 → Field16)
    (rows : RS8Group Field16) (j : Fin 8) :
    rs8Bits points rows j ≠ 0 ↔ j ∈ rs8Support points rows := by
  simp [rs8Bits, rs8Support]

/-- Packet count after the eight independent sampled column adjoints,
before the within-group packet shuffle. -/
def rs8SampledPacketWeight (points : Fin 8 → Field16)
    (rows : RS8Group Field16)
    (coeff : Fin 8 → {c : Field32 // c ≠ 0}) : ℕ :=
  ∑ j : Fin 8, (packetColumnSupport (bits32PacketEquiv
    (coordinateAdjoint field32Coordinates (coeff j).val
      (rs8Bits points rows j)))).card

private theorem coordinateAdjoint_zero (c : Field32) :
    coordinateAdjoint field32Coordinates c 0 = 0 := by
  simp [coordinateAdjoint]

/-- A fixed RS8 group with `h` active columns has the sampled packet
moment `P(z)^h/(2^32-1)^h`. No independence across different messages is
assumed; only the eight setup coefficients within this group are sampled
independently. -/
theorem rs8_fixed_group_packet_gf
    (points : Fin 8 → Field16) (rows : RS8Group Field16)
    (z : ℝ) :
    (piPMF (fun _ : Fin 8 =>
      FinPMF.uniform {c : Field32 // c ≠ 0})).expect
      (fun coeff => z ^ rs8SampledPacketWeight points rows coeff) =
      (((1 + 15 * z) ^ 8 - 1) / (2 ^ 32 - 1 : ℝ)) ^
        (rs8Support points rows).card := by
  let q : ℝ := ((1 + 15 * z) ^ 8 - 1) / (2 ^ 32 - 1)
  have hpow (coeff : Fin 8 → {c : Field32 // c ≠ 0}) :
      z ^ rs8SampledPacketWeight points rows coeff =
      ∏ j : Fin 8, z ^ (packetColumnSupport (bits32PacketEquiv
        (coordinateAdjoint field32Coordinates (coeff j).val
          (rs8Bits points rows j)))).card := by
    exact (Finset.prod_pow_eq_pow_sum Finset.univ
      (fun j : Fin 8 => (packetColumnSupport (bits32PacketEquiv
        (coordinateAdjoint field32Coordinates (coeff j).val
          (rs8Bits points rows j)))).card) z).symm
  simp_rw [hpow]
  calc
    _ = ∏ j : Fin 8,
        (FinPMF.uniform {c : Field32 // c ≠ 0}).expect
          (fun c => z ^ (packetColumnSupport (bits32PacketEquiv
            (coordinateAdjoint field32Coordinates c.val
              (rs8Bits points rows j)))).card) :=
      piPMF_expect_prod
        (fun _ : Fin 8 => FinPMF.uniform {c : Field32 // c ≠ 0})
        (fun j c => z ^ (packetColumnSupport (bits32PacketEquiv
          (coordinateAdjoint field32Coordinates c.val
            (rs8Bits points rows j)))).card)
    _ = ∏ j : Fin 8, if j ∈ rs8Support points rows then q else 1 := by
      apply Finset.prod_congr rfl
      intro j _
      by_cases hj : j ∈ rs8Support points rows
      · have hx := (rs8Bits_ne_zero_iff points rows j).mpr hj
        simpa [hj, q] using
          field32_adjoint_packet_gf (rs8Bits points rows j) hx z
      · have hz : rs8Bits points rows j = 0 := by
          by_contra hx
          exact hj ((rs8Bits_ne_zero_iff points rows j).mp hx)
        have hconst (c : {c : Field32 // c ≠ 0}) :
            z ^ (packetColumnSupport (bits32PacketEquiv
              (coordinateAdjoint field32Coordinates c.val
                (rs8Bits points rows j)))).card = 1 := by
          rw [hz, coordinateAdjoint_zero, bits32PacketEquiv_zero]
          simp [packetColumnSupport]
        calc
          _ = (FinPMF.uniform {c : Field32 // c ≠ 0}).expect
              (fun _ => (1 : ℝ)) := by
              unfold FinPMF.expect
              apply Finset.sum_congr rfl
              intro c _
              dsimp only
              rw [hconst c]
          _ = 1 := FinPMF.expect_const _ _
          _ = (if j ∈ rs8Support points rows then q else 1) := by
            simp only [hj, ↓reduceIte]
    _ = q ^ (rs8Support points rows).card := by
      rw [Finset.prod_ite_mem_eq, Finset.prod_const]
    _ = _ := rfl

/-- Group the sampled packet moment by the actual number of nonzero RS8
columns. This is an exact sum over messages, including the zero message;
the latter will be removed when forming the expected nonzero-message count. -/
theorem rs8_group_packet_gf_by_column_weight
    (points : Fin 8 → Field16) (z : ℝ) :
    (∑ rows : RS8Group Field16,
      (piPMF (fun _ : Fin 8 =>
        FinPMF.uniform {c : Field32 // c ≠ 0})).expect
        (fun coeff => z ^ rs8SampledPacketWeight points rows coeff)) =
      ∑ h ∈ Finset.range 9,
        ((Finset.univ.filter (fun rows : RS8Group Field16 =>
          (rs8Support points rows).card = h)).card : ℝ) *
          (((1 + 15 * z) ^ 8 - 1) / (2 ^ 32 - 1 : ℝ)) ^ h := by
  let q : ℝ := ((1 + 15 * z) ^ 8 - 1) / (2 ^ 32 - 1)
  have hmaps (rows : RS8Group Field16) :
      (rs8Support points rows).card ∈ Finset.range 9 := by
    simp only [Finset.mem_range]
    have hle := Finset.card_le_card
      (Finset.filter_subset (p := fun j => rs8Column points rows j ≠ 0)
        (Finset.univ : Finset (Fin 8)))
    simp only [Finset.card_univ, Fintype.card_fin] at hle
    change (rs8Support points rows).card ≤ 8 at hle
    omega
  calc
    _ = ∑ rows : RS8Group Field16,
        q ^ (rs8Support points rows).card := by
      apply Finset.sum_congr rfl
      intro rows _
      exact rs8_fixed_group_packet_gf points rows z
    _ = ∑ h ∈ Finset.range 9,
        ∑ rows ∈ (Finset.univ.filter (fun rows : RS8Group Field16 =>
          (rs8Support points rows).card = h)),
          q ^ (rs8Support points rows).card :=
      (Finset.sum_fiberwise_of_maps_to
        (fun rows _ => hmaps rows) _).symm
    _ = _ := by
      apply Finset.sum_congr rfl
      intro h _
      rw [Finset.sum_congr rfl (fun rows hr => by
        rw [(Finset.mem_filter.mp hr).2]), Finset.sum_const, nsmul_eq_mul]

private theorem rs8Support_zero (points : Fin 8 → Field16) :
    rs8Support points (0 : RS8Group Field16) = ∅ := by
  have hcol (j : Fin 8) :
      rs8Column points (0 : RS8Group Field16) j = 0 := by
    funext r
    change (0 : Polynomial Field16).eval (points j) = 0
    simp
  ext j
  simp [rs8Support, hcol]

/-- The nonzero-message expected packet generating function is the sum
over the actual RS8 column-support classes. The separate checked `A_h`
theorem gives each class count its closed MDS expression. -/
theorem rs8_nonzero_group_packet_gf
    (points : Fin 8 → Field16) (hpoints : Function.Injective points)
    (z : ℝ) :
    (∑ rows ∈ nonzeroMsgs (RS8Group Field16),
      (piPMF (fun _ : Fin 8 =>
        FinPMF.uniform {c : Field32 // c ≠ 0})).expect
        (fun coeff => z ^ rs8SampledPacketWeight points rows coeff)) =
      ∑ h ∈ Finset.Icc 5 8,
        (rs8ColumnWeightCount points h : ℝ) *
          (((1 + 15 * z) ^ 8 - 1) / (2 ^ 32 - 1 : ℝ)) ^ h := by
  let q : ℝ := ((1 + 15 * z) ^ 8 - 1) / (2 ^ 32 - 1)
  have hupper (rows : RS8Group Field16) :
      (rs8Support points rows).card ≤ 8 := by
    have hle := Finset.card_le_card
      (Finset.filter_subset (p := fun j => rs8Column points rows j ≠ 0)
        (Finset.univ : Finset (Fin 8)))
    simpa only [rs8Support, Finset.card_univ, Fintype.card_fin] using hle
  have hmaps (rows : RS8Group Field16)
      (hr : rows ∈ nonzeroMsgs (RS8Group Field16)) :
      (rs8Support points rows).card ∈ Finset.Icc 5 8 := by
    have hn : rows ≠ 0 := (Finset.mem_filter.mp hr).2
    exact Finset.mem_Icc.mpr
      ⟨rs8Support_ge_five points hpoints rows hn, hupper rows⟩
  have hclass (h : ℕ) (hh : h ∈ Finset.Icc 5 8) :
      (nonzeroMsgs (RS8Group Field16)).filter
        (fun rows => (rs8Support points rows).card = h) =
      Finset.univ.filter
        (fun rows : RS8Group Field16 =>
          (rs8Support points rows).card = h) := by
    ext rows
    simp only [nonzeroMsgs, Finset.mem_filter, Finset.mem_univ,
      true_and]
    constructor
    · exact And.right
    · intro hr
      refine ⟨?_, hr⟩
      intro hzero
      subst rows
      have hz : (rs8Support points (0 : RS8Group Field16)).card = 0 := by
        rw [rs8Support_zero]
        rfl
      have hlo := (Finset.mem_Icc.mp hh).1
      omega
  calc
    _ = ∑ rows ∈ nonzeroMsgs (RS8Group Field16),
        q ^ (rs8Support points rows).card := by
      apply Finset.sum_congr rfl
      intro rows _
      exact rs8_fixed_group_packet_gf points rows z
    _ = ∑ h ∈ Finset.Icc 5 8,
        ∑ rows ∈ (nonzeroMsgs (RS8Group Field16)).filter
          (fun rows => (rs8Support points rows).card = h),
          q ^ (rs8Support points rows).card :=
      (Finset.sum_fiberwise_of_maps_to hmaps _).symm
    _ = ∑ h ∈ Finset.Icc 5 8,
        (((nonzeroMsgs (RS8Group Field16)).filter
          (fun rows => (rs8Support points rows).card = h)).card : ℝ) *
          q ^ h := by
      apply Finset.sum_congr rfl
      intro h _
      rw [Finset.sum_congr rfl (fun rows hr => by
        rw [(Finset.mem_filter.mp hr).2]), Finset.sum_const, nsmul_eq_mul]
    _ = _ := by
      apply Finset.sum_congr rfl
      intro h hh
      rw [hclass h hh]
      have hc : ((Finset.univ.filter
          (fun rows : RS8Group Field16 =>
            (rs8Support points rows).card = h)).card : ℝ) =
          (rs8ColumnWeightCount points h : ℝ) := by
        exact_mod_cast (rs8ColumnWeightCount_eq_card points h).symm
      rw [hc]

/-- The paper's expected outer packet polynomial, with its exact MDS
coefficient formula and the actual independently sampled adjoint maps. -/
theorem rs8_nonzero_group_packet_gf_closed
    (points : Fin 8 → Field16) (hpoints : Function.Injective points)
    (z : ℝ) :
    (∑ rows ∈ nonzeroMsgs (RS8Group Field16),
      (piPMF (fun _ : Fin 8 =>
        FinPMF.uniform {c : Field32 // c ≠ 0})).expect
        (fun coeff => z ^ rs8SampledPacketWeight points rows coeff)) =
      ∑ h ∈ Finset.Icc 5 8,
        (((Nat.choose 8 h : ℤ) *
          ∑ j ∈ Finset.range (h-4),
            (-1 : ℤ) ^ j * (h.choose j : ℤ) *
              (((2^32 : ℕ) ^ (h-4-j) - 1 : ℕ) : ℤ) : ℤ) : ℝ) *
          (((1 + 15 * z) ^ 8 - 1) / (2 ^ 32 - 1 : ℝ)) ^ h := by
  rw [rs8_nonzero_group_packet_gf points hpoints z]
  apply Finset.sum_congr rfl
  intro h hh
  rw [rs8ColumnWeightCount_MDS field16_card points hpoints h
    (Finset.mem_Icc.mp hh).1]

private theorem packet_support_pos {p : PacketColumn} (hp : p ≠ 0) :
    1 ≤ (packetColumnSupport p).card := by
  have hex : ∃ i : Fin 8, p i ≠ 0 := by
    by_contra hn
    apply hp
    funext i
    by_contra hi
    exact hn ⟨i, hi⟩
  obtain ⟨i, hi⟩ := hex
  exact Finset.card_pos.mpr
    ⟨i, Finset.mem_filter.mpr ⟨Finset.mem_univ _, hi⟩⟩

/-- Every active RS output column remains nonzero after its sampled
adjoint, and therefore contributes at least one active packet. -/
theorem rs8_sampled_packet_weight_ge_support
    (points : Fin 8 → Field16) (rows : RS8Group Field16)
    (coeff : Fin 8 → {c : Field32 // c ≠ 0}) :
    (rs8Support points rows).card ≤
      rs8SampledPacketWeight points rows coeff := by
  have hlocal (j : Fin 8) :
      (if j ∈ rs8Support points rows then 1 else 0) ≤
        (packetColumnSupport (bits32PacketEquiv
          (coordinateAdjoint field32Coordinates (coeff j).val
            (rs8Bits points rows j)))).card := by
    by_cases hj : j ∈ rs8Support points rows
    · have hx := (rs8Bits_ne_zero_iff points rows j).mpr hj
      have ha : coordinateAdjoint field32Coordinates (coeff j).val
          (rs8Bits points rows j) ≠ 0 :=
        adjoint_ne_zero field32Coordinates.toAddEquiv binaryPair
          binaryPair_nondegenerate binaryPair_zero_left
          (coordinateAdjoint field32Coordinates)
          (coordinateAdjoint_pairing field32Coordinates) (coeff j).2 hx
      have hp : bits32PacketEquiv
          (coordinateAdjoint field32Coordinates (coeff j).val
            (rs8Bits points rows j)) ≠ 0 := by
        intro hz
        apply ha
        exact bits32PacketEquiv.injective
          (hz.trans bits32PacketEquiv_zero.symm)
      simpa only [hj, ↓reduceIte] using packet_support_pos hp
    · simp only [hj, ↓reduceIte, Nat.zero_le]
  have hsum := Finset.sum_le_sum
    (fun j (_ : j ∈ (Finset.univ : Finset (Fin 8))) => hlocal j)
  simpa only [Finset.sum_ite, Finset.sum_const_zero,
    Finset.card_eq_sum_ones, rs8SampledPacketWeight,
    Nat.add_zero, Finset.filter_univ_mem] using hsum

/-- The sampled packet count of a nonzero group lies in the claimed
range `5..64`, independently of the sampled nonzero coefficients. -/
theorem rs8_sampled_packet_weight_range
    (points : Fin 8 → Field16) (hpoints : Function.Injective points)
    (rows : RS8Group Field16) (hrows : rows ≠ 0)
    (coeff : Fin 8 → {c : Field32 // c ≠ 0}) :
    5 ≤ rs8SampledPacketWeight points rows coeff ∧
      rs8SampledPacketWeight points rows coeff ≤ 64 := by
  constructor
  · exact (rs8Support_ge_five points hpoints rows hrows).trans
      (rs8_sampled_packet_weight_ge_support points rows coeff)
  · unfold rs8SampledPacketWeight
    calc
      _ ≤ ∑ _j : Fin 8, (8 : ℕ) := by
        apply Finset.sum_le_sum
        intro j _
        exact Finset.card_le_card (Finset.filter_subset
          (p := fun i => bits32PacketEquiv
            (coordinateAdjoint field32Coordinates (coeff j).val
              (rs8Bits points rows j)) i ≠ 0)
            (Finset.univ : Finset (Fin 8)))
      _ = 64 := by norm_num

/-- `C_v`: the expected number of nonzero RS8 group messages producing
exactly `v` active packets under the ideal independently sampled column
adjoints. Different messages share the same sampled setup. -/
def rs8ExpectedPacketCount (points : Fin 8 → Field16) (v : ℕ) : ℝ :=
  ∑ rows ∈ nonzeroMsgs (RS8Group Field16),
    (piPMF (fun _ : Fin 8 =>
      FinPMF.uniform {c : Field32 // c ≠ 0})).prob
      (fun coeff => rs8SampledPacketWeight points rows coeff = v)

/-- The coefficientwise expected counts `C_v` have exactly the paper's
generating polynomial. The left side is the expected message count, not
the packet-weight law of one fixed message or a realized code. -/
theorem rs8_expected_packet_count_gf
    (points : Fin 8 → Field16) (hpoints : Function.Injective points)
    (z : ℝ) :
    (∑ v ∈ Finset.Icc 5 64, rs8ExpectedPacketCount points v * z ^ v) =
      ∑ h ∈ Finset.Icc 5 8,
        (((Nat.choose 8 h : ℤ) *
          ∑ j ∈ Finset.range (h-4),
            (-1 : ℤ) ^ j * (h.choose j : ℤ) *
              (((2^32 : ℕ) ^ (h-4-j) - 1 : ℕ) : ℤ) : ℤ) : ℝ) *
          (((1 + 15 * z) ^ 8 - 1) / (2 ^ 32 - 1 : ℝ)) ^ h := by
  rw [← rs8_nonzero_group_packet_gf_closed points hpoints z]
  unfold rs8ExpectedPacketCount
  calc
    _ = ∑ v ∈ Finset.Icc 5 64,
        ∑ rows ∈ nonzeroMsgs (RS8Group Field16),
          (piPMF (fun _ : Fin 8 =>
            FinPMF.uniform {c : Field32 // c ≠ 0})).prob
              (fun coeff => rs8SampledPacketWeight points rows coeff = v) *
                z ^ v := by
      apply Finset.sum_congr rfl
      intro v _
      rw [Finset.sum_mul]
    _ = ∑ rows ∈ nonzeroMsgs (RS8Group Field16),
        ∑ v ∈ Finset.Icc 5 64,
          (piPMF (fun _ : Fin 8 =>
            FinPMF.uniform {c : Field32 // c ≠ 0})).prob
              (fun coeff => rs8SampledPacketWeight points rows coeff = v) *
                z ^ v := by rw [Finset.sum_comm]
    _ = _ := by
      apply Finset.sum_congr rfl
      intro rows hr
      have hn : rows ≠ 0 := (Finset.mem_filter.mp hr).2
      exact (FinPMF.expect_comp
        (piPMF (fun _ : Fin 8 =>
          FinPMF.uniform {c : Field32 // c ≠ 0}))
        (rs8SampledPacketWeight points rows) (Finset.Icc 5 64)
        (fun coeff => Finset.mem_Icc.mpr
          (rs8_sampled_packet_weight_range points hpoints rows hn coeff))
        (fun v => z ^ v)).symm

/-- Summing the expected packet counts returns exactly one contribution
from each nonzero RS8 message. -/
theorem rs8_expected_packet_count_total
    (points : Fin 8 → Field16) (hpoints : Function.Injective points) :
    (∑ v ∈ Finset.Icc 5 64, rs8ExpectedPacketCount points v) =
      ((nonzeroMsgs (RS8Group Field16)).card : ℝ) := by
  unfold rs8ExpectedPacketCount
  calc
    _ = ∑ rows ∈ nonzeroMsgs (RS8Group Field16),
        ∑ v ∈ Finset.Icc 5 64,
          (piPMF (fun _ : Fin 8 =>
            FinPMF.uniform {c : Field32 // c ≠ 0})).prob
              (fun coeff => rs8SampledPacketWeight points rows coeff = v) := by
      rw [Finset.sum_comm]
    _ = ∑ _rows ∈ nonzeroMsgs (RS8Group Field16), (1 : ℝ) := by
      apply Finset.sum_congr rfl
      intro rows hr
      have hn : rows ≠ 0 := (Finset.mem_filter.mp hr).2
      have h := FinPMF.expect_comp
        (piPMF (fun _ : Fin 8 =>
          FinPMF.uniform {c : Field32 // c ≠ 0}))
        (rs8SampledPacketWeight points rows) (Finset.Icc 5 64)
        (fun coeff => Finset.mem_Icc.mpr
          (rs8_sampled_packet_weight_range points hpoints rows hn coeff))
        (fun _ => (1 : ℝ))
      simpa only [FinPMF.expect_const, mul_one] using h.symm
    _ = _ := by simp

/-- The normalization in the paper: one group has `2^128-1` nonzero
messages, so the expected packet counts sum to that integer. -/
theorem rs8_expected_packet_count_total_numeric
    (points : Fin 8 → Field16) (hpoints : Function.Injective points) :
    (∑ v ∈ Finset.Icc 5 64, rs8ExpectedPacketCount points v) =
      (2 ^ 128 - 1 : ℝ) := by
  rw [rs8_expected_packet_count_total points hpoints]
  have hnonzero : nonzeroMsgs (RS8Group Field16) =
      (Finset.univ : Finset (RS8Group Field16)).erase 0 := by
    ext rows
    simp [nonzeroMsgs]
  rw [hnonzero, Finset.card_erase_of_mem (Finset.mem_univ _),
    Finset.card_univ, rs8Group_card, field16_card]
  norm_num

#print axioms rs8ColumnBits_ne_zero
#print axioms rs8Bits_ne_zero_iff
#print axioms rs8_fixed_group_packet_gf
#print axioms rs8_group_packet_gf_by_column_weight
#print axioms rs8_nonzero_group_packet_gf
#print axioms rs8_nonzero_group_packet_gf_closed
#print axioms rs8_sampled_packet_weight_range
#print axioms rs8_expected_packet_count_gf
#print axioms rs8_expected_packet_count_total_numeric

end Spin.Turbo
