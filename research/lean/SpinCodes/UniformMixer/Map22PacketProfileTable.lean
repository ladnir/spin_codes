import SpinCodes.UniformMixer.Map22PacketProfileCensus
noncomputable section
namespace Spin.UniformMixer
set_option maxHeartbeats 1000000
set_option maxRecDepth 10000
def productionPacketProfileTable : List (Nat × (Fin 5 → Nat)) := [
  (1, ![16,0,0,0,0]),
  (1024, ![0,16,0,0,0]),
  (1440, ![8,0,8,0,0]),
  (107520, ![4,0,12,0,0]),
  (190464, ![0,0,16,0,0]),
  (143360, ![0,12,0,4,0]),
  (458752, ![0,10,0,6,0]),
  (890880, ![0,8,0,8,0]),
  (458752, ![0,6,0,10,0]),
  (143360, ![0,4,0,12,0]),
  (1024, ![0,0,0,16,0]),
  (430080, ![3,0,12,0,1]),
  (40320, ![6,0,8,0,2]),
  (645120, ![2,0,12,0,2]),
  (430080, ![1,0,12,0,3]),
  (140, ![12,0,0,0,4]),
  (100800, ![4,0,8,0,4]),
  (107520, ![0,0,12,0,4]),
  (448, ![10,0,0,0,6]),
  (40320, ![2,0,8,0,6]),
  (870, ![8,0,0,0,8]),
  (1440, ![0,0,8,0,8]),
  (448, ![6,0,0,0,10]),
  (140, ![4,0,0,0,12]),
  (1, ![0,0,0,0,16])
]
theorem productionPacketProfile_table_sum {M : Type*} [AddCommMonoid M]
    (f : (Fin 5 → Nat) → M) :
    (∑ q : State22, f (packetProfile64 (fun j => coordMap j q))) =
      (productionPacketProfileTable.map (fun e => e.1 • f e.2)).sum := by
  rw [productionPacketProfile_census]
  simp only [profileEvenClassSum, Fin.sum_univ_succ]
  simp [profileQuad4WeightCounts, profileQuad3WeightCounts, profileQuad2WeightCounts,
    productionPacketProfileTable, profileOddShape, profileEvenShape,
    smul_smul]
  <;> abel
theorem productionPacketProfileTable_total :
    (productionPacketProfileTable.map Prod.fst).sum = 2^22 := by decide
theorem productionPacketProfileTable_length : productionPacketProfileTable.length = 25 := by decide
theorem productionPacketProfileTable_counts_checked : ∀ i : Fin 25,
    (productionPacketProfileTable.map (fun e =>
      if e.2 = (productionPacketProfileTable.getD i.val (0,0)).2 then e.1 else 0)).sum =
      (productionPacketProfileTable.getD i.val (0,0)).1 := by decide
theorem productionPacketProfile_table_count (i : Fin 25) :
    (Finset.univ.filter (fun q : State22 => packetProfile64 (fun j => coordMap j q) =
      (productionPacketProfileTable.getD i.val (0,0)).2)).card =
      (productionPacketProfileTable.getD i.val (0,0)).1 := by
  rw [Finset.card_filter, productionPacketProfile_table_sum (fun p =>
    if p = (productionPacketProfileTable.getD i.val (0,0)).2 then (1 : Nat) else 0)]
  simpa only [nsmul_eq_mul, mul_ite, Nat.mul_one, Nat.mul_zero, Nat.cast_id] using
    productionPacketProfileTable_counts_checked i
#print axioms productionPacketProfile_table_sum
#print axioms productionPacketProfile_table_count
end Spin.UniformMixer
