import SpinCodes.Turbo.RS8ReferenceRoute22
import SpinCodes.UniformMixer.ScalableMapLaw

/-!
# Cell law of the Turbo uniform-value reference

The old Structured proof's `ConcreteRouteTranspose.transpose_density`
reorders an independent bit-cell product. Here each active group has
independent uniform four-bit cells; inactive groups have zero cells.
The transpose therefore changes only the indexing of the product.
-/

noncomputable section

namespace Spin.Turbo

def zeroPacketLaw22 : FinPMF Packet4 :=
  (FinPMF.uniform Unit).map (fun _ => (0 : Packet4))

theorem zeroPacketLaw22_apply (x : Packet4) :
    zeroPacketLaw22.p x = if x = 0 then 1 else 0 := by
  exact finite_map_const_apply (FinPMF.uniform Unit) 0 x

def turboReferenceCellLaw22 {R : ℕ} (S : Finset (Fin (R * 16)))
    (i : Fin (R * 16)) : FinPMF Packet4 :=
  if i ∈ S then FinPMF.uniform Packet4 else zeroPacketLaw22

theorem pi_zeroPacketLaw22_apply (word : PacketWord64) :
    (Spin.piPMF (fun _ : Fin 64 => zeroPacketLaw22)).p word =
      if word = 0 then 1 else 0 := by
  rw [Spin.piPMF_apply]
  simp only [zeroPacketLaw22_apply]
  by_cases h : word = 0
  · subst word
    simp
  · have hj : ∃ j : Fin 64, word j ≠ 0 := by
      by_contra hn
      apply h
      funext j
      by_contra hz
      exact hn ⟨j, hz⟩
    obtain ⟨j, hj⟩ := hj
    rw [Finset.prod_eq_zero (Finset.mem_univ j)]
    · simp [h]
    · simp [hj]

theorem turboReferenceRowLaw22_eq_cell_product
    (points : Fin 8 → Field16) (active : Prop)
    [Decidable active] :
    turboReferenceRowLaw22 points active =
      Spin.piPMF (fun _ : Fin 64 =>
        if active then FinPMF.uniform Packet4 else zeroPacketLaw22) := by
  by_cases h : active
  · simp only [turboReferenceRowLaw22, h, ite_true]
    exact Spin.UniformMixer.pi_uniform_eq_uniform.symm
  · apply FinPMF.ext
    funext word
    simp only [turboReferenceRowLaw22, h, ite_false]
    rw [shuffledRS8RowLaw_zero_apply, pi_zeroPacketLaw22_apply]

/-- Fixed wiring from group-major packet words to region-major values. -/
def packetTransposeEquiv22 (R : ℕ) :
    TurboOuterWord22 R ≃
      (Fin 64 → Fin (R * 16) → Packet4) where
  toFun rows j i := rows i j
  invFun regions i j := regions j i
  left_inv := by intro rows; funext i j; rfl
  right_inv := by intro regions; funext j i; rfl

def turboReferenceRegionPreLaw22 (R : ℕ)
    (S : Finset (Fin (R * 16))) :
    FinPMF (Fin 64 → Fin (R * 16) → Packet4) :=
  Spin.piPMF (fun _ : Fin 64 =>
    Spin.piPMF (fun i : Fin (R * 16) => turboReferenceCellLaw22 S i))

/-- Under the normalized reference, transposition gives independent
complete regional value rows. Each selected group position has one
uniform four-bit value in each region; all other positions are zero. -/
theorem turboReferenceArrayLaw22_transpose
    (R : ℕ) (points : Fin 8 → Field16)
    (S : Finset (Fin (R * 16))) :
    (turboReferenceArrayLaw22 R points S).map
      (packetTransposeEquiv22 R) =
      turboReferenceRegionPreLaw22 R S := by
  apply FinPMF.ext
  funext regions
  rw [FinPMF.map_equiv_apply]
  simp only [turboReferenceArrayLaw22, turboReferenceRegionPreLaw22,
    Spin.piPMF_apply, turboReferenceRowLaw22_eq_cell_product,
    turboReferenceCellLaw22]
  exact Finset.prod_comm

def turboReferenceRoutedLaw22 (R : ℕ)
    (points : Fin 8 → Field16)
    (S : Finset (Fin (R * 16))) :
    FinPMF (Fin 64 → Fin (R * 16) → Packet4) :=
  ((turboReferenceArrayLaw22 R points S).prod
      (turboRegionPermLaw22 R)).map
    (fun seed => turboRouteRegions22 seed.1 seed.2)

def turboReferenceOneRegionLaw22 (R : ℕ)
    (S : Finset (Fin (R * 16))) :
    FinPMF (Fin (R * 16) → Packet4) :=
  ((Spin.piPMF (fun i : Fin (R * 16) =>
      turboReferenceCellLaw22 S i)).prod
      (FinPMF.uniform (Equiv.Perm (Fin (R * 16))))).map
    (fun pair i => pair.1 (pair.2 i))

def turboReferenceRegionalProduct22 (R : ℕ)
    (S : Finset (Fin (R * 16))) :
    FinPMF (Fin 64 → Fin (R * 16) → Packet4) :=
  Spin.piPMF (fun _ : Fin 64 => turboReferenceOneRegionLaw22 R S)

/-- The actual independent regional permutations preserve the product
structure across regions. Each region has one shared permutation of all
group positions, while the four-bit values remain complete. -/
theorem turboReferenceRoutedLaw22_product
    (R : ℕ) (points : Fin 8 → Field16)
    (S : Finset (Fin (R * 16))) :
    turboReferenceRoutedLaw22 R points S =
      turboReferenceRegionalProduct22 R S := by
  let P : Fin 64 → FinPMF (Fin (R * 16) → Packet4) :=
    fun _ => Spin.piPMF (fun i : Fin (R * 16) =>
      turboReferenceCellLaw22 S i)
  let Q : Fin 64 → FinPMF (Equiv.Perm (Fin (R * 16))) :=
    fun _ => FinPMF.uniform (Equiv.Perm (Fin (R * 16)))
  let f : ∀ _ : Fin 64,
      ((Fin (R * 16) → Packet4) × Equiv.Perm (Fin (R * 16))) →
        (Fin (R * 16) → Packet4) :=
    fun _ pair i => pair.1 (pair.2 i)
  have hfirst :
      (turboReferenceArrayLaw22 R points S).map
        (packetTransposeEquiv22 R) = Spin.piPMF P :=
    turboReferenceArrayLaw22_transpose R points S
  calc
    turboReferenceRoutedLaw22 R points S =
        (((turboReferenceArrayLaw22 R points S).map
            (packetTransposeEquiv22 R)).prod
            (turboRegionPermLaw22 R)).map
          (fun seed j i => seed.1 j (seed.2 j i)) := by
      apply FinPMF.eq_of_expect_eq
      intro g
      simp only [turboReferenceRoutedLaw22, FinPMF.expect_map,
        FinPMF.expect_prod]
      rfl
    _ = (((Spin.piPMF P).prod (Spin.piPMF Q)).map
          piPairEquiv).map
          (fun pairs j => f j (pairs j)) := by
      rw [hfirst]
      change ((Spin.piPMF P).prod (Spin.piPMF Q)).map _ = _
      rw [FinPMF.map_comp]
      rfl
    _ = (Spin.piPMF (fun j => (P j).prod (Q j))).map
          (fun pairs j => f j (pairs j)) := by
      rw [piPMF_prod_pair]
    _ = Spin.piPMF (fun j => ((P j).prod (Q j)).map (f j)) :=
      Spin.piPMF_map _ _
    _ = turboReferenceRegionalProduct22 R S := by
      rfl

#print axioms turboReferenceRowLaw22_eq_cell_product
#print axioms turboReferenceArrayLaw22_transpose
#print axioms turboReferenceRoutedLaw22_product

end Spin.Turbo
