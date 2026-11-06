import SpinCodes.UniformMixer.RegionalEncoder

/-!
# Chronological regional kernel

The encoder processes one 64-bit input group per step. This file relates
the product of those step kernels to the empty-gap and singleton-input
decomposition used by the regional proof.
-/

noncomputable section

namespace Spin.UniformMixer

def inputListKernel (z : ℝ) :
    List (Fin 64 → ZMod 2) → State22 → State22 → ℝ
  | [] => identity
  | x :: xs => compose (stepKernel z x) (inputListKernel z xs)

/-- The 64-bit input groups in their actual encoder order. -/
def chronologicalInputList {R : ℕ}
    (T : Finset (Fin (64 * R))) : List (Fin 64 → ZMod 2) :=
  List.ofFn (fun i : Fin R => regionInput T i)

theorem chronologicalInputList_length {R : ℕ}
    (T : Finset (Fin (64 * R))) :
    (chronologicalInputList T).length = R := by
  simp [chronologicalInputList]

theorem chronologicalInputList_empty (R : ℕ) :
    chronologicalInputList (∅ : Finset (Fin (64 * R))) =
      List.replicate R 0 := by
  have hz : ∀ i : Fin R,
      regionInput (∅ : Finset (Fin (64 * R))) i = 0 := by
    intro i
    funext j
    simp [regionInput]
  simp [chronologicalInputList, hz]

/-- Full-state kernel induced by the routed support and the encoder's
chronological step law. -/
def chronologicalRegionKernel {R : ℕ} (z : ℝ)
    (T : Finset (Fin (64 * R))) : State22 → State22 → ℝ :=
  inputListKernel z (chronologicalInputList T)

theorem inputListKernel_append (z : ℝ)
    (xs ys : List (Fin 64 → ZMod 2)) :
    inputListKernel z (xs ++ ys) =
      compose (inputListKernel z xs) (inputListKernel z ys) := by
  induction xs with
  | nil =>
      funext q r
      exact (compose_identity_left (inputListKernel z ys) q r).symm
  | cons x tail ih =>
      funext q r
      change compose (stepKernel z x) (inputListKernel z (tail ++ ys)) q r =
        compose (compose (stepKernel z x) (inputListKernel z tail))
          (inputListKernel z ys) q r
      rw [ih]
      exact (compose_assoc (stepKernel z x)
        (inputListKernel z tail) (inputListKernel z ys) q r).symm

theorem inputListKernel_empty_run (z : ℝ) (g : ℕ) :
    inputListKernel z (List.replicate g (0 : Fin 64 → ZMod 2)) =
      emptyGap z (coordinateWeight coordMap) g := by
  induction g with
  | zero => rfl
  | succ g ih =>
      change compose (stepKernel z 0)
        (inputListKernel z (List.replicate g 0)) =
        emptyGap z (coordinateWeight coordMap) (g + 1)
      rw [stepKernel_empty, ih]
      have hOne : emptyGap z (coordinateWeight coordMap) 1 =
          emptyStep z (coordinateWeight coordMap) := by
        funext q r
        exact emptyGap_one z (coordinateWeight coordMap) q r
      calc
        compose (emptyStep z (coordinateWeight coordMap))
            (emptyGap z (coordinateWeight coordMap) g) =
          compose (emptyGap z (coordinateWeight coordMap) 1)
            (emptyGap z (coordinateWeight coordMap) g) := by rw [hOne]
        _ = emptyGap z (coordinateWeight coordMap) (g + 1) := by
          rw [← emptyGap_add]
          simp [Nat.add_comm]

def markedInputList :
    List (ℕ × Fin 64) → ℕ → List (Fin 64 → ZMod 2)
  | [], last => List.replicate last 0
  | (g, j) :: tail, last =>
      List.replicate g 0 ++ singletonInput j :: markedInputList tail last

theorem inputListKernel_markedInputList (z : ℝ)
    (marks : List (ℕ × Fin 64)) (last : ℕ) :
    inputListKernel z (markedInputList marks last) =
      markedRegionKernel z marks last := by
  induction marks with
  | nil => exact inputListKernel_empty_run z last
  | cons p tail ih =>
      rcases p with ⟨g, j⟩
      change inputListKernel z
        (List.replicate g 0 ++ singletonInput j ::
          markedInputList tail last) =
        compose (emptyGap z (coordinateWeight coordMap) g)
          (compose (productionImpulse z j)
            (markedRegionKernel z tail last))
      rw [inputListKernel_append, inputListKernel_empty_run]
      change compose (emptyGap z (coordinateWeight coordMap) g)
        (compose (stepKernel z (singletonInput j))
          (inputListKernel z (markedInputList tail last))) = _
      rw [stepKernel_singleton, ih]

theorem markedInputList_support_length {R prev : ℕ}
    {xs : List (Fin (64 * R))} (h : ValidSupport R prev xs) :
    (markedInputList (supportMarks R xs prev)
      (supportLastGap R xs prev)).length = R - prev := by
  induction h with
  | nil hprev => simp [supportMarks, supportLastGap, markedInputList]
  | cons hgap hinterior htail ih =>
      simp only [supportMarks, supportLastGap, markedInputList,
        List.length_append, List.length_replicate, List.length_cons, ih]
      omega

theorem markedInputList_goodSupport_length {R : ℕ} (hR : 0 < R)
    (T : Finset (Fin (64 * R))) (hgood : ¬ badRegionalSites T) :
    (markedInputList (supportMarks R (T.sort (· ≤ ·)) 0)
      (supportLastGap R (T.sort (· ≤ ·)) 0)).length = R := by
  simpa using markedInputList_support_length
    (goodRegionalSites_validSupport hR T hgood)

theorem empty_support_trace (R : ℕ) :
    chronologicalInputList (∅ : Finset (Fin (64 * R))) =
      markedInputList
        (supportMarks R ((∅ : Finset (Fin (64 * R))).sort (· ≤ ·)) 0)
        (supportLastGap R ((∅ : Finset (Fin (64 * R))).sort (· ≤ ·)) 0) := by
  simp [chronologicalInputList_empty, supportMarks, supportLastGap,
    markedInputList]

/-- The regional upper bound for the actual chronological encoder,
conditional on the remaining input-list identification. -/
theorem good_chronologicalRegion_coarse_upper {R : ℕ} (hR : 0 < R)
    (T : Finset (Fin (64 * R))) (hgood : ¬ badRegionalSites T)
    (htrace : chronologicalInputList T =
      markedInputList (supportMarks R (T.sort (· ≤ ·)) 0)
        (supportLastGap R (T.sort (· ≤ ·)) 0))
    (θ : ℝ) (hθ : 0 ≤ θ) (q r : State22) :
    let L : ℝ := 64 * R
    let C := ((64 * θ + 64 * θ ^ 2) / L) / Real.exp (-θ)
    chronologicalRegionKernel (Real.exp (-(θ / L))) T q r ≤
      (1 + C) ^ (T.card + 1) *
        (liftM * coarseRegionProduct
          (fun g => Real.exp (-((θ / L) * liveWeightAverage coordMap)) ^ g)
          (supportMarks R (T.sort (· ≤ ·)) 0)
          (supportLastGap R (T.sort (· ≤ ·)) 0) * projectM) q r := by
  simpa [chronologicalRegionKernel, supportRegionKernel, htrace,
    inputListKernel_markedInputList] using
    good_supportRegion_coarse_upper hR T hgood θ hθ q r

end Spin.UniformMixer
