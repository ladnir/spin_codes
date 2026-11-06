import SpinCodes.UniformMixer.RegionalEncoder

/-!
# Chronological regional kernel

The encoder processes one 64-bit input group per step. This file relates
the product of those step kernels to the empty-gap and singleton-input
decomposition used by the regional proof.
-/

noncomputable section

namespace Spin.UniformMixer

/-- Input group obtained by reading one step from an ordered support.
The definition works for any list; ordering is used only in the trace proof. -/
def stepFromSupport {R : ℕ} (xs : List (Fin (64 * R)))
    (k : ℕ) : Fin 64 → ZMod 2 :=
  fun j => if ∃ x ∈ xs, x.val = 64 * k + j.val then 1 else 0

theorem validSupport_steps_gt {R prev : ℕ}
    {xs : List (Fin (64 * R))} (h : ValidSupport R prev xs) :
    ∀ x ∈ xs, prev < x.val / 64 := by
  induction h with
  | nil _ => simp
  | cons hgap _ htail ih =>
      intro x hx
      simp only [List.mem_cons] at hx
      rcases hx with rfl | hx
      · exact hgap
      · have ht := ih x hx
        omega

theorem validSupport_unique_step {R prev : ℕ}
    {xs : List (Fin (64 * R))} (h : ValidSupport R prev xs)
    {x y : Fin (64 * R)} (hx : x ∈ xs) (hy : y ∈ xs)
    (hstep : x.val / 64 = y.val / 64) : x = y := by
  induction h with
  | nil _ => simp at hx
  | cons hgap _ htail ih =>
      simp only [List.mem_cons] at hx hy
      rcases hx with rfl | hx
      · rcases hy with rfl | hy
        · rfl
        · have ht := validSupport_steps_gt htail y hy
          omega
      · rcases hy with rfl | hy
        · have ht := validSupport_steps_gt htail x hx
          omega
        · exact ih hx hy

theorem stepFromSupport_zero_of_no_step {R : ℕ}
    (xs : List (Fin (64 * R))) (k : ℕ)
    (hnone : ∀ x ∈ xs, x.val / 64 ≠ k) :
    stepFromSupport xs k = 0 := by
  funext j
  have hnoBit : ¬ ∃ x ∈ xs, x.val = 64 * k + j.val := by
    rintro ⟨x, hx, heq⟩
    have hdiv : x.val / 64 = k := by
      rw [heq]
      have hj := j.isLt
      omega
    exact (hnone x hx) hdiv
  simp [stepFromSupport, hnoBit]

theorem stepFromSupport_active {R prev : ℕ}
    {xs : List (Fin (64 * R))} (h : ValidSupport R prev xs)
    (x : Fin (64 * R)) (hx : x ∈ xs) :
    stepFromSupport xs (x.val / 64) =
      singletonInput (regionLocal x) := by
  funext j
  by_cases hj : j = regionLocal x
  · subst j
    have hval : x.val = 64 * (x.val / 64) + (regionLocal x).val := by
      simp only [regionLocal]
      exact (Nat.div_add_mod x.val 64).symm
    have hwitness : ∃ y ∈ xs,
        y.val = 64 * (x.val / 64) + (regionLocal x).val :=
      ⟨x, hx, hval⟩
    simp [stepFromSupport, singletonInput, hwitness]
  · have hnoBit : ¬ ∃ y ∈ xs,
        y.val = 64 * (x.val / 64) + j.val := by
      rintro ⟨y, hy, hval⟩
      have hydiv : y.val / 64 = x.val / 64 := by
        rw [hval]
        have hjlt := j.isLt
        omega
      have hxy := validSupport_unique_step h hx hy hydiv.symm
      subst y
      have hmod := Nat.mod_lt x.val (by omega : 0 < 64)
      have hrem : j.val = x.val % 64 := by
        have hxdiv := Nat.div_add_mod x.val 64
        omega
      exact hj (Fin.ext hrem)
    simp [stepFromSupport, singletonInput, hj, hnoBit]

theorem stepFromSupport_cons_ne {R : ℕ}
    (x : Fin (64 * R)) (xs : List (Fin (64 * R))) (k : ℕ)
    (hne : x.val / 64 ≠ k) :
    stepFromSupport (x :: xs) k = stepFromSupport xs k := by
  funext j
  have hx : x.val ≠ 64 * k + j.val := by
    intro heq
    have hdiv : x.val / 64 = k := by
      rw [heq]
      have hj := j.isLt
      omega
    exact hne hdiv
  have heq :
      (∃ y ∈ x :: xs, y.val = 64 * k + j.val) ↔
        (∃ y ∈ xs, y.val = 64 * k + j.val) := by
    constructor
    · rintro ⟨y, hy, hval⟩
      simp only [List.mem_cons] at hy
      rcases hy with rfl | hy
      · exact (hx hval).elim
      · exact ⟨y, hy, hval⟩
    · rintro ⟨y, hy, hval⟩
      exact ⟨y, List.mem_cons_of_mem x hy, hval⟩
  change (if ∃ y ∈ x :: xs, y.val = 64 * k + j.val then
      (1 : ZMod 2) else 0) =
    (if ∃ y ∈ xs, y.val = 64 * k + j.val then 1 else 0)
  simp only [heq]

theorem stepFromSupport_zero_of_before {R : ℕ}
    (xs : List (Fin (64 * R))) (k : ℕ)
    (hbefore : ∀ x ∈ xs, k < x.val / 64) :
    stepFromSupport xs k = 0 := by
  funext j
  have hnone : ¬ ∃ x ∈ xs, x.val = 64 * k + j.val := by
    rintro ⟨x, hx, heq⟩
    have h := hbefore x hx
    have hj := j.isLt
    have hdiv : x.val / 64 = k := by
      rw [heq]
      omega
    omega
  simp [stepFromSupport, hnone]

theorem regionInput_sorted_support {R : ℕ}
    (T : Finset (Fin (64 * R))) (i : Fin R) :
    regionInput T i = stepFromSupport (T.sort (· ≤ ·)) i.val := by
  funext j
  have hmem :
      (∃ x ∈ T.sort (· ≤ ·), x.val = 64 * i.val + j.val) ↔
        regionBit i j ∈ T := by
    constructor
    · rintro ⟨x, hx, hv⟩
      have heq : x = regionBit i j := Fin.ext hv
      rw [heq] at hx
      exact (Finset.mem_sort (· ≤ ·)).mp hx
    · intro hx
      exact ⟨regionBit i j, (Finset.mem_sort (· ≤ ·)).mpr hx, rfl⟩
  change (if regionBit i j ∈ T then (1 : ZMod 2) else 0) =
    (if ∃ x ∈ T.sort (· ≤ ·), x.val = 64 * i.val + j.val then 1 else 0)
  simp only [hmem]

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

theorem markedInputList_getElem_support {R prev : ℕ}
    {xs : List (Fin (64 * R))} (h : ValidSupport R prev xs) :
    ∀ k (hk : k < R - prev),
      (markedInputList (supportMarks R xs prev)
        (supportLastGap R xs prev))[k]'(by
          rw [markedInputList_support_length h]
          exact hk) = stepFromSupport xs (prev + k) := by
  induction h with
  | nil hprev =>
      intro k hk
      simp only [supportMarks, supportLastGap, markedInputList,
        List.getElem_replicate]
      funext j
      simp [stepFromSupport]
  | @cons prev x xs hgap hinterior htail ih =>
      intro k hk
      let s := x.val / 64
      let g := s - prev
      by_cases hprefix : k < g
      · have hbefore : ∀ y ∈ x :: xs,
            prev + k < y.val / 64 := by
          intro y hy
          simp only [List.mem_cons] at hy
          rcases hy with rfl | hy
          · dsimp [g, s] at hprefix
            omega
          · have ht := validSupport_steps_gt htail y hy
            dsimp [g, s] at hprefix
            omega
        have hz := stepFromSupport_zero_of_before (x :: xs)
          (prev + k) hbefore
        have hn : ¬ x.val / 64 ≤ k + prev := by
          dsimp [g, s] at hprefix
          omega
        simp [supportMarks, supportLastGap, markedInputList,
          List.getElem_append, hn, hz]
      · by_cases hactive : k = g
        · subst k
          have hstep : prev + g = x.val / 64 := by
            dsimp [g, s]
            omega
          have hact : stepFromSupport (x :: xs) (prev + g) =
              singletonInput (regionLocal x) := by
            rw [hstep]
            exact stepFromSupport_active
              (ValidSupport.cons hgap hinterior htail) x (by simp)
          simp [supportMarks, supportLastGap, markedInputList,
            List.getElem_append, g, hact]
          simp [s]
          rfl
        · have hafter : g < k := by omega
          let k' := k - (g + 1)
          have hk' : k' < R - (x.val / 64 + 1) := by
            dsimp [k', g, s]
            omega
          have hnot : ¬ k < x.val / 64 - prev := by
            dsimp [g, s] at hafter
            omega
          have hne0 : k - (x.val / 64 - prev) ≠ 0 := by
            dsimp [g, s] at hafter
            omega
          have hidx : k - (x.val / 64 - prev) - 1 = k' := by
            dsimp [k', g, s]
            omega
          have hkfull : k < (markedInputList
              (supportMarks R (x :: xs) prev)
              (supportLastGap R (x :: xs) prev)).length := by
            rw [markedInputList_support_length
              (ValidSupport.cons hgap hinterior htail)]
            exact hk
          have hktail : k' < (markedInputList
              (supportMarks R xs (x.val / 64 + 1))
              (supportLastGap R xs (x.val / 64 + 1))).length := by
            rw [markedInputList_support_length htail]
            exact hk'
          have hget :
              (markedInputList (supportMarks R (x :: xs) prev)
                (supportLastGap R (x :: xs) prev))[k]'hkfull =
              (markedInputList (supportMarks R xs (x.val / 64 + 1))
                (supportLastGap R xs (x.val / 64 + 1)))[k']'hktail := by
            simp [supportMarks, supportLastGap, markedInputList,
              List.getElem_append, List.getElem_cons, hnot, hne0, hidx]
          have hstepeq : x.val / 64 + 1 + k' = prev + k := by
            dsimp [k', g, s]
            omega
          have hih := ih k' hk'
          rw [hstepeq] at hih
          have hne : x.val / 64 ≠ prev + k := by
            dsimp [g, s] at hafter
            omega
          rw [hget, hih, stepFromSupport_cons_ne x xs (prev + k) hne]

/-- A good routed support and its marked gaps induce the same ordered
64-bit inputs. This closes the combinatorial trace identification. -/
theorem good_support_trace {R : ℕ} (hR : 0 < R)
    (T : Finset (Fin (64 * R))) (hgood : ¬ badRegionalSites T) :
    chronologicalInputList T =
      markedInputList (supportMarks R (T.sort (· ≤ ·)) 0)
        (supportLastGap R (T.sort (· ≤ ·)) 0) := by
  have hvalid := goodRegionalSites_validSupport hR T hgood
  apply List.ext_getElem
  · rw [chronologicalInputList_length,
      markedInputList_goodSupport_length hR T hgood]
  · intro k hkleft hkright
    have hkR : k < R := by
      simpa only [chronologicalInputList_length] using hkleft
    have hmark := markedInputList_getElem_support hvalid k (by simpa using hkR)
    simp only [Nat.zero_add] at hmark
    rw [hmark]
    simpa [chronologicalInputList] using
      regionInput_sorted_support T ⟨k, hkR⟩

theorem empty_support_trace (R : ℕ) :
    chronologicalInputList (∅ : Finset (Fin (64 * R))) =
      markedInputList
        (supportMarks R ((∅ : Finset (Fin (64 * R))).sort (· ≤ ·)) 0)
        (supportLastGap R ((∅ : Finset (Fin (64 * R))).sort (· ≤ ·)) 0) := by
  simp [chronologicalInputList_empty, supportMarks, supportLastGap,
    markedInputList]

/-- The regional upper bound for the actual chronological encoder on
every good routed support. -/
theorem good_chronologicalRegion_coarse_upper {R : ℕ} (hR : 0 < R)
    (T : Finset (Fin (64 * R))) (hgood : ¬ badRegionalSites T)
    (θ : ℝ) (hθ : 0 ≤ θ) (q r : State22) :
    let L : ℝ := 64 * R
    let C := ((64 * θ + 64 * θ ^ 2) / L) / Real.exp (-θ)
    chronologicalRegionKernel (Real.exp (-(θ / L))) T q r ≤
      (1 + C) ^ (T.card + 1) *
        (liftM * coarseRegionProduct
          (fun g => Real.exp (-((θ / L) * liveWeightAverage coordMap)) ^ g)
          (supportMarks R (T.sort (· ≤ ·)) 0)
          (supportLastGap R (T.sort (· ≤ ·)) 0) * projectM) q r := by
  simpa [chronologicalRegionKernel, supportRegionKernel,
    good_support_trace hR T hgood,
    inputListKernel_markedInputList] using
    good_supportRegion_coarse_upper hR T hgood θ hθ q r

end Spin.UniformMixer
