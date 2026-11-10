import SpinCodes.Turbo.RoutingPairedProfiles22
import SpinCodes.Turbo.PhysicalProduct22

/-!
# Chronological count pairing and complete physical envelope products

The pair indexing is the original stream order, not a permutation of steps.
Complete physical rows then transfer the entire continuous-state endpoint
before the fractional matrix-path bound is applied.
-/
noncomputable section
namespace Spin.Turbo
open Finset Spin.UniformMixer

theorem pairedCountEquiv22_list (n : ℕ) (counts : Fin (n * 2) → Fin 17) :
    List.ofFn counts =
      (List.ofFn (fun b : Fin n =>
        [(pairedCountEquiv22 n counts b).1, (pairedCountEquiv22 n counts b).2])).flatten := by
  rw [List.ofFn_mul]
  apply congrArg List.flatten
  apply congrArg List.ofFn
  funext b
  change List.ofFn (fun t : Fin 2 => counts ⟨b.val * 2 + t.val, _⟩) =
    List.ofFn (fun t : Fin 2 => counts (finProdFinEquiv (b, t)))
  apply congrArg List.ofFn
  funext t
  apply congrArg counts
  ext
  change b.val * 2 + t.val = t.val + 2 * b.val
  omega

def turboCountStream22 {R : ℕ} (counts : TurboPotentialCounts22 R) : List (Fin 17) :=
  (List.ofFn (fun region : Fin 64 => List.ofFn (counts region))).flatten

def turboCountPairs22 {n : ℕ} (counts : TurboPotentialCounts22 (n * 2)) :
    List (Fin 17 × Fin 17) :=
  (List.ofFn (fun region : Fin 64 =>
    List.ofFn (pairedCountEquiv22 n (counts region)))).flatten

theorem turboCountStream22_eq_pairs {n : ℕ} (counts : TurboPotentialCounts22 (n * 2)) :
    turboCountStream22 counts = (turboCountPairs22 counts).flatMap (fun ab => [ab.1, ab.2]) := by
  unfold turboCountStream22 turboCountPairs22
  simp_rw [pairedCountEquiv22_list]
  simp only [List.flatMap_def, List.map_flatten, List.flatten_flatten,
    List.map_map, List.map_ofFn, Function.comp_def]

/-- Full-state kernel factors follow exactly the chronological count stream. -/
theorem turboCountKernel22_eq_stream {R : ℕ} (counts : TurboPotentialCounts22 R) (z : ℝ) :
    turboCountKernel22 counts z =
      ((turboCountStream22 counts).map
        (fun j => Matrix.of (turboSelectedKernel22 j z))).prod := by
  simp only [turboCountKernel22, turboCountStream22, List.map_flatten,
    List.map_ofFn, List.prod_flatten, Function.comp_def]

/-- The scalar full-state endpoint is bounded by paired complete rows, with
chronological pair products taken before entrywise fractional powers. -/
theorem turboCountKernel22_fractional_le_pairs {n : ℕ}
    (counts : TurboPotentialCounts22 (n * 2)) (z : ℝ)
    (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hs : ∀ j, physicalRowSound22 (Matrix.of (turboSelectedKernel22 j z)) (B j))
    (hB : ∀ j a b, 0 ≤ B j a b)
    {α : ℝ} (hα0 : 0 < α) (hα1 : α ≤ 1) :
    (∑ q : State22, turboCountKernel22 counts z 0 q) ^ α ≤
      ∑ k : Fin 11, (List.ofFn (fun region : Fin 64 =>
        regionalPairedMatrix22 B α n (counts region))).prod 0 k := by
  let pairs : List (PhysicalPair22 × PhysicalPair22) :=
    (turboCountPairs22 counts).map (fun ab =>
      ((Matrix.of (turboSelectedKernel22 ab.1 z), B ab.1),
        (Matrix.of (turboSelectedKernel22 ab.2 z), B ab.2)))
  have h := physical_pairs_product_fractional_le pairs (by
    intro p hp
    obtain ⟨ab, _, rfl⟩ := List.mem_map.mp hp
    exact ⟨hs ab.1, hs ab.2⟩) (by
      intro p hp
      obtain ⟨ab, _, rfl⟩ := List.mem_map.mp hp
      exact ⟨hB ab.1, hB ab.2⟩) hα0 hα1
  have hleft : ((pairs.flatMap (fun p => [p.1, p.2])).map Prod.fst).prod =
      turboCountKernel22 counts z := by
    rw [turboCountKernel22_eq_stream, turboCountStream22_eq_pairs]
    simp only [pairs, List.flatMap_map, List.map_flatMap, List.map_cons, List.map_nil]
  have hright : (pairs.map (fun p => matrixEntrywiseRpow (p.1.2 * p.2.2) α)).prod =
      (List.ofFn (fun region : Fin 64 => regionalPairedMatrix22 B α n (counts region))).prod := by
    simp only [pairs, turboCountPairs22, List.map_map, Function.comp_def,
      List.map_flatten, List.map_ofFn, List.prod_flatten, regionalPairedMatrix22]
  rw [hleft, hright] at h
  exact h

#print axioms turboCountKernel22_fractional_le_pairs
end Spin.Turbo

