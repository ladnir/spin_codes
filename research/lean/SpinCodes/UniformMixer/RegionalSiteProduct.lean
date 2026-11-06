import SpinCodes.UniformMixer.RegionalProductBridge

/-!
# Ordered bit sites and the continuum product

The ideal product built recursively from a sorted list of bit positions
is the continuum product evaluated at those positions after normalization.
-/

noncomputable section

namespace Spin.UniformMixer

def positionGapsFrom {a : ℕ} (first : ℝ) (x : Fin a → ℝ)
    (i : Fin (a + 1)) : ℝ :=
  (Fin.snoc x (1 : ℝ) : Fin (a + 1) → ℝ) i -
    (Fin.cons first x : Fin (a + 1) → ℝ) i

def siteProductFrom (γ first : ℝ) {a : ℕ} (x : Fin a → ℝ) :
    Matrix Bool Bool ℝ :=
  timeRegionProduct γ (fun i => positionGapsFrom first x i.castSucc)
    (positionGapsFrom first x (Fin.last a))

theorem siteProductFrom_zero (γ first : ℝ) :
    siteProductFrom γ first (Fin.elim0 : Fin 0 → ℝ) =
      timeRegionGap γ (1 - first) := by
  simp only [siteProductFrom, timeRegionProduct, positionGapsFrom,
    Fin.snoc_last]
  rfl

theorem siteProductFrom_cons (γ first head : ℝ) {a : ℕ}
    (tail : Fin a → ℝ) :
    siteProductFrom γ first (Fin.cons head tail) =
      timeRegionGap γ (head - first) * coarseP *
        siteProductFrom γ head tail := by
  have hhead : positionGapsFrom first (Fin.cons head tail)
      (0 : Fin (a + 1 + 1)) = head - first := by
    simp [positionGapsFrom]
  have htail :
      (Fin.tail (fun i : Fin (a + 1) =>
        positionGapsFrom first (Fin.cons head tail) i.castSucc)) =
      (fun i : Fin a => positionGapsFrom head tail i.castSucc) := by
    funext i
    have hidx : i.castSucc.succ = i.succ.castSucc := by ext; rfl
    simp only [positionGapsFrom, Fin.tail, Fin.snoc_castSucc,
      Fin.cons_succ]
    rw [show i.succ.castSucc = i.castSucc.succ from hidx.symm,
      Fin.cons_succ]
  have hlast : positionGapsFrom first (Fin.cons head tail)
      (Fin.last (a + 1)) = positionGapsFrom head tail (Fin.last a) := by
    simp [positionGapsFrom]
  simp only [siteProductFrom, timeRegionProduct]
  change timeRegionGap γ (positionGapsFrom first (Fin.cons head tail) 0) *
    coarseP * timeRegionProduct γ
      (Fin.tail fun i => positionGapsFrom first (Fin.cons head tail) i.castSucc)
      (positionGapsFrom first (Fin.cons head tail) (Fin.last (a + 1))) = _
  rw [hhead, htail, hlast]

theorem idealSupportProduct_ofFn (γ : ℝ) (R a : ℕ)
    (x : Fin a → Fin (64 * R)) (prevBit : ℕ) :
    idealSupportProduct γ R (List.ofFn x) prevBit =
      siteProductFrom γ ((prevBit : ℝ) / (64 * R))
        (fun i => ((x i).val : ℝ) / (64 * R)) := by
  induction a generalizing prevBit with
  | zero =>
      simp only [List.ofFn_zero, idealSupportProduct]
      have hx : (fun i : Fin 0 => ((x i).val : ℝ) / (64 * R)) = Fin.elim0 :=
        Subsingleton.elim _ _
      rw [hx, siteProductFrom_zero]
  | succ a ih =>
      rw [List.ofFn_succ]
      change timeRegionGap γ
        ((((x 0).val : ℝ) - prevBit) / (64 * R)) * coarseP *
          idealSupportProduct γ R (List.ofFn (Fin.tail x)) (x 0).val = _
      have hx : (fun i : Fin (a + 1) => ((x i).val : ℝ) / (64 * R)) =
          Fin.cons (((x 0).val : ℝ) / (64 * R))
            (fun i : Fin a => ((Fin.tail x i).val : ℝ) / (64 * R)) := by
        funext i
        refine Fin.cases ?_ ?_ i
        · rfl
        · intro j
          rfl
      rw [hx, siteProductFrom_cons, ih]
      congr 2
      ring

theorem siteProductFrom_zero_origin (γ : ℝ) {a : ℕ}
    (x : Fin a → ℝ) :
    siteProductFrom γ 0 x = continuumSiteProduct γ x := by
  unfold siteProductFrom continuumSiteProduct positionGapsFrom
    Spin.Structured.Placement.positionGaps
  rfl

theorem sorted_sites_ofFn {R a : ℕ}
    (B : Spin.Structured.Placement.BlockSubset (64 * R) a) :
    B.val.sort (· ≤ ·) =
      List.ofFn (Spin.Structured.Placement.orderedBlocks B) := by
  rw [List.ofFn_eq_map]
  exact (B.val.listMap_orderEmbOfFin_finRange B.property).symm

theorem idealSupportProduct_sorted (γ : ℝ) {R a : ℕ}
    (B : Spin.Structured.Placement.BlockSubset (64 * R) a) :
    idealSupportProduct γ R (B.val.sort (· ≤ ·)) 0 =
      continuumSiteProduct γ (Spin.Structured.Placement.normalizedSites B) := by
  rw [sorted_sites_ofFn, idealSupportProduct_ofFn]
  simp only [Nat.cast_zero, zero_div, siteProductFrom_zero_origin]
  congr 1
  funext i
  change ((Spin.Structured.Placement.orderedBlocks B i).val : ℝ) /
      (64 * (R : ℝ)) =
    ((Spin.Structured.Placement.orderedBlocks B i).val : ℝ) /
      (((64 * R : ℕ) : ℝ))
  push_cast
  rfl

end Spin.UniformMixer
