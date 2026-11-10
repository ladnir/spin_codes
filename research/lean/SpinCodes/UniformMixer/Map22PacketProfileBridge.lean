import SpinCodes.UniformMixer.Map22PacketProfileDefs
import SpinCodes.UniformMixer.QuadraticCoefficientTriple

/-!
# The production map's low-two-variable decomposition

Reuse the checked identification of the production image with all quadratic
Boolean functions on six variables. Reindex its 22 coefficients as one
quadratic coefficient in the low variables, two affine functions of the four
high variables, and one quadratic function of those high variables. This
is a bijection of coefficient spaces, not a new distributional assumption.
-/
noncomputable section
namespace Spin.UniformMixer

def profileVectorToSplit (v : State22) : ProfileSplit6 :=
  (v 7, ![v 1, v 8, v 9, v 10, v 11],
    ![v 2, v 12, v 13, v 14, v 15],
    ![v 0, v 3, v 4, v 5, v 6, v 16, v 17, v 18, v 19, v 20, v 21])

def profileSplitToVector (t : ProfileSplit6) : State22 := ![
  t.2.2.2 0, t.2.1 0, t.2.2.1 0,
  t.2.2.2 1, t.2.2.2 2, t.2.2.2 3, t.2.2.2 4, t.1,
  t.2.1 1, t.2.1 2, t.2.1 3, t.2.1 4,
  t.2.2.1 1, t.2.2.1 2, t.2.2.1 3, t.2.2.1 4,
  t.2.2.2 5, t.2.2.2 6, t.2.2.2 7, t.2.2.2 8, t.2.2.2 9, t.2.2.2 10]

theorem profileSplitToVector_toSplit (v : State22) :
    profileSplitToVector (profileVectorToSplit v) = v := by
  funext i
  fin_cases i <;> rfl

theorem profileVectorToSplit_toVector (t : ProfileSplit6) :
    profileVectorToSplit (profileSplitToVector t) = t := by
  apply Prod.ext
  · rfl
  · apply Prod.ext
    · funext i; fin_cases i <;> rfl
    · apply Prod.ext
      · funext i; fin_cases i <;> rfl
      · funext i; fin_cases i <;> rfl

def profileVectorSplitEquiv : State22 ≃ ProfileSplit6 where
  toFun := profileVectorToSplit
  invFun := profileSplitToVector
  left_inv := profileSplitToVector_toSplit
  right_inv := profileVectorToSplit_toVector

def profileCanonicalVectorEquiv : QuadraticCoefficients ≃ State22 where
  toFun c i := quadraticCoefficientRestriction c (canonicalMonomialEquiv i)
  invFun v := quadraticCoefficientRestriction.symm
    (fun m => v (canonicalMonomialEquiv.symm m))
  left_inv c := by
    apply quadraticCoefficientRestriction.injective
    funext m
    simp
  right_inv v := by
    funext i
    simp

def productionProfileSplitEquiv : State22 ≃ ProfileSplit6 :=
  (Equiv.ofBijective productionQuadraticCoefficients
    productionQuadraticCoefficients_bijective).trans
    (profileCanonicalVectorEquiv.trans profileVectorSplitEquiv)

/- Check four positions per lemma, following the historical split-check
pattern. Each coefficient remains symbolic; no state space is enumerated. -/
private theorem canonicalProfileSplit_evaluation_0 (v : State22) (l : Fin 4) :
    (∑ i : Fin 22, if isSubmask (canonicalMonomial i) (finProdFinEquiv ((0 : Fin 16),l))
      then v i else 0) = profileSplit6Eval (profileVectorToSplit v) (0,l) := by
  fin_cases l <;>
    simp [Fin.sum_univ_succ, canonicalMonomial, isSubmask,
      profileSplit6Eval, profileNibbleEval, profileVectorToSplit,
      profileAffine4Eval, profileQuad4Eval, profileQuad4Rows, bitScalar,
      finProdFinEquiv, Nat.testBit_succ] <;> ring

private theorem canonicalProfileSplit_evaluation_1 (v : State22) (l : Fin 4) :
    (∑ i : Fin 22, if isSubmask (canonicalMonomial i) (finProdFinEquiv ((1 : Fin 16),l))
      then v i else 0) = profileSplit6Eval (profileVectorToSplit v) (1,l) := by
  fin_cases l <;>
    simp [Fin.sum_univ_succ, canonicalMonomial, isSubmask,
      profileSplit6Eval, profileNibbleEval, profileVectorToSplit,
      profileAffine4Eval, profileQuad4Eval, profileQuad4Rows, bitScalar,
      finProdFinEquiv, Nat.testBit_succ] <;> ring

private theorem canonicalProfileSplit_evaluation_2 (v : State22) (l : Fin 4) :
    (∑ i : Fin 22, if isSubmask (canonicalMonomial i) (finProdFinEquiv ((2 : Fin 16),l))
      then v i else 0) = profileSplit6Eval (profileVectorToSplit v) (2,l) := by
  fin_cases l <;>
    simp [Fin.sum_univ_succ, canonicalMonomial, isSubmask,
      profileSplit6Eval, profileNibbleEval, profileVectorToSplit,
      profileAffine4Eval, profileQuad4Eval, profileQuad4Rows, bitScalar,
      finProdFinEquiv, Nat.testBit_succ] <;> ring

private theorem canonicalProfileSplit_evaluation_3 (v : State22) (l : Fin 4) :
    (∑ i : Fin 22, if isSubmask (canonicalMonomial i) (finProdFinEquiv ((3 : Fin 16),l))
      then v i else 0) = profileSplit6Eval (profileVectorToSplit v) (3,l) := by
  fin_cases l <;>
    simp [Fin.sum_univ_succ, canonicalMonomial, isSubmask,
      profileSplit6Eval, profileNibbleEval, profileVectorToSplit,
      profileAffine4Eval, profileQuad4Eval, profileQuad4Rows, bitScalar,
      finProdFinEquiv, Nat.testBit_succ] <;> ring

private theorem canonicalProfileSplit_evaluation_4 (v : State22) (l : Fin 4) :
    (∑ i : Fin 22, if isSubmask (canonicalMonomial i) (finProdFinEquiv ((4 : Fin 16),l))
      then v i else 0) = profileSplit6Eval (profileVectorToSplit v) (4,l) := by
  fin_cases l <;>
    simp [Fin.sum_univ_succ, canonicalMonomial, isSubmask,
      profileSplit6Eval, profileNibbleEval, profileVectorToSplit,
      profileAffine4Eval, profileQuad4Eval, profileQuad4Rows, bitScalar,
      finProdFinEquiv, Nat.testBit_succ] <;> ring

private theorem canonicalProfileSplit_evaluation_5 (v : State22) (l : Fin 4) :
    (∑ i : Fin 22, if isSubmask (canonicalMonomial i) (finProdFinEquiv ((5 : Fin 16),l))
      then v i else 0) = profileSplit6Eval (profileVectorToSplit v) (5,l) := by
  fin_cases l <;>
    simp [Fin.sum_univ_succ, canonicalMonomial, isSubmask,
      profileSplit6Eval, profileNibbleEval, profileVectorToSplit,
      profileAffine4Eval, profileQuad4Eval, profileQuad4Rows, bitScalar,
      finProdFinEquiv, Nat.testBit_succ] <;> ring

private theorem canonicalProfileSplit_evaluation_6 (v : State22) (l : Fin 4) :
    (∑ i : Fin 22, if isSubmask (canonicalMonomial i) (finProdFinEquiv ((6 : Fin 16),l))
      then v i else 0) = profileSplit6Eval (profileVectorToSplit v) (6,l) := by
  fin_cases l <;>
    simp [Fin.sum_univ_succ, canonicalMonomial, isSubmask,
      profileSplit6Eval, profileNibbleEval, profileVectorToSplit,
      profileAffine4Eval, profileQuad4Eval, profileQuad4Rows, bitScalar,
      finProdFinEquiv, Nat.testBit_succ] <;> ring

private theorem canonicalProfileSplit_evaluation_7 (v : State22) (l : Fin 4) :
    (∑ i : Fin 22, if isSubmask (canonicalMonomial i) (finProdFinEquiv ((7 : Fin 16),l))
      then v i else 0) = profileSplit6Eval (profileVectorToSplit v) (7,l) := by
  fin_cases l <;>
    simp [Fin.sum_univ_succ, canonicalMonomial, isSubmask,
      profileSplit6Eval, profileNibbleEval, profileVectorToSplit,
      profileAffine4Eval, profileQuad4Eval, profileQuad4Rows, bitScalar,
      finProdFinEquiv, Nat.testBit_succ] <;> ring

private theorem canonicalProfileSplit_evaluation_8 (v : State22) (l : Fin 4) :
    (∑ i : Fin 22, if isSubmask (canonicalMonomial i) (finProdFinEquiv ((8 : Fin 16),l))
      then v i else 0) = profileSplit6Eval (profileVectorToSplit v) (8,l) := by
  fin_cases l <;>
    simp [Fin.sum_univ_succ, canonicalMonomial, isSubmask,
      profileSplit6Eval, profileNibbleEval, profileVectorToSplit,
      profileAffine4Eval, profileQuad4Eval, profileQuad4Rows, bitScalar,
      finProdFinEquiv, Nat.testBit_succ] <;> ring

private theorem canonicalProfileSplit_evaluation_9 (v : State22) (l : Fin 4) :
    (∑ i : Fin 22, if isSubmask (canonicalMonomial i) (finProdFinEquiv ((9 : Fin 16),l))
      then v i else 0) = profileSplit6Eval (profileVectorToSplit v) (9,l) := by
  fin_cases l <;>
    simp [Fin.sum_univ_succ, canonicalMonomial, isSubmask,
      profileSplit6Eval, profileNibbleEval, profileVectorToSplit,
      profileAffine4Eval, profileQuad4Eval, profileQuad4Rows, bitScalar,
      finProdFinEquiv, Nat.testBit_succ] <;> ring

private theorem canonicalProfileSplit_evaluation_10 (v : State22) (l : Fin 4) :
    (∑ i : Fin 22, if isSubmask (canonicalMonomial i) (finProdFinEquiv ((10 : Fin 16),l))
      then v i else 0) = profileSplit6Eval (profileVectorToSplit v) (10,l) := by
  fin_cases l <;>
    simp [Fin.sum_univ_succ, canonicalMonomial, isSubmask,
      profileSplit6Eval, profileNibbleEval, profileVectorToSplit,
      profileAffine4Eval, profileQuad4Eval, profileQuad4Rows, bitScalar,
      finProdFinEquiv, Nat.testBit_succ] <;> ring

private theorem canonicalProfileSplit_evaluation_11 (v : State22) (l : Fin 4) :
    (∑ i : Fin 22, if isSubmask (canonicalMonomial i) (finProdFinEquiv ((11 : Fin 16),l))
      then v i else 0) = profileSplit6Eval (profileVectorToSplit v) (11,l) := by
  fin_cases l <;>
    simp [Fin.sum_univ_succ, canonicalMonomial, isSubmask,
      profileSplit6Eval, profileNibbleEval, profileVectorToSplit,
      profileAffine4Eval, profileQuad4Eval, profileQuad4Rows, bitScalar,
      finProdFinEquiv, Nat.testBit_succ] <;> ring

private theorem canonicalProfileSplit_evaluation_12 (v : State22) (l : Fin 4) :
    (∑ i : Fin 22, if isSubmask (canonicalMonomial i) (finProdFinEquiv ((12 : Fin 16),l))
      then v i else 0) = profileSplit6Eval (profileVectorToSplit v) (12,l) := by
  fin_cases l <;>
    simp [Fin.sum_univ_succ, canonicalMonomial, isSubmask,
      profileSplit6Eval, profileNibbleEval, profileVectorToSplit,
      profileAffine4Eval, profileQuad4Eval, profileQuad4Rows, bitScalar,
      finProdFinEquiv, Nat.testBit_succ] <;> ring

private theorem canonicalProfileSplit_evaluation_13 (v : State22) (l : Fin 4) :
    (∑ i : Fin 22, if isSubmask (canonicalMonomial i) (finProdFinEquiv ((13 : Fin 16),l))
      then v i else 0) = profileSplit6Eval (profileVectorToSplit v) (13,l) := by
  fin_cases l <;>
    simp [Fin.sum_univ_succ, canonicalMonomial, isSubmask,
      profileSplit6Eval, profileNibbleEval, profileVectorToSplit,
      profileAffine4Eval, profileQuad4Eval, profileQuad4Rows, bitScalar,
      finProdFinEquiv, Nat.testBit_succ] <;> ring

private theorem canonicalProfileSplit_evaluation_14 (v : State22) (l : Fin 4) :
    (∑ i : Fin 22, if isSubmask (canonicalMonomial i) (finProdFinEquiv ((14 : Fin 16),l))
      then v i else 0) = profileSplit6Eval (profileVectorToSplit v) (14,l) := by
  fin_cases l <;>
    simp [Fin.sum_univ_succ, canonicalMonomial, isSubmask,
      profileSplit6Eval, profileNibbleEval, profileVectorToSplit,
      profileAffine4Eval, profileQuad4Eval, profileQuad4Rows, bitScalar,
      finProdFinEquiv, Nat.testBit_succ] <;> ring

private theorem canonicalProfileSplit_evaluation_15 (v : State22) (l : Fin 4) :
    (∑ i : Fin 22, if isSubmask (canonicalMonomial i) (finProdFinEquiv ((15 : Fin 16),l))
      then v i else 0) = profileSplit6Eval (profileVectorToSplit v) (15,l) := by
  fin_cases l <;>
    simp [Fin.sum_univ_succ, canonicalMonomial, isSubmask,
      profileSplit6Eval, profileNibbleEval, profileVectorToSplit,
      profileAffine4Eval, profileQuad4Eval, profileQuad4Rows, bitScalar,
      finProdFinEquiv, Nat.testBit_succ] <;> ring

theorem canonicalProfileSplit_evaluation (v : State22) (h : Fin 16) (l : Fin 4) :
    (∑ i : Fin 22, if isSubmask (canonicalMonomial i) (finProdFinEquiv (h,l))
      then v i else 0) = profileSplit6Eval (profileVectorToSplit v) (h,l) := by
  fin_cases h
  · exact canonicalProfileSplit_evaluation_0 v l
  · exact canonicalProfileSplit_evaluation_1 v l
  · exact canonicalProfileSplit_evaluation_2 v l
  · exact canonicalProfileSplit_evaluation_3 v l
  · exact canonicalProfileSplit_evaluation_4 v l
  · exact canonicalProfileSplit_evaluation_5 v l
  · exact canonicalProfileSplit_evaluation_6 v l
  · exact canonicalProfileSplit_evaluation_7 v l
  · exact canonicalProfileSplit_evaluation_8 v l
  · exact canonicalProfileSplit_evaluation_9 v l
  · exact canonicalProfileSplit_evaluation_10 v l
  · exact canonicalProfileSplit_evaluation_11 v l
  · exact canonicalProfileSplit_evaluation_12 v l
  · exact canonicalProfileSplit_evaluation_13 v l
  · exact canonicalProfileSplit_evaluation_14 v l
  · exact canonicalProfileSplit_evaluation_15 v l

/-- The exact packet ordering used by the production map. -/
theorem productionProfileSplit_evaluation (q : State22) (vu : Fin 16 × Fin 4) :
    coordMap (finProdFinEquiv vu) q =
      profileSplit6Eval (productionProfileSplitEquiv q) vu := by
  rw [production_quadratic_evaluation, quadraticEvaluation_canonical]
  exact canonicalProfileSplit_evaluation
    (profileCanonicalVectorEquiv (productionQuadraticCoefficients q)) vu.1 vu.2

theorem productionProfileSplit_profile (q : State22) :
    packetProfile64 (fun j => coordMap j q) =
      profileOfNibbles (profileSplit6Eval (productionProfileSplitEquiv q)) := by
  unfold packetProfile64
  have he : (fun vu => coordMap (finProdFinEquiv vu) q) =
      profileSplit6Eval (productionProfileSplitEquiv q) := by
    funext vu
    exact productionProfileSplit_evaluation q vu
  exact congrArg profileOfNibbles he

#print axioms productionProfileSplit_evaluation
#print axioms productionProfileSplit_profile
end Spin.UniformMixer
