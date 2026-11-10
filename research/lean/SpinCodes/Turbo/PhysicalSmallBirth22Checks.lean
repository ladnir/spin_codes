import SpinCodes.Turbo.PhysicalSmallBirth22PairChecks
import Mathlib.Data.Fin.VecNotation
import Mathlib.Algebra.BigOperators.Fin
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.NormNum

namespace Spin.Turbo
set_option maxRecDepth 100000
set_option maxHeartbeats 2000000

theorem smallBirthPairHistograms22_checked (a b : Fin 16) (h : a < b) :
    smallBirthCodeHistogram22 (smallBirthPairCodes22 a b) =
      smallBirthPairHistograms22 a b := by
  fin_cases a
  · fin_cases b
    · norm_num at h
    · exact smallBirthPair0_1Histogram22_checked
    · exact smallBirthPair0_2Histogram22_checked
    · exact smallBirthPair0_3Histogram22_checked
    · exact smallBirthPair0_4Histogram22_checked
    · exact smallBirthPair0_5Histogram22_checked
    · exact smallBirthPair0_6Histogram22_checked
    · exact smallBirthPair0_7Histogram22_checked
    · exact smallBirthPair0_8Histogram22_checked
    · exact smallBirthPair0_9Histogram22_checked
    · exact smallBirthPair0_10Histogram22_checked
    · exact smallBirthPair0_11Histogram22_checked
    · exact smallBirthPair0_12Histogram22_checked
    · exact smallBirthPair0_13Histogram22_checked
    · exact smallBirthPair0_14Histogram22_checked
    · exact smallBirthPair0_15Histogram22_checked
  · fin_cases b
    · norm_num at h
    · norm_num at h
    · exact smallBirthPair1_2Histogram22_checked
    · exact smallBirthPair1_3Histogram22_checked
    · exact smallBirthPair1_4Histogram22_checked
    · exact smallBirthPair1_5Histogram22_checked
    · exact smallBirthPair1_6Histogram22_checked
    · exact smallBirthPair1_7Histogram22_checked
    · exact smallBirthPair1_8Histogram22_checked
    · exact smallBirthPair1_9Histogram22_checked
    · exact smallBirthPair1_10Histogram22_checked
    · exact smallBirthPair1_11Histogram22_checked
    · exact smallBirthPair1_12Histogram22_checked
    · exact smallBirthPair1_13Histogram22_checked
    · exact smallBirthPair1_14Histogram22_checked
    · exact smallBirthPair1_15Histogram22_checked
  · fin_cases b
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · exact smallBirthPair2_3Histogram22_checked
    · exact smallBirthPair2_4Histogram22_checked
    · exact smallBirthPair2_5Histogram22_checked
    · exact smallBirthPair2_6Histogram22_checked
    · exact smallBirthPair2_7Histogram22_checked
    · exact smallBirthPair2_8Histogram22_checked
    · exact smallBirthPair2_9Histogram22_checked
    · exact smallBirthPair2_10Histogram22_checked
    · exact smallBirthPair2_11Histogram22_checked
    · exact smallBirthPair2_12Histogram22_checked
    · exact smallBirthPair2_13Histogram22_checked
    · exact smallBirthPair2_14Histogram22_checked
    · exact smallBirthPair2_15Histogram22_checked
  · fin_cases b
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · exact smallBirthPair3_4Histogram22_checked
    · exact smallBirthPair3_5Histogram22_checked
    · exact smallBirthPair3_6Histogram22_checked
    · exact smallBirthPair3_7Histogram22_checked
    · exact smallBirthPair3_8Histogram22_checked
    · exact smallBirthPair3_9Histogram22_checked
    · exact smallBirthPair3_10Histogram22_checked
    · exact smallBirthPair3_11Histogram22_checked
    · exact smallBirthPair3_12Histogram22_checked
    · exact smallBirthPair3_13Histogram22_checked
    · exact smallBirthPair3_14Histogram22_checked
    · exact smallBirthPair3_15Histogram22_checked
  · fin_cases b
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · exact smallBirthPair4_5Histogram22_checked
    · exact smallBirthPair4_6Histogram22_checked
    · exact smallBirthPair4_7Histogram22_checked
    · exact smallBirthPair4_8Histogram22_checked
    · exact smallBirthPair4_9Histogram22_checked
    · exact smallBirthPair4_10Histogram22_checked
    · exact smallBirthPair4_11Histogram22_checked
    · exact smallBirthPair4_12Histogram22_checked
    · exact smallBirthPair4_13Histogram22_checked
    · exact smallBirthPair4_14Histogram22_checked
    · exact smallBirthPair4_15Histogram22_checked
  · fin_cases b
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · exact smallBirthPair5_6Histogram22_checked
    · exact smallBirthPair5_7Histogram22_checked
    · exact smallBirthPair5_8Histogram22_checked
    · exact smallBirthPair5_9Histogram22_checked
    · exact smallBirthPair5_10Histogram22_checked
    · exact smallBirthPair5_11Histogram22_checked
    · exact smallBirthPair5_12Histogram22_checked
    · exact smallBirthPair5_13Histogram22_checked
    · exact smallBirthPair5_14Histogram22_checked
    · exact smallBirthPair5_15Histogram22_checked
  · fin_cases b
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · exact smallBirthPair6_7Histogram22_checked
    · exact smallBirthPair6_8Histogram22_checked
    · exact smallBirthPair6_9Histogram22_checked
    · exact smallBirthPair6_10Histogram22_checked
    · exact smallBirthPair6_11Histogram22_checked
    · exact smallBirthPair6_12Histogram22_checked
    · exact smallBirthPair6_13Histogram22_checked
    · exact smallBirthPair6_14Histogram22_checked
    · exact smallBirthPair6_15Histogram22_checked
  · fin_cases b
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · exact smallBirthPair7_8Histogram22_checked
    · exact smallBirthPair7_9Histogram22_checked
    · exact smallBirthPair7_10Histogram22_checked
    · exact smallBirthPair7_11Histogram22_checked
    · exact smallBirthPair7_12Histogram22_checked
    · exact smallBirthPair7_13Histogram22_checked
    · exact smallBirthPair7_14Histogram22_checked
    · exact smallBirthPair7_15Histogram22_checked
  · fin_cases b
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · exact smallBirthPair8_9Histogram22_checked
    · exact smallBirthPair8_10Histogram22_checked
    · exact smallBirthPair8_11Histogram22_checked
    · exact smallBirthPair8_12Histogram22_checked
    · exact smallBirthPair8_13Histogram22_checked
    · exact smallBirthPair8_14Histogram22_checked
    · exact smallBirthPair8_15Histogram22_checked
  · fin_cases b
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · exact smallBirthPair9_10Histogram22_checked
    · exact smallBirthPair9_11Histogram22_checked
    · exact smallBirthPair9_12Histogram22_checked
    · exact smallBirthPair9_13Histogram22_checked
    · exact smallBirthPair9_14Histogram22_checked
    · exact smallBirthPair9_15Histogram22_checked
  · fin_cases b
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · exact smallBirthPair10_11Histogram22_checked
    · exact smallBirthPair10_12Histogram22_checked
    · exact smallBirthPair10_13Histogram22_checked
    · exact smallBirthPair10_14Histogram22_checked
    · exact smallBirthPair10_15Histogram22_checked
  · fin_cases b
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · exact smallBirthPair11_12Histogram22_checked
    · exact smallBirthPair11_13Histogram22_checked
    · exact smallBirthPair11_14Histogram22_checked
    · exact smallBirthPair11_15Histogram22_checked
  · fin_cases b
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · exact smallBirthPair12_13Histogram22_checked
    · exact smallBirthPair12_14Histogram22_checked
    · exact smallBirthPair12_15Histogram22_checked
  · fin_cases b
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · exact smallBirthPair13_14Histogram22_checked
    · exact smallBirthPair13_15Histogram22_checked
  · fin_cases b
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · exact smallBirthPair14_15Histogram22_checked
  · fin_cases b
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
    · norm_num at h
def smallBirthTwoShardHistograms22 : Fin 16 → List Nat := ![
smallBirthTwoShard0Histogram22, smallBirthTwoShard1Histogram22, smallBirthTwoShard2Histogram22, smallBirthTwoShard3Histogram22, smallBirthTwoShard4Histogram22, smallBirthTwoShard5Histogram22, smallBirthTwoShard6Histogram22, smallBirthTwoShard7Histogram22, smallBirthTwoShard8Histogram22, smallBirthTwoShard9Histogram22, smallBirthTwoShard10Histogram22, smallBirthTwoShard11Histogram22, smallBirthTwoShard12Histogram22, smallBirthTwoShard13Histogram22, smallBirthTwoShard14Histogram22, smallBirthTwoShard15Histogram22]

theorem smallBirthTwoShard0Sum22_checked : ∀ k : Fin 99,
    smallBirthTwoShard0Histogram22.getD k.val 0 =
      ∑ b : Fin 16, if (0 : Fin 16) < b then
        (smallBirthPairHistograms22 0 b).getD k.val 0 else 0 := by decide +kernel

theorem smallBirthTwoShard1Sum22_checked : ∀ k : Fin 99,
    smallBirthTwoShard1Histogram22.getD k.val 0 =
      ∑ b : Fin 16, if (1 : Fin 16) < b then
        (smallBirthPairHistograms22 1 b).getD k.val 0 else 0 := by decide +kernel

theorem smallBirthTwoShard2Sum22_checked : ∀ k : Fin 99,
    smallBirthTwoShard2Histogram22.getD k.val 0 =
      ∑ b : Fin 16, if (2 : Fin 16) < b then
        (smallBirthPairHistograms22 2 b).getD k.val 0 else 0 := by decide +kernel

theorem smallBirthTwoShard3Sum22_checked : ∀ k : Fin 99,
    smallBirthTwoShard3Histogram22.getD k.val 0 =
      ∑ b : Fin 16, if (3 : Fin 16) < b then
        (smallBirthPairHistograms22 3 b).getD k.val 0 else 0 := by decide +kernel

theorem smallBirthTwoShard4Sum22_checked : ∀ k : Fin 99,
    smallBirthTwoShard4Histogram22.getD k.val 0 =
      ∑ b : Fin 16, if (4 : Fin 16) < b then
        (smallBirthPairHistograms22 4 b).getD k.val 0 else 0 := by decide +kernel

theorem smallBirthTwoShard5Sum22_checked : ∀ k : Fin 99,
    smallBirthTwoShard5Histogram22.getD k.val 0 =
      ∑ b : Fin 16, if (5 : Fin 16) < b then
        (smallBirthPairHistograms22 5 b).getD k.val 0 else 0 := by decide +kernel

theorem smallBirthTwoShard6Sum22_checked : ∀ k : Fin 99,
    smallBirthTwoShard6Histogram22.getD k.val 0 =
      ∑ b : Fin 16, if (6 : Fin 16) < b then
        (smallBirthPairHistograms22 6 b).getD k.val 0 else 0 := by decide +kernel

theorem smallBirthTwoShard7Sum22_checked : ∀ k : Fin 99,
    smallBirthTwoShard7Histogram22.getD k.val 0 =
      ∑ b : Fin 16, if (7 : Fin 16) < b then
        (smallBirthPairHistograms22 7 b).getD k.val 0 else 0 := by decide +kernel

theorem smallBirthTwoShard8Sum22_checked : ∀ k : Fin 99,
    smallBirthTwoShard8Histogram22.getD k.val 0 =
      ∑ b : Fin 16, if (8 : Fin 16) < b then
        (smallBirthPairHistograms22 8 b).getD k.val 0 else 0 := by decide +kernel

theorem smallBirthTwoShard9Sum22_checked : ∀ k : Fin 99,
    smallBirthTwoShard9Histogram22.getD k.val 0 =
      ∑ b : Fin 16, if (9 : Fin 16) < b then
        (smallBirthPairHistograms22 9 b).getD k.val 0 else 0 := by decide +kernel

theorem smallBirthTwoShard10Sum22_checked : ∀ k : Fin 99,
    smallBirthTwoShard10Histogram22.getD k.val 0 =
      ∑ b : Fin 16, if (10 : Fin 16) < b then
        (smallBirthPairHistograms22 10 b).getD k.val 0 else 0 := by decide +kernel

theorem smallBirthTwoShard11Sum22_checked : ∀ k : Fin 99,
    smallBirthTwoShard11Histogram22.getD k.val 0 =
      ∑ b : Fin 16, if (11 : Fin 16) < b then
        (smallBirthPairHistograms22 11 b).getD k.val 0 else 0 := by decide +kernel

theorem smallBirthTwoShard12Sum22_checked : ∀ k : Fin 99,
    smallBirthTwoShard12Histogram22.getD k.val 0 =
      ∑ b : Fin 16, if (12 : Fin 16) < b then
        (smallBirthPairHistograms22 12 b).getD k.val 0 else 0 := by decide +kernel

theorem smallBirthTwoShard13Sum22_checked : ∀ k : Fin 99,
    smallBirthTwoShard13Histogram22.getD k.val 0 =
      ∑ b : Fin 16, if (13 : Fin 16) < b then
        (smallBirthPairHistograms22 13 b).getD k.val 0 else 0 := by decide +kernel

theorem smallBirthTwoShard14Sum22_checked : ∀ k : Fin 99,
    smallBirthTwoShard14Histogram22.getD k.val 0 =
      ∑ b : Fin 16, if (14 : Fin 16) < b then
        (smallBirthPairHistograms22 14 b).getD k.val 0 else 0 := by decide +kernel

theorem smallBirthTwoShard15Sum22_checked : ∀ k : Fin 99,
    smallBirthTwoShard15Histogram22.getD k.val 0 =
      ∑ b : Fin 16, if (15 : Fin 16) < b then
        (smallBirthPairHistograms22 15 b).getD k.val 0 else 0 := by decide +kernel

theorem smallBirthTwoShardHistograms22_sum_checked :
    ∀ (a : Fin 16) (k : Fin 99),
      (smallBirthTwoShardHistograms22 a).getD k.val 0 =
        ∑ b : Fin 16, if a < b then (smallBirthPairHistograms22 a b).getD k.val 0 else 0 := by
  intro a k
  fin_cases a
  · exact smallBirthTwoShard0Sum22_checked k
  · exact smallBirthTwoShard1Sum22_checked k
  · exact smallBirthTwoShard2Sum22_checked k
  · exact smallBirthTwoShard3Sum22_checked k
  · exact smallBirthTwoShard4Sum22_checked k
  · exact smallBirthTwoShard5Sum22_checked k
  · exact smallBirthTwoShard6Sum22_checked k
  · exact smallBirthTwoShard7Sum22_checked k
  · exact smallBirthTwoShard8Sum22_checked k
  · exact smallBirthTwoShard9Sum22_checked k
  · exact smallBirthTwoShard10Sum22_checked k
  · exact smallBirthTwoShard11Sum22_checked k
  · exact smallBirthTwoShard12Sum22_checked k
  · exact smallBirthTwoShard13Sum22_checked k
  · exact smallBirthTwoShard14Sum22_checked k
  · exact smallBirthTwoShard15Sum22_checked k

theorem smallBirthTwoHistogram22_sum_checked :
    (List.ofFn (fun b : Fin 99 => ∑ a : Fin 16,
      (smallBirthTwoShardHistograms22 a).getD b.val 0)) =
        smallBirthTwoHistogram22 := by decide +kernel

end Spin.Turbo
