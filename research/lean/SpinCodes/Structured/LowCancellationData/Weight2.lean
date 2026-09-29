import SpinCodes.Structured.LowCancellationData.Weight2Block0
import SpinCodes.Structured.LowCancellationData.Weight2Block1
import SpinCodes.Structured.LowCancellationData.Weight2Block2
import SpinCodes.Structured.LowCancellationData.Weight2Block3
import SpinCodes.Structured.LowCancellationData.Weight2Block4
import SpinCodes.Structured.LowCancellationData.Weight2Block5
import SpinCodes.Structured.LowCancellationData.Weight2Block6
import SpinCodes.Structured.LowCancellationData.Weight2Block7
import SpinCodes.Structured.LowCancellationData.Weight2Block8
import SpinCodes.Structured.LowCancellationData.Weight2Block9
import SpinCodes.Structured.LowCancellationData.Weight2Block10
import SpinCodes.Structured.LowCancellationData.Weight2Block11
import SpinCodes.Structured.LowCancellationData.Weight2Block12
import SpinCodes.Structured.LowCancellationData.Weight2Block13
import SpinCodes.Structured.LowCancellationData.Weight2Block14
import SpinCodes.Structured.LowCancellationData.Weight2Block15
import SpinCodes.Structured.LowCancellationData.Weight2Block16
import SpinCodes.Structured.LowCancellationData.Weight2Block17
import SpinCodes.Structured.LowCancellationData.Weight2Block18
import SpinCodes.Structured.LowCancellationData.Weight2Block19
import SpinCodes.Structured.LowCancellationData.Weight2Block20
import SpinCodes.Structured.LowCancellationData.Weight2Block21
import SpinCodes.Structured.LowCancellationData.Weight2Block22
import SpinCodes.Structured.LowCancellationData.Weight2Block23
import SpinCodes.Structured.LowCancellationData.Weight2Block24
import SpinCodes.Structured.LowCancellationData.Weight2Block25
import SpinCodes.Structured.LowCancellationData.Weight2Block26
import SpinCodes.Structured.LowCancellationData.Weight2Block27
import SpinCodes.Structured.LowCancellationData.Weight2Block28
import SpinCodes.Structured.LowCancellationData.Weight2Block29
import SpinCodes.Structured.LowCancellationData.Weight2Block30
import SpinCodes.Structured.LowCancellationData.Weight2Block31
import SpinCodes.Structured.LowCancellationData.Weight2Block32
import SpinCodes.Structured.LowCancellationData.Weight2Block33
import SpinCodes.Structured.LowCancellationData.Weight2Block34
import SpinCodes.Structured.LowCancellationData.Weight2Block35
import SpinCodes.Structured.LowCancellationData.Weight2Block36
import SpinCodes.Structured.LowCancellationData.Weight2Block37
import SpinCodes.Structured.LowCancellationData.Weight2Block38
import SpinCodes.Structured.LowCancellationData.Weight2Block39
import SpinCodes.Structured.LowCancellationData.Weight2Block40
import SpinCodes.Structured.LowCancellationData.Weight2Block41
import SpinCodes.Structured.LowCancellationData.Weight2Block42
import SpinCodes.Structured.LowCancellationData.Weight2Block43
import SpinCodes.Structured.LowCancellationData.Weight2Block44
import SpinCodes.Structured.LowCancellationData.Weight2Block45
import SpinCodes.Structured.LowCancellationData.Weight2Block46
import SpinCodes.Structured.LowCancellationData.Weight2Block47
import SpinCodes.Structured.LowCancellationData.Weight2Block48
import SpinCodes.Structured.LowCancellationData.Weight2Block49
import SpinCodes.Structured.LowCancellationData.Weight2Block50
import SpinCodes.Structured.LowCancellationData.Weight2Block51
import SpinCodes.Structured.LowCancellationData.Weight2Block52
import SpinCodes.Structured.LowCancellationData.Weight2Block53
import SpinCodes.Structured.LowCancellationData.Weight2Block54
import SpinCodes.Structured.LowCancellationData.Weight2Block55
import SpinCodes.Structured.LowCancellationData.Weight2Block56
import SpinCodes.Structured.LowCancellationData.Weight2Block57
import SpinCodes.Structured.LowCancellationData.Weight2Block58
import SpinCodes.Structured.LowCancellationData.Weight2Block59
import SpinCodes.Structured.LowCancellationData.Weight2Block60
import SpinCodes.Structured.LowCancellationData.Weight2Block61
import SpinCodes.Structured.LowCancellationBridge
namespace Spin.Structured.LowCancellation.Data
set_option maxRecDepth 100000
set_option maxHeartbeats 0
def groups2 : List Group := Weight2Block0 ++ Weight2Block1 ++ Weight2Block2 ++ Weight2Block3 ++ Weight2Block4 ++ Weight2Block5 ++ Weight2Block6 ++ Weight2Block7 ++ Weight2Block8 ++ Weight2Block9 ++ Weight2Block10 ++ Weight2Block11 ++ Weight2Block12 ++ Weight2Block13 ++ Weight2Block14 ++ Weight2Block15 ++ Weight2Block16 ++ Weight2Block17 ++ Weight2Block18 ++ Weight2Block19 ++ Weight2Block20 ++ Weight2Block21 ++ Weight2Block22 ++ Weight2Block23 ++ Weight2Block24 ++ Weight2Block25 ++ Weight2Block26 ++ Weight2Block27 ++ Weight2Block28 ++ Weight2Block29 ++ Weight2Block30 ++ Weight2Block31 ++ Weight2Block32 ++ Weight2Block33 ++ Weight2Block34 ++ Weight2Block35 ++ Weight2Block36 ++ Weight2Block37 ++ Weight2Block38 ++ Weight2Block39 ++ Weight2Block40 ++ Weight2Block41 ++ Weight2Block42 ++ Weight2Block43 ++ Weight2Block44 ++ Weight2Block45 ++ Weight2Block46 ++ Weight2Block47 ++ Weight2Block48 ++ Weight2Block49 ++ Weight2Block50 ++ Weight2Block51 ++ Weight2Block52 ++ Weight2Block53 ++ Weight2Block54 ++ Weight2Block55 ++ Weight2Block56 ++ Weight2Block57 ++ Weight2Block58 ++ Weight2Block59 ++ Weight2Block60 ++ Weight2Block61
theorem groups2_valid : ∀ g ∈ groups2, valid 2 g := by
  simp only [groups2, List.mem_append, or_imp, forall_and, and_assoc]
  exact ⟨Weight2Block0_valid, Weight2Block1_valid, Weight2Block2_valid, Weight2Block3_valid, Weight2Block4_valid, Weight2Block5_valid, Weight2Block6_valid, Weight2Block7_valid, Weight2Block8_valid, Weight2Block9_valid, Weight2Block10_valid, Weight2Block11_valid, Weight2Block12_valid, Weight2Block13_valid, Weight2Block14_valid, Weight2Block15_valid, Weight2Block16_valid, Weight2Block17_valid, Weight2Block18_valid, Weight2Block19_valid, Weight2Block20_valid, Weight2Block21_valid, Weight2Block22_valid, Weight2Block23_valid, Weight2Block24_valid, Weight2Block25_valid, Weight2Block26_valid, Weight2Block27_valid, Weight2Block28_valid, Weight2Block29_valid, Weight2Block30_valid, Weight2Block31_valid, Weight2Block32_valid, Weight2Block33_valid, Weight2Block34_valid, Weight2Block35_valid, Weight2Block36_valid, Weight2Block37_valid, Weight2Block38_valid, Weight2Block39_valid, Weight2Block40_valid, Weight2Block41_valid, Weight2Block42_valid, Weight2Block43_valid, Weight2Block44_valid, Weight2Block45_valid, Weight2Block46_valid, Weight2Block47_valid, Weight2Block48_valid, Weight2Block49_valid, Weight2Block50_valid, Weight2Block51_valid, Weight2Block52_valid, Weight2Block53_valid, Weight2Block54_valid, Weight2Block55_valid, Weight2Block56_valid, Weight2Block57_valid, Weight2Block58_valid, Weight2Block59_valid, Weight2Block60_valid, Weight2Block61_valid⟩
theorem groups2_patterns : ∀ g ∈ groups2, g.entries.map Prod.snd ∈
    ((SparsePolynomial.Data.weight2).lowPatterns.getD []).map expand := by
  simp only [groups2, List.mem_append, or_imp, forall_and, and_assoc]
  exact ⟨Weight2Block0_patterns, Weight2Block1_patterns, Weight2Block2_patterns, Weight2Block3_patterns, Weight2Block4_patterns, Weight2Block5_patterns, Weight2Block6_patterns, Weight2Block7_patterns, Weight2Block8_patterns, Weight2Block9_patterns, Weight2Block10_patterns, Weight2Block11_patterns, Weight2Block12_patterns, Weight2Block13_patterns, Weight2Block14_patterns, Weight2Block15_patterns, Weight2Block16_patterns, Weight2Block17_patterns, Weight2Block18_patterns, Weight2Block19_patterns, Weight2Block20_patterns, Weight2Block21_patterns, Weight2Block22_patterns, Weight2Block23_patterns, Weight2Block24_patterns, Weight2Block25_patterns, Weight2Block26_patterns, Weight2Block27_patterns, Weight2Block28_patterns, Weight2Block29_patterns, Weight2Block30_patterns, Weight2Block31_patterns, Weight2Block32_patterns, Weight2Block33_patterns, Weight2Block34_patterns, Weight2Block35_patterns, Weight2Block36_patterns, Weight2Block37_patterns, Weight2Block38_patterns, Weight2Block39_patterns, Weight2Block40_patterns, Weight2Block41_patterns, Weight2Block42_patterns, Weight2Block43_patterns, Weight2Block44_patterns, Weight2Block45_patterns, Weight2Block46_patterns, Weight2Block47_patterns, Weight2Block48_patterns, Weight2Block49_patterns, Weight2Block50_patterns, Weight2Block51_patterns, Weight2Block52_patterns, Weight2Block53_patterns, Weight2Block54_patterns, Weight2Block55_patterns, Weight2Block56_patterns, Weight2Block57_patterns, Weight2Block58_patterns, Weight2Block59_patterns, Weight2Block60_patterns, Weight2Block61_patterns⟩
theorem groups2_length : (inputs groups2).length = Nat.choose 128 2 := by decide +kernel
theorem groups2_inputs : (inputs groups2).Nodup :=
  sorted_nodup _ (by decide +kernel)
theorem groups2_keys : (groups2.map Group.syndrome).Nodup := by
  apply List.Pairwise.imp (fun h => Nat.ne_of_lt h)
  apply List.isChain_iff_pairwise.mp
  decide +kernel
theorem groups2_shell0 : certSort (shellExponents groups2 48) =
    expand ((SparsePolynomial.Data.weight2).lowShells.getD 0 []) := by decide +kernel
theorem groups2_shell1 : certSort (shellExponents groups2 56) =
    expand ((SparsePolynomial.Data.weight2).lowShells.getD 1 []) := by decide +kernel
theorem groups2_shell2 : certSort (shellExponents groups2 64) =
    expand ((SparsePolynomial.Data.weight2).lowShells.getD 2 []) := by decide +kernel
theorem groups2_shell3 : certSort (shellExponents groups2 72) =
    expand ((SparsePolynomial.Data.weight2).lowShells.getD 3 []) := by decide +kernel
theorem groups2_shell4 : certSort (shellExponents groups2 80) =
    expand ((SparsePolynomial.Data.weight2).lowShells.getD 4 []) := by decide +kernel
end Spin.Structured.LowCancellation.Data
