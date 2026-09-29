import SpinCodes.Structured.LowCancellationData.Weight1Block0
import SpinCodes.Structured.LowCancellationBridge
namespace Spin.Structured.LowCancellation.Data
set_option maxRecDepth 100000
set_option maxHeartbeats 0
def groups1 : List Group := Weight1Block0
theorem groups1_valid : ∀ g ∈ groups1, valid 1 g := by
  simp only [groups1, List.mem_append, or_imp, forall_and, and_assoc]
  exact Weight1Block0_valid
theorem groups1_patterns : ∀ g ∈ groups1, g.entries.map Prod.snd ∈
    ((SparsePolynomial.Data.weight1).lowPatterns.getD []).map expand := by
  simp only [groups1, List.mem_append, or_imp, forall_and, and_assoc]
  exact Weight1Block0_patterns
theorem groups1_length : (inputs groups1).length = Nat.choose 128 1 := by decide +kernel
theorem groups1_inputs : (inputs groups1).Nodup :=
  sorted_nodup _ (by decide +kernel)
theorem groups1_keys : (groups1.map Group.syndrome).Nodup := by
  apply List.Pairwise.imp (fun h => Nat.ne_of_lt h)
  apply List.isChain_iff_pairwise.mp
  decide +kernel
theorem groups1_shell0 : certSort (shellExponents groups1 48) =
    expand ((SparsePolynomial.Data.weight1).lowShells.getD 0 []) := by decide +kernel
theorem groups1_shell1 : certSort (shellExponents groups1 56) =
    expand ((SparsePolynomial.Data.weight1).lowShells.getD 1 []) := by decide +kernel
theorem groups1_shell2 : certSort (shellExponents groups1 64) =
    expand ((SparsePolynomial.Data.weight1).lowShells.getD 2 []) := by decide +kernel
theorem groups1_shell3 : certSort (shellExponents groups1 72) =
    expand ((SparsePolynomial.Data.weight1).lowShells.getD 3 []) := by decide +kernel
theorem groups1_shell4 : certSort (shellExponents groups1 80) =
    expand ((SparsePolynomial.Data.weight1).lowShells.getD 4 []) := by decide +kernel
end Spin.Structured.LowCancellation.Data
