import SpinCodes.Majorant.Segment18.Part0
import SpinCodes.Majorant.Segment18.Part1
import SpinCodes.Numeric.MajorantEval
set_option maxRecDepth 8000000
set_option maxHeartbeats 0
namespace Spin.Majorant.S18
open Spin.Numeric Spin.Numeric.Fix
theorem part0 {s i : ℝ} (hs : Mem P0.S s) (hi : Mem P0.I i)
    (a b w : ℝ) (ha : Mem P0.A a) (hb : Mem P0.B b) (hw : Mem P0.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P0.logs_ok) hs hi P0.tree P0.tree_ok ha hb hw
theorem part1 {s i : ℝ} (hs : Mem P1.S s) (hi : Mem P1.I i)
    (a b w : ℝ) (ha : Mem P1.A a) (hb : Mem P1.B b) (hw : Mem P1.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P1.logs_ok) hs hi P1.tree P1.tree_ok ha hb hw
theorem cover {s i : ℝ} (hs : Mem ⟨10000000000000000000000000000, 10000000000000000000000000000⟩ s) (hi : Mem ⟨341596090227889668594869802323, 341596090227889668594869802324⟩ i)
    (a b w : ℝ) (ha : Mem ⟨0, 1000000000000000000000000000000⟩ a) (hb : Mem ⟨0, 1000000000000000000000000000000⟩ b) (hw : Mem ⟨495000208320695308716712518240, 497750005208303140513019767608⟩ w) :
    MajorantClaim s i a b w :=
  (claim_split_fst (A := ⟨0, 1000000000000000000000000000000⟩) (B := ⟨0, 1000000000000000000000000000000⟩) (W := ⟨495000208320695308716712518240, 497750005208303140513019767608⟩) (m := 500000000000000000000000000000) (part0 hs hi) (part1 hs hi)) a b w ha hb hw
end Spin.Majorant.S18
