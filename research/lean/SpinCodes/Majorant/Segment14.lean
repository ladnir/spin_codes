import SpinCodes.Majorant.Segment14.Part0
import SpinCodes.Majorant.Segment14.Part1
import SpinCodes.Majorant.Segment14.Part2
import SpinCodes.Majorant.Segment14.Part3
import SpinCodes.Numeric.MajorantEval
set_option maxRecDepth 8000000
set_option maxHeartbeats 0
namespace Spin.Majorant.S14
open Spin.Numeric Spin.Numeric.Fix
theorem part0 {s i : ℝ} (hs : Mem P0.S s) (hi : Mem P0.I i)
    (a b w : ℝ) (ha : Mem P0.A a) (hb : Mem P0.B b) (hw : Mem P0.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P0.logs_ok) hs hi P0.tree P0.tree_ok ha hb hw
theorem part1 {s i : ℝ} (hs : Mem P1.S s) (hi : Mem P1.I i)
    (a b w : ℝ) (ha : Mem P1.A a) (hb : Mem P1.B b) (hw : Mem P1.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P1.logs_ok) hs hi P1.tree P1.tree_ok ha hb hw
theorem part2 {s i : ℝ} (hs : Mem P2.S s) (hi : Mem P2.I i)
    (a b w : ℝ) (ha : Mem P2.A a) (hb : Mem P2.B b) (hw : Mem P2.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P2.logs_ok) hs hi P2.tree P2.tree_ok ha hb hw
theorem part3 {s i : ℝ} (hs : Mem P3.S s) (hi : Mem P3.I i)
    (a b w : ℝ) (ha : Mem P3.A a) (hb : Mem P3.B b) (hw : Mem P3.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P3.logs_ok) hs hi P3.tree P3.tree_ok ha hb hw
theorem cover {s i : ℝ} (hs : Mem ⟨90000000000000000000000000000, 90000000000000000000000000000⟩ s) (hi : Mem ⟨302595748745637418755008535073, 302595748745637418755008535074⟩ i)
    (a b w : ℝ) (ha : Mem ⟨0, 1000000000000000000000000000000⟩ a) (hb : Mem ⟨0, 1000000000000000000000000000000⟩ b) (hw : Mem ⟨475460124185441875500853507329, 480010826154353232496514123203⟩ w) :
    MajorantClaim s i a b w :=
  (claim_split_fst (A := ⟨0, 1000000000000000000000000000000⟩) (B := ⟨0, 1000000000000000000000000000000⟩) (W := ⟨475460124185441875500853507329, 480010826154353232496514123203⟩) (m := 500000000000000000000000000000) (claim_split_snd (A := ⟨0, 500000000000000000000000000000⟩) (B := ⟨0, 1000000000000000000000000000000⟩) (W := ⟨475460124185441875500853507329, 480010826154353232496514123203⟩) (m := 500000000000000000000000000000) (part0 hs hi) (part1 hs hi)) (claim_split_snd (A := ⟨500000000000000000000000000000, 1000000000000000000000000000000⟩) (B := ⟨0, 1000000000000000000000000000000⟩) (W := ⟨475460124185441875500853507329, 480010826154353232496514123203⟩) (m := 500000000000000000000000000000) (part2 hs hi) (part3 hs hi))) a b w ha hb hw
end Spin.Majorant.S14
