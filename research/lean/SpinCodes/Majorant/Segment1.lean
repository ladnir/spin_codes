import SpinCodes.Majorant.Segment1.Part0
import SpinCodes.Majorant.Segment1.Part1
import SpinCodes.Majorant.Segment1.Part2
import SpinCodes.Majorant.Segment1.Part3
import SpinCodes.Majorant.Segment1.Part4
import SpinCodes.Majorant.Segment1.Part5
import SpinCodes.Majorant.Segment1.Part6
import SpinCodes.Majorant.Segment1.Part7
import SpinCodes.Majorant.Segment1.Part8
import SpinCodes.Majorant.Segment1.Part9
import SpinCodes.Majorant.Segment1.Part10
import SpinCodes.Majorant.Segment1.Part11
import SpinCodes.Majorant.Segment1.Part12
import SpinCodes.Majorant.Segment1.Part13
import SpinCodes.Majorant.Segment1.Part14
import SpinCodes.Majorant.Segment1.Part15
import SpinCodes.Numeric.MajorantEval
set_option maxRecDepth 8000000
set_option maxHeartbeats 0
namespace Spin.Majorant.S1
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
theorem part4 {s i : ℝ} (hs : Mem P4.S s) (hi : Mem P4.I i)
    (a b w : ℝ) (ha : Mem P4.A a) (hb : Mem P4.B b) (hw : Mem P4.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P4.logs_ok) hs hi P4.tree P4.tree_ok ha hb hw
theorem part5 {s i : ℝ} (hs : Mem P5.S s) (hi : Mem P5.I i)
    (a b w : ℝ) (ha : Mem P5.A a) (hb : Mem P5.B b) (hw : Mem P5.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P5.logs_ok) hs hi P5.tree P5.tree_ok ha hb hw
theorem part6 {s i : ℝ} (hs : Mem P6.S s) (hi : Mem P6.I i)
    (a b w : ℝ) (ha : Mem P6.A a) (hb : Mem P6.B b) (hw : Mem P6.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P6.logs_ok) hs hi P6.tree P6.tree_ok ha hb hw
theorem part7 {s i : ℝ} (hs : Mem P7.S s) (hi : Mem P7.I i)
    (a b w : ℝ) (ha : Mem P7.A a) (hb : Mem P7.B b) (hw : Mem P7.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P7.logs_ok) hs hi P7.tree P7.tree_ok ha hb hw
theorem part8 {s i : ℝ} (hs : Mem P8.S s) (hi : Mem P8.I i)
    (a b w : ℝ) (ha : Mem P8.A a) (hb : Mem P8.B b) (hw : Mem P8.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P8.logs_ok) hs hi P8.tree P8.tree_ok ha hb hw
theorem part9 {s i : ℝ} (hs : Mem P9.S s) (hi : Mem P9.I i)
    (a b w : ℝ) (ha : Mem P9.A a) (hb : Mem P9.B b) (hw : Mem P9.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P9.logs_ok) hs hi P9.tree P9.tree_ok ha hb hw
theorem part10 {s i : ℝ} (hs : Mem P10.S s) (hi : Mem P10.I i)
    (a b w : ℝ) (ha : Mem P10.A a) (hb : Mem P10.B b) (hw : Mem P10.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P10.logs_ok) hs hi P10.tree P10.tree_ok ha hb hw
theorem part11 {s i : ℝ} (hs : Mem P11.S s) (hi : Mem P11.I i)
    (a b w : ℝ) (ha : Mem P11.A a) (hb : Mem P11.B b) (hw : Mem P11.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P11.logs_ok) hs hi P11.tree P11.tree_ok ha hb hw
theorem part12 {s i : ℝ} (hs : Mem P12.S s) (hi : Mem P12.I i)
    (a b w : ℝ) (ha : Mem P12.A a) (hb : Mem P12.B b) (hw : Mem P12.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P12.logs_ok) hs hi P12.tree P12.tree_ok ha hb hw
theorem part13 {s i : ℝ} (hs : Mem P13.S s) (hi : Mem P13.I i)
    (a b w : ℝ) (ha : Mem P13.A a) (hb : Mem P13.B b) (hw : Mem P13.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P13.logs_ok) hs hi P13.tree P13.tree_ok ha hb hw
theorem part14 {s i : ℝ} (hs : Mem P14.S s) (hi : Mem P14.I i)
    (a b w : ℝ) (ha : Mem P14.A a) (hb : Mem P14.B b) (hw : Mem P14.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P14.logs_ok) hs hi P14.tree P14.tree_ok ha hb hw
theorem part15 {s i : ℝ} (hs : Mem P15.S s) (hi : Mem P15.I i)
    (a b w : ℝ) (ha : Mem P15.A a) (hb : Mem P15.B b) (hw : Mem P15.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P15.logs_ok) hs hi P15.tree P15.tree_ok ha hb hw
theorem cover {s i : ℝ} (hs : Mem ⟨1660000000000000000000000000000, 1660000000000000000000000000000⟩ s) (hi : Mem ⟨-172329198693950000000000000000, -172329198693950000000000000000⟩ i)
    (a b w : ℝ) (ha : Mem ⟨0, 1000000000000000000000000000000⟩ a) (hb : Mem ⟨0, 1000000000000000000000000000000⟩ b) (hw : Mem ⟨152466884341666666666666666666, 158319841724600000000000000000⟩ w) :
    MajorantClaim s i a b w :=
  (claim_split_fst (A := ⟨0, 1000000000000000000000000000000⟩) (B := ⟨0, 1000000000000000000000000000000⟩) (W := ⟨152466884341666666666666666666, 158319841724600000000000000000⟩) (m := 500000000000000000000000000000) (claim_split_snd (A := ⟨0, 500000000000000000000000000000⟩) (B := ⟨0, 1000000000000000000000000000000⟩) (W := ⟨152466884341666666666666666666, 158319841724600000000000000000⟩) (m := 500000000000000000000000000000) (claim_split_fst (A := ⟨0, 500000000000000000000000000000⟩) (B := ⟨0, 500000000000000000000000000000⟩) (W := ⟨152466884341666666666666666666, 158319841724600000000000000000⟩) (m := 250000000000000000000000000000) (part0 hs hi) (claim_split_snd (A := ⟨250000000000000000000000000000, 500000000000000000000000000000⟩) (B := ⟨0, 500000000000000000000000000000⟩) (W := ⟨152466884341666666666666666666, 158319841724600000000000000000⟩) (m := 250000000000000000000000000000) (part1 hs hi) (claim_split_fst (A := ⟨250000000000000000000000000000, 500000000000000000000000000000⟩) (B := ⟨250000000000000000000000000000, 500000000000000000000000000000⟩) (W := ⟨152466884341666666666666666666, 158319841724600000000000000000⟩) (m := 375000000000000000000000000000) (part2 hs hi) (claim_split_snd (A := ⟨375000000000000000000000000000, 500000000000000000000000000000⟩) (B := ⟨250000000000000000000000000000, 500000000000000000000000000000⟩) (W := ⟨152466884341666666666666666666, 158319841724600000000000000000⟩) (m := 375000000000000000000000000000) (claim_split_fst (A := ⟨375000000000000000000000000000, 500000000000000000000000000000⟩) (B := ⟨250000000000000000000000000000, 375000000000000000000000000000⟩) (W := ⟨152466884341666666666666666666, 158319841724600000000000000000⟩) (m := 437500000000000000000000000000) (claim_split_snd (A := ⟨375000000000000000000000000000, 437500000000000000000000000000⟩) (B := ⟨250000000000000000000000000000, 375000000000000000000000000000⟩) (W := ⟨152466884341666666666666666666, 158319841724600000000000000000⟩) (m := 312500000000000000000000000000) (claim_split_fst (A := ⟨375000000000000000000000000000, 437500000000000000000000000000⟩) (B := ⟨250000000000000000000000000000, 312500000000000000000000000000⟩) (W := ⟨152466884341666666666666666666, 158319841724600000000000000000⟩) (m := 406250000000000000000000000000) (claim_split_snd (A := ⟨375000000000000000000000000000, 406250000000000000000000000000⟩) (B := ⟨250000000000000000000000000000, 312500000000000000000000000000⟩) (W := ⟨152466884341666666666666666666, 158319841724600000000000000000⟩) (m := 281250000000000000000000000000) (claim_split_fst (A := ⟨375000000000000000000000000000, 406250000000000000000000000000⟩) (B := ⟨250000000000000000000000000000, 281250000000000000000000000000⟩) (W := ⟨152466884341666666666666666666, 158319841724600000000000000000⟩) (m := 390625000000000000000000000000) (claim_split_snd (A := ⟨375000000000000000000000000000, 390625000000000000000000000000⟩) (B := ⟨250000000000000000000000000000, 281250000000000000000000000000⟩) (W := ⟨152466884341666666666666666666, 158319841724600000000000000000⟩) (m := 265625000000000000000000000000) (claim_split_fst (A := ⟨375000000000000000000000000000, 390625000000000000000000000000⟩) (B := ⟨250000000000000000000000000000, 265625000000000000000000000000⟩) (W := ⟨152466884341666666666666666666, 158319841724600000000000000000⟩) (m := 382812500000000000000000000000) (part3 hs hi) (claim_split_snd (A := ⟨382812500000000000000000000000, 390625000000000000000000000000⟩) (B := ⟨250000000000000000000000000000, 265625000000000000000000000000⟩) (W := ⟨152466884341666666666666666666, 158319841724600000000000000000⟩) (m := 257812500000000000000000000000) (part4 hs hi) (claim_split_thd (A := ⟨382812500000000000000000000000, 390625000000000000000000000000⟩) (B := ⟨257812500000000000000000000000, 265625000000000000000000000000⟩) (W := ⟨152466884341666666666666666666, 158319841724600000000000000000⟩) (m := 155393363033133333333333333333) (part5 hs hi) (part6 hs hi)))) (part7 hs hi)) (part8 hs hi)) (part9 hs hi)) (part10 hs hi)) (part11 hs hi)) (part12 hs hi)) (part13 hs hi))))) (part14 hs hi)) (part15 hs hi)) a b w ha hb hw
end Spin.Majorant.S1
