import SpinCodes.Majorant.Segment10.Part0
import SpinCodes.Majorant.Segment10.Part1
import SpinCodes.Majorant.Segment10.Part2
import SpinCodes.Majorant.Segment10.Part3
import SpinCodes.Majorant.Segment10.Part4
import SpinCodes.Majorant.Segment10.Part5
import SpinCodes.Majorant.Segment10.Part6
import SpinCodes.Majorant.Segment10.Part7
import SpinCodes.Majorant.Segment10.Part8
import SpinCodes.Majorant.Segment10.Part9
import SpinCodes.Majorant.Segment10.Part10
import SpinCodes.Majorant.Segment10.Part11
import SpinCodes.Majorant.Segment10.Part12
import SpinCodes.Majorant.Segment10.Part13
import SpinCodes.Majorant.Segment10.Part14
import SpinCodes.Majorant.Segment10.Part15
import SpinCodes.Majorant.Segment10.Part16
import SpinCodes.Majorant.Segment10.Part17
import SpinCodes.Majorant.Segment10.Part18
import SpinCodes.Majorant.Segment10.Part19
import SpinCodes.Majorant.Segment10.Part20
import SpinCodes.Majorant.Segment10.Part21
import SpinCodes.Majorant.Segment10.Part22
import SpinCodes.Majorant.Segment10.Part23
import SpinCodes.Majorant.Segment10.Part24
import SpinCodes.Majorant.Segment10.Part25
import SpinCodes.Majorant.Segment10.Part26
import SpinCodes.Majorant.Segment10.Part27
import SpinCodes.Numeric.MajorantEval
set_option maxRecDepth 8000000
set_option maxHeartbeats 0
namespace Spin.Majorant.S10
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
theorem part16 {s i : ℝ} (hs : Mem P16.S s) (hi : Mem P16.I i)
    (a b w : ℝ) (ha : Mem P16.A a) (hb : Mem P16.B b) (hw : Mem P16.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P16.logs_ok) hs hi P16.tree P16.tree_ok ha hb hw
theorem part17 {s i : ℝ} (hs : Mem P17.S s) (hi : Mem P17.I i)
    (a b w : ℝ) (ha : Mem P17.A a) (hb : Mem P17.B b) (hw : Mem P17.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P17.logs_ok) hs hi P17.tree P17.tree_ok ha hb hw
theorem part18 {s i : ℝ} (hs : Mem P18.S s) (hi : Mem P18.I i)
    (a b w : ℝ) (ha : Mem P18.A a) (hb : Mem P18.B b) (hw : Mem P18.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P18.logs_ok) hs hi P18.tree P18.tree_ok ha hb hw
theorem part19 {s i : ℝ} (hs : Mem P19.S s) (hi : Mem P19.I i)
    (a b w : ℝ) (ha : Mem P19.A a) (hb : Mem P19.B b) (hw : Mem P19.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P19.logs_ok) hs hi P19.tree P19.tree_ok ha hb hw
theorem part20 {s i : ℝ} (hs : Mem P20.S s) (hi : Mem P20.I i)
    (a b w : ℝ) (ha : Mem P20.A a) (hb : Mem P20.B b) (hw : Mem P20.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P20.logs_ok) hs hi P20.tree P20.tree_ok ha hb hw
theorem part21 {s i : ℝ} (hs : Mem P21.S s) (hi : Mem P21.I i)
    (a b w : ℝ) (ha : Mem P21.A a) (hb : Mem P21.B b) (hw : Mem P21.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P21.logs_ok) hs hi P21.tree P21.tree_ok ha hb hw
theorem part22 {s i : ℝ} (hs : Mem P22.S s) (hi : Mem P22.I i)
    (a b w : ℝ) (ha : Mem P22.A a) (hb : Mem P22.B b) (hw : Mem P22.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P22.logs_ok) hs hi P22.tree P22.tree_ok ha hb hw
theorem part23 {s i : ℝ} (hs : Mem P23.S s) (hi : Mem P23.I i)
    (a b w : ℝ) (ha : Mem P23.A a) (hb : Mem P23.B b) (hw : Mem P23.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P23.logs_ok) hs hi P23.tree P23.tree_ok ha hb hw
theorem part24 {s i : ℝ} (hs : Mem P24.S s) (hi : Mem P24.I i)
    (a b w : ℝ) (ha : Mem P24.A a) (hb : Mem P24.B b) (hw : Mem P24.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P24.logs_ok) hs hi P24.tree P24.tree_ok ha hb hw
theorem part25 {s i : ℝ} (hs : Mem P25.S s) (hi : Mem P25.I i)
    (a b w : ℝ) (ha : Mem P25.A a) (hb : Mem P25.B b) (hw : Mem P25.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P25.logs_ok) hs hi P25.tree P25.tree_ok ha hb hw
theorem part26 {s i : ℝ} (hs : Mem P26.S s) (hi : Mem P26.I i)
    (a b w : ℝ) (ha : Mem P26.A a) (hb : Mem P26.B b) (hw : Mem P26.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P26.logs_ok) hs hi P26.tree P26.tree_ok ha hb hw
theorem part27 {s i : ℝ} (hs : Mem P27.S s) (hi : Mem P27.I i)
    (a b w : ℝ) (ha : Mem P27.A a) (hb : Mem P27.B b) (hw : Mem P27.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P27.logs_ok) hs hi P27.tree P27.tree_ok ha hb hw
theorem cover {s i : ℝ} (hs : Mem ⟨600000000000000000000000000000, 600000000000000000000000000000⟩ s) (hi : Mem ⟨90934102009248100000000000000, 90934102009248100000000000000⟩ i)
    (a b w : ℝ) (ha : Mem ⟨0, 1000000000000000000000000000000⟩ a) (hb : Mem ⟨0, 1000000000000000000000000000000⟩ b) (hw : Mem ⟨331935139312516500000000000000, 377632881027959500000000000000⟩ w) :
    MajorantClaim s i a b w :=
  (claim_split_fst (A := ⟨0, 1000000000000000000000000000000⟩) (B := ⟨0, 1000000000000000000000000000000⟩) (W := ⟨331935139312516500000000000000, 377632881027959500000000000000⟩) (m := 500000000000000000000000000000) (claim_split_snd (A := ⟨0, 500000000000000000000000000000⟩) (B := ⟨0, 1000000000000000000000000000000⟩) (W := ⟨331935139312516500000000000000, 377632881027959500000000000000⟩) (m := 500000000000000000000000000000) (claim_split_fst (A := ⟨0, 500000000000000000000000000000⟩) (B := ⟨0, 500000000000000000000000000000⟩) (W := ⟨331935139312516500000000000000, 377632881027959500000000000000⟩) (m := 250000000000000000000000000000) (part0 hs hi) (claim_split_snd (A := ⟨250000000000000000000000000000, 500000000000000000000000000000⟩) (B := ⟨0, 500000000000000000000000000000⟩) (W := ⟨331935139312516500000000000000, 377632881027959500000000000000⟩) (m := 250000000000000000000000000000) (part1 hs hi) (claim_split_fst (A := ⟨250000000000000000000000000000, 500000000000000000000000000000⟩) (B := ⟨250000000000000000000000000000, 500000000000000000000000000000⟩) (W := ⟨331935139312516500000000000000, 377632881027959500000000000000⟩) (m := 375000000000000000000000000000) (part2 hs hi) (claim_split_snd (A := ⟨375000000000000000000000000000, 500000000000000000000000000000⟩) (B := ⟨250000000000000000000000000000, 500000000000000000000000000000⟩) (W := ⟨331935139312516500000000000000, 377632881027959500000000000000⟩) (m := 375000000000000000000000000000) (part3 hs hi) (claim_split_fst (A := ⟨375000000000000000000000000000, 500000000000000000000000000000⟩) (B := ⟨375000000000000000000000000000, 500000000000000000000000000000⟩) (W := ⟨331935139312516500000000000000, 377632881027959500000000000000⟩) (m := 437500000000000000000000000000) (part4 hs hi) (claim_split_snd (A := ⟨437500000000000000000000000000, 500000000000000000000000000000⟩) (B := ⟨375000000000000000000000000000, 500000000000000000000000000000⟩) (W := ⟨331935139312516500000000000000, 377632881027959500000000000000⟩) (m := 437500000000000000000000000000) (part5 hs hi) (claim_split_thd (A := ⟨437500000000000000000000000000, 500000000000000000000000000000⟩) (B := ⟨437500000000000000000000000000, 500000000000000000000000000000⟩) (W := ⟨331935139312516500000000000000, 377632881027959500000000000000⟩) (m := 354784010170238000000000000000) (claim_split_fst (A := ⟨437500000000000000000000000000, 500000000000000000000000000000⟩) (B := ⟨437500000000000000000000000000, 500000000000000000000000000000⟩) (W := ⟨331935139312516500000000000000, 354784010170238000000000000000⟩) (m := 468750000000000000000000000000) (part6 hs hi) (claim_split_snd (A := ⟨468750000000000000000000000000, 500000000000000000000000000000⟩) (B := ⟨437500000000000000000000000000, 500000000000000000000000000000⟩) (W := ⟨331935139312516500000000000000, 354784010170238000000000000000⟩) (m := 468750000000000000000000000000) (claim_split_thd (A := ⟨468750000000000000000000000000, 500000000000000000000000000000⟩) (B := ⟨437500000000000000000000000000, 468750000000000000000000000000⟩) (W := ⟨331935139312516500000000000000, 354784010170238000000000000000⟩) (m := 343359574741377250000000000000) (part7 hs hi) (claim_split_fst (A := ⟨468750000000000000000000000000, 500000000000000000000000000000⟩) (B := ⟨437500000000000000000000000000, 468750000000000000000000000000⟩) (W := ⟨343359574741377250000000000000, 354784010170238000000000000000⟩) (m := 484375000000000000000000000000) (part8 hs hi) (claim_split_snd (A := ⟨484375000000000000000000000000, 500000000000000000000000000000⟩) (B := ⟨437500000000000000000000000000, 468750000000000000000000000000⟩) (W := ⟨343359574741377250000000000000, 354784010170238000000000000000⟩) (m := 453125000000000000000000000000) (part9 hs hi) (part10 hs hi)))) (part11 hs hi))) (claim_split_fst (A := ⟨437500000000000000000000000000, 500000000000000000000000000000⟩) (B := ⟨437500000000000000000000000000, 500000000000000000000000000000⟩) (W := ⟨354784010170238000000000000000, 377632881027959500000000000000⟩) (m := 468750000000000000000000000000) (part12 hs hi) (claim_split_snd (A := ⟨468750000000000000000000000000, 500000000000000000000000000000⟩) (B := ⟨437500000000000000000000000000, 500000000000000000000000000000⟩) (W := ⟨354784010170238000000000000000, 377632881027959500000000000000⟩) (m := 468750000000000000000000000000) (claim_split_thd (A := ⟨468750000000000000000000000000, 500000000000000000000000000000⟩) (B := ⟨437500000000000000000000000000, 468750000000000000000000000000⟩) (W := ⟨354784010170238000000000000000, 377632881027959500000000000000⟩) (m := 366208445599098750000000000000) (claim_split_fst (A := ⟨468750000000000000000000000000, 500000000000000000000000000000⟩) (B := ⟨437500000000000000000000000000, 468750000000000000000000000000⟩) (W := ⟨354784010170238000000000000000, 366208445599098750000000000000⟩) (m := 484375000000000000000000000000) (part13 hs hi) (claim_split_snd (A := ⟨484375000000000000000000000000, 500000000000000000000000000000⟩) (B := ⟨437500000000000000000000000000, 468750000000000000000000000000⟩) (W := ⟨354784010170238000000000000000, 366208445599098750000000000000⟩) (m := 453125000000000000000000000000) (part14 hs hi) (part15 hs hi))) (part16 hs hi)) (part17 hs hi)))))))))) (part18 hs hi)) (claim_split_snd (A := ⟨500000000000000000000000000000, 1000000000000000000000000000000⟩) (B := ⟨0, 1000000000000000000000000000000⟩) (W := ⟨331935139312516500000000000000, 377632881027959500000000000000⟩) (m := 500000000000000000000000000000) (claim_split_fst (A := ⟨500000000000000000000000000000, 1000000000000000000000000000000⟩) (B := ⟨0, 500000000000000000000000000000⟩) (W := ⟨331935139312516500000000000000, 377632881027959500000000000000⟩) (m := 750000000000000000000000000000) (claim_split_snd (A := ⟨500000000000000000000000000000, 750000000000000000000000000000⟩) (B := ⟨0, 500000000000000000000000000000⟩) (W := ⟨331935139312516500000000000000, 377632881027959500000000000000⟩) (m := 250000000000000000000000000000) (part19 hs hi) (claim_split_fst (A := ⟨500000000000000000000000000000, 750000000000000000000000000000⟩) (B := ⟨250000000000000000000000000000, 500000000000000000000000000000⟩) (W := ⟨331935139312516500000000000000, 377632881027959500000000000000⟩) (m := 625000000000000000000000000000) (claim_split_snd (A := ⟨500000000000000000000000000000, 625000000000000000000000000000⟩) (B := ⟨250000000000000000000000000000, 500000000000000000000000000000⟩) (W := ⟨331935139312516500000000000000, 377632881027959500000000000000⟩) (m := 375000000000000000000000000000) (part20 hs hi) (claim_split_fst (A := ⟨500000000000000000000000000000, 625000000000000000000000000000⟩) (B := ⟨375000000000000000000000000000, 500000000000000000000000000000⟩) (W := ⟨331935139312516500000000000000, 377632881027959500000000000000⟩) (m := 562500000000000000000000000000) (claim_split_snd (A := ⟨500000000000000000000000000000, 562500000000000000000000000000⟩) (B := ⟨375000000000000000000000000000, 500000000000000000000000000000⟩) (W := ⟨331935139312516500000000000000, 377632881027959500000000000000⟩) (m := 437500000000000000000000000000) (part21 hs hi) (claim_split_thd (A := ⟨500000000000000000000000000000, 562500000000000000000000000000⟩) (B := ⟨437500000000000000000000000000, 500000000000000000000000000000⟩) (W := ⟨331935139312516500000000000000, 377632881027959500000000000000⟩) (m := 354784010170238000000000000000) (part22 hs hi) (part23 hs hi))) (part24 hs hi))) (part25 hs hi))) (part26 hs hi)) (part27 hs hi))) a b w ha hb hw
end Spin.Majorant.S10
