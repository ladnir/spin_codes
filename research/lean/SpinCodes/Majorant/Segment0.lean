import SpinCodes.Majorant.Segment0.Part0
import SpinCodes.Majorant.Segment0.Part1
import SpinCodes.Majorant.Segment0.Part2
import SpinCodes.Majorant.Segment0.Part3
import SpinCodes.Majorant.Segment0.Part4
import SpinCodes.Majorant.Segment0.Part5
import SpinCodes.Majorant.Segment0.Part6
import SpinCodes.Majorant.Segment0.Part7
import SpinCodes.Majorant.Segment0.Part8
import SpinCodes.Majorant.Segment0.Part9
import SpinCodes.Majorant.Segment0.Part10
import SpinCodes.Majorant.Segment0.Part11
import SpinCodes.Majorant.Segment0.Part12
import SpinCodes.Majorant.Segment0.Part13
import SpinCodes.Majorant.Segment0.Part14
import SpinCodes.Majorant.Segment0.Part15
import SpinCodes.Majorant.Segment0.Part16
import SpinCodes.Majorant.Segment0.Part17
import SpinCodes.Majorant.Segment0.Part18
import SpinCodes.Majorant.Segment0.Part19
import SpinCodes.Majorant.Segment0.Part20
import SpinCodes.Majorant.Segment0.Part21
import SpinCodes.Majorant.Segment0.Part22
import SpinCodes.Majorant.Segment0.Part23
import SpinCodes.Majorant.Segment0.Part24
import SpinCodes.Majorant.Segment0.Part25
import SpinCodes.Majorant.Segment0.Part26
import SpinCodes.Majorant.Segment0.Part27
import SpinCodes.Majorant.Segment0.Part28
import SpinCodes.Majorant.Segment0.Part29
import SpinCodes.Majorant.Segment0.Part30
import SpinCodes.Majorant.Segment0.Part31
import SpinCodes.Majorant.Segment0.Part32
import SpinCodes.Majorant.Segment0.Part33
import SpinCodes.Majorant.Segment0.Part34
import SpinCodes.Majorant.Segment0.Part35
import SpinCodes.Majorant.Segment0.Part36
import SpinCodes.Majorant.Segment0.Part37
import SpinCodes.Majorant.Segment0.Part38
import SpinCodes.Majorant.Segment0.Part39
import SpinCodes.Majorant.Segment0.Part40
import SpinCodes.Majorant.Segment0.Part41
import SpinCodes.Majorant.Segment0.Part42
import SpinCodes.Majorant.Segment0.Part43
import SpinCodes.Majorant.Segment0.Part44
import SpinCodes.Majorant.Segment0.Part45
import SpinCodes.Majorant.Segment0.Part46
import SpinCodes.Majorant.Segment0.Part47
import SpinCodes.Majorant.Segment0.Part48
import SpinCodes.Majorant.Segment0.Part49
import SpinCodes.Majorant.Segment0.Part50
import SpinCodes.Majorant.Segment0.Part51
import SpinCodes.Majorant.Segment0.Part52
import SpinCodes.Majorant.Segment0.Part53
import SpinCodes.Majorant.Segment0.Part54
import SpinCodes.Majorant.Segment0.Part55
import SpinCodes.Majorant.Segment0.Part56
import SpinCodes.Majorant.Segment0.Part57
import SpinCodes.Majorant.Segment0.Part58
import SpinCodes.Majorant.Segment0.Part59
import SpinCodes.Majorant.Segment0.Part60
import SpinCodes.Majorant.Segment0.Part61
import SpinCodes.Majorant.Segment0.Part62
import SpinCodes.Majorant.Segment0.Part63
import SpinCodes.Majorant.Segment0.Part64
import SpinCodes.Majorant.Segment0.Part65
import SpinCodes.Majorant.Segment0.Part66
import SpinCodes.Majorant.Segment0.Part67
import SpinCodes.Majorant.Segment0.Part68
import SpinCodes.Majorant.Segment0.Part69
import SpinCodes.Majorant.Segment0.Part70
import SpinCodes.Majorant.Segment0.Part71
import SpinCodes.Majorant.Segment0.Part72
import SpinCodes.Majorant.Segment0.Part73
import SpinCodes.Majorant.Segment0.Part74
import SpinCodes.Majorant.Segment0.Part75
import SpinCodes.Majorant.Segment0.Part76
import SpinCodes.Majorant.Segment0.Part77
import SpinCodes.Majorant.Segment0.Part78
import SpinCodes.Majorant.Segment0.Part79
import SpinCodes.Majorant.Segment0.Part80
import SpinCodes.Majorant.Segment0.Part81
import SpinCodes.Majorant.Segment0.Part82
import SpinCodes.Majorant.Segment0.Part83
import SpinCodes.Majorant.Segment0.Part84
import SpinCodes.Majorant.Segment0.Part85
import SpinCodes.Majorant.Segment0.Part86
import SpinCodes.Majorant.Segment0.Part87
import SpinCodes.Majorant.Segment0.Part88
import SpinCodes.Majorant.Segment0.Part89
import SpinCodes.Numeric.MajorantEval
set_option maxRecDepth 8000000
set_option maxHeartbeats 0
namespace Spin.Majorant.S0
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
theorem part28 {s i : ℝ} (hs : Mem P28.S s) (hi : Mem P28.I i)
    (a b w : ℝ) (ha : Mem P28.A a) (hb : Mem P28.B b) (hw : Mem P28.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P28.logs_ok) hs hi P28.tree P28.tree_ok ha hb hw
theorem part29 {s i : ℝ} (hs : Mem P29.S s) (hi : Mem P29.I i)
    (a b w : ℝ) (ha : Mem P29.A a) (hb : Mem P29.B b) (hw : Mem P29.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P29.logs_ok) hs hi P29.tree P29.tree_ok ha hb hw
theorem part30 {s i : ℝ} (hs : Mem P30.S s) (hi : Mem P30.I i)
    (a b w : ℝ) (ha : Mem P30.A a) (hb : Mem P30.B b) (hw : Mem P30.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P30.logs_ok) hs hi P30.tree P30.tree_ok ha hb hw
theorem part31 {s i : ℝ} (hs : Mem P31.S s) (hi : Mem P31.I i)
    (a b w : ℝ) (ha : Mem P31.A a) (hb : Mem P31.B b) (hw : Mem P31.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P31.logs_ok) hs hi P31.tree P31.tree_ok ha hb hw
theorem part32 {s i : ℝ} (hs : Mem P32.S s) (hi : Mem P32.I i)
    (a b w : ℝ) (ha : Mem P32.A a) (hb : Mem P32.B b) (hw : Mem P32.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P32.logs_ok) hs hi P32.tree P32.tree_ok ha hb hw
theorem part33 {s i : ℝ} (hs : Mem P33.S s) (hi : Mem P33.I i)
    (a b w : ℝ) (ha : Mem P33.A a) (hb : Mem P33.B b) (hw : Mem P33.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P33.logs_ok) hs hi P33.tree P33.tree_ok ha hb hw
theorem part34 {s i : ℝ} (hs : Mem P34.S s) (hi : Mem P34.I i)
    (a b w : ℝ) (ha : Mem P34.A a) (hb : Mem P34.B b) (hw : Mem P34.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P34.logs_ok) hs hi P34.tree P34.tree_ok ha hb hw
theorem part35 {s i : ℝ} (hs : Mem P35.S s) (hi : Mem P35.I i)
    (a b w : ℝ) (ha : Mem P35.A a) (hb : Mem P35.B b) (hw : Mem P35.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P35.logs_ok) hs hi P35.tree P35.tree_ok ha hb hw
theorem part36 {s i : ℝ} (hs : Mem P36.S s) (hi : Mem P36.I i)
    (a b w : ℝ) (ha : Mem P36.A a) (hb : Mem P36.B b) (hw : Mem P36.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P36.logs_ok) hs hi P36.tree P36.tree_ok ha hb hw
theorem part37 {s i : ℝ} (hs : Mem P37.S s) (hi : Mem P37.I i)
    (a b w : ℝ) (ha : Mem P37.A a) (hb : Mem P37.B b) (hw : Mem P37.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P37.logs_ok) hs hi P37.tree P37.tree_ok ha hb hw
theorem part38 {s i : ℝ} (hs : Mem P38.S s) (hi : Mem P38.I i)
    (a b w : ℝ) (ha : Mem P38.A a) (hb : Mem P38.B b) (hw : Mem P38.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P38.logs_ok) hs hi P38.tree P38.tree_ok ha hb hw
theorem part39 {s i : ℝ} (hs : Mem P39.S s) (hi : Mem P39.I i)
    (a b w : ℝ) (ha : Mem P39.A a) (hb : Mem P39.B b) (hw : Mem P39.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P39.logs_ok) hs hi P39.tree P39.tree_ok ha hb hw
theorem part40 {s i : ℝ} (hs : Mem P40.S s) (hi : Mem P40.I i)
    (a b w : ℝ) (ha : Mem P40.A a) (hb : Mem P40.B b) (hw : Mem P40.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P40.logs_ok) hs hi P40.tree P40.tree_ok ha hb hw
theorem part41 {s i : ℝ} (hs : Mem P41.S s) (hi : Mem P41.I i)
    (a b w : ℝ) (ha : Mem P41.A a) (hb : Mem P41.B b) (hw : Mem P41.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P41.logs_ok) hs hi P41.tree P41.tree_ok ha hb hw
theorem part42 {s i : ℝ} (hs : Mem P42.S s) (hi : Mem P42.I i)
    (a b w : ℝ) (ha : Mem P42.A a) (hb : Mem P42.B b) (hw : Mem P42.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P42.logs_ok) hs hi P42.tree P42.tree_ok ha hb hw
theorem part43 {s i : ℝ} (hs : Mem P43.S s) (hi : Mem P43.I i)
    (a b w : ℝ) (ha : Mem P43.A a) (hb : Mem P43.B b) (hw : Mem P43.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P43.logs_ok) hs hi P43.tree P43.tree_ok ha hb hw
theorem part44 {s i : ℝ} (hs : Mem P44.S s) (hi : Mem P44.I i)
    (a b w : ℝ) (ha : Mem P44.A a) (hb : Mem P44.B b) (hw : Mem P44.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P44.logs_ok) hs hi P44.tree P44.tree_ok ha hb hw
theorem part45 {s i : ℝ} (hs : Mem P45.S s) (hi : Mem P45.I i)
    (a b w : ℝ) (ha : Mem P45.A a) (hb : Mem P45.B b) (hw : Mem P45.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P45.logs_ok) hs hi P45.tree P45.tree_ok ha hb hw
theorem part46 {s i : ℝ} (hs : Mem P46.S s) (hi : Mem P46.I i)
    (a b w : ℝ) (ha : Mem P46.A a) (hb : Mem P46.B b) (hw : Mem P46.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P46.logs_ok) hs hi P46.tree P46.tree_ok ha hb hw
theorem part47 {s i : ℝ} (hs : Mem P47.S s) (hi : Mem P47.I i)
    (a b w : ℝ) (ha : Mem P47.A a) (hb : Mem P47.B b) (hw : Mem P47.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P47.logs_ok) hs hi P47.tree P47.tree_ok ha hb hw
theorem part48 {s i : ℝ} (hs : Mem P48.S s) (hi : Mem P48.I i)
    (a b w : ℝ) (ha : Mem P48.A a) (hb : Mem P48.B b) (hw : Mem P48.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P48.logs_ok) hs hi P48.tree P48.tree_ok ha hb hw
theorem part49 {s i : ℝ} (hs : Mem P49.S s) (hi : Mem P49.I i)
    (a b w : ℝ) (ha : Mem P49.A a) (hb : Mem P49.B b) (hw : Mem P49.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P49.logs_ok) hs hi P49.tree P49.tree_ok ha hb hw
theorem part50 {s i : ℝ} (hs : Mem P50.S s) (hi : Mem P50.I i)
    (a b w : ℝ) (ha : Mem P50.A a) (hb : Mem P50.B b) (hw : Mem P50.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P50.logs_ok) hs hi P50.tree P50.tree_ok ha hb hw
theorem part51 {s i : ℝ} (hs : Mem P51.S s) (hi : Mem P51.I i)
    (a b w : ℝ) (ha : Mem P51.A a) (hb : Mem P51.B b) (hw : Mem P51.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P51.logs_ok) hs hi P51.tree P51.tree_ok ha hb hw
theorem part52 {s i : ℝ} (hs : Mem P52.S s) (hi : Mem P52.I i)
    (a b w : ℝ) (ha : Mem P52.A a) (hb : Mem P52.B b) (hw : Mem P52.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P52.logs_ok) hs hi P52.tree P52.tree_ok ha hb hw
theorem part53 {s i : ℝ} (hs : Mem P53.S s) (hi : Mem P53.I i)
    (a b w : ℝ) (ha : Mem P53.A a) (hb : Mem P53.B b) (hw : Mem P53.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P53.logs_ok) hs hi P53.tree P53.tree_ok ha hb hw
theorem part54 {s i : ℝ} (hs : Mem P54.S s) (hi : Mem P54.I i)
    (a b w : ℝ) (ha : Mem P54.A a) (hb : Mem P54.B b) (hw : Mem P54.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P54.logs_ok) hs hi P54.tree P54.tree_ok ha hb hw
theorem part55 {s i : ℝ} (hs : Mem P55.S s) (hi : Mem P55.I i)
    (a b w : ℝ) (ha : Mem P55.A a) (hb : Mem P55.B b) (hw : Mem P55.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P55.logs_ok) hs hi P55.tree P55.tree_ok ha hb hw
theorem part56 {s i : ℝ} (hs : Mem P56.S s) (hi : Mem P56.I i)
    (a b w : ℝ) (ha : Mem P56.A a) (hb : Mem P56.B b) (hw : Mem P56.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P56.logs_ok) hs hi P56.tree P56.tree_ok ha hb hw
theorem part57 {s i : ℝ} (hs : Mem P57.S s) (hi : Mem P57.I i)
    (a b w : ℝ) (ha : Mem P57.A a) (hb : Mem P57.B b) (hw : Mem P57.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P57.logs_ok) hs hi P57.tree P57.tree_ok ha hb hw
theorem part58 {s i : ℝ} (hs : Mem P58.S s) (hi : Mem P58.I i)
    (a b w : ℝ) (ha : Mem P58.A a) (hb : Mem P58.B b) (hw : Mem P58.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P58.logs_ok) hs hi P58.tree P58.tree_ok ha hb hw
theorem part59 {s i : ℝ} (hs : Mem P59.S s) (hi : Mem P59.I i)
    (a b w : ℝ) (ha : Mem P59.A a) (hb : Mem P59.B b) (hw : Mem P59.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P59.logs_ok) hs hi P59.tree P59.tree_ok ha hb hw
theorem part60 {s i : ℝ} (hs : Mem P60.S s) (hi : Mem P60.I i)
    (a b w : ℝ) (ha : Mem P60.A a) (hb : Mem P60.B b) (hw : Mem P60.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P60.logs_ok) hs hi P60.tree P60.tree_ok ha hb hw
theorem part61 {s i : ℝ} (hs : Mem P61.S s) (hi : Mem P61.I i)
    (a b w : ℝ) (ha : Mem P61.A a) (hb : Mem P61.B b) (hw : Mem P61.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P61.logs_ok) hs hi P61.tree P61.tree_ok ha hb hw
theorem part62 {s i : ℝ} (hs : Mem P62.S s) (hi : Mem P62.I i)
    (a b w : ℝ) (ha : Mem P62.A a) (hb : Mem P62.B b) (hw : Mem P62.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P62.logs_ok) hs hi P62.tree P62.tree_ok ha hb hw
theorem part63 {s i : ℝ} (hs : Mem P63.S s) (hi : Mem P63.I i)
    (a b w : ℝ) (ha : Mem P63.A a) (hb : Mem P63.B b) (hw : Mem P63.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P63.logs_ok) hs hi P63.tree P63.tree_ok ha hb hw
theorem part64 {s i : ℝ} (hs : Mem P64.S s) (hi : Mem P64.I i)
    (a b w : ℝ) (ha : Mem P64.A a) (hb : Mem P64.B b) (hw : Mem P64.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P64.logs_ok) hs hi P64.tree P64.tree_ok ha hb hw
theorem part65 {s i : ℝ} (hs : Mem P65.S s) (hi : Mem P65.I i)
    (a b w : ℝ) (ha : Mem P65.A a) (hb : Mem P65.B b) (hw : Mem P65.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P65.logs_ok) hs hi P65.tree P65.tree_ok ha hb hw
theorem part66 {s i : ℝ} (hs : Mem P66.S s) (hi : Mem P66.I i)
    (a b w : ℝ) (ha : Mem P66.A a) (hb : Mem P66.B b) (hw : Mem P66.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P66.logs_ok) hs hi P66.tree P66.tree_ok ha hb hw
theorem part67 {s i : ℝ} (hs : Mem P67.S s) (hi : Mem P67.I i)
    (a b w : ℝ) (ha : Mem P67.A a) (hb : Mem P67.B b) (hw : Mem P67.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P67.logs_ok) hs hi P67.tree P67.tree_ok ha hb hw
theorem part68 {s i : ℝ} (hs : Mem P68.S s) (hi : Mem P68.I i)
    (a b w : ℝ) (ha : Mem P68.A a) (hb : Mem P68.B b) (hw : Mem P68.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P68.logs_ok) hs hi P68.tree P68.tree_ok ha hb hw
theorem part69 {s i : ℝ} (hs : Mem P69.S s) (hi : Mem P69.I i)
    (a b w : ℝ) (ha : Mem P69.A a) (hb : Mem P69.B b) (hw : Mem P69.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P69.logs_ok) hs hi P69.tree P69.tree_ok ha hb hw
theorem part70 {s i : ℝ} (hs : Mem P70.S s) (hi : Mem P70.I i)
    (a b w : ℝ) (ha : Mem P70.A a) (hb : Mem P70.B b) (hw : Mem P70.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P70.logs_ok) hs hi P70.tree P70.tree_ok ha hb hw
theorem part71 {s i : ℝ} (hs : Mem P71.S s) (hi : Mem P71.I i)
    (a b w : ℝ) (ha : Mem P71.A a) (hb : Mem P71.B b) (hw : Mem P71.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P71.logs_ok) hs hi P71.tree P71.tree_ok ha hb hw
theorem part72 {s i : ℝ} (hs : Mem P72.S s) (hi : Mem P72.I i)
    (a b w : ℝ) (ha : Mem P72.A a) (hb : Mem P72.B b) (hw : Mem P72.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P72.logs_ok) hs hi P72.tree P72.tree_ok ha hb hw
theorem part73 {s i : ℝ} (hs : Mem P73.S s) (hi : Mem P73.I i)
    (a b w : ℝ) (ha : Mem P73.A a) (hb : Mem P73.B b) (hw : Mem P73.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P73.logs_ok) hs hi P73.tree P73.tree_ok ha hb hw
theorem part74 {s i : ℝ} (hs : Mem P74.S s) (hi : Mem P74.I i)
    (a b w : ℝ) (ha : Mem P74.A a) (hb : Mem P74.B b) (hw : Mem P74.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P74.logs_ok) hs hi P74.tree P74.tree_ok ha hb hw
theorem part75 {s i : ℝ} (hs : Mem P75.S s) (hi : Mem P75.I i)
    (a b w : ℝ) (ha : Mem P75.A a) (hb : Mem P75.B b) (hw : Mem P75.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P75.logs_ok) hs hi P75.tree P75.tree_ok ha hb hw
theorem part76 {s i : ℝ} (hs : Mem P76.S s) (hi : Mem P76.I i)
    (a b w : ℝ) (ha : Mem P76.A a) (hb : Mem P76.B b) (hw : Mem P76.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P76.logs_ok) hs hi P76.tree P76.tree_ok ha hb hw
theorem part77 {s i : ℝ} (hs : Mem P77.S s) (hi : Mem P77.I i)
    (a b w : ℝ) (ha : Mem P77.A a) (hb : Mem P77.B b) (hw : Mem P77.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P77.logs_ok) hs hi P77.tree P77.tree_ok ha hb hw
theorem part78 {s i : ℝ} (hs : Mem P78.S s) (hi : Mem P78.I i)
    (a b w : ℝ) (ha : Mem P78.A a) (hb : Mem P78.B b) (hw : Mem P78.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P78.logs_ok) hs hi P78.tree P78.tree_ok ha hb hw
theorem part79 {s i : ℝ} (hs : Mem P79.S s) (hi : Mem P79.I i)
    (a b w : ℝ) (ha : Mem P79.A a) (hb : Mem P79.B b) (hw : Mem P79.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P79.logs_ok) hs hi P79.tree P79.tree_ok ha hb hw
theorem part80 {s i : ℝ} (hs : Mem P80.S s) (hi : Mem P80.I i)
    (a b w : ℝ) (ha : Mem P80.A a) (hb : Mem P80.B b) (hw : Mem P80.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P80.logs_ok) hs hi P80.tree P80.tree_ok ha hb hw
theorem part81 {s i : ℝ} (hs : Mem P81.S s) (hi : Mem P81.I i)
    (a b w : ℝ) (ha : Mem P81.A a) (hb : Mem P81.B b) (hw : Mem P81.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P81.logs_ok) hs hi P81.tree P81.tree_ok ha hb hw
theorem part82 {s i : ℝ} (hs : Mem P82.S s) (hi : Mem P82.I i)
    (a b w : ℝ) (ha : Mem P82.A a) (hb : Mem P82.B b) (hw : Mem P82.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P82.logs_ok) hs hi P82.tree P82.tree_ok ha hb hw
theorem part83 {s i : ℝ} (hs : Mem P83.S s) (hi : Mem P83.I i)
    (a b w : ℝ) (ha : Mem P83.A a) (hb : Mem P83.B b) (hw : Mem P83.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P83.logs_ok) hs hi P83.tree P83.tree_ok ha hb hw
theorem part84 {s i : ℝ} (hs : Mem P84.S s) (hi : Mem P84.I i)
    (a b w : ℝ) (ha : Mem P84.A a) (hb : Mem P84.B b) (hw : Mem P84.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P84.logs_ok) hs hi P84.tree P84.tree_ok ha hb hw
theorem part85 {s i : ℝ} (hs : Mem P85.S s) (hi : Mem P85.I i)
    (a b w : ℝ) (ha : Mem P85.A a) (hb : Mem P85.B b) (hw : Mem P85.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P85.logs_ok) hs hi P85.tree P85.tree_ok ha hb hw
theorem part86 {s i : ℝ} (hs : Mem P86.S s) (hi : Mem P86.I i)
    (a b w : ℝ) (ha : Mem P86.A a) (hb : Mem P86.B b) (hw : Mem P86.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P86.logs_ok) hs hi P86.tree P86.tree_ok ha hb hw
theorem part87 {s i : ℝ} (hs : Mem P87.S s) (hi : Mem P87.I i)
    (a b w : ℝ) (ha : Mem P87.A a) (hb : Mem P87.B b) (hw : Mem P87.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P87.logs_ok) hs hi P87.tree P87.tree_ok ha hb hw
theorem part88 {s i : ℝ} (hs : Mem P88.S s) (hi : Mem P88.I i)
    (a b w : ℝ) (ha : Mem P88.A a) (hb : Mem P88.B b) (hw : Mem P88.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P88.logs_ok) hs hi P88.tree P88.tree_ok ha hb hw
theorem part89 {s i : ℝ} (hs : Mem P89.S s) (hi : Mem P89.I i)
    (a b w : ℝ) (ha : Mem P89.A a) (hb : Mem P89.B b) (hw : Mem P89.W w) :
    MajorantClaim s i a b w :=
  majorant_of_checkTree (treeLog_oracle P89.logs_ok) hs hi P89.tree P89.tree_ok ha hb hw
theorem cover {s i : ℝ} (hs : Mem ⟨1666000000000000000000000000000, 1666000000000000000000000000000⟩ s) (hi : Mem ⟨-173244000000000000000000000000, -173244000000000000000000000000⟩ i)
    (a b w : ℝ) (ha : Mem ⟨0, 1000000000000000000000000000000⟩ a) (hb : Mem ⟨0, 1000000000000000000000000000000⟩ b) (hw : Mem ⟨104000000000000000000000000000, 152466884341666666666666666667⟩ w) :
    MajorantClaim s i a b w :=
  (claim_split_fst (A := ⟨0, 1000000000000000000000000000000⟩) (B := ⟨0, 1000000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 500000000000000000000000000000) (claim_split_snd (A := ⟨0, 500000000000000000000000000000⟩) (B := ⟨0, 1000000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 500000000000000000000000000000) (claim_split_fst (A := ⟨0, 500000000000000000000000000000⟩) (B := ⟨0, 500000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 250000000000000000000000000000) (claim_split_snd (A := ⟨0, 250000000000000000000000000000⟩) (B := ⟨0, 500000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 250000000000000000000000000000) (claim_split_fst (A := ⟨0, 250000000000000000000000000000⟩) (B := ⟨0, 250000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 125000000000000000000000000000) (claim_split_snd (A := ⟨0, 125000000000000000000000000000⟩) (B := ⟨0, 250000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 125000000000000000000000000000) (claim_split_fst (A := ⟨0, 125000000000000000000000000000⟩) (B := ⟨0, 125000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 62500000000000000000000000000) (claim_split_snd (A := ⟨0, 62500000000000000000000000000⟩) (B := ⟨0, 125000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 62500000000000000000000000000) (claim_split_thd (A := ⟨0, 62500000000000000000000000000⟩) (B := ⟨0, 62500000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 128233442170833333333333333333) (claim_split_fst (A := ⟨0, 62500000000000000000000000000⟩) (B := ⟨0, 62500000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 128233442170833333333333333333⟩) (m := 31250000000000000000000000000) (claim_split_snd (A := ⟨0, 31250000000000000000000000000⟩) (B := ⟨0, 62500000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 128233442170833333333333333333⟩) (m := 31250000000000000000000000000) (claim_split_thd (A := ⟨0, 31250000000000000000000000000⟩) (B := ⟨0, 31250000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 128233442170833333333333333333⟩) (m := 116116721085416666666666666666) (claim_split_fst (A := ⟨0, 31250000000000000000000000000⟩) (B := ⟨0, 31250000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 116116721085416666666666666666⟩) (m := 15625000000000000000000000000) (claim_split_snd (A := ⟨0, 15625000000000000000000000000⟩) (B := ⟨0, 31250000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 116116721085416666666666666666⟩) (m := 15625000000000000000000000000) (claim_split_thd (A := ⟨0, 15625000000000000000000000000⟩) (B := ⟨0, 15625000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 116116721085416666666666666666⟩) (m := 110058360542708333333333333333) (claim_split_fst (A := ⟨0, 15625000000000000000000000000⟩) (B := ⟨0, 15625000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 110058360542708333333333333333⟩) (m := 7812500000000000000000000000) (claim_split_snd (A := ⟨0, 7812500000000000000000000000⟩) (B := ⟨0, 15625000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 110058360542708333333333333333⟩) (m := 7812500000000000000000000000) (claim_split_thd (A := ⟨0, 7812500000000000000000000000⟩) (B := ⟨0, 7812500000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 110058360542708333333333333333⟩) (m := 107029180271354166666666666666) (claim_split_fst (A := ⟨0, 7812500000000000000000000000⟩) (B := ⟨0, 7812500000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 107029180271354166666666666666⟩) (m := 3906250000000000000000000000) (claim_split_snd (A := ⟨0, 3906250000000000000000000000⟩) (B := ⟨0, 7812500000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 107029180271354166666666666666⟩) (m := 3906250000000000000000000000) (claim_split_thd (A := ⟨0, 3906250000000000000000000000⟩) (B := ⟨0, 3906250000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 107029180271354166666666666666⟩) (m := 105514590135677083333333333333) (claim_split_fst (A := ⟨0, 3906250000000000000000000000⟩) (B := ⟨0, 3906250000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 105514590135677083333333333333⟩) (m := 1953125000000000000000000000) (claim_split_snd (A := ⟨0, 1953125000000000000000000000⟩) (B := ⟨0, 3906250000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 105514590135677083333333333333⟩) (m := 1953125000000000000000000000) (claim_split_thd (A := ⟨0, 1953125000000000000000000000⟩) (B := ⟨0, 1953125000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 105514590135677083333333333333⟩) (m := 104757295067838541666666666666) (claim_split_fst (A := ⟨0, 1953125000000000000000000000⟩) (B := ⟨0, 1953125000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 104757295067838541666666666666⟩) (m := 976562500000000000000000000) (claim_split_snd (A := ⟨0, 976562500000000000000000000⟩) (B := ⟨0, 1953125000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 104757295067838541666666666666⟩) (m := 976562500000000000000000000) (claim_split_thd (A := ⟨0, 976562500000000000000000000⟩) (B := ⟨0, 976562500000000000000000000⟩) (W := ⟨104000000000000000000000000000, 104757295067838541666666666666⟩) (m := 104378647533919270833333333333) (claim_split_fst (A := ⟨0, 976562500000000000000000000⟩) (B := ⟨0, 976562500000000000000000000⟩) (W := ⟨104000000000000000000000000000, 104378647533919270833333333333⟩) (m := 488281250000000000000000000) (claim_split_snd (A := ⟨0, 488281250000000000000000000⟩) (B := ⟨0, 976562500000000000000000000⟩) (W := ⟨104000000000000000000000000000, 104378647533919270833333333333⟩) (m := 488281250000000000000000000) (claim_split_thd (A := ⟨0, 488281250000000000000000000⟩) (B := ⟨0, 488281250000000000000000000⟩) (W := ⟨104000000000000000000000000000, 104378647533919270833333333333⟩) (m := 104189323766959635416666666666) (claim_split_fst (A := ⟨0, 488281250000000000000000000⟩) (B := ⟨0, 488281250000000000000000000⟩) (W := ⟨104000000000000000000000000000, 104189323766959635416666666666⟩) (m := 244140625000000000000000000) (claim_split_snd (A := ⟨0, 244140625000000000000000000⟩) (B := ⟨0, 488281250000000000000000000⟩) (W := ⟨104000000000000000000000000000, 104189323766959635416666666666⟩) (m := 244140625000000000000000000) (claim_split_thd (A := ⟨0, 244140625000000000000000000⟩) (B := ⟨0, 244140625000000000000000000⟩) (W := ⟨104000000000000000000000000000, 104189323766959635416666666666⟩) (m := 104094661883479817708333333333) (claim_split_fst (A := ⟨0, 244140625000000000000000000⟩) (B := ⟨0, 244140625000000000000000000⟩) (W := ⟨104000000000000000000000000000, 104094661883479817708333333333⟩) (m := 122070312500000000000000000) (claim_split_snd (A := ⟨0, 122070312500000000000000000⟩) (B := ⟨0, 244140625000000000000000000⟩) (W := ⟨104000000000000000000000000000, 104094661883479817708333333333⟩) (m := 122070312500000000000000000) (claim_split_thd (A := ⟨0, 122070312500000000000000000⟩) (B := ⟨0, 122070312500000000000000000⟩) (W := ⟨104000000000000000000000000000, 104094661883479817708333333333⟩) (m := 104047330941739908854166666666) (claim_split_fst (A := ⟨0, 122070312500000000000000000⟩) (B := ⟨0, 122070312500000000000000000⟩) (W := ⟨104000000000000000000000000000, 104047330941739908854166666666⟩) (m := 61035156250000000000000000) (claim_split_snd (A := ⟨0, 61035156250000000000000000⟩) (B := ⟨0, 122070312500000000000000000⟩) (W := ⟨104000000000000000000000000000, 104047330941739908854166666666⟩) (m := 61035156250000000000000000) (claim_split_thd (A := ⟨0, 61035156250000000000000000⟩) (B := ⟨0, 61035156250000000000000000⟩) (W := ⟨104000000000000000000000000000, 104047330941739908854166666666⟩) (m := 104023665470869954427083333333) (claim_split_fst (A := ⟨0, 61035156250000000000000000⟩) (B := ⟨0, 61035156250000000000000000⟩) (W := ⟨104000000000000000000000000000, 104023665470869954427083333333⟩) (m := 30517578125000000000000000) (part0 hs hi) (part1 hs hi)) (part2 hs hi)) (part3 hs hi)) (part4 hs hi)) (part5 hs hi)) (part6 hs hi)) (part7 hs hi)) (part8 hs hi)) (part9 hs hi)) (part10 hs hi)) (part11 hs hi)) (part12 hs hi)) (part13 hs hi)) (part14 hs hi)) (part15 hs hi)) (part16 hs hi)) (part17 hs hi)) (part18 hs hi)) (part19 hs hi)) (part20 hs hi)) (part21 hs hi)) (part22 hs hi)) (part23 hs hi)) (part24 hs hi)) (part25 hs hi)) (part26 hs hi)) (part27 hs hi)) (part28 hs hi)) (part29 hs hi)) (part30 hs hi)) (part31 hs hi)) (part32 hs hi)) (part33 hs hi)) (part34 hs hi)) (part35 hs hi)) (claim_split_snd (A := ⟨125000000000000000000000000000, 250000000000000000000000000000⟩) (B := ⟨0, 250000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 125000000000000000000000000000) (claim_split_fst (A := ⟨125000000000000000000000000000, 250000000000000000000000000000⟩) (B := ⟨0, 125000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 187500000000000000000000000000) (part36 hs hi) (part37 hs hi)) (claim_split_fst (A := ⟨125000000000000000000000000000, 250000000000000000000000000000⟩) (B := ⟨125000000000000000000000000000, 250000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 187500000000000000000000000000) (part38 hs hi) (claim_split_snd (A := ⟨187500000000000000000000000000, 250000000000000000000000000000⟩) (B := ⟨125000000000000000000000000000, 250000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 187500000000000000000000000000) (part39 hs hi) (part40 hs hi))))) (part41 hs hi)) (claim_split_snd (A := ⟨250000000000000000000000000000, 500000000000000000000000000000⟩) (B := ⟨0, 500000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 250000000000000000000000000000) (claim_split_fst (A := ⟨250000000000000000000000000000, 500000000000000000000000000000⟩) (B := ⟨0, 250000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 375000000000000000000000000000) (claim_split_snd (A := ⟨250000000000000000000000000000, 375000000000000000000000000000⟩) (B := ⟨0, 250000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 125000000000000000000000000000) (part42 hs hi) (claim_split_fst (A := ⟨250000000000000000000000000000, 375000000000000000000000000000⟩) (B := ⟨125000000000000000000000000000, 250000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 312500000000000000000000000000) (claim_split_snd (A := ⟨250000000000000000000000000000, 312500000000000000000000000000⟩) (B := ⟨125000000000000000000000000000, 250000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 187500000000000000000000000000) (claim_split_thd (A := ⟨250000000000000000000000000000, 312500000000000000000000000000⟩) (B := ⟨125000000000000000000000000000, 187500000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 128233442170833333333333333333) (claim_split_fst (A := ⟨250000000000000000000000000000, 312500000000000000000000000000⟩) (B := ⟨125000000000000000000000000000, 187500000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 128233442170833333333333333333⟩) (m := 281250000000000000000000000000) (part43 hs hi) (part44 hs hi)) (part45 hs hi)) (claim_split_thd (A := ⟨250000000000000000000000000000, 312500000000000000000000000000⟩) (B := ⟨187500000000000000000000000000, 250000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 128233442170833333333333333333) (claim_split_fst (A := ⟨250000000000000000000000000000, 312500000000000000000000000000⟩) (B := ⟨187500000000000000000000000000, 250000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 128233442170833333333333333333⟩) (m := 281250000000000000000000000000) (part46 hs hi) (claim_split_snd (A := ⟨281250000000000000000000000000, 312500000000000000000000000000⟩) (B := ⟨187500000000000000000000000000, 250000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 128233442170833333333333333333⟩) (m := 218750000000000000000000000000) (part47 hs hi) (part48 hs hi))) (part49 hs hi))) (claim_split_snd (A := ⟨312500000000000000000000000000, 375000000000000000000000000000⟩) (B := ⟨125000000000000000000000000000, 250000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 187500000000000000000000000000) (part50 hs hi) (claim_split_thd (A := ⟨312500000000000000000000000000, 375000000000000000000000000000⟩) (B := ⟨187500000000000000000000000000, 250000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 128233442170833333333333333333) (claim_split_fst (A := ⟨312500000000000000000000000000, 375000000000000000000000000000⟩) (B := ⟨187500000000000000000000000000, 250000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 128233442170833333333333333333⟩) (m := 343750000000000000000000000000) (claim_split_snd (A := ⟨312500000000000000000000000000, 343750000000000000000000000000⟩) (B := ⟨187500000000000000000000000000, 250000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 128233442170833333333333333333⟩) (m := 218750000000000000000000000000) (part51 hs hi) (claim_split_thd (A := ⟨312500000000000000000000000000, 343750000000000000000000000000⟩) (B := ⟨218750000000000000000000000000, 250000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 128233442170833333333333333333⟩) (m := 116116721085416666666666666666) (part52 hs hi) (part53 hs hi))) (claim_split_snd (A := ⟨343750000000000000000000000000, 375000000000000000000000000000⟩) (B := ⟨187500000000000000000000000000, 250000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 128233442170833333333333333333⟩) (m := 218750000000000000000000000000) (part54 hs hi) (claim_split_thd (A := ⟨343750000000000000000000000000, 375000000000000000000000000000⟩) (B := ⟨218750000000000000000000000000, 250000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 128233442170833333333333333333⟩) (m := 116116721085416666666666666666) (part55 hs hi) (part56 hs hi)))) (claim_split_fst (A := ⟨312500000000000000000000000000, 375000000000000000000000000000⟩) (B := ⟨187500000000000000000000000000, 250000000000000000000000000000⟩) (W := ⟨128233442170833333333333333333, 152466884341666666666666666667⟩) (m := 343750000000000000000000000000) (part57 hs hi) (claim_split_snd (A := ⟨343750000000000000000000000000, 375000000000000000000000000000⟩) (B := ⟨187500000000000000000000000000, 250000000000000000000000000000⟩) (W := ⟨128233442170833333333333333333, 152466884341666666666666666667⟩) (m := 218750000000000000000000000000) (part58 hs hi) (claim_split_thd (A := ⟨343750000000000000000000000000, 375000000000000000000000000000⟩) (B := ⟨218750000000000000000000000000, 250000000000000000000000000000⟩) (W := ⟨128233442170833333333333333333, 152466884341666666666666666667⟩) (m := 140350163256250000000000000000) (part59 hs hi) (part60 hs hi)))))))) (claim_split_snd (A := ⟨375000000000000000000000000000, 500000000000000000000000000000⟩) (B := ⟨0, 250000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 125000000000000000000000000000) (part61 hs hi) (claim_split_fst (A := ⟨375000000000000000000000000000, 500000000000000000000000000000⟩) (B := ⟨125000000000000000000000000000, 250000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 437500000000000000000000000000) (claim_split_snd (A := ⟨375000000000000000000000000000, 437500000000000000000000000000⟩) (B := ⟨125000000000000000000000000000, 250000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 187500000000000000000000000000) (part62 hs hi) (claim_split_thd (A := ⟨375000000000000000000000000000, 437500000000000000000000000000⟩) (B := ⟨187500000000000000000000000000, 250000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 128233442170833333333333333333) (claim_split_fst (A := ⟨375000000000000000000000000000, 437500000000000000000000000000⟩) (B := ⟨187500000000000000000000000000, 250000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 128233442170833333333333333333⟩) (m := 406250000000000000000000000000) (claim_split_snd (A := ⟨375000000000000000000000000000, 406250000000000000000000000000⟩) (B := ⟨187500000000000000000000000000, 250000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 128233442170833333333333333333⟩) (m := 218750000000000000000000000000) (part63 hs hi) (part64 hs hi)) (part65 hs hi)) (part66 hs hi))) (part67 hs hi)))) (claim_split_fst (A := ⟨250000000000000000000000000000, 500000000000000000000000000000⟩) (B := ⟨250000000000000000000000000000, 500000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 375000000000000000000000000000) (claim_split_snd (A := ⟨250000000000000000000000000000, 375000000000000000000000000000⟩) (B := ⟨250000000000000000000000000000, 500000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 375000000000000000000000000000) (claim_split_fst (A := ⟨250000000000000000000000000000, 375000000000000000000000000000⟩) (B := ⟨250000000000000000000000000000, 375000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 312500000000000000000000000000) (part68 hs hi) (claim_split_snd (A := ⟨312500000000000000000000000000, 375000000000000000000000000000⟩) (B := ⟨250000000000000000000000000000, 375000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 312500000000000000000000000000) (claim_split_thd (A := ⟨312500000000000000000000000000, 375000000000000000000000000000⟩) (B := ⟨250000000000000000000000000000, 312500000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 128233442170833333333333333333) (part69 hs hi) (claim_split_fst (A := ⟨312500000000000000000000000000, 375000000000000000000000000000⟩) (B := ⟨250000000000000000000000000000, 312500000000000000000000000000⟩) (W := ⟨128233442170833333333333333333, 152466884341666666666666666667⟩) (m := 343750000000000000000000000000) (part70 hs hi) (claim_split_snd (A := ⟨343750000000000000000000000000, 375000000000000000000000000000⟩) (B := ⟨250000000000000000000000000000, 312500000000000000000000000000⟩) (W := ⟨128233442170833333333333333333, 152466884341666666666666666667⟩) (m := 281250000000000000000000000000) (claim_split_thd (A := ⟨343750000000000000000000000000, 375000000000000000000000000000⟩) (B := ⟨250000000000000000000000000000, 281250000000000000000000000000⟩) (W := ⟨128233442170833333333333333333, 152466884341666666666666666667⟩) (m := 140350163256250000000000000000) (part71 hs hi) (part72 hs hi)) (part73 hs hi)))) (part74 hs hi))) (part75 hs hi)) (claim_split_snd (A := ⟨375000000000000000000000000000, 500000000000000000000000000000⟩) (B := ⟨250000000000000000000000000000, 500000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 375000000000000000000000000000) (claim_split_fst (A := ⟨375000000000000000000000000000, 500000000000000000000000000000⟩) (B := ⟨250000000000000000000000000000, 375000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 437500000000000000000000000000) (claim_split_snd (A := ⟨375000000000000000000000000000, 437500000000000000000000000000⟩) (B := ⟨250000000000000000000000000000, 375000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 312500000000000000000000000000) (claim_split_thd (A := ⟨375000000000000000000000000000, 437500000000000000000000000000⟩) (B := ⟨250000000000000000000000000000, 312500000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 128233442170833333333333333333) (part76 hs hi) (claim_split_fst (A := ⟨375000000000000000000000000000, 437500000000000000000000000000⟩) (B := ⟨250000000000000000000000000000, 312500000000000000000000000000⟩) (W := ⟨128233442170833333333333333333, 152466884341666666666666666667⟩) (m := 406250000000000000000000000000) (claim_split_snd (A := ⟨375000000000000000000000000000, 406250000000000000000000000000⟩) (B := ⟨250000000000000000000000000000, 312500000000000000000000000000⟩) (W := ⟨128233442170833333333333333333, 152466884341666666666666666667⟩) (m := 281250000000000000000000000000) (claim_split_thd (A := ⟨375000000000000000000000000000, 406250000000000000000000000000⟩) (B := ⟨250000000000000000000000000000, 281250000000000000000000000000⟩) (W := ⟨128233442170833333333333333333, 152466884341666666666666666667⟩) (m := 140350163256250000000000000000) (part77 hs hi) (part78 hs hi)) (part79 hs hi)) (claim_split_snd (A := ⟨406250000000000000000000000000, 437500000000000000000000000000⟩) (B := ⟨250000000000000000000000000000, 312500000000000000000000000000⟩) (W := ⟨128233442170833333333333333333, 152466884341666666666666666667⟩) (m := 281250000000000000000000000000) (part80 hs hi) (part81 hs hi)))) (part82 hs hi)) (claim_split_snd (A := ⟨437500000000000000000000000000, 500000000000000000000000000000⟩) (B := ⟨250000000000000000000000000000, 375000000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 312500000000000000000000000000) (claim_split_thd (A := ⟨437500000000000000000000000000, 500000000000000000000000000000⟩) (B := ⟨250000000000000000000000000000, 312500000000000000000000000000⟩) (W := ⟨104000000000000000000000000000, 152466884341666666666666666667⟩) (m := 128233442170833333333333333333) (part83 hs hi) (claim_split_fst (A := ⟨437500000000000000000000000000, 500000000000000000000000000000⟩) (B := ⟨250000000000000000000000000000, 312500000000000000000000000000⟩) (W := ⟨128233442170833333333333333333, 152466884341666666666666666667⟩) (m := 468750000000000000000000000000) (part84 hs hi) (part85 hs hi))) (part86 hs hi))) (part87 hs hi))))) (part88 hs hi)) (part89 hs hi)) a b w ha hb hw
end Spin.Majorant.S0
