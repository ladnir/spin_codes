import SpinCodes.Structured.ConcreteOuterDefectNumeric
import SpinCodes.Majorant.RefinedData

set_option maxRecDepth 100000
set_option maxHeartbeats 0
noncomputable section
namespace Spin.Structured.ConcreteOuter
open Spin.Numeric Spin.Numeric.Fix

lemma defect_endpoint_0_lo :
    ((833:ℝ)/500)*((13:ℝ)/125)+((-43311:ℝ)/250000)-hEnt ((13:ℝ)/125)+Real.log 2/2 ≤ 1281/100000 := by
  have hc : defectOK (ofFrac (13) 125) (ofFrac (833) 500) (ofFrac (-43311) 250000) = true := by decide +kernel
  have hh := defectOK_sound (ofFrac_mem (p:=13) (q:=125) (by norm_num))
    (ofFrac_mem (p:=833) (q:=500) (by norm_num))
    (ofFrac_mem (p:=-43311) (q:=250000) (by norm_num)) hc
  norm_num at hh ⊢
  exact hh

lemma defect_endpoint_0_hi :
    ((833:ℝ)/500)*((3:ℝ)/20)+((-43311:ℝ)/250000)-hEnt ((3:ℝ)/20)+Real.log 2/2 ≤ 1281/100000 := by
  have hc : defectOK (ofFrac (3) 20) (ofFrac (833) 500) (ofFrac (-43311) 250000) = true := by decide +kernel
  have hh := defectOK_sound (ofFrac_mem (p:=3) (q:=20) (by norm_num))
    (ofFrac_mem (p:=833) (q:=500) (by norm_num))
    (ofFrac_mem (p:=-43311) (q:=250000) (by norm_num)) hc
  norm_num at hh ⊢
  exact hh

lemma defect_interval_0 {x : ℝ} (hx : x ∈ Set.Icc ((13:ℝ)/125) ((3:ℝ)/20)) :
    Spin.Majorant.refined.toFun x-hEnt x ≤ -Real.log 2/2+1281/100000 := by
  have hl := defect_endpoint_0_lo
  have hh := defect_endpoint_0_hi
  have hbound := entropy_defect_interval (d := -Real.log 2/2+1281/100000) (s:=((833:ℝ)/500)) (i:=((-43311:ℝ)/250000))
    (by norm_num) (by norm_num) hx (by linarith) (by linarith)
  have hm := Spin.Majorant.refined.toFun_le_support Spin.Majorant.member0 x
  norm_num at hm
  linarith

lemma defect_endpoint_1_lo :
    ((31:ℝ)/20)*((3:ℝ)/20)+((-77022601619127:ℝ)/500000000000000)-hEnt ((3:ℝ)/20)+Real.log 2/2 ≤ 1281/100000 := by
  have hc : defectOK (ofFrac (3) 20) (ofFrac (31) 20) (ofFrac (-77022601619127) 500000000000000) = true := by decide +kernel
  have hh := defectOK_sound (ofFrac_mem (p:=3) (q:=20) (by norm_num))
    (ofFrac_mem (p:=31) (q:=20) (by norm_num))
    (ofFrac_mem (p:=-77022601619127) (q:=500000000000000) (by norm_num)) hc
  norm_num at hh ⊢
  exact hh

lemma defect_endpoint_1_hi :
    ((31:ℝ)/20)*((1:ℝ)/5)+((-77022601619127:ℝ)/500000000000000)-hEnt ((1:ℝ)/5)+Real.log 2/2 ≤ 1281/100000 := by
  have hc : defectOK (ofFrac (1) 5) (ofFrac (31) 20) (ofFrac (-77022601619127) 500000000000000) = true := by decide +kernel
  have hh := defectOK_sound (ofFrac_mem (p:=1) (q:=5) (by norm_num))
    (ofFrac_mem (p:=31) (q:=20) (by norm_num))
    (ofFrac_mem (p:=-77022601619127) (q:=500000000000000) (by norm_num)) hc
  norm_num at hh ⊢
  exact hh

lemma defect_interval_1 {x : ℝ} (hx : x ∈ Set.Icc ((3:ℝ)/20) ((1:ℝ)/5)) :
    Spin.Majorant.refined.toFun x-hEnt x ≤ -Real.log 2/2+1281/100000 := by
  have hl := defect_endpoint_1_lo
  have hh := defect_endpoint_1_hi
  have hbound := entropy_defect_interval (d := -Real.log 2/2+1281/100000) (s:=((31:ℝ)/20)) (i:=((-77022601619127:ℝ)/500000000000000))
    (by norm_num) (by norm_num) hx (by linarith) (by linarith)
  have hm := Spin.Majorant.refined.toFun_le_support Spin.Majorant.member5 x
  norm_num at hm
  linarith

lemma defect_endpoint_2_lo :
    ((6:ℝ)/5)*((1:ℝ)/5)+((-832718284462267:ℝ)/10000000000000000)-hEnt ((1:ℝ)/5)+Real.log 2/2 ≤ 1281/100000 := by
  have hc : defectOK (ofFrac (1) 5) (ofFrac (6) 5) (ofFrac (-832718284462267) 10000000000000000) = true := by decide +kernel
  have hh := defectOK_sound (ofFrac_mem (p:=1) (q:=5) (by norm_num))
    (ofFrac_mem (p:=6) (q:=5) (by norm_num))
    (ofFrac_mem (p:=-832718284462267) (q:=10000000000000000) (by norm_num)) hc
  norm_num at hh ⊢
  exact hh

lemma defect_endpoint_2_hi :
    ((6:ℝ)/5)*((1:ℝ)/4)+((-832718284462267:ℝ)/10000000000000000)-hEnt ((1:ℝ)/4)+Real.log 2/2 ≤ 1281/100000 := by
  have hc : defectOK (ofFrac (1) 4) (ofFrac (6) 5) (ofFrac (-832718284462267) 10000000000000000) = true := by decide +kernel
  have hh := defectOK_sound (ofFrac_mem (p:=1) (q:=4) (by norm_num))
    (ofFrac_mem (p:=6) (q:=5) (by norm_num))
    (ofFrac_mem (p:=-832718284462267) (q:=10000000000000000) (by norm_num)) hc
  norm_num at hh ⊢
  exact hh

lemma defect_interval_2 {x : ℝ} (hx : x ∈ Set.Icc ((1:ℝ)/5) ((1:ℝ)/4)) :
    Spin.Majorant.refined.toFun x-hEnt x ≤ -Real.log 2/2+1281/100000 := by
  have hl := defect_endpoint_2_lo
  have hh := defect_endpoint_2_hi
  have hbound := entropy_defect_interval (d := -Real.log 2/2+1281/100000) (s:=((6:ℝ)/5)) (i:=((-832718284462267:ℝ)/10000000000000000))
    (by norm_num) (by norm_num) hx (by linarith) (by linarith)
  have hm := Spin.Majorant.refined.toFun_le_support Spin.Majorant.member7 x
  norm_num at hm
  linarith

lemma defect_endpoint_3_lo :
    ((1:ℝ)/1)*((1:ℝ)/4)+((-16647078810249:ℝ)/500000000000000)-hEnt ((1:ℝ)/4)+Real.log 2/2 ≤ 1281/100000 := by
  have hc : defectOK (ofFrac (1) 4) (ofFrac (1) 1) (ofFrac (-16647078810249) 500000000000000) = true := by decide +kernel
  have hh := defectOK_sound (ofFrac_mem (p:=1) (q:=4) (by norm_num))
    (ofFrac_mem (p:=1) (q:=1) (by norm_num))
    (ofFrac_mem (p:=-16647078810249) (q:=500000000000000) (by norm_num)) hc
  norm_num at hh ⊢
  exact hh

lemma defect_endpoint_3_hi :
    ((1:ℝ)/1)*((3:ℝ)/10)+((-16647078810249:ℝ)/500000000000000)-hEnt ((3:ℝ)/10)+Real.log 2/2 ≤ 1281/100000 := by
  have hc : defectOK (ofFrac (3) 10) (ofFrac (1) 1) (ofFrac (-16647078810249) 500000000000000) = true := by decide +kernel
  have hh := defectOK_sound (ofFrac_mem (p:=3) (q:=10) (by norm_num))
    (ofFrac_mem (p:=1) (q:=1) (by norm_num))
    (ofFrac_mem (p:=-16647078810249) (q:=500000000000000) (by norm_num)) hc
  norm_num at hh ⊢
  exact hh

lemma defect_interval_3 {x : ℝ} (hx : x ∈ Set.Icc ((1:ℝ)/4) ((3:ℝ)/10)) :
    Spin.Majorant.refined.toFun x-hEnt x ≤ -Real.log 2/2+1281/100000 := by
  have hl := defect_endpoint_3_lo
  have hh := defect_endpoint_3_hi
  have hbound := entropy_defect_interval (d := -Real.log 2/2+1281/100000) (s:=((1:ℝ)/1)) (i:=((-16647078810249:ℝ)/500000000000000))
    (by norm_num) (by norm_num) hx (by linarith) (by linarith)
  have hm := Spin.Majorant.refined.toFun_le_support Spin.Majorant.member8 x
  norm_num at hm
  linarith

lemma defect_endpoint_4_lo :
    ((4:ℝ)/5)*((3:ℝ)/10)+((30683842683431:ℝ)/1250000000000000)-hEnt ((3:ℝ)/10)+Real.log 2/2 ≤ 1281/100000 := by
  have hc : defectOK (ofFrac (3) 10) (ofFrac (4) 5) (ofFrac (30683842683431) 1250000000000000) = true := by decide +kernel
  have hh := defectOK_sound (ofFrac_mem (p:=3) (q:=10) (by norm_num))
    (ofFrac_mem (p:=4) (q:=5) (by norm_num))
    (ofFrac_mem (p:=30683842683431) (q:=1250000000000000) (by norm_num)) hc
  norm_num at hh ⊢
  exact hh

lemma defect_endpoint_4_hi :
    ((4:ℝ)/5)*((7:ℝ)/20)+((30683842683431:ℝ)/1250000000000000)-hEnt ((7:ℝ)/20)+Real.log 2/2 ≤ 1281/100000 := by
  have hc : defectOK (ofFrac (7) 20) (ofFrac (4) 5) (ofFrac (30683842683431) 1250000000000000) = true := by decide +kernel
  have hh := defectOK_sound (ofFrac_mem (p:=7) (q:=20) (by norm_num))
    (ofFrac_mem (p:=4) (q:=5) (by norm_num))
    (ofFrac_mem (p:=30683842683431) (q:=1250000000000000) (by norm_num)) hc
  norm_num at hh ⊢
  exact hh

lemma defect_interval_4 {x : ℝ} (hx : x ∈ Set.Icc ((3:ℝ)/10) ((7:ℝ)/20)) :
    Spin.Majorant.refined.toFun x-hEnt x ≤ -Real.log 2/2+1281/100000 := by
  have hl := defect_endpoint_4_lo
  have hh := defect_endpoint_4_hi
  have hbound := entropy_defect_interval (d := -Real.log 2/2+1281/100000) (s:=((4:ℝ)/5)) (i:=((30683842683431:ℝ)/1250000000000000))
    (by norm_num) (by norm_num) hx (by linarith) (by linarith)
  have hm := Spin.Majorant.refined.toFun_le_support Spin.Majorant.member9 x
  norm_num at hm
  linarith

lemma defect_endpoint_5_lo :
    ((3:ℝ)/5)*((7:ℝ)/20)+((909341020092481:ℝ)/10000000000000000)-hEnt ((7:ℝ)/20)+Real.log 2/2 ≤ 1281/100000 := by
  have hc : defectOK (ofFrac (7) 20) (ofFrac (3) 5) (ofFrac (909341020092481) 10000000000000000) = true := by decide +kernel
  have hh := defectOK_sound (ofFrac_mem (p:=7) (q:=20) (by norm_num))
    (ofFrac_mem (p:=3) (q:=5) (by norm_num))
    (ofFrac_mem (p:=909341020092481) (q:=10000000000000000) (by norm_num)) hc
  norm_num at hh ⊢
  exact hh

lemma defect_endpoint_5_hi :
    ((3:ℝ)/5)*((2:ℝ)/5)+((909341020092481:ℝ)/10000000000000000)-hEnt ((2:ℝ)/5)+Real.log 2/2 ≤ 1281/100000 := by
  have hc : defectOK (ofFrac (2) 5) (ofFrac (3) 5) (ofFrac (909341020092481) 10000000000000000) = true := by decide +kernel
  have hh := defectOK_sound (ofFrac_mem (p:=2) (q:=5) (by norm_num))
    (ofFrac_mem (p:=3) (q:=5) (by norm_num))
    (ofFrac_mem (p:=909341020092481) (q:=10000000000000000) (by norm_num)) hc
  norm_num at hh ⊢
  exact hh

lemma defect_interval_5 {x : ℝ} (hx : x ∈ Set.Icc ((7:ℝ)/20) ((2:ℝ)/5)) :
    Spin.Majorant.refined.toFun x-hEnt x ≤ -Real.log 2/2+1281/100000 := by
  have hl := defect_endpoint_5_lo
  have hh := defect_endpoint_5_hi
  have hbound := entropy_defect_interval (d := -Real.log 2/2+1281/100000) (s:=((3:ℝ)/5)) (i:=((909341020092481:ℝ)/10000000000000000))
    (by norm_num) (by norm_num) hx (by linarith) (by linarith)
  have hm := Spin.Majorant.refined.toFun_le_support Spin.Majorant.member10 x
  norm_num at hm
  linarith

lemma defect_endpoint_6_lo :
    ((2:ℝ)/5)*((2:ℝ)/5)+((4161516955371:ℝ)/25000000000000)-hEnt ((2:ℝ)/5)+Real.log 2/2 ≤ 1281/100000 := by
  have hc : defectOK (ofFrac (2) 5) (ofFrac (2) 5) (ofFrac (4161516955371) 25000000000000) = true := by decide +kernel
  have hh := defectOK_sound (ofFrac_mem (p:=2) (q:=5) (by norm_num))
    (ofFrac_mem (p:=2) (q:=5) (by norm_num))
    (ofFrac_mem (p:=4161516955371) (q:=25000000000000) (by norm_num)) hc
  norm_num at hh ⊢
  exact hh

lemma defect_endpoint_6_hi :
    ((2:ℝ)/5)*((9:ℝ)/20)+((4161516955371:ℝ)/25000000000000)-hEnt ((9:ℝ)/20)+Real.log 2/2 ≤ 1281/100000 := by
  have hc : defectOK (ofFrac (9) 20) (ofFrac (2) 5) (ofFrac (4161516955371) 25000000000000) = true := by decide +kernel
  have hh := defectOK_sound (ofFrac_mem (p:=9) (q:=20) (by norm_num))
    (ofFrac_mem (p:=2) (q:=5) (by norm_num))
    (ofFrac_mem (p:=4161516955371) (q:=25000000000000) (by norm_num)) hc
  norm_num at hh ⊢
  exact hh

lemma defect_interval_6 {x : ℝ} (hx : x ∈ Set.Icc ((2:ℝ)/5) ((9:ℝ)/20)) :
    Spin.Majorant.refined.toFun x-hEnt x ≤ -Real.log 2/2+1281/100000 := by
  have hl := defect_endpoint_6_lo
  have hh := defect_endpoint_6_hi
  have hbound := entropy_defect_interval (d := -Real.log 2/2+1281/100000) (s:=((2:ℝ)/5)) (i:=((4161516955371:ℝ)/25000000000000))
    (by norm_num) (by norm_num) hx (by linarith) (by linarith)
  have hm := Spin.Majorant.refined.toFun_le_support Spin.Majorant.member11 x
  norm_num at hm
  linarith

lemma defect_endpoint_7_lo :
    ((1:ℝ)/10)*((9:ℝ)/20)+((297841147503783:ℝ)/1000000000000000)-hEnt ((9:ℝ)/20)+Real.log 2/2 ≤ 1281/100000 := by
  have hc : defectOK (ofFrac (9) 20) (ofFrac (1) 10) (ofFrac (297841147503783) 1000000000000000) = true := by decide +kernel
  have hh := defectOK_sound (ofFrac_mem (p:=9) (q:=20) (by norm_num))
    (ofFrac_mem (p:=1) (q:=10) (by norm_num))
    (ofFrac_mem (p:=297841147503783) (q:=1000000000000000) (by norm_num)) hc
  norm_num at hh ⊢
  exact hh

lemma defect_endpoint_7_hi :
    ((1:ℝ)/10)*((1:ℝ)/2)+((297841147503783:ℝ)/1000000000000000)-hEnt ((1:ℝ)/2)+Real.log 2/2 ≤ 1281/100000 := by
  have hc : defectOK (ofFrac (1) 2) (ofFrac (1) 10) (ofFrac (297841147503783) 1000000000000000) = true := by decide +kernel
  have hh := defectOK_sound (ofFrac_mem (p:=1) (q:=2) (by norm_num))
    (ofFrac_mem (p:=1) (q:=10) (by norm_num))
    (ofFrac_mem (p:=297841147503783) (q:=1000000000000000) (by norm_num)) hc
  norm_num at hh ⊢
  exact hh

lemma defect_interval_7 {x : ℝ} (hx : x ∈ Set.Icc ((9:ℝ)/20) ((1:ℝ)/2)) :
    Spin.Majorant.refined.toFun x-hEnt x ≤ -Real.log 2/2+1281/100000 := by
  have hl := defect_endpoint_7_lo
  have hh := defect_endpoint_7_hi
  have hbound := entropy_defect_interval (d := -Real.log 2/2+1281/100000) (s:=((1:ℝ)/10)) (i:=((297841147503783:ℝ)/1000000000000000))
    (by norm_num) (by norm_num) hx (by linarith) (by linarith)
  have hm := Spin.Majorant.refined.toFun_le_support Spin.Majorant.member13 x
  norm_num at hm
  linarith

theorem refined_defect_left {x : ℝ} (hx : x ∈ Set.Icc ((13:ℝ)/125) (1/2)) :
    Spin.Majorant.refined.toFun x-hEnt x ≤ -Real.log 2/2+1281/100000 := by
  by_cases h0 : x ≤ ((3:ℝ)/20)
  · exact defect_interval_0 ⟨by linarith [hx.1], h0⟩
  by_cases h1 : x ≤ ((1:ℝ)/5)
  · exact defect_interval_1 ⟨by linarith [hx.1], h1⟩
  by_cases h2 : x ≤ ((1:ℝ)/4)
  · exact defect_interval_2 ⟨by linarith [hx.1], h2⟩
  by_cases h3 : x ≤ ((3:ℝ)/10)
  · exact defect_interval_3 ⟨by linarith [hx.1], h3⟩
  by_cases h4 : x ≤ ((7:ℝ)/20)
  · exact defect_interval_4 ⟨by linarith [hx.1], h4⟩
  by_cases h5 : x ≤ ((2:ℝ)/5)
  · exact defect_interval_5 ⟨by linarith [hx.1], h5⟩
  by_cases h6 : x ≤ ((9:ℝ)/20)
  · exact defect_interval_6 ⟨by linarith [hx.1], h6⟩
  exact defect_interval_7 ⟨by linarith, hx.2⟩

theorem refined_defect {x : ℝ} (hx : x ∈ Set.Icc ((13:ℝ)/125) (112/125)) :
    Spin.Majorant.refined.toFun x-hEnt x ≤ -Real.log 2/2+1281/100000 := by
  by_cases h : x ≤ (1:ℝ)/2
  · exact refined_defect_left ⟨hx.1,h⟩
  · have hh := refined_defect_left (x:=1-x) ⟨by linarith [hx.2], by linarith⟩
    rw [Spin.Majorant.refined_reflect] at hh
    have he : hEnt (1-x) = hEnt x := by rw [hEnt_eq_binEntropy, Real.binEntropy_one_sub, ←hEnt_eq_binEntropy]
    rwa [he] at hh

end Spin.Structured.ConcreteOuter

