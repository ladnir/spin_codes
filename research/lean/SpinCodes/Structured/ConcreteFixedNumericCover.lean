import SpinCodes.Structured.ConcreteFixedNumericSegments

noncomputable section
namespace Spin.Structured.ConcreteFixedNumeric
open Spin.Numeric DenseOccupationFixed Set
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def fugacities : Finset ℚ := {(96838234476759233/100000000000000000:ℚ), (96838234476759233/100000000000000000:ℚ), (96838234476759233/100000000000000000:ℚ), (96838234476759233/100000000000000000:ℚ), (96838234476759233/100000000000000000:ℚ), (96838234476759233/100000000000000000:ℚ), (96838234476759233/100000000000000000:ℚ), (96838234476759233/100000000000000000:ℚ), (96838234476759233/100000000000000000:ℚ), (96838234476759233/100000000000000000:ℚ), (96838234476759233/100000000000000000:ℚ), (96838234476759233/100000000000000000:ℚ), (96838234476759233/100000000000000000:ℚ), (96838234476759233/100000000000000000:ℚ), (96838234476759233/100000000000000000:ℚ), (49900171601753829/50000000000000000:ℚ), (10499206513492543/10000000000000000:ℚ), (11167521538448069/10000000000000000:ℚ), (12492915353143481/10000000000000000:ℚ), (1486892840550039/1000000000000000:ℚ), (907492662852861/500000000000000:ℚ), (4430894293586297/2000000000000000:ℚ), (6761285149803551/2500000000000000:ℚ), (8254819521857547/2500000000000000:ℚ), (7974382895176137/2000000000000000:ℚ), (918788069754107/200000000000000:ℚ), (49510688747396809/10000000000000000:ℚ), (51330699179736179/10000000000000000:ℚ), (52546977145860421/10000000000000000:ℚ), (10848858740933581/2000000000000000:ℚ), (16961761337730057/2500000000000000:ℚ)}
theorem rate_uniform {x : ℝ} (hx : x∈Icc (13/125) (112/125)) :
    ∃ u : ℚ, u∈fugacities ∧ 0<(u:ℝ) ∧ rate x u ≤ -(8679/10000000) := by
  have hl00 : Segment00.lo.real≤x := by norm_num [Segment00.lo,QInput.real]; exact hx.1
  by_cases hh00 : x≤Segment00.hi.real
  · refine ⟨(96838234476759233/100000000000000000:ℚ), by decide +kernel, by norm_num, ?_⟩
    have h := Segment00.interval_bound ⟨hl00,hh00⟩
    norm_num [Segment00.u,QInput.real,Rat.cast_div,Rat.cast_natCast,Rat.cast_intCast] at h ⊢; exact h
  have hl01 : Segment01.lo.real≤x := by
    have h := (not_le.mp hh00).le
    norm_num [Segment00.hi,Segment01.lo,QInput.real] at h ⊢; exact h
  by_cases hh01 : x≤Segment01.hi.real
  · refine ⟨(96838234476759233/100000000000000000:ℚ), by decide +kernel, by norm_num, ?_⟩
    have h := Segment01.interval_bound ⟨hl01,hh01⟩
    norm_num [Segment01.u,QInput.real,Rat.cast_div,Rat.cast_natCast,Rat.cast_intCast] at h ⊢; exact h
  have hl02 : Segment02.lo.real≤x := by
    have h := (not_le.mp hh01).le
    norm_num [Segment01.hi,Segment02.lo,QInput.real] at h ⊢; exact h
  by_cases hh02 : x≤Segment02.hi.real
  · refine ⟨(96838234476759233/100000000000000000:ℚ), by decide +kernel, by norm_num, ?_⟩
    have h := Segment02.interval_bound ⟨hl02,hh02⟩
    norm_num [Segment02.u,QInput.real,Rat.cast_div,Rat.cast_natCast,Rat.cast_intCast] at h ⊢; exact h
  have hl03 : Segment03.lo.real≤x := by
    have h := (not_le.mp hh02).le
    norm_num [Segment02.hi,Segment03.lo,QInput.real] at h ⊢; exact h
  by_cases hh03 : x≤Segment03.hi.real
  · refine ⟨(96838234476759233/100000000000000000:ℚ), by decide +kernel, by norm_num, ?_⟩
    have h := Segment03.interval_bound ⟨hl03,hh03⟩
    norm_num [Segment03.u,QInput.real,Rat.cast_div,Rat.cast_natCast,Rat.cast_intCast] at h ⊢; exact h
  have hl04 : Segment04.lo.real≤x := by
    have h := (not_le.mp hh03).le
    norm_num [Segment03.hi,Segment04.lo,QInput.real] at h ⊢; exact h
  by_cases hh04 : x≤Segment04.hi.real
  · refine ⟨(96838234476759233/100000000000000000:ℚ), by decide +kernel, by norm_num, ?_⟩
    have h := Segment04.interval_bound ⟨hl04,hh04⟩
    norm_num [Segment04.u,QInput.real,Rat.cast_div,Rat.cast_natCast,Rat.cast_intCast] at h ⊢; exact h
  have hl05 : Segment05.lo.real≤x := by
    have h := (not_le.mp hh04).le
    norm_num [Segment04.hi,Segment05.lo,QInput.real] at h ⊢; exact h
  by_cases hh05 : x≤Segment05.hi.real
  · refine ⟨(96838234476759233/100000000000000000:ℚ), by decide +kernel, by norm_num, ?_⟩
    have h := Segment05.interval_bound ⟨hl05,hh05⟩
    norm_num [Segment05.u,QInput.real,Rat.cast_div,Rat.cast_natCast,Rat.cast_intCast] at h ⊢; exact h
  have hl06 : Segment06.lo.real≤x := by
    have h := (not_le.mp hh05).le
    norm_num [Segment05.hi,Segment06.lo,QInput.real] at h ⊢; exact h
  by_cases hh06 : x≤Segment06.hi.real
  · refine ⟨(96838234476759233/100000000000000000:ℚ), by decide +kernel, by norm_num, ?_⟩
    have h := Segment06.interval_bound ⟨hl06,hh06⟩
    norm_num [Segment06.u,QInput.real,Rat.cast_div,Rat.cast_natCast,Rat.cast_intCast] at h ⊢; exact h
  have hl07 : Segment07.lo.real≤x := by
    have h := (not_le.mp hh06).le
    norm_num [Segment06.hi,Segment07.lo,QInput.real] at h ⊢; exact h
  by_cases hh07 : x≤Segment07.hi.real
  · refine ⟨(96838234476759233/100000000000000000:ℚ), by decide +kernel, by norm_num, ?_⟩
    have h := Segment07.interval_bound ⟨hl07,hh07⟩
    norm_num [Segment07.u,QInput.real,Rat.cast_div,Rat.cast_natCast,Rat.cast_intCast] at h ⊢; exact h
  have hl08 : Segment08.lo.real≤x := by
    have h := (not_le.mp hh07).le
    norm_num [Segment07.hi,Segment08.lo,QInput.real] at h ⊢; exact h
  by_cases hh08 : x≤Segment08.hi.real
  · refine ⟨(96838234476759233/100000000000000000:ℚ), by decide +kernel, by norm_num, ?_⟩
    have h := Segment08.interval_bound ⟨hl08,hh08⟩
    norm_num [Segment08.u,QInput.real,Rat.cast_div,Rat.cast_natCast,Rat.cast_intCast] at h ⊢; exact h
  have hl09 : Segment09.lo.real≤x := by
    have h := (not_le.mp hh08).le
    norm_num [Segment08.hi,Segment09.lo,QInput.real] at h ⊢; exact h
  by_cases hh09 : x≤Segment09.hi.real
  · refine ⟨(96838234476759233/100000000000000000:ℚ), by decide +kernel, by norm_num, ?_⟩
    have h := Segment09.interval_bound ⟨hl09,hh09⟩
    norm_num [Segment09.u,QInput.real,Rat.cast_div,Rat.cast_natCast,Rat.cast_intCast] at h ⊢; exact h
  have hl10 : Segment10.lo.real≤x := by
    have h := (not_le.mp hh09).le
    norm_num [Segment09.hi,Segment10.lo,QInput.real] at h ⊢; exact h
  by_cases hh10 : x≤Segment10.hi.real
  · refine ⟨(96838234476759233/100000000000000000:ℚ), by decide +kernel, by norm_num, ?_⟩
    have h := Segment10.interval_bound ⟨hl10,hh10⟩
    norm_num [Segment10.u,QInput.real,Rat.cast_div,Rat.cast_natCast,Rat.cast_intCast] at h ⊢; exact h
  have hl11 : Segment11.lo.real≤x := by
    have h := (not_le.mp hh10).le
    norm_num [Segment10.hi,Segment11.lo,QInput.real] at h ⊢; exact h
  by_cases hh11 : x≤Segment11.hi.real
  · refine ⟨(96838234476759233/100000000000000000:ℚ), by decide +kernel, by norm_num, ?_⟩
    have h := Segment11.interval_bound ⟨hl11,hh11⟩
    norm_num [Segment11.u,QInput.real,Rat.cast_div,Rat.cast_natCast,Rat.cast_intCast] at h ⊢; exact h
  have hl12 : Segment12.lo.real≤x := by
    have h := (not_le.mp hh11).le
    norm_num [Segment11.hi,Segment12.lo,QInput.real] at h ⊢; exact h
  by_cases hh12 : x≤Segment12.hi.real
  · refine ⟨(96838234476759233/100000000000000000:ℚ), by decide +kernel, by norm_num, ?_⟩
    have h := Segment12.interval_bound ⟨hl12,hh12⟩
    norm_num [Segment12.u,QInput.real,Rat.cast_div,Rat.cast_natCast,Rat.cast_intCast] at h ⊢; exact h
  have hl13 : Segment13.lo.real≤x := by
    have h := (not_le.mp hh12).le
    norm_num [Segment12.hi,Segment13.lo,QInput.real] at h ⊢; exact h
  by_cases hh13 : x≤Segment13.hi.real
  · refine ⟨(96838234476759233/100000000000000000:ℚ), by decide +kernel, by norm_num, ?_⟩
    have h := Segment13.interval_bound ⟨hl13,hh13⟩
    norm_num [Segment13.u,QInput.real,Rat.cast_div,Rat.cast_natCast,Rat.cast_intCast] at h ⊢; exact h
  have hl14 : Segment14.lo.real≤x := by
    have h := (not_le.mp hh13).le
    norm_num [Segment13.hi,Segment14.lo,QInput.real] at h ⊢; exact h
  by_cases hh14 : x≤Segment14.hi.real
  · refine ⟨(96838234476759233/100000000000000000:ℚ), by decide +kernel, by norm_num, ?_⟩
    have h := Segment14.interval_bound ⟨hl14,hh14⟩
    norm_num [Segment14.u,QInput.real,Rat.cast_div,Rat.cast_natCast,Rat.cast_intCast] at h ⊢; exact h
  have hl15 : Segment15.lo.real≤x := by
    have h := (not_le.mp hh14).le
    norm_num [Segment14.hi,Segment15.lo,QInput.real] at h ⊢; exact h
  by_cases hh15 : x≤Segment15.hi.real
  · refine ⟨(49900171601753829/50000000000000000:ℚ), by decide +kernel, by norm_num, ?_⟩
    have h := Segment15.interval_bound ⟨hl15,hh15⟩
    norm_num [Segment15.u,QInput.real,Rat.cast_div,Rat.cast_natCast,Rat.cast_intCast] at h ⊢; exact h
  have hl16 : Segment16.lo.real≤x := by
    have h := (not_le.mp hh15).le
    norm_num [Segment15.hi,Segment16.lo,QInput.real] at h ⊢; exact h
  by_cases hh16 : x≤Segment16.hi.real
  · refine ⟨(10499206513492543/10000000000000000:ℚ), by decide +kernel, by norm_num, ?_⟩
    have h := Segment16.interval_bound ⟨hl16,hh16⟩
    norm_num [Segment16.u,QInput.real,Rat.cast_div,Rat.cast_natCast,Rat.cast_intCast] at h ⊢; exact h
  have hl17 : Segment17.lo.real≤x := by
    have h := (not_le.mp hh16).le
    norm_num [Segment16.hi,Segment17.lo,QInput.real] at h ⊢; exact h
  by_cases hh17 : x≤Segment17.hi.real
  · refine ⟨(11167521538448069/10000000000000000:ℚ), by decide +kernel, by norm_num, ?_⟩
    have h := Segment17.interval_bound ⟨hl17,hh17⟩
    norm_num [Segment17.u,QInput.real,Rat.cast_div,Rat.cast_natCast,Rat.cast_intCast] at h ⊢; exact h
  have hl18 : Segment18.lo.real≤x := by
    have h := (not_le.mp hh17).le
    norm_num [Segment17.hi,Segment18.lo,QInput.real] at h ⊢; exact h
  by_cases hh18 : x≤Segment18.hi.real
  · refine ⟨(12492915353143481/10000000000000000:ℚ), by decide +kernel, by norm_num, ?_⟩
    have h := Segment18.interval_bound ⟨hl18,hh18⟩
    norm_num [Segment18.u,QInput.real,Rat.cast_div,Rat.cast_natCast,Rat.cast_intCast] at h ⊢; exact h
  have hl19 : Segment19.lo.real≤x := by
    have h := (not_le.mp hh18).le
    norm_num [Segment18.hi,Segment19.lo,QInput.real] at h ⊢; exact h
  by_cases hh19 : x≤Segment19.hi.real
  · refine ⟨(1486892840550039/1000000000000000:ℚ), by decide +kernel, by norm_num, ?_⟩
    have h := Segment19.interval_bound ⟨hl19,hh19⟩
    norm_num [Segment19.u,QInput.real,Rat.cast_div,Rat.cast_natCast,Rat.cast_intCast] at h ⊢; exact h
  have hl20 : Segment20.lo.real≤x := by
    have h := (not_le.mp hh19).le
    norm_num [Segment19.hi,Segment20.lo,QInput.real] at h ⊢; exact h
  by_cases hh20 : x≤Segment20.hi.real
  · refine ⟨(907492662852861/500000000000000:ℚ), by decide +kernel, by norm_num, ?_⟩
    have h := Segment20.interval_bound ⟨hl20,hh20⟩
    norm_num [Segment20.u,QInput.real,Rat.cast_div,Rat.cast_natCast,Rat.cast_intCast] at h ⊢; exact h
  have hl21 : Segment21.lo.real≤x := by
    have h := (not_le.mp hh20).le
    norm_num [Segment20.hi,Segment21.lo,QInput.real] at h ⊢; exact h
  by_cases hh21 : x≤Segment21.hi.real
  · refine ⟨(4430894293586297/2000000000000000:ℚ), by decide +kernel, by norm_num, ?_⟩
    have h := Segment21.interval_bound ⟨hl21,hh21⟩
    norm_num [Segment21.u,QInput.real,Rat.cast_div,Rat.cast_natCast,Rat.cast_intCast] at h ⊢; exact h
  have hl22 : Segment22.lo.real≤x := by
    have h := (not_le.mp hh21).le
    norm_num [Segment21.hi,Segment22.lo,QInput.real] at h ⊢; exact h
  by_cases hh22 : x≤Segment22.hi.real
  · refine ⟨(6761285149803551/2500000000000000:ℚ), by decide +kernel, by norm_num, ?_⟩
    have h := Segment22.interval_bound ⟨hl22,hh22⟩
    norm_num [Segment22.u,QInput.real,Rat.cast_div,Rat.cast_natCast,Rat.cast_intCast] at h ⊢; exact h
  have hl23 : Segment23.lo.real≤x := by
    have h := (not_le.mp hh22).le
    norm_num [Segment22.hi,Segment23.lo,QInput.real] at h ⊢; exact h
  by_cases hh23 : x≤Segment23.hi.real
  · refine ⟨(8254819521857547/2500000000000000:ℚ), by decide +kernel, by norm_num, ?_⟩
    have h := Segment23.interval_bound ⟨hl23,hh23⟩
    norm_num [Segment23.u,QInput.real,Rat.cast_div,Rat.cast_natCast,Rat.cast_intCast] at h ⊢; exact h
  have hl24 : Segment24.lo.real≤x := by
    have h := (not_le.mp hh23).le
    norm_num [Segment23.hi,Segment24.lo,QInput.real] at h ⊢; exact h
  by_cases hh24 : x≤Segment24.hi.real
  · refine ⟨(7974382895176137/2000000000000000:ℚ), by decide +kernel, by norm_num, ?_⟩
    have h := Segment24.interval_bound ⟨hl24,hh24⟩
    norm_num [Segment24.u,QInput.real,Rat.cast_div,Rat.cast_natCast,Rat.cast_intCast] at h ⊢; exact h
  have hl25 : Segment25.lo.real≤x := by
    have h := (not_le.mp hh24).le
    norm_num [Segment24.hi,Segment25.lo,QInput.real] at h ⊢; exact h
  by_cases hh25 : x≤Segment25.hi.real
  · refine ⟨(918788069754107/200000000000000:ℚ), by decide +kernel, by norm_num, ?_⟩
    have h := Segment25.interval_bound ⟨hl25,hh25⟩
    norm_num [Segment25.u,QInput.real,Rat.cast_div,Rat.cast_natCast,Rat.cast_intCast] at h ⊢; exact h
  have hl26 : Segment26.lo.real≤x := by
    have h := (not_le.mp hh25).le
    norm_num [Segment25.hi,Segment26.lo,QInput.real] at h ⊢; exact h
  by_cases hh26 : x≤Segment26.hi.real
  · refine ⟨(49510688747396809/10000000000000000:ℚ), by decide +kernel, by norm_num, ?_⟩
    have h := Segment26.interval_bound ⟨hl26,hh26⟩
    norm_num [Segment26.u,QInput.real,Rat.cast_div,Rat.cast_natCast,Rat.cast_intCast] at h ⊢; exact h
  have hl27 : Segment27.lo.real≤x := by
    have h := (not_le.mp hh26).le
    norm_num [Segment26.hi,Segment27.lo,QInput.real] at h ⊢; exact h
  by_cases hh27 : x≤Segment27.hi.real
  · refine ⟨(51330699179736179/10000000000000000:ℚ), by decide +kernel, by norm_num, ?_⟩
    have h := Segment27.interval_bound ⟨hl27,hh27⟩
    norm_num [Segment27.u,QInput.real,Rat.cast_div,Rat.cast_natCast,Rat.cast_intCast] at h ⊢; exact h
  have hl28 : Segment28.lo.real≤x := by
    have h := (not_le.mp hh27).le
    norm_num [Segment27.hi,Segment28.lo,QInput.real] at h ⊢; exact h
  by_cases hh28 : x≤Segment28.hi.real
  · refine ⟨(52546977145860421/10000000000000000:ℚ), by decide +kernel, by norm_num, ?_⟩
    have h := Segment28.interval_bound ⟨hl28,hh28⟩
    norm_num [Segment28.u,QInput.real,Rat.cast_div,Rat.cast_natCast,Rat.cast_intCast] at h ⊢; exact h
  have hl29 : Segment29.lo.real≤x := by
    have h := (not_le.mp hh28).le
    norm_num [Segment28.hi,Segment29.lo,QInput.real] at h ⊢; exact h
  by_cases hh29 : x≤Segment29.hi.real
  · refine ⟨(10848858740933581/2000000000000000:ℚ), by decide +kernel, by norm_num, ?_⟩
    have h := Segment29.interval_bound ⟨hl29,hh29⟩
    norm_num [Segment29.u,QInput.real,Rat.cast_div,Rat.cast_natCast,Rat.cast_intCast] at h ⊢; exact h
  have hl30 : Segment30.lo.real≤x := by
    have h := (not_le.mp hh29).le
    norm_num [Segment29.hi,Segment30.lo,QInput.real] at h ⊢; exact h
  have hh30 : x≤Segment30.hi.real := by norm_num [Segment30.hi,QInput.real]; exact hx.2
  refine ⟨(16961761337730057/2500000000000000:ℚ), by decide +kernel, by norm_num, ?_⟩
  have h := Segment30.interval_bound ⟨hl30,hh30⟩
  norm_num [Segment30.u,QInput.real,Rat.cast_div,Rat.cast_natCast,Rat.cast_intCast] at h ⊢; exact h
#print axioms rate_uniform
end Spin.Structured.ConcreteFixedNumeric
