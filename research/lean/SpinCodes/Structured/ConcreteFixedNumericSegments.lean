import SpinCodes.Structured.ConcreteFixedNumericData
import SpinCodes.Structured.ConcreteFixedNumericInterval

noncomputable section
namespace Spin.Structured.ConcreteFixedNumeric
open Spin.Numeric DenseOccupationFixed Set
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace Segment00
theorem bound_lo : rowRate lo.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_lo
theorem bound_hi : rowRate hi.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_hi
def supportLine : ℚ×ℚ := ((s.num:ℚ)/s.den,(c.num:ℚ)/c.den)
theorem support_mem : supportLine∈Spin.Majorant.refined.supports := by decide +kernel
theorem v_eq : v.real=Spin.mv (127/250) (3/1600) (1/524287) u.real := by norm_num [v,u,QInput.real,Spin.mv]
theorem u_pos : 0<u.real := by norm_num [u,QInput.real]
theorem interval_bound {x : ℝ} (hx : x∈Icc lo.real hi.real) : rate x u.real≤-(8679/10000000) := by
  have hh := rowRate_interval (by norm_num [lo,QInput.real]) (by norm_num [hi,QInput.real]) hx bound_lo bound_hi
  have hm := rate_le_rowRate support_mem x u.real v.real v_eq
  have hm2 : rate x u.real ≤ rowRate x s.real c.real u.real v.real := by
    simpa only [supportLine,QInput.real,Rat.cast_div,Rat.cast_intCast,Rat.cast_natCast] using hm
  exact hm2.trans hh
#print axioms interval_bound
end Segment00
namespace Segment01
theorem bound_lo : rowRate lo.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_lo
theorem bound_hi : rowRate hi.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_hi
def supportLine : ℚ×ℚ := ((s.num:ℚ)/s.den,(c.num:ℚ)/c.den)
theorem support_mem : supportLine∈Spin.Majorant.refined.supports := by decide +kernel
theorem v_eq : v.real=Spin.mv (127/250) (3/1600) (1/524287) u.real := by norm_num [v,u,QInput.real,Spin.mv]
theorem u_pos : 0<u.real := by norm_num [u,QInput.real]
theorem interval_bound {x : ℝ} (hx : x∈Icc lo.real hi.real) : rate x u.real≤-(8679/10000000) := by
  have hh := rowRate_interval (by norm_num [lo,QInput.real]) (by norm_num [hi,QInput.real]) hx bound_lo bound_hi
  have hm := rate_le_rowRate support_mem x u.real v.real v_eq
  have hm2 : rate x u.real ≤ rowRate x s.real c.real u.real v.real := by
    simpa only [supportLine,QInput.real,Rat.cast_div,Rat.cast_intCast,Rat.cast_natCast] using hm
  exact hm2.trans hh
#print axioms interval_bound
end Segment01
namespace Segment02
theorem bound_lo : rowRate lo.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_lo
theorem bound_hi : rowRate hi.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_hi
def supportLine : ℚ×ℚ := ((s.num:ℚ)/s.den,(c.num:ℚ)/c.den)
theorem support_mem : supportLine∈Spin.Majorant.refined.supports := by decide +kernel
theorem v_eq : v.real=Spin.mv (127/250) (3/1600) (1/524287) u.real := by norm_num [v,u,QInput.real,Spin.mv]
theorem u_pos : 0<u.real := by norm_num [u,QInput.real]
theorem interval_bound {x : ℝ} (hx : x∈Icc lo.real hi.real) : rate x u.real≤-(8679/10000000) := by
  have hh := rowRate_interval (by norm_num [lo,QInput.real]) (by norm_num [hi,QInput.real]) hx bound_lo bound_hi
  have hm := rate_le_rowRate support_mem x u.real v.real v_eq
  have hm2 : rate x u.real ≤ rowRate x s.real c.real u.real v.real := by
    simpa only [supportLine,QInput.real,Rat.cast_div,Rat.cast_intCast,Rat.cast_natCast] using hm
  exact hm2.trans hh
#print axioms interval_bound
end Segment02
namespace Segment03
theorem bound_lo : rowRate lo.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_lo
theorem bound_hi : rowRate hi.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_hi
def supportLine : ℚ×ℚ := ((s.num:ℚ)/s.den,(c.num:ℚ)/c.den)
theorem support_mem : supportLine∈Spin.Majorant.refined.supports := by decide +kernel
theorem v_eq : v.real=Spin.mv (127/250) (3/1600) (1/524287) u.real := by norm_num [v,u,QInput.real,Spin.mv]
theorem u_pos : 0<u.real := by norm_num [u,QInput.real]
theorem interval_bound {x : ℝ} (hx : x∈Icc lo.real hi.real) : rate x u.real≤-(8679/10000000) := by
  have hh := rowRate_interval (by norm_num [lo,QInput.real]) (by norm_num [hi,QInput.real]) hx bound_lo bound_hi
  have hm := rate_le_rowRate support_mem x u.real v.real v_eq
  have hm2 : rate x u.real ≤ rowRate x s.real c.real u.real v.real := by
    simpa only [supportLine,QInput.real,Rat.cast_div,Rat.cast_intCast,Rat.cast_natCast] using hm
  exact hm2.trans hh
#print axioms interval_bound
end Segment03
namespace Segment04
theorem bound_lo : rowRate lo.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_lo
theorem bound_hi : rowRate hi.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_hi
def supportLine : ℚ×ℚ := ((s.num:ℚ)/s.den,(c.num:ℚ)/c.den)
theorem support_mem : supportLine∈Spin.Majorant.refined.supports := by decide +kernel
theorem v_eq : v.real=Spin.mv (127/250) (3/1600) (1/524287) u.real := by norm_num [v,u,QInput.real,Spin.mv]
theorem u_pos : 0<u.real := by norm_num [u,QInput.real]
theorem interval_bound {x : ℝ} (hx : x∈Icc lo.real hi.real) : rate x u.real≤-(8679/10000000) := by
  have hh := rowRate_interval (by norm_num [lo,QInput.real]) (by norm_num [hi,QInput.real]) hx bound_lo bound_hi
  have hm := rate_le_rowRate support_mem x u.real v.real v_eq
  have hm2 : rate x u.real ≤ rowRate x s.real c.real u.real v.real := by
    simpa only [supportLine,QInput.real,Rat.cast_div,Rat.cast_intCast,Rat.cast_natCast] using hm
  exact hm2.trans hh
#print axioms interval_bound
end Segment04
namespace Segment05
theorem bound_lo : rowRate lo.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_lo
theorem bound_hi : rowRate hi.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_hi
def supportLine : ℚ×ℚ := ((s.num:ℚ)/s.den,(c.num:ℚ)/c.den)
theorem support_mem : supportLine∈Spin.Majorant.refined.supports := by decide +kernel
theorem v_eq : v.real=Spin.mv (127/250) (3/1600) (1/524287) u.real := by norm_num [v,u,QInput.real,Spin.mv]
theorem u_pos : 0<u.real := by norm_num [u,QInput.real]
theorem interval_bound {x : ℝ} (hx : x∈Icc lo.real hi.real) : rate x u.real≤-(8679/10000000) := by
  have hh := rowRate_interval (by norm_num [lo,QInput.real]) (by norm_num [hi,QInput.real]) hx bound_lo bound_hi
  have hm := rate_le_rowRate support_mem x u.real v.real v_eq
  have hm2 : rate x u.real ≤ rowRate x s.real c.real u.real v.real := by
    simpa only [supportLine,QInput.real,Rat.cast_div,Rat.cast_intCast,Rat.cast_natCast] using hm
  exact hm2.trans hh
#print axioms interval_bound
end Segment05
namespace Segment06
theorem bound_lo : rowRate lo.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_lo
theorem bound_hi : rowRate hi.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_hi
def supportLine : ℚ×ℚ := ((s.num:ℚ)/s.den,(c.num:ℚ)/c.den)
theorem support_mem : supportLine∈Spin.Majorant.refined.supports := by decide +kernel
theorem v_eq : v.real=Spin.mv (127/250) (3/1600) (1/524287) u.real := by norm_num [v,u,QInput.real,Spin.mv]
theorem u_pos : 0<u.real := by norm_num [u,QInput.real]
theorem interval_bound {x : ℝ} (hx : x∈Icc lo.real hi.real) : rate x u.real≤-(8679/10000000) := by
  have hh := rowRate_interval (by norm_num [lo,QInput.real]) (by norm_num [hi,QInput.real]) hx bound_lo bound_hi
  have hm := rate_le_rowRate support_mem x u.real v.real v_eq
  have hm2 : rate x u.real ≤ rowRate x s.real c.real u.real v.real := by
    simpa only [supportLine,QInput.real,Rat.cast_div,Rat.cast_intCast,Rat.cast_natCast] using hm
  exact hm2.trans hh
#print axioms interval_bound
end Segment06
namespace Segment07
theorem bound_lo : rowRate lo.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_lo
theorem bound_hi : rowRate hi.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_hi
def supportLine : ℚ×ℚ := ((s.num:ℚ)/s.den,(c.num:ℚ)/c.den)
theorem support_mem : supportLine∈Spin.Majorant.refined.supports := by decide +kernel
theorem v_eq : v.real=Spin.mv (127/250) (3/1600) (1/524287) u.real := by norm_num [v,u,QInput.real,Spin.mv]
theorem u_pos : 0<u.real := by norm_num [u,QInput.real]
theorem interval_bound {x : ℝ} (hx : x∈Icc lo.real hi.real) : rate x u.real≤-(8679/10000000) := by
  have hh := rowRate_interval (by norm_num [lo,QInput.real]) (by norm_num [hi,QInput.real]) hx bound_lo bound_hi
  have hm := rate_le_rowRate support_mem x u.real v.real v_eq
  have hm2 : rate x u.real ≤ rowRate x s.real c.real u.real v.real := by
    simpa only [supportLine,QInput.real,Rat.cast_div,Rat.cast_intCast,Rat.cast_natCast] using hm
  exact hm2.trans hh
#print axioms interval_bound
end Segment07
namespace Segment08
theorem bound_lo : rowRate lo.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_lo
theorem bound_hi : rowRate hi.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_hi
def supportLine : ℚ×ℚ := ((s.num:ℚ)/s.den,(c.num:ℚ)/c.den)
theorem support_mem : supportLine∈Spin.Majorant.refined.supports := by decide +kernel
theorem v_eq : v.real=Spin.mv (127/250) (3/1600) (1/524287) u.real := by norm_num [v,u,QInput.real,Spin.mv]
theorem u_pos : 0<u.real := by norm_num [u,QInput.real]
theorem interval_bound {x : ℝ} (hx : x∈Icc lo.real hi.real) : rate x u.real≤-(8679/10000000) := by
  have hh := rowRate_interval (by norm_num [lo,QInput.real]) (by norm_num [hi,QInput.real]) hx bound_lo bound_hi
  have hm := rate_le_rowRate support_mem x u.real v.real v_eq
  have hm2 : rate x u.real ≤ rowRate x s.real c.real u.real v.real := by
    simpa only [supportLine,QInput.real,Rat.cast_div,Rat.cast_intCast,Rat.cast_natCast] using hm
  exact hm2.trans hh
#print axioms interval_bound
end Segment08
namespace Segment09
theorem bound_lo : rowRate lo.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_lo
theorem bound_hi : rowRate hi.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_hi
def supportLine : ℚ×ℚ := ((s.num:ℚ)/s.den,(c.num:ℚ)/c.den)
theorem support_mem : supportLine∈Spin.Majorant.refined.supports := by decide +kernel
theorem v_eq : v.real=Spin.mv (127/250) (3/1600) (1/524287) u.real := by norm_num [v,u,QInput.real,Spin.mv]
theorem u_pos : 0<u.real := by norm_num [u,QInput.real]
theorem interval_bound {x : ℝ} (hx : x∈Icc lo.real hi.real) : rate x u.real≤-(8679/10000000) := by
  have hh := rowRate_interval (by norm_num [lo,QInput.real]) (by norm_num [hi,QInput.real]) hx bound_lo bound_hi
  have hm := rate_le_rowRate support_mem x u.real v.real v_eq
  have hm2 : rate x u.real ≤ rowRate x s.real c.real u.real v.real := by
    simpa only [supportLine,QInput.real,Rat.cast_div,Rat.cast_intCast,Rat.cast_natCast] using hm
  exact hm2.trans hh
#print axioms interval_bound
end Segment09
namespace Segment10
theorem bound_lo : rowRate lo.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_lo
theorem bound_hi : rowRate hi.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_hi
def supportLine : ℚ×ℚ := ((s.num:ℚ)/s.den,(c.num:ℚ)/c.den)
theorem support_mem : supportLine∈Spin.Majorant.refined.supports := by decide +kernel
theorem v_eq : v.real=Spin.mv (127/250) (3/1600) (1/524287) u.real := by norm_num [v,u,QInput.real,Spin.mv]
theorem u_pos : 0<u.real := by norm_num [u,QInput.real]
theorem interval_bound {x : ℝ} (hx : x∈Icc lo.real hi.real) : rate x u.real≤-(8679/10000000) := by
  have hh := rowRate_interval (by norm_num [lo,QInput.real]) (by norm_num [hi,QInput.real]) hx bound_lo bound_hi
  have hm := rate_le_rowRate support_mem x u.real v.real v_eq
  have hm2 : rate x u.real ≤ rowRate x s.real c.real u.real v.real := by
    simpa only [supportLine,QInput.real,Rat.cast_div,Rat.cast_intCast,Rat.cast_natCast] using hm
  exact hm2.trans hh
#print axioms interval_bound
end Segment10
namespace Segment11
theorem bound_lo : rowRate lo.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_lo
theorem bound_hi : rowRate hi.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_hi
def supportLine : ℚ×ℚ := ((s.num:ℚ)/s.den,(c.num:ℚ)/c.den)
theorem support_mem : supportLine∈Spin.Majorant.refined.supports := by decide +kernel
theorem v_eq : v.real=Spin.mv (127/250) (3/1600) (1/524287) u.real := by norm_num [v,u,QInput.real,Spin.mv]
theorem u_pos : 0<u.real := by norm_num [u,QInput.real]
theorem interval_bound {x : ℝ} (hx : x∈Icc lo.real hi.real) : rate x u.real≤-(8679/10000000) := by
  have hh := rowRate_interval (by norm_num [lo,QInput.real]) (by norm_num [hi,QInput.real]) hx bound_lo bound_hi
  have hm := rate_le_rowRate support_mem x u.real v.real v_eq
  have hm2 : rate x u.real ≤ rowRate x s.real c.real u.real v.real := by
    simpa only [supportLine,QInput.real,Rat.cast_div,Rat.cast_intCast,Rat.cast_natCast] using hm
  exact hm2.trans hh
#print axioms interval_bound
end Segment11
namespace Segment12
theorem bound_lo : rowRate lo.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_lo
theorem bound_hi : rowRate hi.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_hi
def supportLine : ℚ×ℚ := ((s.num:ℚ)/s.den,(c.num:ℚ)/c.den)
theorem support_mem : supportLine∈Spin.Majorant.refined.supports := by decide +kernel
theorem v_eq : v.real=Spin.mv (127/250) (3/1600) (1/524287) u.real := by norm_num [v,u,QInput.real,Spin.mv]
theorem u_pos : 0<u.real := by norm_num [u,QInput.real]
theorem interval_bound {x : ℝ} (hx : x∈Icc lo.real hi.real) : rate x u.real≤-(8679/10000000) := by
  have hh := rowRate_interval (by norm_num [lo,QInput.real]) (by norm_num [hi,QInput.real]) hx bound_lo bound_hi
  have hm := rate_le_rowRate support_mem x u.real v.real v_eq
  have hm2 : rate x u.real ≤ rowRate x s.real c.real u.real v.real := by
    simpa only [supportLine,QInput.real,Rat.cast_div,Rat.cast_intCast,Rat.cast_natCast] using hm
  exact hm2.trans hh
#print axioms interval_bound
end Segment12
namespace Segment13
theorem bound_lo : rowRate lo.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_lo
theorem bound_hi : rowRate hi.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_hi
def supportLine : ℚ×ℚ := ((s.num:ℚ)/s.den,(c.num:ℚ)/c.den)
theorem support_mem : supportLine∈Spin.Majorant.refined.supports := by decide +kernel
theorem v_eq : v.real=Spin.mv (127/250) (3/1600) (1/524287) u.real := by norm_num [v,u,QInput.real,Spin.mv]
theorem u_pos : 0<u.real := by norm_num [u,QInput.real]
theorem interval_bound {x : ℝ} (hx : x∈Icc lo.real hi.real) : rate x u.real≤-(8679/10000000) := by
  have hh := rowRate_interval (by norm_num [lo,QInput.real]) (by norm_num [hi,QInput.real]) hx bound_lo bound_hi
  have hm := rate_le_rowRate support_mem x u.real v.real v_eq
  have hm2 : rate x u.real ≤ rowRate x s.real c.real u.real v.real := by
    simpa only [supportLine,QInput.real,Rat.cast_div,Rat.cast_intCast,Rat.cast_natCast] using hm
  exact hm2.trans hh
#print axioms interval_bound
end Segment13
namespace Segment14
theorem bound_lo : rowRate lo.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_lo
theorem bound_hi : rowRate hi.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_hi
def supportLine : ℚ×ℚ := ((s.num:ℚ)/s.den,(c.num:ℚ)/c.den)
theorem support_mem : supportLine∈Spin.Majorant.refined.supports := by decide +kernel
theorem v_eq : v.real=Spin.mv (127/250) (3/1600) (1/524287) u.real := by norm_num [v,u,QInput.real,Spin.mv]
theorem u_pos : 0<u.real := by norm_num [u,QInput.real]
theorem interval_bound {x : ℝ} (hx : x∈Icc lo.real hi.real) : rate x u.real≤-(8679/10000000) := by
  have hh := rowRate_interval (by norm_num [lo,QInput.real]) (by norm_num [hi,QInput.real]) hx bound_lo bound_hi
  have hm := rate_le_rowRate support_mem x u.real v.real v_eq
  have hm2 : rate x u.real ≤ rowRate x s.real c.real u.real v.real := by
    simpa only [supportLine,QInput.real,Rat.cast_div,Rat.cast_intCast,Rat.cast_natCast] using hm
  exact hm2.trans hh
#print axioms interval_bound
end Segment14
namespace Segment15
theorem bound_lo : rowRate lo.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_lo
theorem bound_hi : rowRate hi.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_hi
def supportLine : ℚ×ℚ := ((s.num:ℚ)/s.den,(c.num:ℚ)/c.den)
theorem support_mem : supportLine∈Spin.Majorant.refined.supports := by decide +kernel
theorem v_eq : v.real=Spin.mv (127/250) (3/1600) (1/524287) u.real := by norm_num [v,u,QInput.real,Spin.mv]
theorem u_pos : 0<u.real := by norm_num [u,QInput.real]
theorem interval_bound {x : ℝ} (hx : x∈Icc lo.real hi.real) : rate x u.real≤-(8679/10000000) := by
  have hh := rowRate_interval (by norm_num [lo,QInput.real]) (by norm_num [hi,QInput.real]) hx bound_lo bound_hi
  have hm := rate_le_rowRate support_mem x u.real v.real v_eq
  have hm2 : rate x u.real ≤ rowRate x s.real c.real u.real v.real := by
    simpa only [supportLine,QInput.real,Rat.cast_div,Rat.cast_intCast,Rat.cast_natCast] using hm
  exact hm2.trans hh
#print axioms interval_bound
end Segment15
namespace Segment16
theorem bound_lo : rowRate lo.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_lo
theorem bound_hi : rowRate hi.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_hi
def supportLine : ℚ×ℚ := ((s.num:ℚ)/s.den,(c.num:ℚ)/c.den)
theorem support_mem : supportLine∈Spin.Majorant.refined.supports := by decide +kernel
theorem v_eq : v.real=Spin.mv (127/250) (3/1600) (1/524287) u.real := by norm_num [v,u,QInput.real,Spin.mv]
theorem u_pos : 0<u.real := by norm_num [u,QInput.real]
theorem interval_bound {x : ℝ} (hx : x∈Icc lo.real hi.real) : rate x u.real≤-(8679/10000000) := by
  have hh := rowRate_interval (by norm_num [lo,QInput.real]) (by norm_num [hi,QInput.real]) hx bound_lo bound_hi
  have hm := rate_le_rowRate support_mem x u.real v.real v_eq
  have hm2 : rate x u.real ≤ rowRate x s.real c.real u.real v.real := by
    simpa only [supportLine,QInput.real,Rat.cast_div,Rat.cast_intCast,Rat.cast_natCast] using hm
  exact hm2.trans hh
#print axioms interval_bound
end Segment16
namespace Segment17
theorem bound_lo : rowRate lo.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_lo
theorem bound_hi : rowRate hi.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_hi
def supportLine : ℚ×ℚ := ((s.num:ℚ)/s.den,(c.num:ℚ)/c.den)
theorem support_mem : supportLine∈Spin.Majorant.refined.supports := by decide +kernel
theorem v_eq : v.real=Spin.mv (127/250) (3/1600) (1/524287) u.real := by norm_num [v,u,QInput.real,Spin.mv]
theorem u_pos : 0<u.real := by norm_num [u,QInput.real]
theorem interval_bound {x : ℝ} (hx : x∈Icc lo.real hi.real) : rate x u.real≤-(8679/10000000) := by
  have hh := rowRate_interval (by norm_num [lo,QInput.real]) (by norm_num [hi,QInput.real]) hx bound_lo bound_hi
  have hm := rate_le_rowRate support_mem x u.real v.real v_eq
  have hm2 : rate x u.real ≤ rowRate x s.real c.real u.real v.real := by
    simpa only [supportLine,QInput.real,Rat.cast_div,Rat.cast_intCast,Rat.cast_natCast] using hm
  exact hm2.trans hh
#print axioms interval_bound
end Segment17
namespace Segment18
theorem bound_lo : rowRate lo.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_lo
theorem bound_hi : rowRate hi.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_hi
def supportLine : ℚ×ℚ := ((s.num:ℚ)/s.den,(c.num:ℚ)/c.den)
theorem support_mem : supportLine∈Spin.Majorant.refined.supports := by decide +kernel
theorem v_eq : v.real=Spin.mv (127/250) (3/1600) (1/524287) u.real := by norm_num [v,u,QInput.real,Spin.mv]
theorem u_pos : 0<u.real := by norm_num [u,QInput.real]
theorem interval_bound {x : ℝ} (hx : x∈Icc lo.real hi.real) : rate x u.real≤-(8679/10000000) := by
  have hh := rowRate_interval (by norm_num [lo,QInput.real]) (by norm_num [hi,QInput.real]) hx bound_lo bound_hi
  have hm := rate_le_rowRate support_mem x u.real v.real v_eq
  have hm2 : rate x u.real ≤ rowRate x s.real c.real u.real v.real := by
    simpa only [supportLine,QInput.real,Rat.cast_div,Rat.cast_intCast,Rat.cast_natCast] using hm
  exact hm2.trans hh
#print axioms interval_bound
end Segment18
namespace Segment19
theorem bound_lo : rowRate lo.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_lo
theorem bound_hi : rowRate hi.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_hi
def supportLine : ℚ×ℚ := ((s.num:ℚ)/s.den,(c.num:ℚ)/c.den)
theorem support_mem : supportLine∈Spin.Majorant.refined.supports := by decide +kernel
theorem v_eq : v.real=Spin.mv (127/250) (3/1600) (1/524287) u.real := by norm_num [v,u,QInput.real,Spin.mv]
theorem u_pos : 0<u.real := by norm_num [u,QInput.real]
theorem interval_bound {x : ℝ} (hx : x∈Icc lo.real hi.real) : rate x u.real≤-(8679/10000000) := by
  have hh := rowRate_interval (by norm_num [lo,QInput.real]) (by norm_num [hi,QInput.real]) hx bound_lo bound_hi
  have hm := rate_le_rowRate support_mem x u.real v.real v_eq
  have hm2 : rate x u.real ≤ rowRate x s.real c.real u.real v.real := by
    simpa only [supportLine,QInput.real,Rat.cast_div,Rat.cast_intCast,Rat.cast_natCast] using hm
  exact hm2.trans hh
#print axioms interval_bound
end Segment19
namespace Segment20
theorem bound_lo : rowRate lo.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_lo
theorem bound_hi : rowRate hi.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_hi
def supportLine : ℚ×ℚ := ((s.num:ℚ)/s.den,(c.num:ℚ)/c.den)
theorem support_mem : supportLine∈Spin.Majorant.refined.supports := by decide +kernel
theorem v_eq : v.real=Spin.mv (127/250) (3/1600) (1/524287) u.real := by norm_num [v,u,QInput.real,Spin.mv]
theorem u_pos : 0<u.real := by norm_num [u,QInput.real]
theorem interval_bound {x : ℝ} (hx : x∈Icc lo.real hi.real) : rate x u.real≤-(8679/10000000) := by
  have hh := rowRate_interval (by norm_num [lo,QInput.real]) (by norm_num [hi,QInput.real]) hx bound_lo bound_hi
  have hm := rate_le_rowRate support_mem x u.real v.real v_eq
  have hm2 : rate x u.real ≤ rowRate x s.real c.real u.real v.real := by
    simpa only [supportLine,QInput.real,Rat.cast_div,Rat.cast_intCast,Rat.cast_natCast] using hm
  exact hm2.trans hh
#print axioms interval_bound
end Segment20
namespace Segment21
theorem bound_lo : rowRate lo.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_lo
theorem bound_hi : rowRate hi.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_hi
def supportLine : ℚ×ℚ := ((s.num:ℚ)/s.den,(c.num:ℚ)/c.den)
theorem support_mem : supportLine∈Spin.Majorant.refined.supports := by decide +kernel
theorem v_eq : v.real=Spin.mv (127/250) (3/1600) (1/524287) u.real := by norm_num [v,u,QInput.real,Spin.mv]
theorem u_pos : 0<u.real := by norm_num [u,QInput.real]
theorem interval_bound {x : ℝ} (hx : x∈Icc lo.real hi.real) : rate x u.real≤-(8679/10000000) := by
  have hh := rowRate_interval (by norm_num [lo,QInput.real]) (by norm_num [hi,QInput.real]) hx bound_lo bound_hi
  have hm := rate_le_rowRate support_mem x u.real v.real v_eq
  have hm2 : rate x u.real ≤ rowRate x s.real c.real u.real v.real := by
    simpa only [supportLine,QInput.real,Rat.cast_div,Rat.cast_intCast,Rat.cast_natCast] using hm
  exact hm2.trans hh
#print axioms interval_bound
end Segment21
namespace Segment22
theorem bound_lo : rowRate lo.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_lo
theorem bound_hi : rowRate hi.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_hi
def supportLine : ℚ×ℚ := ((s.num:ℚ)/s.den,(c.num:ℚ)/c.den)
theorem support_mem : supportLine∈Spin.Majorant.refined.supports := by decide +kernel
theorem v_eq : v.real=Spin.mv (127/250) (3/1600) (1/524287) u.real := by norm_num [v,u,QInput.real,Spin.mv]
theorem u_pos : 0<u.real := by norm_num [u,QInput.real]
theorem interval_bound {x : ℝ} (hx : x∈Icc lo.real hi.real) : rate x u.real≤-(8679/10000000) := by
  have hh := rowRate_interval (by norm_num [lo,QInput.real]) (by norm_num [hi,QInput.real]) hx bound_lo bound_hi
  have hm := rate_le_rowRate support_mem x u.real v.real v_eq
  have hm2 : rate x u.real ≤ rowRate x s.real c.real u.real v.real := by
    simpa only [supportLine,QInput.real,Rat.cast_div,Rat.cast_intCast,Rat.cast_natCast] using hm
  exact hm2.trans hh
#print axioms interval_bound
end Segment22
namespace Segment23
theorem bound_lo : rowRate lo.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_lo
theorem bound_hi : rowRate hi.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_hi
def supportLine : ℚ×ℚ := ((s.num:ℚ)/s.den,(c.num:ℚ)/c.den)
theorem support_mem : supportLine∈Spin.Majorant.refined.supports := by decide +kernel
theorem v_eq : v.real=Spin.mv (127/250) (3/1600) (1/524287) u.real := by norm_num [v,u,QInput.real,Spin.mv]
theorem u_pos : 0<u.real := by norm_num [u,QInput.real]
theorem interval_bound {x : ℝ} (hx : x∈Icc lo.real hi.real) : rate x u.real≤-(8679/10000000) := by
  have hh := rowRate_interval (by norm_num [lo,QInput.real]) (by norm_num [hi,QInput.real]) hx bound_lo bound_hi
  have hm := rate_le_rowRate support_mem x u.real v.real v_eq
  have hm2 : rate x u.real ≤ rowRate x s.real c.real u.real v.real := by
    simpa only [supportLine,QInput.real,Rat.cast_div,Rat.cast_intCast,Rat.cast_natCast] using hm
  exact hm2.trans hh
#print axioms interval_bound
end Segment23
namespace Segment24
theorem bound_lo : rowRate lo.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_lo
theorem bound_hi : rowRate hi.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_hi
def supportLine : ℚ×ℚ := ((s.num:ℚ)/s.den,(c.num:ℚ)/c.den)
theorem support_mem : supportLine∈Spin.Majorant.refined.supports := by decide +kernel
theorem v_eq : v.real=Spin.mv (127/250) (3/1600) (1/524287) u.real := by norm_num [v,u,QInput.real,Spin.mv]
theorem u_pos : 0<u.real := by norm_num [u,QInput.real]
theorem interval_bound {x : ℝ} (hx : x∈Icc lo.real hi.real) : rate x u.real≤-(8679/10000000) := by
  have hh := rowRate_interval (by norm_num [lo,QInput.real]) (by norm_num [hi,QInput.real]) hx bound_lo bound_hi
  have hm := rate_le_rowRate support_mem x u.real v.real v_eq
  have hm2 : rate x u.real ≤ rowRate x s.real c.real u.real v.real := by
    simpa only [supportLine,QInput.real,Rat.cast_div,Rat.cast_intCast,Rat.cast_natCast] using hm
  exact hm2.trans hh
#print axioms interval_bound
end Segment24
namespace Segment25
theorem bound_lo : rowRate lo.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_lo
theorem bound_hi : rowRate hi.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_hi
def supportLine : ℚ×ℚ := ((s.num:ℚ)/s.den,(c.num:ℚ)/c.den)
theorem support_mem : supportLine∈Spin.Majorant.refined.supports := by decide +kernel
theorem v_eq : v.real=Spin.mv (127/250) (3/1600) (1/524287) u.real := by norm_num [v,u,QInput.real,Spin.mv]
theorem u_pos : 0<u.real := by norm_num [u,QInput.real]
theorem interval_bound {x : ℝ} (hx : x∈Icc lo.real hi.real) : rate x u.real≤-(8679/10000000) := by
  have hh := rowRate_interval (by norm_num [lo,QInput.real]) (by norm_num [hi,QInput.real]) hx bound_lo bound_hi
  have hm := rate_le_rowRate support_mem x u.real v.real v_eq
  have hm2 : rate x u.real ≤ rowRate x s.real c.real u.real v.real := by
    simpa only [supportLine,QInput.real,Rat.cast_div,Rat.cast_intCast,Rat.cast_natCast] using hm
  exact hm2.trans hh
#print axioms interval_bound
end Segment25
namespace Segment26
theorem bound_lo : rowRate lo.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_lo
theorem bound_hi : rowRate hi.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_hi
def supportLine : ℚ×ℚ := ((s.num:ℚ)/s.den,(c.num:ℚ)/c.den)
theorem support_mem : supportLine∈Spin.Majorant.refined.supports := by decide +kernel
theorem v_eq : v.real=Spin.mv (127/250) (3/1600) (1/524287) u.real := by norm_num [v,u,QInput.real,Spin.mv]
theorem u_pos : 0<u.real := by norm_num [u,QInput.real]
theorem interval_bound {x : ℝ} (hx : x∈Icc lo.real hi.real) : rate x u.real≤-(8679/10000000) := by
  have hh := rowRate_interval (by norm_num [lo,QInput.real]) (by norm_num [hi,QInput.real]) hx bound_lo bound_hi
  have hm := rate_le_rowRate support_mem x u.real v.real v_eq
  have hm2 : rate x u.real ≤ rowRate x s.real c.real u.real v.real := by
    simpa only [supportLine,QInput.real,Rat.cast_div,Rat.cast_intCast,Rat.cast_natCast] using hm
  exact hm2.trans hh
#print axioms interval_bound
end Segment26
namespace Segment27
theorem bound_lo : rowRate lo.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_lo
theorem bound_hi : rowRate hi.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_hi
def supportLine : ℚ×ℚ := ((s.num:ℚ)/s.den,(c.num:ℚ)/c.den)
theorem support_mem : supportLine∈Spin.Majorant.refined.supports := by decide +kernel
theorem v_eq : v.real=Spin.mv (127/250) (3/1600) (1/524287) u.real := by norm_num [v,u,QInput.real,Spin.mv]
theorem u_pos : 0<u.real := by norm_num [u,QInput.real]
theorem interval_bound {x : ℝ} (hx : x∈Icc lo.real hi.real) : rate x u.real≤-(8679/10000000) := by
  have hh := rowRate_interval (by norm_num [lo,QInput.real]) (by norm_num [hi,QInput.real]) hx bound_lo bound_hi
  have hm := rate_le_rowRate support_mem x u.real v.real v_eq
  have hm2 : rate x u.real ≤ rowRate x s.real c.real u.real v.real := by
    simpa only [supportLine,QInput.real,Rat.cast_div,Rat.cast_intCast,Rat.cast_natCast] using hm
  exact hm2.trans hh
#print axioms interval_bound
end Segment27
namespace Segment28
theorem bound_lo : rowRate lo.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_lo
theorem bound_hi : rowRate hi.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_hi
def supportLine : ℚ×ℚ := ((s.num:ℚ)/s.den,(c.num:ℚ)/c.den)
theorem support_mem : supportLine∈Spin.Majorant.refined.supports := by decide +kernel
theorem v_eq : v.real=Spin.mv (127/250) (3/1600) (1/524287) u.real := by norm_num [v,u,QInput.real,Spin.mv]
theorem u_pos : 0<u.real := by norm_num [u,QInput.real]
theorem interval_bound {x : ℝ} (hx : x∈Icc lo.real hi.real) : rate x u.real≤-(8679/10000000) := by
  have hh := rowRate_interval (by norm_num [lo,QInput.real]) (by norm_num [hi,QInput.real]) hx bound_lo bound_hi
  have hm := rate_le_rowRate support_mem x u.real v.real v_eq
  have hm2 : rate x u.real ≤ rowRate x s.real c.real u.real v.real := by
    simpa only [supportLine,QInput.real,Rat.cast_div,Rat.cast_intCast,Rat.cast_natCast] using hm
  exact hm2.trans hh
#print axioms interval_bound
end Segment28
namespace Segment29
theorem bound_lo : rowRate lo.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_lo
theorem bound_hi : rowRate hi.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_hi
def supportLine : ℚ×ℚ := ((s.num:ℚ)/s.den,(c.num:ℚ)/c.den)
theorem support_mem : supportLine∈Spin.Majorant.refined.supports := by decide +kernel
theorem v_eq : v.real=Spin.mv (127/250) (3/1600) (1/524287) u.real := by norm_num [v,u,QInput.real,Spin.mv]
theorem u_pos : 0<u.real := by norm_num [u,QInput.real]
theorem interval_bound {x : ℝ} (hx : x∈Icc lo.real hi.real) : rate x u.real≤-(8679/10000000) := by
  have hh := rowRate_interval (by norm_num [lo,QInput.real]) (by norm_num [hi,QInput.real]) hx bound_lo bound_hi
  have hm := rate_le_rowRate support_mem x u.real v.real v_eq
  have hm2 : rate x u.real ≤ rowRate x s.real c.real u.real v.real := by
    simpa only [supportLine,QInput.real,Rat.cast_div,Rat.cast_intCast,Rat.cast_natCast] using hm
  exact hm2.trans hh
#print axioms interval_bound
end Segment29
namespace Segment30
theorem bound_lo : rowRate lo.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_lo
theorem bound_hi : rowRate hi.real s.real c.real u.real v.real ≤ -(8679/10000000) :=
  rowCheck_sound (by decide) (by decide) (by decide) check_hi
def supportLine : ℚ×ℚ := ((s.num:ℚ)/s.den,(c.num:ℚ)/c.den)
theorem support_mem : supportLine∈Spin.Majorant.refined.supports := by decide +kernel
theorem v_eq : v.real=Spin.mv (127/250) (3/1600) (1/524287) u.real := by norm_num [v,u,QInput.real,Spin.mv]
theorem u_pos : 0<u.real := by norm_num [u,QInput.real]
theorem interval_bound {x : ℝ} (hx : x∈Icc lo.real hi.real) : rate x u.real≤-(8679/10000000) := by
  have hh := rowRate_interval (by norm_num [lo,QInput.real]) (by norm_num [hi,QInput.real]) hx bound_lo bound_hi
  have hm := rate_le_rowRate support_mem x u.real v.real v_eq
  have hm2 : rate x u.real ≤ rowRate x s.real c.real u.real v.real := by
    simpa only [supportLine,QInput.real,Rat.cast_div,Rat.cast_intCast,Rat.cast_natCast] using hm
  exact hm2.trans hh
#print axioms interval_bound
end Segment30
end Spin.Structured.ConcreteFixedNumeric
