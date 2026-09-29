import SpinCodes.Structured.SparsePowers
import SpinCodes.Structured.PolyIdentity
namespace Spin.Structured.SparsePowers
set_option maxRecDepth 100000
set_option maxHeartbeats 0
theorem z0_eval (x : ℝ) : z0.eval x = (1 - x / 6250) ^ 0 := by
  norm_num [z0, RatPoly.eval, listEval]
theorem z1_eval (x : ℝ) : z1.eval x = (1 - x / 6250) := by
  simp only [z1, RatPoly.eval, listEval]
  push_cast
  ring
theorem z2_eval (x : ℝ) : z2.eval x = (1 - x / 6250) ^ 2 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step1_checked,
    z1_eval, pow_two]
theorem z3_eval (x : ℝ) : z3.eval x = (1 - x / 6250) ^ 3 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step2_checked,
    z2_eval, z1_eval, ← pow_succ]
theorem z4_eval (x : ℝ) : z4.eval x = (1 - x / 6250) ^ 4 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step3_checked,
    z3_eval, z1_eval, ← pow_succ]
theorem z5_eval (x : ℝ) : z5.eval x = (1 - x / 6250) ^ 5 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step4_checked,
    z4_eval, z1_eval, ← pow_succ]
theorem z6_eval (x : ℝ) : z6.eval x = (1 - x / 6250) ^ 6 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step5_checked,
    z5_eval, z1_eval, ← pow_succ]
theorem z7_eval (x : ℝ) : z7.eval x = (1 - x / 6250) ^ 7 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step6_checked,
    z6_eval, z1_eval, ← pow_succ]
theorem z8_eval (x : ℝ) : z8.eval x = (1 - x / 6250) ^ 8 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step7_checked,
    z7_eval, z1_eval, ← pow_succ]
theorem z9_eval (x : ℝ) : z9.eval x = (1 - x / 6250) ^ 9 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step8_checked,
    z8_eval, z1_eval, ← pow_succ]
theorem z10_eval (x : ℝ) : z10.eval x = (1 - x / 6250) ^ 10 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step9_checked,
    z9_eval, z1_eval, ← pow_succ]
theorem z11_eval (x : ℝ) : z11.eval x = (1 - x / 6250) ^ 11 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step10_checked,
    z10_eval, z1_eval, ← pow_succ]
theorem z12_eval (x : ℝ) : z12.eval x = (1 - x / 6250) ^ 12 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step11_checked,
    z11_eval, z1_eval, ← pow_succ]
theorem z13_eval (x : ℝ) : z13.eval x = (1 - x / 6250) ^ 13 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step12_checked,
    z12_eval, z1_eval, ← pow_succ]
theorem z14_eval (x : ℝ) : z14.eval x = (1 - x / 6250) ^ 14 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step13_checked,
    z13_eval, z1_eval, ← pow_succ]
theorem z15_eval (x : ℝ) : z15.eval x = (1 - x / 6250) ^ 15 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step14_checked,
    z14_eval, z1_eval, ← pow_succ]
theorem z16_eval (x : ℝ) : z16.eval x = (1 - x / 6250) ^ 16 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step15_checked,
    z15_eval, z1_eval, ← pow_succ]
theorem z17_eval (x : ℝ) : z17.eval x = (1 - x / 6250) ^ 17 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step16_checked,
    z16_eval, z1_eval, ← pow_succ]
theorem z18_eval (x : ℝ) : z18.eval x = (1 - x / 6250) ^ 18 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step17_checked,
    z17_eval, z1_eval, ← pow_succ]
theorem z19_eval (x : ℝ) : z19.eval x = (1 - x / 6250) ^ 19 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step18_checked,
    z18_eval, z1_eval, ← pow_succ]
theorem z20_eval (x : ℝ) : z20.eval x = (1 - x / 6250) ^ 20 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step19_checked,
    z19_eval, z1_eval, ← pow_succ]
theorem z21_eval (x : ℝ) : z21.eval x = (1 - x / 6250) ^ 21 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step20_checked,
    z20_eval, z1_eval, ← pow_succ]
theorem z22_eval (x : ℝ) : z22.eval x = (1 - x / 6250) ^ 22 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step21_checked,
    z21_eval, z1_eval, ← pow_succ]
theorem z23_eval (x : ℝ) : z23.eval x = (1 - x / 6250) ^ 23 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step22_checked,
    z22_eval, z1_eval, ← pow_succ]
theorem z24_eval (x : ℝ) : z24.eval x = (1 - x / 6250) ^ 24 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step23_checked,
    z23_eval, z1_eval, ← pow_succ]
theorem z25_eval (x : ℝ) : z25.eval x = (1 - x / 6250) ^ 25 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step24_checked,
    z24_eval, z1_eval, ← pow_succ]
theorem z26_eval (x : ℝ) : z26.eval x = (1 - x / 6250) ^ 26 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step25_checked,
    z25_eval, z1_eval, ← pow_succ]
theorem z27_eval (x : ℝ) : z27.eval x = (1 - x / 6250) ^ 27 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step26_checked,
    z26_eval, z1_eval, ← pow_succ]
theorem z28_eval (x : ℝ) : z28.eval x = (1 - x / 6250) ^ 28 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step27_checked,
    z27_eval, z1_eval, ← pow_succ]
theorem z29_eval (x : ℝ) : z29.eval x = (1 - x / 6250) ^ 29 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step28_checked,
    z28_eval, z1_eval, ← pow_succ]
theorem z30_eval (x : ℝ) : z30.eval x = (1 - x / 6250) ^ 30 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step29_checked,
    z29_eval, z1_eval, ← pow_succ]
theorem z31_eval (x : ℝ) : z31.eval x = (1 - x / 6250) ^ 31 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step30_checked,
    z30_eval, z1_eval, ← pow_succ]
theorem z32_eval (x : ℝ) : z32.eval x = (1 - x / 6250) ^ 32 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step31_checked,
    z31_eval, z1_eval, ← pow_succ]
theorem z33_eval (x : ℝ) : z33.eval x = (1 - x / 6250) ^ 33 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step32_checked,
    z32_eval, z1_eval, ← pow_succ]
theorem z34_eval (x : ℝ) : z34.eval x = (1 - x / 6250) ^ 34 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step33_checked,
    z33_eval, z1_eval, ← pow_succ]
theorem z35_eval (x : ℝ) : z35.eval x = (1 - x / 6250) ^ 35 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step34_checked,
    z34_eval, z1_eval, ← pow_succ]
theorem z36_eval (x : ℝ) : z36.eval x = (1 - x / 6250) ^ 36 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step35_checked,
    z35_eval, z1_eval, ← pow_succ]
theorem z37_eval (x : ℝ) : z37.eval x = (1 - x / 6250) ^ 37 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step36_checked,
    z36_eval, z1_eval, ← pow_succ]
theorem z38_eval (x : ℝ) : z38.eval x = (1 - x / 6250) ^ 38 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step37_checked,
    z37_eval, z1_eval, ← pow_succ]
theorem z39_eval (x : ℝ) : z39.eval x = (1 - x / 6250) ^ 39 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step38_checked,
    z38_eval, z1_eval, ← pow_succ]
theorem z40_eval (x : ℝ) : z40.eval x = (1 - x / 6250) ^ 40 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step39_checked,
    z39_eval, z1_eval, ← pow_succ]
theorem z41_eval (x : ℝ) : z41.eval x = (1 - x / 6250) ^ 41 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step40_checked,
    z40_eval, z1_eval, ← pow_succ]
theorem z42_eval (x : ℝ) : z42.eval x = (1 - x / 6250) ^ 42 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step41_checked,
    z41_eval, z1_eval, ← pow_succ]
theorem z43_eval (x : ℝ) : z43.eval x = (1 - x / 6250) ^ 43 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step42_checked,
    z42_eval, z1_eval, ← pow_succ]
theorem z44_eval (x : ℝ) : z44.eval x = (1 - x / 6250) ^ 44 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step43_checked,
    z43_eval, z1_eval, ← pow_succ]
theorem z45_eval (x : ℝ) : z45.eval x = (1 - x / 6250) ^ 45 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step44_checked,
    z44_eval, z1_eval, ← pow_succ]
theorem z46_eval (x : ℝ) : z46.eval x = (1 - x / 6250) ^ 46 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step45_checked,
    z45_eval, z1_eval, ← pow_succ]
theorem z47_eval (x : ℝ) : z47.eval x = (1 - x / 6250) ^ 47 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step46_checked,
    z46_eval, z1_eval, ← pow_succ]
theorem z48_eval (x : ℝ) : z48.eval x = (1 - x / 6250) ^ 48 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step47_checked,
    z47_eval, z1_eval, ← pow_succ]
theorem z49_eval (x : ℝ) : z49.eval x = (1 - x / 6250) ^ 49 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step48_checked,
    z48_eval, z1_eval, ← pow_succ]
theorem z50_eval (x : ℝ) : z50.eval x = (1 - x / 6250) ^ 50 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step49_checked,
    z49_eval, z1_eval, ← pow_succ]
theorem z51_eval (x : ℝ) : z51.eval x = (1 - x / 6250) ^ 51 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step50_checked,
    z50_eval, z1_eval, ← pow_succ]
theorem z52_eval (x : ℝ) : z52.eval x = (1 - x / 6250) ^ 52 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step51_checked,
    z51_eval, z1_eval, ← pow_succ]
theorem z53_eval (x : ℝ) : z53.eval x = (1 - x / 6250) ^ 53 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step52_checked,
    z52_eval, z1_eval, ← pow_succ]
theorem z54_eval (x : ℝ) : z54.eval x = (1 - x / 6250) ^ 54 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step53_checked,
    z53_eval, z1_eval, ← pow_succ]
theorem z55_eval (x : ℝ) : z55.eval x = (1 - x / 6250) ^ 55 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step54_checked,
    z54_eval, z1_eval, ← pow_succ]
theorem z56_eval (x : ℝ) : z56.eval x = (1 - x / 6250) ^ 56 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step55_checked,
    z55_eval, z1_eval, ← pow_succ]
theorem z57_eval (x : ℝ) : z57.eval x = (1 - x / 6250) ^ 57 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step56_checked,
    z56_eval, z1_eval, ← pow_succ]
theorem z58_eval (x : ℝ) : z58.eval x = (1 - x / 6250) ^ 58 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step57_checked,
    z57_eval, z1_eval, ← pow_succ]
theorem z59_eval (x : ℝ) : z59.eval x = (1 - x / 6250) ^ 59 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step58_checked,
    z58_eval, z1_eval, ← pow_succ]
theorem z60_eval (x : ℝ) : z60.eval x = (1 - x / 6250) ^ 60 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step59_checked,
    z59_eval, z1_eval, ← pow_succ]
theorem z61_eval (x : ℝ) : z61.eval x = (1 - x / 6250) ^ 61 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step60_checked,
    z60_eval, z1_eval, ← pow_succ]
theorem z62_eval (x : ℝ) : z62.eval x = (1 - x / 6250) ^ 62 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step61_checked,
    z61_eval, z1_eval, ← pow_succ]
theorem z63_eval (x : ℝ) : z63.eval x = (1 - x / 6250) ^ 63 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step62_checked,
    z62_eval, z1_eval, ← pow_succ]
theorem z64_eval (x : ℝ) : z64.eval x = (1 - x / 6250) ^ 64 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step63_checked,
    z63_eval, z1_eval, ← pow_succ]
theorem z65_eval (x : ℝ) : z65.eval x = (1 - x / 6250) ^ 65 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step64_checked,
    z64_eval, z1_eval, ← pow_succ]
theorem z66_eval (x : ℝ) : z66.eval x = (1 - x / 6250) ^ 66 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step65_checked,
    z65_eval, z1_eval, ← pow_succ]
theorem z67_eval (x : ℝ) : z67.eval x = (1 - x / 6250) ^ 67 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step66_checked,
    z66_eval, z1_eval, ← pow_succ]
theorem z68_eval (x : ℝ) : z68.eval x = (1 - x / 6250) ^ 68 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step67_checked,
    z67_eval, z1_eval, ← pow_succ]
theorem z69_eval (x : ℝ) : z69.eval x = (1 - x / 6250) ^ 69 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step68_checked,
    z68_eval, z1_eval, ← pow_succ]
theorem z70_eval (x : ℝ) : z70.eval x = (1 - x / 6250) ^ 70 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step69_checked,
    z69_eval, z1_eval, ← pow_succ]
theorem z71_eval (x : ℝ) : z71.eval x = (1 - x / 6250) ^ 71 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step70_checked,
    z70_eval, z1_eval, ← pow_succ]
theorem z72_eval (x : ℝ) : z72.eval x = (1 - x / 6250) ^ 72 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step71_checked,
    z71_eval, z1_eval, ← pow_succ]
theorem z73_eval (x : ℝ) : z73.eval x = (1 - x / 6250) ^ 73 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step72_checked,
    z72_eval, z1_eval, ← pow_succ]
theorem z74_eval (x : ℝ) : z74.eval x = (1 - x / 6250) ^ 74 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step73_checked,
    z73_eval, z1_eval, ← pow_succ]
theorem z75_eval (x : ℝ) : z75.eval x = (1 - x / 6250) ^ 75 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step74_checked,
    z74_eval, z1_eval, ← pow_succ]
theorem z76_eval (x : ℝ) : z76.eval x = (1 - x / 6250) ^ 76 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step75_checked,
    z75_eval, z1_eval, ← pow_succ]
theorem z77_eval (x : ℝ) : z77.eval x = (1 - x / 6250) ^ 77 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step76_checked,
    z76_eval, z1_eval, ← pow_succ]
theorem z78_eval (x : ℝ) : z78.eval x = (1 - x / 6250) ^ 78 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step77_checked,
    z77_eval, z1_eval, ← pow_succ]
theorem z79_eval (x : ℝ) : z79.eval x = (1 - x / 6250) ^ 79 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step78_checked,
    z78_eval, z1_eval, ← pow_succ]
theorem z80_eval (x : ℝ) : z80.eval x = (1 - x / 6250) ^ 80 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step79_checked,
    z79_eval, z1_eval, ← pow_succ]
theorem z81_eval (x : ℝ) : z81.eval x = (1 - x / 6250) ^ 81 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step80_checked,
    z80_eval, z1_eval, ← pow_succ]
theorem z82_eval (x : ℝ) : z82.eval x = (1 - x / 6250) ^ 82 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step81_checked,
    z81_eval, z1_eval, ← pow_succ]
theorem z83_eval (x : ℝ) : z83.eval x = (1 - x / 6250) ^ 83 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step82_checked,
    z82_eval, z1_eval, ← pow_succ]
theorem z84_eval (x : ℝ) : z84.eval x = (1 - x / 6250) ^ 84 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step83_checked,
    z83_eval, z1_eval, ← pow_succ]
theorem z85_eval (x : ℝ) : z85.eval x = (1 - x / 6250) ^ 85 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step84_checked,
    z84_eval, z1_eval, ← pow_succ]
theorem z86_eval (x : ℝ) : z86.eval x = (1 - x / 6250) ^ 86 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step85_checked,
    z85_eval, z1_eval, ← pow_succ]
theorem z87_eval (x : ℝ) : z87.eval x = (1 - x / 6250) ^ 87 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step86_checked,
    z86_eval, z1_eval, ← pow_succ]
theorem z88_eval (x : ℝ) : z88.eval x = (1 - x / 6250) ^ 88 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step87_checked,
    z87_eval, z1_eval, ← pow_succ]
theorem z89_eval (x : ℝ) : z89.eval x = (1 - x / 6250) ^ 89 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step88_checked,
    z88_eval, z1_eval, ← pow_succ]
theorem z90_eval (x : ℝ) : z90.eval x = (1 - x / 6250) ^ 90 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step89_checked,
    z89_eval, z1_eval, ← pow_succ]
theorem z91_eval (x : ℝ) : z91.eval x = (1 - x / 6250) ^ 91 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step90_checked,
    z90_eval, z1_eval, ← pow_succ]
theorem z92_eval (x : ℝ) : z92.eval x = (1 - x / 6250) ^ 92 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step91_checked,
    z91_eval, z1_eval, ← pow_succ]
theorem z93_eval (x : ℝ) : z93.eval x = (1 - x / 6250) ^ 93 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step92_checked,
    z92_eval, z1_eval, ← pow_succ]
theorem z94_eval (x : ℝ) : z94.eval x = (1 - x / 6250) ^ 94 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step93_checked,
    z93_eval, z1_eval, ← pow_succ]
theorem z95_eval (x : ℝ) : z95.eval x = (1 - x / 6250) ^ 95 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step94_checked,
    z94_eval, z1_eval, ← pow_succ]
theorem z96_eval (x : ℝ) : z96.eval x = (1 - x / 6250) ^ 96 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step95_checked,
    z95_eval, z1_eval, ← pow_succ]
theorem z97_eval (x : ℝ) : z97.eval x = (1 - x / 6250) ^ 97 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step96_checked,
    z96_eval, z1_eval, ← pow_succ]
theorem z98_eval (x : ℝ) : z98.eval x = (1 - x / 6250) ^ 98 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step97_checked,
    z97_eval, z1_eval, ← pow_succ]
theorem z99_eval (x : ℝ) : z99.eval x = (1 - x / 6250) ^ 99 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step98_checked,
    z98_eval, z1_eval, ← pow_succ]
theorem z100_eval (x : ℝ) : z100.eval x = (1 - x / 6250) ^ 100 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step99_checked,
    z99_eval, z1_eval, ← pow_succ]
theorem z101_eval (x : ℝ) : z101.eval x = (1 - x / 6250) ^ 101 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step100_checked,
    z100_eval, z1_eval, ← pow_succ]
theorem z102_eval (x : ℝ) : z102.eval x = (1 - x / 6250) ^ 102 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step101_checked,
    z101_eval, z1_eval, ← pow_succ]
theorem z103_eval (x : ℝ) : z103.eval x = (1 - x / 6250) ^ 103 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step102_checked,
    z102_eval, z1_eval, ← pow_succ]
theorem z104_eval (x : ℝ) : z104.eval x = (1 - x / 6250) ^ 104 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step103_checked,
    z103_eval, z1_eval, ← pow_succ]
theorem z105_eval (x : ℝ) : z105.eval x = (1 - x / 6250) ^ 105 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step104_checked,
    z104_eval, z1_eval, ← pow_succ]
theorem z106_eval (x : ℝ) : z106.eval x = (1 - x / 6250) ^ 106 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step105_checked,
    z105_eval, z1_eval, ← pow_succ]
theorem z107_eval (x : ℝ) : z107.eval x = (1 - x / 6250) ^ 107 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step106_checked,
    z106_eval, z1_eval, ← pow_succ]
theorem z108_eval (x : ℝ) : z108.eval x = (1 - x / 6250) ^ 108 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step107_checked,
    z107_eval, z1_eval, ← pow_succ]
theorem z109_eval (x : ℝ) : z109.eval x = (1 - x / 6250) ^ 109 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step108_checked,
    z108_eval, z1_eval, ← pow_succ]
theorem z110_eval (x : ℝ) : z110.eval x = (1 - x / 6250) ^ 110 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step109_checked,
    z109_eval, z1_eval, ← pow_succ]
theorem z111_eval (x : ℝ) : z111.eval x = (1 - x / 6250) ^ 111 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step110_checked,
    z110_eval, z1_eval, ← pow_succ]
theorem z112_eval (x : ℝ) : z112.eval x = (1 - x / 6250) ^ 112 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step111_checked,
    z111_eval, z1_eval, ← pow_succ]
theorem z113_eval (x : ℝ) : z113.eval x = (1 - x / 6250) ^ 113 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step112_checked,
    z112_eval, z1_eval, ← pow_succ]
theorem z114_eval (x : ℝ) : z114.eval x = (1 - x / 6250) ^ 114 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step113_checked,
    z113_eval, z1_eval, ← pow_succ]
theorem z115_eval (x : ℝ) : z115.eval x = (1 - x / 6250) ^ 115 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step114_checked,
    z114_eval, z1_eval, ← pow_succ]
theorem z116_eval (x : ℝ) : z116.eval x = (1 - x / 6250) ^ 116 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step115_checked,
    z115_eval, z1_eval, ← pow_succ]
theorem z117_eval (x : ℝ) : z117.eval x = (1 - x / 6250) ^ 117 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step116_checked,
    z116_eval, z1_eval, ← pow_succ]
theorem z118_eval (x : ℝ) : z118.eval x = (1 - x / 6250) ^ 118 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step117_checked,
    z117_eval, z1_eval, ← pow_succ]
theorem z119_eval (x : ℝ) : z119.eval x = (1 - x / 6250) ^ 119 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step118_checked,
    z118_eval, z1_eval, ← pow_succ]
theorem z120_eval (x : ℝ) : z120.eval x = (1 - x / 6250) ^ 120 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step119_checked,
    z119_eval, z1_eval, ← pow_succ]
theorem z121_eval (x : ℝ) : z121.eval x = (1 - x / 6250) ^ 121 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step120_checked,
    z120_eval, z1_eval, ← pow_succ]
theorem z122_eval (x : ℝ) : z122.eval x = (1 - x / 6250) ^ 122 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step121_checked,
    z121_eval, z1_eval, ← pow_succ]
theorem z123_eval (x : ℝ) : z123.eval x = (1 - x / 6250) ^ 123 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step122_checked,
    z122_eval, z1_eval, ← pow_succ]
theorem z124_eval (x : ℝ) : z124.eval x = (1 - x / 6250) ^ 124 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step123_checked,
    z123_eval, z1_eval, ← pow_succ]
theorem z125_eval (x : ℝ) : z125.eval x = (1 - x / 6250) ^ 125 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step124_checked,
    z124_eval, z1_eval, ← pow_succ]
theorem z126_eval (x : ℝ) : z126.eval x = (1 - x / 6250) ^ 126 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step125_checked,
    z125_eval, z1_eval, ← pow_succ]
theorem z127_eval (x : ℝ) : z127.eval x = (1 - x / 6250) ^ 127 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step126_checked,
    z126_eval, z1_eval, ← pow_succ]
theorem z128_eval (x : ℝ) : z128.eval x = (1 - x / 6250) ^ 128 := by
  rw [RatPoly.checkMul_sound _ _ _ z_step127_checked,
    z127_eval, z1_eval, ← pow_succ]
theorem zPow_eval (n : ℕ) (x : ℝ) : (zPow n).eval x = (1 - x / 6250) ^ n := by
  by_cases h : n < 129
  · interval_cases n
    · exact z0_eval x
    · rw [pow_one]; exact z1_eval x
    · exact z2_eval x
    · exact z3_eval x
    · exact z4_eval x
    · exact z5_eval x
    · exact z6_eval x
    · exact z7_eval x
    · exact z8_eval x
    · exact z9_eval x
    · exact z10_eval x
    · exact z11_eval x
    · exact z12_eval x
    · exact z13_eval x
    · exact z14_eval x
    · exact z15_eval x
    · exact z16_eval x
    · exact z17_eval x
    · exact z18_eval x
    · exact z19_eval x
    · exact z20_eval x
    · exact z21_eval x
    · exact z22_eval x
    · exact z23_eval x
    · exact z24_eval x
    · exact z25_eval x
    · exact z26_eval x
    · exact z27_eval x
    · exact z28_eval x
    · exact z29_eval x
    · exact z30_eval x
    · exact z31_eval x
    · exact z32_eval x
    · exact z33_eval x
    · exact z34_eval x
    · exact z35_eval x
    · exact z36_eval x
    · exact z37_eval x
    · exact z38_eval x
    · exact z39_eval x
    · exact z40_eval x
    · exact z41_eval x
    · exact z42_eval x
    · exact z43_eval x
    · exact z44_eval x
    · exact z45_eval x
    · exact z46_eval x
    · exact z47_eval x
    · exact z48_eval x
    · exact z49_eval x
    · exact z50_eval x
    · exact z51_eval x
    · exact z52_eval x
    · exact z53_eval x
    · exact z54_eval x
    · exact z55_eval x
    · exact z56_eval x
    · exact z57_eval x
    · exact z58_eval x
    · exact z59_eval x
    · exact z60_eval x
    · exact z61_eval x
    · exact z62_eval x
    · exact z63_eval x
    · exact z64_eval x
    · exact z65_eval x
    · exact z66_eval x
    · exact z67_eval x
    · exact z68_eval x
    · exact z69_eval x
    · exact z70_eval x
    · exact z71_eval x
    · exact z72_eval x
    · exact z73_eval x
    · exact z74_eval x
    · exact z75_eval x
    · exact z76_eval x
    · exact z77_eval x
    · exact z78_eval x
    · exact z79_eval x
    · exact z80_eval x
    · exact z81_eval x
    · exact z82_eval x
    · exact z83_eval x
    · exact z84_eval x
    · exact z85_eval x
    · exact z86_eval x
    · exact z87_eval x
    · exact z88_eval x
    · exact z89_eval x
    · exact z90_eval x
    · exact z91_eval x
    · exact z92_eval x
    · exact z93_eval x
    · exact z94_eval x
    · exact z95_eval x
    · exact z96_eval x
    · exact z97_eval x
    · exact z98_eval x
    · exact z99_eval x
    · exact z100_eval x
    · exact z101_eval x
    · exact z102_eval x
    · exact z103_eval x
    · exact z104_eval x
    · exact z105_eval x
    · exact z106_eval x
    · exact z107_eval x
    · exact z108_eval x
    · exact z109_eval x
    · exact z110_eval x
    · exact z111_eval x
    · exact z112_eval x
    · exact z113_eval x
    · exact z114_eval x
    · exact z115_eval x
    · exact z116_eval x
    · exact z117_eval x
    · exact z118_eval x
    · exact z119_eval x
    · exact z120_eval x
    · exact z121_eval x
    · exact z122_eval x
    · exact z123_eval x
    · exact z124_eval x
    · exact z125_eval x
    · exact z126_eval x
    · exact z127_eval x
    · exact z128_eval x
  · rw [zPow, List.getD_eq_default _ _ (by simpa only [zPowers, List.length_cons, List.length_nil] using Nat.le_of_not_gt h)]
    rw [RatPoly.eval_pow, z1_eval]
theorem zPow_den_pos (n : ℕ) : 0 < (zPow n).den := by
  by_cases h : n < 129
  · interval_cases n <;> decide
  · rw [zPow, List.getD_eq_default _ _ (by simpa only [zPowers, List.length_cons, List.length_nil] using Nat.le_of_not_gt h)]
    exact RatPoly.pow_den_pos _ (by decide) n
theorem beta0_eval (x : ℝ) : beta0.eval x = (x / 12500) ^ 0 := by
  norm_num [beta0, RatPoly.eval, listEval]
theorem beta1_eval (x : ℝ) : beta1.eval x = (x / 12500) := by
  simp only [beta1, RatPoly.eval, listEval]
  push_cast
  ring
theorem beta2_eval (x : ℝ) : beta2.eval x = (x / 12500) ^ 2 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step1_checked,
    beta1_eval, pow_two]
theorem beta3_eval (x : ℝ) : beta3.eval x = (x / 12500) ^ 3 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step2_checked,
    beta2_eval, beta1_eval, ← pow_succ]
theorem beta4_eval (x : ℝ) : beta4.eval x = (x / 12500) ^ 4 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step3_checked,
    beta3_eval, beta1_eval, ← pow_succ]
theorem beta5_eval (x : ℝ) : beta5.eval x = (x / 12500) ^ 5 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step4_checked,
    beta4_eval, beta1_eval, ← pow_succ]
theorem beta6_eval (x : ℝ) : beta6.eval x = (x / 12500) ^ 6 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step5_checked,
    beta5_eval, beta1_eval, ← pow_succ]
theorem beta7_eval (x : ℝ) : beta7.eval x = (x / 12500) ^ 7 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step6_checked,
    beta6_eval, beta1_eval, ← pow_succ]
theorem beta8_eval (x : ℝ) : beta8.eval x = (x / 12500) ^ 8 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step7_checked,
    beta7_eval, beta1_eval, ← pow_succ]
theorem beta9_eval (x : ℝ) : beta9.eval x = (x / 12500) ^ 9 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step8_checked,
    beta8_eval, beta1_eval, ← pow_succ]
theorem beta10_eval (x : ℝ) : beta10.eval x = (x / 12500) ^ 10 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step9_checked,
    beta9_eval, beta1_eval, ← pow_succ]
theorem beta11_eval (x : ℝ) : beta11.eval x = (x / 12500) ^ 11 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step10_checked,
    beta10_eval, beta1_eval, ← pow_succ]
theorem beta12_eval (x : ℝ) : beta12.eval x = (x / 12500) ^ 12 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step11_checked,
    beta11_eval, beta1_eval, ← pow_succ]
theorem beta13_eval (x : ℝ) : beta13.eval x = (x / 12500) ^ 13 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step12_checked,
    beta12_eval, beta1_eval, ← pow_succ]
theorem beta14_eval (x : ℝ) : beta14.eval x = (x / 12500) ^ 14 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step13_checked,
    beta13_eval, beta1_eval, ← pow_succ]
theorem beta15_eval (x : ℝ) : beta15.eval x = (x / 12500) ^ 15 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step14_checked,
    beta14_eval, beta1_eval, ← pow_succ]
theorem beta16_eval (x : ℝ) : beta16.eval x = (x / 12500) ^ 16 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step15_checked,
    beta15_eval, beta1_eval, ← pow_succ]
theorem beta17_eval (x : ℝ) : beta17.eval x = (x / 12500) ^ 17 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step16_checked,
    beta16_eval, beta1_eval, ← pow_succ]
theorem beta18_eval (x : ℝ) : beta18.eval x = (x / 12500) ^ 18 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step17_checked,
    beta17_eval, beta1_eval, ← pow_succ]
theorem beta19_eval (x : ℝ) : beta19.eval x = (x / 12500) ^ 19 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step18_checked,
    beta18_eval, beta1_eval, ← pow_succ]
theorem beta20_eval (x : ℝ) : beta20.eval x = (x / 12500) ^ 20 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step19_checked,
    beta19_eval, beta1_eval, ← pow_succ]
theorem beta21_eval (x : ℝ) : beta21.eval x = (x / 12500) ^ 21 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step20_checked,
    beta20_eval, beta1_eval, ← pow_succ]
theorem beta22_eval (x : ℝ) : beta22.eval x = (x / 12500) ^ 22 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step21_checked,
    beta21_eval, beta1_eval, ← pow_succ]
theorem beta23_eval (x : ℝ) : beta23.eval x = (x / 12500) ^ 23 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step22_checked,
    beta22_eval, beta1_eval, ← pow_succ]
theorem beta24_eval (x : ℝ) : beta24.eval x = (x / 12500) ^ 24 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step23_checked,
    beta23_eval, beta1_eval, ← pow_succ]
theorem beta25_eval (x : ℝ) : beta25.eval x = (x / 12500) ^ 25 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step24_checked,
    beta24_eval, beta1_eval, ← pow_succ]
theorem beta26_eval (x : ℝ) : beta26.eval x = (x / 12500) ^ 26 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step25_checked,
    beta25_eval, beta1_eval, ← pow_succ]
theorem beta27_eval (x : ℝ) : beta27.eval x = (x / 12500) ^ 27 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step26_checked,
    beta26_eval, beta1_eval, ← pow_succ]
theorem beta28_eval (x : ℝ) : beta28.eval x = (x / 12500) ^ 28 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step27_checked,
    beta27_eval, beta1_eval, ← pow_succ]
theorem beta29_eval (x : ℝ) : beta29.eval x = (x / 12500) ^ 29 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step28_checked,
    beta28_eval, beta1_eval, ← pow_succ]
theorem beta30_eval (x : ℝ) : beta30.eval x = (x / 12500) ^ 30 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step29_checked,
    beta29_eval, beta1_eval, ← pow_succ]
theorem beta31_eval (x : ℝ) : beta31.eval x = (x / 12500) ^ 31 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step30_checked,
    beta30_eval, beta1_eval, ← pow_succ]
theorem beta32_eval (x : ℝ) : beta32.eval x = (x / 12500) ^ 32 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step31_checked,
    beta31_eval, beta1_eval, ← pow_succ]
theorem beta33_eval (x : ℝ) : beta33.eval x = (x / 12500) ^ 33 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step32_checked,
    beta32_eval, beta1_eval, ← pow_succ]
theorem beta34_eval (x : ℝ) : beta34.eval x = (x / 12500) ^ 34 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step33_checked,
    beta33_eval, beta1_eval, ← pow_succ]
theorem beta35_eval (x : ℝ) : beta35.eval x = (x / 12500) ^ 35 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step34_checked,
    beta34_eval, beta1_eval, ← pow_succ]
theorem beta36_eval (x : ℝ) : beta36.eval x = (x / 12500) ^ 36 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step35_checked,
    beta35_eval, beta1_eval, ← pow_succ]
theorem beta37_eval (x : ℝ) : beta37.eval x = (x / 12500) ^ 37 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step36_checked,
    beta36_eval, beta1_eval, ← pow_succ]
theorem beta38_eval (x : ℝ) : beta38.eval x = (x / 12500) ^ 38 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step37_checked,
    beta37_eval, beta1_eval, ← pow_succ]
theorem beta39_eval (x : ℝ) : beta39.eval x = (x / 12500) ^ 39 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step38_checked,
    beta38_eval, beta1_eval, ← pow_succ]
theorem beta40_eval (x : ℝ) : beta40.eval x = (x / 12500) ^ 40 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step39_checked,
    beta39_eval, beta1_eval, ← pow_succ]
theorem beta41_eval (x : ℝ) : beta41.eval x = (x / 12500) ^ 41 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step40_checked,
    beta40_eval, beta1_eval, ← pow_succ]
theorem beta42_eval (x : ℝ) : beta42.eval x = (x / 12500) ^ 42 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step41_checked,
    beta41_eval, beta1_eval, ← pow_succ]
theorem beta43_eval (x : ℝ) : beta43.eval x = (x / 12500) ^ 43 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step42_checked,
    beta42_eval, beta1_eval, ← pow_succ]
theorem beta44_eval (x : ℝ) : beta44.eval x = (x / 12500) ^ 44 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step43_checked,
    beta43_eval, beta1_eval, ← pow_succ]
theorem beta45_eval (x : ℝ) : beta45.eval x = (x / 12500) ^ 45 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step44_checked,
    beta44_eval, beta1_eval, ← pow_succ]
theorem beta46_eval (x : ℝ) : beta46.eval x = (x / 12500) ^ 46 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step45_checked,
    beta45_eval, beta1_eval, ← pow_succ]
theorem beta47_eval (x : ℝ) : beta47.eval x = (x / 12500) ^ 47 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step46_checked,
    beta46_eval, beta1_eval, ← pow_succ]
theorem beta48_eval (x : ℝ) : beta48.eval x = (x / 12500) ^ 48 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step47_checked,
    beta47_eval, beta1_eval, ← pow_succ]
theorem beta49_eval (x : ℝ) : beta49.eval x = (x / 12500) ^ 49 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step48_checked,
    beta48_eval, beta1_eval, ← pow_succ]
theorem beta50_eval (x : ℝ) : beta50.eval x = (x / 12500) ^ 50 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step49_checked,
    beta49_eval, beta1_eval, ← pow_succ]
theorem beta51_eval (x : ℝ) : beta51.eval x = (x / 12500) ^ 51 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step50_checked,
    beta50_eval, beta1_eval, ← pow_succ]
theorem beta52_eval (x : ℝ) : beta52.eval x = (x / 12500) ^ 52 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step51_checked,
    beta51_eval, beta1_eval, ← pow_succ]
theorem beta53_eval (x : ℝ) : beta53.eval x = (x / 12500) ^ 53 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step52_checked,
    beta52_eval, beta1_eval, ← pow_succ]
theorem beta54_eval (x : ℝ) : beta54.eval x = (x / 12500) ^ 54 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step53_checked,
    beta53_eval, beta1_eval, ← pow_succ]
theorem beta55_eval (x : ℝ) : beta55.eval x = (x / 12500) ^ 55 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step54_checked,
    beta54_eval, beta1_eval, ← pow_succ]
theorem beta56_eval (x : ℝ) : beta56.eval x = (x / 12500) ^ 56 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step55_checked,
    beta55_eval, beta1_eval, ← pow_succ]
theorem beta57_eval (x : ℝ) : beta57.eval x = (x / 12500) ^ 57 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step56_checked,
    beta56_eval, beta1_eval, ← pow_succ]
theorem beta58_eval (x : ℝ) : beta58.eval x = (x / 12500) ^ 58 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step57_checked,
    beta57_eval, beta1_eval, ← pow_succ]
theorem beta59_eval (x : ℝ) : beta59.eval x = (x / 12500) ^ 59 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step58_checked,
    beta58_eval, beta1_eval, ← pow_succ]
theorem beta60_eval (x : ℝ) : beta60.eval x = (x / 12500) ^ 60 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step59_checked,
    beta59_eval, beta1_eval, ← pow_succ]
theorem beta61_eval (x : ℝ) : beta61.eval x = (x / 12500) ^ 61 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step60_checked,
    beta60_eval, beta1_eval, ← pow_succ]
theorem beta62_eval (x : ℝ) : beta62.eval x = (x / 12500) ^ 62 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step61_checked,
    beta61_eval, beta1_eval, ← pow_succ]
theorem beta63_eval (x : ℝ) : beta63.eval x = (x / 12500) ^ 63 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step62_checked,
    beta62_eval, beta1_eval, ← pow_succ]
theorem beta64_eval (x : ℝ) : beta64.eval x = (x / 12500) ^ 64 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step63_checked,
    beta63_eval, beta1_eval, ← pow_succ]
theorem beta65_eval (x : ℝ) : beta65.eval x = (x / 12500) ^ 65 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step64_checked,
    beta64_eval, beta1_eval, ← pow_succ]
theorem beta66_eval (x : ℝ) : beta66.eval x = (x / 12500) ^ 66 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step65_checked,
    beta65_eval, beta1_eval, ← pow_succ]
theorem beta67_eval (x : ℝ) : beta67.eval x = (x / 12500) ^ 67 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step66_checked,
    beta66_eval, beta1_eval, ← pow_succ]
theorem beta68_eval (x : ℝ) : beta68.eval x = (x / 12500) ^ 68 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step67_checked,
    beta67_eval, beta1_eval, ← pow_succ]
theorem beta69_eval (x : ℝ) : beta69.eval x = (x / 12500) ^ 69 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step68_checked,
    beta68_eval, beta1_eval, ← pow_succ]
theorem beta70_eval (x : ℝ) : beta70.eval x = (x / 12500) ^ 70 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step69_checked,
    beta69_eval, beta1_eval, ← pow_succ]
theorem beta71_eval (x : ℝ) : beta71.eval x = (x / 12500) ^ 71 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step70_checked,
    beta70_eval, beta1_eval, ← pow_succ]
theorem beta72_eval (x : ℝ) : beta72.eval x = (x / 12500) ^ 72 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step71_checked,
    beta71_eval, beta1_eval, ← pow_succ]
theorem beta73_eval (x : ℝ) : beta73.eval x = (x / 12500) ^ 73 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step72_checked,
    beta72_eval, beta1_eval, ← pow_succ]
theorem beta74_eval (x : ℝ) : beta74.eval x = (x / 12500) ^ 74 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step73_checked,
    beta73_eval, beta1_eval, ← pow_succ]
theorem beta75_eval (x : ℝ) : beta75.eval x = (x / 12500) ^ 75 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step74_checked,
    beta74_eval, beta1_eval, ← pow_succ]
theorem beta76_eval (x : ℝ) : beta76.eval x = (x / 12500) ^ 76 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step75_checked,
    beta75_eval, beta1_eval, ← pow_succ]
theorem beta77_eval (x : ℝ) : beta77.eval x = (x / 12500) ^ 77 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step76_checked,
    beta76_eval, beta1_eval, ← pow_succ]
theorem beta78_eval (x : ℝ) : beta78.eval x = (x / 12500) ^ 78 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step77_checked,
    beta77_eval, beta1_eval, ← pow_succ]
theorem beta79_eval (x : ℝ) : beta79.eval x = (x / 12500) ^ 79 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step78_checked,
    beta78_eval, beta1_eval, ← pow_succ]
theorem beta80_eval (x : ℝ) : beta80.eval x = (x / 12500) ^ 80 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step79_checked,
    beta79_eval, beta1_eval, ← pow_succ]
theorem beta81_eval (x : ℝ) : beta81.eval x = (x / 12500) ^ 81 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step80_checked,
    beta80_eval, beta1_eval, ← pow_succ]
theorem beta82_eval (x : ℝ) : beta82.eval x = (x / 12500) ^ 82 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step81_checked,
    beta81_eval, beta1_eval, ← pow_succ]
theorem beta83_eval (x : ℝ) : beta83.eval x = (x / 12500) ^ 83 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step82_checked,
    beta82_eval, beta1_eval, ← pow_succ]
theorem beta84_eval (x : ℝ) : beta84.eval x = (x / 12500) ^ 84 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step83_checked,
    beta83_eval, beta1_eval, ← pow_succ]
theorem beta85_eval (x : ℝ) : beta85.eval x = (x / 12500) ^ 85 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step84_checked,
    beta84_eval, beta1_eval, ← pow_succ]
theorem beta86_eval (x : ℝ) : beta86.eval x = (x / 12500) ^ 86 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step85_checked,
    beta85_eval, beta1_eval, ← pow_succ]
theorem beta87_eval (x : ℝ) : beta87.eval x = (x / 12500) ^ 87 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step86_checked,
    beta86_eval, beta1_eval, ← pow_succ]
theorem beta88_eval (x : ℝ) : beta88.eval x = (x / 12500) ^ 88 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step87_checked,
    beta87_eval, beta1_eval, ← pow_succ]
theorem beta89_eval (x : ℝ) : beta89.eval x = (x / 12500) ^ 89 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step88_checked,
    beta88_eval, beta1_eval, ← pow_succ]
theorem beta90_eval (x : ℝ) : beta90.eval x = (x / 12500) ^ 90 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step89_checked,
    beta89_eval, beta1_eval, ← pow_succ]
theorem beta91_eval (x : ℝ) : beta91.eval x = (x / 12500) ^ 91 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step90_checked,
    beta90_eval, beta1_eval, ← pow_succ]
theorem beta92_eval (x : ℝ) : beta92.eval x = (x / 12500) ^ 92 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step91_checked,
    beta91_eval, beta1_eval, ← pow_succ]
theorem beta93_eval (x : ℝ) : beta93.eval x = (x / 12500) ^ 93 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step92_checked,
    beta92_eval, beta1_eval, ← pow_succ]
theorem beta94_eval (x : ℝ) : beta94.eval x = (x / 12500) ^ 94 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step93_checked,
    beta93_eval, beta1_eval, ← pow_succ]
theorem beta95_eval (x : ℝ) : beta95.eval x = (x / 12500) ^ 95 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step94_checked,
    beta94_eval, beta1_eval, ← pow_succ]
theorem beta96_eval (x : ℝ) : beta96.eval x = (x / 12500) ^ 96 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step95_checked,
    beta95_eval, beta1_eval, ← pow_succ]
theorem beta97_eval (x : ℝ) : beta97.eval x = (x / 12500) ^ 97 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step96_checked,
    beta96_eval, beta1_eval, ← pow_succ]
theorem beta98_eval (x : ℝ) : beta98.eval x = (x / 12500) ^ 98 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step97_checked,
    beta97_eval, beta1_eval, ← pow_succ]
theorem beta99_eval (x : ℝ) : beta99.eval x = (x / 12500) ^ 99 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step98_checked,
    beta98_eval, beta1_eval, ← pow_succ]
theorem beta100_eval (x : ℝ) : beta100.eval x = (x / 12500) ^ 100 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step99_checked,
    beta99_eval, beta1_eval, ← pow_succ]
theorem beta101_eval (x : ℝ) : beta101.eval x = (x / 12500) ^ 101 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step100_checked,
    beta100_eval, beta1_eval, ← pow_succ]
theorem beta102_eval (x : ℝ) : beta102.eval x = (x / 12500) ^ 102 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step101_checked,
    beta101_eval, beta1_eval, ← pow_succ]
theorem beta103_eval (x : ℝ) : beta103.eval x = (x / 12500) ^ 103 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step102_checked,
    beta102_eval, beta1_eval, ← pow_succ]
theorem beta104_eval (x : ℝ) : beta104.eval x = (x / 12500) ^ 104 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step103_checked,
    beta103_eval, beta1_eval, ← pow_succ]
theorem beta105_eval (x : ℝ) : beta105.eval x = (x / 12500) ^ 105 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step104_checked,
    beta104_eval, beta1_eval, ← pow_succ]
theorem beta106_eval (x : ℝ) : beta106.eval x = (x / 12500) ^ 106 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step105_checked,
    beta105_eval, beta1_eval, ← pow_succ]
theorem beta107_eval (x : ℝ) : beta107.eval x = (x / 12500) ^ 107 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step106_checked,
    beta106_eval, beta1_eval, ← pow_succ]
theorem beta108_eval (x : ℝ) : beta108.eval x = (x / 12500) ^ 108 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step107_checked,
    beta107_eval, beta1_eval, ← pow_succ]
theorem beta109_eval (x : ℝ) : beta109.eval x = (x / 12500) ^ 109 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step108_checked,
    beta108_eval, beta1_eval, ← pow_succ]
theorem beta110_eval (x : ℝ) : beta110.eval x = (x / 12500) ^ 110 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step109_checked,
    beta109_eval, beta1_eval, ← pow_succ]
theorem beta111_eval (x : ℝ) : beta111.eval x = (x / 12500) ^ 111 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step110_checked,
    beta110_eval, beta1_eval, ← pow_succ]
theorem beta112_eval (x : ℝ) : beta112.eval x = (x / 12500) ^ 112 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step111_checked,
    beta111_eval, beta1_eval, ← pow_succ]
theorem beta113_eval (x : ℝ) : beta113.eval x = (x / 12500) ^ 113 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step112_checked,
    beta112_eval, beta1_eval, ← pow_succ]
theorem beta114_eval (x : ℝ) : beta114.eval x = (x / 12500) ^ 114 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step113_checked,
    beta113_eval, beta1_eval, ← pow_succ]
theorem beta115_eval (x : ℝ) : beta115.eval x = (x / 12500) ^ 115 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step114_checked,
    beta114_eval, beta1_eval, ← pow_succ]
theorem beta116_eval (x : ℝ) : beta116.eval x = (x / 12500) ^ 116 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step115_checked,
    beta115_eval, beta1_eval, ← pow_succ]
theorem beta117_eval (x : ℝ) : beta117.eval x = (x / 12500) ^ 117 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step116_checked,
    beta116_eval, beta1_eval, ← pow_succ]
theorem beta118_eval (x : ℝ) : beta118.eval x = (x / 12500) ^ 118 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step117_checked,
    beta117_eval, beta1_eval, ← pow_succ]
theorem beta119_eval (x : ℝ) : beta119.eval x = (x / 12500) ^ 119 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step118_checked,
    beta118_eval, beta1_eval, ← pow_succ]
theorem beta120_eval (x : ℝ) : beta120.eval x = (x / 12500) ^ 120 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step119_checked,
    beta119_eval, beta1_eval, ← pow_succ]
theorem beta121_eval (x : ℝ) : beta121.eval x = (x / 12500) ^ 121 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step120_checked,
    beta120_eval, beta1_eval, ← pow_succ]
theorem beta122_eval (x : ℝ) : beta122.eval x = (x / 12500) ^ 122 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step121_checked,
    beta121_eval, beta1_eval, ← pow_succ]
theorem beta123_eval (x : ℝ) : beta123.eval x = (x / 12500) ^ 123 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step122_checked,
    beta122_eval, beta1_eval, ← pow_succ]
theorem beta124_eval (x : ℝ) : beta124.eval x = (x / 12500) ^ 124 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step123_checked,
    beta123_eval, beta1_eval, ← pow_succ]
theorem beta125_eval (x : ℝ) : beta125.eval x = (x / 12500) ^ 125 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step124_checked,
    beta124_eval, beta1_eval, ← pow_succ]
theorem beta126_eval (x : ℝ) : beta126.eval x = (x / 12500) ^ 126 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step125_checked,
    beta125_eval, beta1_eval, ← pow_succ]
theorem beta127_eval (x : ℝ) : beta127.eval x = (x / 12500) ^ 127 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step126_checked,
    beta126_eval, beta1_eval, ← pow_succ]
theorem beta128_eval (x : ℝ) : beta128.eval x = (x / 12500) ^ 128 := by
  rw [RatPoly.checkMul_sound _ _ _ beta_step127_checked,
    beta127_eval, beta1_eval, ← pow_succ]
theorem betaPow_eval (n : ℕ) (x : ℝ) : (betaPow n).eval x = (x / 12500) ^ n := by
  by_cases h : n < 129
  · interval_cases n
    · exact beta0_eval x
    · rw [pow_one]; exact beta1_eval x
    · exact beta2_eval x
    · exact beta3_eval x
    · exact beta4_eval x
    · exact beta5_eval x
    · exact beta6_eval x
    · exact beta7_eval x
    · exact beta8_eval x
    · exact beta9_eval x
    · exact beta10_eval x
    · exact beta11_eval x
    · exact beta12_eval x
    · exact beta13_eval x
    · exact beta14_eval x
    · exact beta15_eval x
    · exact beta16_eval x
    · exact beta17_eval x
    · exact beta18_eval x
    · exact beta19_eval x
    · exact beta20_eval x
    · exact beta21_eval x
    · exact beta22_eval x
    · exact beta23_eval x
    · exact beta24_eval x
    · exact beta25_eval x
    · exact beta26_eval x
    · exact beta27_eval x
    · exact beta28_eval x
    · exact beta29_eval x
    · exact beta30_eval x
    · exact beta31_eval x
    · exact beta32_eval x
    · exact beta33_eval x
    · exact beta34_eval x
    · exact beta35_eval x
    · exact beta36_eval x
    · exact beta37_eval x
    · exact beta38_eval x
    · exact beta39_eval x
    · exact beta40_eval x
    · exact beta41_eval x
    · exact beta42_eval x
    · exact beta43_eval x
    · exact beta44_eval x
    · exact beta45_eval x
    · exact beta46_eval x
    · exact beta47_eval x
    · exact beta48_eval x
    · exact beta49_eval x
    · exact beta50_eval x
    · exact beta51_eval x
    · exact beta52_eval x
    · exact beta53_eval x
    · exact beta54_eval x
    · exact beta55_eval x
    · exact beta56_eval x
    · exact beta57_eval x
    · exact beta58_eval x
    · exact beta59_eval x
    · exact beta60_eval x
    · exact beta61_eval x
    · exact beta62_eval x
    · exact beta63_eval x
    · exact beta64_eval x
    · exact beta65_eval x
    · exact beta66_eval x
    · exact beta67_eval x
    · exact beta68_eval x
    · exact beta69_eval x
    · exact beta70_eval x
    · exact beta71_eval x
    · exact beta72_eval x
    · exact beta73_eval x
    · exact beta74_eval x
    · exact beta75_eval x
    · exact beta76_eval x
    · exact beta77_eval x
    · exact beta78_eval x
    · exact beta79_eval x
    · exact beta80_eval x
    · exact beta81_eval x
    · exact beta82_eval x
    · exact beta83_eval x
    · exact beta84_eval x
    · exact beta85_eval x
    · exact beta86_eval x
    · exact beta87_eval x
    · exact beta88_eval x
    · exact beta89_eval x
    · exact beta90_eval x
    · exact beta91_eval x
    · exact beta92_eval x
    · exact beta93_eval x
    · exact beta94_eval x
    · exact beta95_eval x
    · exact beta96_eval x
    · exact beta97_eval x
    · exact beta98_eval x
    · exact beta99_eval x
    · exact beta100_eval x
    · exact beta101_eval x
    · exact beta102_eval x
    · exact beta103_eval x
    · exact beta104_eval x
    · exact beta105_eval x
    · exact beta106_eval x
    · exact beta107_eval x
    · exact beta108_eval x
    · exact beta109_eval x
    · exact beta110_eval x
    · exact beta111_eval x
    · exact beta112_eval x
    · exact beta113_eval x
    · exact beta114_eval x
    · exact beta115_eval x
    · exact beta116_eval x
    · exact beta117_eval x
    · exact beta118_eval x
    · exact beta119_eval x
    · exact beta120_eval x
    · exact beta121_eval x
    · exact beta122_eval x
    · exact beta123_eval x
    · exact beta124_eval x
    · exact beta125_eval x
    · exact beta126_eval x
    · exact beta127_eval x
    · exact beta128_eval x
  · rw [betaPow, List.getD_eq_default _ _ (by simpa only [betaPowers, List.length_cons, List.length_nil] using Nat.le_of_not_gt h)]
    rw [RatPoly.eval_pow, beta1_eval]
theorem betaPow_den_pos (n : ℕ) : 0 < (betaPow n).den := by
  by_cases h : n < 129
  · interval_cases n <;> decide
  · rw [betaPow, List.getD_eq_default _ _ (by simpa only [betaPowers, List.length_cons, List.length_nil] using Nat.le_of_not_gt h)]
    exact RatPoly.pow_den_pos _ (by decide) n
theorem complement0_eval (x : ℝ) : complement0.eval x = (1 - x / 12500) ^ 0 := by
  norm_num [complement0, RatPoly.eval, listEval]
theorem complement1_eval (x : ℝ) : complement1.eval x = (1 - x / 12500) := by
  simp only [complement1, RatPoly.eval, listEval]
  push_cast
  ring
theorem complement2_eval (x : ℝ) : complement2.eval x = (1 - x / 12500) ^ 2 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step1_checked,
    complement1_eval, pow_two]
theorem complement3_eval (x : ℝ) : complement3.eval x = (1 - x / 12500) ^ 3 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step2_checked,
    complement2_eval, complement1_eval, ← pow_succ]
theorem complement4_eval (x : ℝ) : complement4.eval x = (1 - x / 12500) ^ 4 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step3_checked,
    complement3_eval, complement1_eval, ← pow_succ]
theorem complement5_eval (x : ℝ) : complement5.eval x = (1 - x / 12500) ^ 5 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step4_checked,
    complement4_eval, complement1_eval, ← pow_succ]
theorem complement6_eval (x : ℝ) : complement6.eval x = (1 - x / 12500) ^ 6 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step5_checked,
    complement5_eval, complement1_eval, ← pow_succ]
theorem complement7_eval (x : ℝ) : complement7.eval x = (1 - x / 12500) ^ 7 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step6_checked,
    complement6_eval, complement1_eval, ← pow_succ]
theorem complement8_eval (x : ℝ) : complement8.eval x = (1 - x / 12500) ^ 8 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step7_checked,
    complement7_eval, complement1_eval, ← pow_succ]
theorem complement9_eval (x : ℝ) : complement9.eval x = (1 - x / 12500) ^ 9 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step8_checked,
    complement8_eval, complement1_eval, ← pow_succ]
theorem complement10_eval (x : ℝ) : complement10.eval x = (1 - x / 12500) ^ 10 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step9_checked,
    complement9_eval, complement1_eval, ← pow_succ]
theorem complement11_eval (x : ℝ) : complement11.eval x = (1 - x / 12500) ^ 11 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step10_checked,
    complement10_eval, complement1_eval, ← pow_succ]
theorem complement12_eval (x : ℝ) : complement12.eval x = (1 - x / 12500) ^ 12 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step11_checked,
    complement11_eval, complement1_eval, ← pow_succ]
theorem complement13_eval (x : ℝ) : complement13.eval x = (1 - x / 12500) ^ 13 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step12_checked,
    complement12_eval, complement1_eval, ← pow_succ]
theorem complement14_eval (x : ℝ) : complement14.eval x = (1 - x / 12500) ^ 14 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step13_checked,
    complement13_eval, complement1_eval, ← pow_succ]
theorem complement15_eval (x : ℝ) : complement15.eval x = (1 - x / 12500) ^ 15 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step14_checked,
    complement14_eval, complement1_eval, ← pow_succ]
theorem complement16_eval (x : ℝ) : complement16.eval x = (1 - x / 12500) ^ 16 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step15_checked,
    complement15_eval, complement1_eval, ← pow_succ]
theorem complement17_eval (x : ℝ) : complement17.eval x = (1 - x / 12500) ^ 17 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step16_checked,
    complement16_eval, complement1_eval, ← pow_succ]
theorem complement18_eval (x : ℝ) : complement18.eval x = (1 - x / 12500) ^ 18 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step17_checked,
    complement17_eval, complement1_eval, ← pow_succ]
theorem complement19_eval (x : ℝ) : complement19.eval x = (1 - x / 12500) ^ 19 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step18_checked,
    complement18_eval, complement1_eval, ← pow_succ]
theorem complement20_eval (x : ℝ) : complement20.eval x = (1 - x / 12500) ^ 20 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step19_checked,
    complement19_eval, complement1_eval, ← pow_succ]
theorem complement21_eval (x : ℝ) : complement21.eval x = (1 - x / 12500) ^ 21 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step20_checked,
    complement20_eval, complement1_eval, ← pow_succ]
theorem complement22_eval (x : ℝ) : complement22.eval x = (1 - x / 12500) ^ 22 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step21_checked,
    complement21_eval, complement1_eval, ← pow_succ]
theorem complement23_eval (x : ℝ) : complement23.eval x = (1 - x / 12500) ^ 23 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step22_checked,
    complement22_eval, complement1_eval, ← pow_succ]
theorem complement24_eval (x : ℝ) : complement24.eval x = (1 - x / 12500) ^ 24 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step23_checked,
    complement23_eval, complement1_eval, ← pow_succ]
theorem complement25_eval (x : ℝ) : complement25.eval x = (1 - x / 12500) ^ 25 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step24_checked,
    complement24_eval, complement1_eval, ← pow_succ]
theorem complement26_eval (x : ℝ) : complement26.eval x = (1 - x / 12500) ^ 26 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step25_checked,
    complement25_eval, complement1_eval, ← pow_succ]
theorem complement27_eval (x : ℝ) : complement27.eval x = (1 - x / 12500) ^ 27 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step26_checked,
    complement26_eval, complement1_eval, ← pow_succ]
theorem complement28_eval (x : ℝ) : complement28.eval x = (1 - x / 12500) ^ 28 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step27_checked,
    complement27_eval, complement1_eval, ← pow_succ]
theorem complement29_eval (x : ℝ) : complement29.eval x = (1 - x / 12500) ^ 29 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step28_checked,
    complement28_eval, complement1_eval, ← pow_succ]
theorem complement30_eval (x : ℝ) : complement30.eval x = (1 - x / 12500) ^ 30 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step29_checked,
    complement29_eval, complement1_eval, ← pow_succ]
theorem complement31_eval (x : ℝ) : complement31.eval x = (1 - x / 12500) ^ 31 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step30_checked,
    complement30_eval, complement1_eval, ← pow_succ]
theorem complement32_eval (x : ℝ) : complement32.eval x = (1 - x / 12500) ^ 32 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step31_checked,
    complement31_eval, complement1_eval, ← pow_succ]
theorem complement33_eval (x : ℝ) : complement33.eval x = (1 - x / 12500) ^ 33 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step32_checked,
    complement32_eval, complement1_eval, ← pow_succ]
theorem complement34_eval (x : ℝ) : complement34.eval x = (1 - x / 12500) ^ 34 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step33_checked,
    complement33_eval, complement1_eval, ← pow_succ]
theorem complement35_eval (x : ℝ) : complement35.eval x = (1 - x / 12500) ^ 35 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step34_checked,
    complement34_eval, complement1_eval, ← pow_succ]
theorem complement36_eval (x : ℝ) : complement36.eval x = (1 - x / 12500) ^ 36 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step35_checked,
    complement35_eval, complement1_eval, ← pow_succ]
theorem complement37_eval (x : ℝ) : complement37.eval x = (1 - x / 12500) ^ 37 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step36_checked,
    complement36_eval, complement1_eval, ← pow_succ]
theorem complement38_eval (x : ℝ) : complement38.eval x = (1 - x / 12500) ^ 38 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step37_checked,
    complement37_eval, complement1_eval, ← pow_succ]
theorem complement39_eval (x : ℝ) : complement39.eval x = (1 - x / 12500) ^ 39 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step38_checked,
    complement38_eval, complement1_eval, ← pow_succ]
theorem complement40_eval (x : ℝ) : complement40.eval x = (1 - x / 12500) ^ 40 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step39_checked,
    complement39_eval, complement1_eval, ← pow_succ]
theorem complement41_eval (x : ℝ) : complement41.eval x = (1 - x / 12500) ^ 41 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step40_checked,
    complement40_eval, complement1_eval, ← pow_succ]
theorem complement42_eval (x : ℝ) : complement42.eval x = (1 - x / 12500) ^ 42 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step41_checked,
    complement41_eval, complement1_eval, ← pow_succ]
theorem complement43_eval (x : ℝ) : complement43.eval x = (1 - x / 12500) ^ 43 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step42_checked,
    complement42_eval, complement1_eval, ← pow_succ]
theorem complement44_eval (x : ℝ) : complement44.eval x = (1 - x / 12500) ^ 44 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step43_checked,
    complement43_eval, complement1_eval, ← pow_succ]
theorem complement45_eval (x : ℝ) : complement45.eval x = (1 - x / 12500) ^ 45 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step44_checked,
    complement44_eval, complement1_eval, ← pow_succ]
theorem complement46_eval (x : ℝ) : complement46.eval x = (1 - x / 12500) ^ 46 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step45_checked,
    complement45_eval, complement1_eval, ← pow_succ]
theorem complement47_eval (x : ℝ) : complement47.eval x = (1 - x / 12500) ^ 47 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step46_checked,
    complement46_eval, complement1_eval, ← pow_succ]
theorem complement48_eval (x : ℝ) : complement48.eval x = (1 - x / 12500) ^ 48 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step47_checked,
    complement47_eval, complement1_eval, ← pow_succ]
theorem complement49_eval (x : ℝ) : complement49.eval x = (1 - x / 12500) ^ 49 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step48_checked,
    complement48_eval, complement1_eval, ← pow_succ]
theorem complement50_eval (x : ℝ) : complement50.eval x = (1 - x / 12500) ^ 50 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step49_checked,
    complement49_eval, complement1_eval, ← pow_succ]
theorem complement51_eval (x : ℝ) : complement51.eval x = (1 - x / 12500) ^ 51 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step50_checked,
    complement50_eval, complement1_eval, ← pow_succ]
theorem complement52_eval (x : ℝ) : complement52.eval x = (1 - x / 12500) ^ 52 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step51_checked,
    complement51_eval, complement1_eval, ← pow_succ]
theorem complement53_eval (x : ℝ) : complement53.eval x = (1 - x / 12500) ^ 53 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step52_checked,
    complement52_eval, complement1_eval, ← pow_succ]
theorem complement54_eval (x : ℝ) : complement54.eval x = (1 - x / 12500) ^ 54 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step53_checked,
    complement53_eval, complement1_eval, ← pow_succ]
theorem complement55_eval (x : ℝ) : complement55.eval x = (1 - x / 12500) ^ 55 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step54_checked,
    complement54_eval, complement1_eval, ← pow_succ]
theorem complement56_eval (x : ℝ) : complement56.eval x = (1 - x / 12500) ^ 56 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step55_checked,
    complement55_eval, complement1_eval, ← pow_succ]
theorem complement57_eval (x : ℝ) : complement57.eval x = (1 - x / 12500) ^ 57 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step56_checked,
    complement56_eval, complement1_eval, ← pow_succ]
theorem complement58_eval (x : ℝ) : complement58.eval x = (1 - x / 12500) ^ 58 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step57_checked,
    complement57_eval, complement1_eval, ← pow_succ]
theorem complement59_eval (x : ℝ) : complement59.eval x = (1 - x / 12500) ^ 59 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step58_checked,
    complement58_eval, complement1_eval, ← pow_succ]
theorem complement60_eval (x : ℝ) : complement60.eval x = (1 - x / 12500) ^ 60 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step59_checked,
    complement59_eval, complement1_eval, ← pow_succ]
theorem complement61_eval (x : ℝ) : complement61.eval x = (1 - x / 12500) ^ 61 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step60_checked,
    complement60_eval, complement1_eval, ← pow_succ]
theorem complement62_eval (x : ℝ) : complement62.eval x = (1 - x / 12500) ^ 62 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step61_checked,
    complement61_eval, complement1_eval, ← pow_succ]
theorem complement63_eval (x : ℝ) : complement63.eval x = (1 - x / 12500) ^ 63 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step62_checked,
    complement62_eval, complement1_eval, ← pow_succ]
theorem complement64_eval (x : ℝ) : complement64.eval x = (1 - x / 12500) ^ 64 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step63_checked,
    complement63_eval, complement1_eval, ← pow_succ]
theorem complement65_eval (x : ℝ) : complement65.eval x = (1 - x / 12500) ^ 65 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step64_checked,
    complement64_eval, complement1_eval, ← pow_succ]
theorem complement66_eval (x : ℝ) : complement66.eval x = (1 - x / 12500) ^ 66 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step65_checked,
    complement65_eval, complement1_eval, ← pow_succ]
theorem complement67_eval (x : ℝ) : complement67.eval x = (1 - x / 12500) ^ 67 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step66_checked,
    complement66_eval, complement1_eval, ← pow_succ]
theorem complement68_eval (x : ℝ) : complement68.eval x = (1 - x / 12500) ^ 68 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step67_checked,
    complement67_eval, complement1_eval, ← pow_succ]
theorem complement69_eval (x : ℝ) : complement69.eval x = (1 - x / 12500) ^ 69 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step68_checked,
    complement68_eval, complement1_eval, ← pow_succ]
theorem complement70_eval (x : ℝ) : complement70.eval x = (1 - x / 12500) ^ 70 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step69_checked,
    complement69_eval, complement1_eval, ← pow_succ]
theorem complement71_eval (x : ℝ) : complement71.eval x = (1 - x / 12500) ^ 71 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step70_checked,
    complement70_eval, complement1_eval, ← pow_succ]
theorem complement72_eval (x : ℝ) : complement72.eval x = (1 - x / 12500) ^ 72 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step71_checked,
    complement71_eval, complement1_eval, ← pow_succ]
theorem complement73_eval (x : ℝ) : complement73.eval x = (1 - x / 12500) ^ 73 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step72_checked,
    complement72_eval, complement1_eval, ← pow_succ]
theorem complement74_eval (x : ℝ) : complement74.eval x = (1 - x / 12500) ^ 74 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step73_checked,
    complement73_eval, complement1_eval, ← pow_succ]
theorem complement75_eval (x : ℝ) : complement75.eval x = (1 - x / 12500) ^ 75 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step74_checked,
    complement74_eval, complement1_eval, ← pow_succ]
theorem complement76_eval (x : ℝ) : complement76.eval x = (1 - x / 12500) ^ 76 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step75_checked,
    complement75_eval, complement1_eval, ← pow_succ]
theorem complement77_eval (x : ℝ) : complement77.eval x = (1 - x / 12500) ^ 77 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step76_checked,
    complement76_eval, complement1_eval, ← pow_succ]
theorem complement78_eval (x : ℝ) : complement78.eval x = (1 - x / 12500) ^ 78 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step77_checked,
    complement77_eval, complement1_eval, ← pow_succ]
theorem complement79_eval (x : ℝ) : complement79.eval x = (1 - x / 12500) ^ 79 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step78_checked,
    complement78_eval, complement1_eval, ← pow_succ]
theorem complement80_eval (x : ℝ) : complement80.eval x = (1 - x / 12500) ^ 80 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step79_checked,
    complement79_eval, complement1_eval, ← pow_succ]
theorem complement81_eval (x : ℝ) : complement81.eval x = (1 - x / 12500) ^ 81 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step80_checked,
    complement80_eval, complement1_eval, ← pow_succ]
theorem complement82_eval (x : ℝ) : complement82.eval x = (1 - x / 12500) ^ 82 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step81_checked,
    complement81_eval, complement1_eval, ← pow_succ]
theorem complement83_eval (x : ℝ) : complement83.eval x = (1 - x / 12500) ^ 83 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step82_checked,
    complement82_eval, complement1_eval, ← pow_succ]
theorem complement84_eval (x : ℝ) : complement84.eval x = (1 - x / 12500) ^ 84 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step83_checked,
    complement83_eval, complement1_eval, ← pow_succ]
theorem complement85_eval (x : ℝ) : complement85.eval x = (1 - x / 12500) ^ 85 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step84_checked,
    complement84_eval, complement1_eval, ← pow_succ]
theorem complement86_eval (x : ℝ) : complement86.eval x = (1 - x / 12500) ^ 86 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step85_checked,
    complement85_eval, complement1_eval, ← pow_succ]
theorem complement87_eval (x : ℝ) : complement87.eval x = (1 - x / 12500) ^ 87 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step86_checked,
    complement86_eval, complement1_eval, ← pow_succ]
theorem complement88_eval (x : ℝ) : complement88.eval x = (1 - x / 12500) ^ 88 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step87_checked,
    complement87_eval, complement1_eval, ← pow_succ]
theorem complement89_eval (x : ℝ) : complement89.eval x = (1 - x / 12500) ^ 89 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step88_checked,
    complement88_eval, complement1_eval, ← pow_succ]
theorem complement90_eval (x : ℝ) : complement90.eval x = (1 - x / 12500) ^ 90 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step89_checked,
    complement89_eval, complement1_eval, ← pow_succ]
theorem complement91_eval (x : ℝ) : complement91.eval x = (1 - x / 12500) ^ 91 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step90_checked,
    complement90_eval, complement1_eval, ← pow_succ]
theorem complement92_eval (x : ℝ) : complement92.eval x = (1 - x / 12500) ^ 92 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step91_checked,
    complement91_eval, complement1_eval, ← pow_succ]
theorem complement93_eval (x : ℝ) : complement93.eval x = (1 - x / 12500) ^ 93 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step92_checked,
    complement92_eval, complement1_eval, ← pow_succ]
theorem complement94_eval (x : ℝ) : complement94.eval x = (1 - x / 12500) ^ 94 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step93_checked,
    complement93_eval, complement1_eval, ← pow_succ]
theorem complement95_eval (x : ℝ) : complement95.eval x = (1 - x / 12500) ^ 95 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step94_checked,
    complement94_eval, complement1_eval, ← pow_succ]
theorem complement96_eval (x : ℝ) : complement96.eval x = (1 - x / 12500) ^ 96 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step95_checked,
    complement95_eval, complement1_eval, ← pow_succ]
theorem complement97_eval (x : ℝ) : complement97.eval x = (1 - x / 12500) ^ 97 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step96_checked,
    complement96_eval, complement1_eval, ← pow_succ]
theorem complement98_eval (x : ℝ) : complement98.eval x = (1 - x / 12500) ^ 98 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step97_checked,
    complement97_eval, complement1_eval, ← pow_succ]
theorem complement99_eval (x : ℝ) : complement99.eval x = (1 - x / 12500) ^ 99 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step98_checked,
    complement98_eval, complement1_eval, ← pow_succ]
theorem complement100_eval (x : ℝ) : complement100.eval x = (1 - x / 12500) ^ 100 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step99_checked,
    complement99_eval, complement1_eval, ← pow_succ]
theorem complement101_eval (x : ℝ) : complement101.eval x = (1 - x / 12500) ^ 101 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step100_checked,
    complement100_eval, complement1_eval, ← pow_succ]
theorem complement102_eval (x : ℝ) : complement102.eval x = (1 - x / 12500) ^ 102 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step101_checked,
    complement101_eval, complement1_eval, ← pow_succ]
theorem complement103_eval (x : ℝ) : complement103.eval x = (1 - x / 12500) ^ 103 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step102_checked,
    complement102_eval, complement1_eval, ← pow_succ]
theorem complement104_eval (x : ℝ) : complement104.eval x = (1 - x / 12500) ^ 104 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step103_checked,
    complement103_eval, complement1_eval, ← pow_succ]
theorem complement105_eval (x : ℝ) : complement105.eval x = (1 - x / 12500) ^ 105 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step104_checked,
    complement104_eval, complement1_eval, ← pow_succ]
theorem complement106_eval (x : ℝ) : complement106.eval x = (1 - x / 12500) ^ 106 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step105_checked,
    complement105_eval, complement1_eval, ← pow_succ]
theorem complement107_eval (x : ℝ) : complement107.eval x = (1 - x / 12500) ^ 107 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step106_checked,
    complement106_eval, complement1_eval, ← pow_succ]
theorem complement108_eval (x : ℝ) : complement108.eval x = (1 - x / 12500) ^ 108 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step107_checked,
    complement107_eval, complement1_eval, ← pow_succ]
theorem complement109_eval (x : ℝ) : complement109.eval x = (1 - x / 12500) ^ 109 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step108_checked,
    complement108_eval, complement1_eval, ← pow_succ]
theorem complement110_eval (x : ℝ) : complement110.eval x = (1 - x / 12500) ^ 110 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step109_checked,
    complement109_eval, complement1_eval, ← pow_succ]
theorem complement111_eval (x : ℝ) : complement111.eval x = (1 - x / 12500) ^ 111 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step110_checked,
    complement110_eval, complement1_eval, ← pow_succ]
theorem complement112_eval (x : ℝ) : complement112.eval x = (1 - x / 12500) ^ 112 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step111_checked,
    complement111_eval, complement1_eval, ← pow_succ]
theorem complement113_eval (x : ℝ) : complement113.eval x = (1 - x / 12500) ^ 113 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step112_checked,
    complement112_eval, complement1_eval, ← pow_succ]
theorem complement114_eval (x : ℝ) : complement114.eval x = (1 - x / 12500) ^ 114 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step113_checked,
    complement113_eval, complement1_eval, ← pow_succ]
theorem complement115_eval (x : ℝ) : complement115.eval x = (1 - x / 12500) ^ 115 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step114_checked,
    complement114_eval, complement1_eval, ← pow_succ]
theorem complement116_eval (x : ℝ) : complement116.eval x = (1 - x / 12500) ^ 116 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step115_checked,
    complement115_eval, complement1_eval, ← pow_succ]
theorem complement117_eval (x : ℝ) : complement117.eval x = (1 - x / 12500) ^ 117 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step116_checked,
    complement116_eval, complement1_eval, ← pow_succ]
theorem complement118_eval (x : ℝ) : complement118.eval x = (1 - x / 12500) ^ 118 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step117_checked,
    complement117_eval, complement1_eval, ← pow_succ]
theorem complement119_eval (x : ℝ) : complement119.eval x = (1 - x / 12500) ^ 119 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step118_checked,
    complement118_eval, complement1_eval, ← pow_succ]
theorem complement120_eval (x : ℝ) : complement120.eval x = (1 - x / 12500) ^ 120 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step119_checked,
    complement119_eval, complement1_eval, ← pow_succ]
theorem complement121_eval (x : ℝ) : complement121.eval x = (1 - x / 12500) ^ 121 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step120_checked,
    complement120_eval, complement1_eval, ← pow_succ]
theorem complement122_eval (x : ℝ) : complement122.eval x = (1 - x / 12500) ^ 122 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step121_checked,
    complement121_eval, complement1_eval, ← pow_succ]
theorem complement123_eval (x : ℝ) : complement123.eval x = (1 - x / 12500) ^ 123 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step122_checked,
    complement122_eval, complement1_eval, ← pow_succ]
theorem complement124_eval (x : ℝ) : complement124.eval x = (1 - x / 12500) ^ 124 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step123_checked,
    complement123_eval, complement1_eval, ← pow_succ]
theorem complement125_eval (x : ℝ) : complement125.eval x = (1 - x / 12500) ^ 125 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step124_checked,
    complement124_eval, complement1_eval, ← pow_succ]
theorem complement126_eval (x : ℝ) : complement126.eval x = (1 - x / 12500) ^ 126 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step125_checked,
    complement125_eval, complement1_eval, ← pow_succ]
theorem complement127_eval (x : ℝ) : complement127.eval x = (1 - x / 12500) ^ 127 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step126_checked,
    complement126_eval, complement1_eval, ← pow_succ]
theorem complement128_eval (x : ℝ) : complement128.eval x = (1 - x / 12500) ^ 128 := by
  rw [RatPoly.checkMul_sound _ _ _ complement_step127_checked,
    complement127_eval, complement1_eval, ← pow_succ]
theorem complementPow_eval (n : ℕ) (x : ℝ) : (complementPow n).eval x = (1 - x / 12500) ^ n := by
  by_cases h : n < 129
  · interval_cases n
    · exact complement0_eval x
    · rw [pow_one]; exact complement1_eval x
    · exact complement2_eval x
    · exact complement3_eval x
    · exact complement4_eval x
    · exact complement5_eval x
    · exact complement6_eval x
    · exact complement7_eval x
    · exact complement8_eval x
    · exact complement9_eval x
    · exact complement10_eval x
    · exact complement11_eval x
    · exact complement12_eval x
    · exact complement13_eval x
    · exact complement14_eval x
    · exact complement15_eval x
    · exact complement16_eval x
    · exact complement17_eval x
    · exact complement18_eval x
    · exact complement19_eval x
    · exact complement20_eval x
    · exact complement21_eval x
    · exact complement22_eval x
    · exact complement23_eval x
    · exact complement24_eval x
    · exact complement25_eval x
    · exact complement26_eval x
    · exact complement27_eval x
    · exact complement28_eval x
    · exact complement29_eval x
    · exact complement30_eval x
    · exact complement31_eval x
    · exact complement32_eval x
    · exact complement33_eval x
    · exact complement34_eval x
    · exact complement35_eval x
    · exact complement36_eval x
    · exact complement37_eval x
    · exact complement38_eval x
    · exact complement39_eval x
    · exact complement40_eval x
    · exact complement41_eval x
    · exact complement42_eval x
    · exact complement43_eval x
    · exact complement44_eval x
    · exact complement45_eval x
    · exact complement46_eval x
    · exact complement47_eval x
    · exact complement48_eval x
    · exact complement49_eval x
    · exact complement50_eval x
    · exact complement51_eval x
    · exact complement52_eval x
    · exact complement53_eval x
    · exact complement54_eval x
    · exact complement55_eval x
    · exact complement56_eval x
    · exact complement57_eval x
    · exact complement58_eval x
    · exact complement59_eval x
    · exact complement60_eval x
    · exact complement61_eval x
    · exact complement62_eval x
    · exact complement63_eval x
    · exact complement64_eval x
    · exact complement65_eval x
    · exact complement66_eval x
    · exact complement67_eval x
    · exact complement68_eval x
    · exact complement69_eval x
    · exact complement70_eval x
    · exact complement71_eval x
    · exact complement72_eval x
    · exact complement73_eval x
    · exact complement74_eval x
    · exact complement75_eval x
    · exact complement76_eval x
    · exact complement77_eval x
    · exact complement78_eval x
    · exact complement79_eval x
    · exact complement80_eval x
    · exact complement81_eval x
    · exact complement82_eval x
    · exact complement83_eval x
    · exact complement84_eval x
    · exact complement85_eval x
    · exact complement86_eval x
    · exact complement87_eval x
    · exact complement88_eval x
    · exact complement89_eval x
    · exact complement90_eval x
    · exact complement91_eval x
    · exact complement92_eval x
    · exact complement93_eval x
    · exact complement94_eval x
    · exact complement95_eval x
    · exact complement96_eval x
    · exact complement97_eval x
    · exact complement98_eval x
    · exact complement99_eval x
    · exact complement100_eval x
    · exact complement101_eval x
    · exact complement102_eval x
    · exact complement103_eval x
    · exact complement104_eval x
    · exact complement105_eval x
    · exact complement106_eval x
    · exact complement107_eval x
    · exact complement108_eval x
    · exact complement109_eval x
    · exact complement110_eval x
    · exact complement111_eval x
    · exact complement112_eval x
    · exact complement113_eval x
    · exact complement114_eval x
    · exact complement115_eval x
    · exact complement116_eval x
    · exact complement117_eval x
    · exact complement118_eval x
    · exact complement119_eval x
    · exact complement120_eval x
    · exact complement121_eval x
    · exact complement122_eval x
    · exact complement123_eval x
    · exact complement124_eval x
    · exact complement125_eval x
    · exact complement126_eval x
    · exact complement127_eval x
    · exact complement128_eval x
  · rw [complementPow, List.getD_eq_default _ _ (by simpa only [complementPowers, List.length_cons, List.length_nil] using Nat.le_of_not_gt h)]
    rw [RatPoly.eval_pow, complement1_eval]
theorem complementPow_den_pos (n : ℕ) : 0 < (complementPow n).den := by
  by_cases h : n < 129
  · interval_cases n <;> decide
  · rw [complementPow, List.getD_eq_default _ _ (by simpa only [complementPowers, List.length_cons, List.length_nil] using Nat.le_of_not_gt h)]
    exact RatPoly.pow_den_pos _ (by decide) n
end Spin.Structured.SparsePowers


