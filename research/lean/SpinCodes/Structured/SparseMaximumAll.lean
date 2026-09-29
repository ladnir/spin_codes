import SpinCodes.Structured.SparseMaxima
import SpinCodes.Structured.SparseWeights
import SpinCodes.Structured.SparseDataChecks
noncomputable section
namespace Spin.Structured.SparsePolynomial
set_option maxRecDepth 100000
set_option maxHeartbeats 0
theorem maximum_all (x : ℝ) (h0 : 0 ≤ x) (h1 : x ≤ 1) (j : Fin 129) :
    Spin.Imt.Occupation.Sparse.maximumMoment j (1 - x / 6250) ≤
      (arbitrary j (Data.weight j)).eval x := by
  fin_cases j
  · exact maximum_0 x h0 h1
  · exact maximum_1 x h0 h1
  · exact maximum_2 x h0 h1
  · exact maximum_3 x h0 h1
  · exact maximum_4 x h0 h1
  · exact maximum_5 x h0 h1
  · exact maximum_6 x h0 h1
  · exact maximum_7 x h0 h1
  · exact maximum_8 x h0 h1
  · exact maximum_9 x h0 h1
  · exact maximum_10 x h0 h1
  · exact maximum_11 x h0 h1
  · exact maximum_12 x h0 h1
  · exact maximum_13 x h0 h1
  · exact maximum_14 x h0 h1
  · exact maximum_15 x h0 h1
  · exact maximum_16 x h0 h1
  · exact maximum_17 x h0 h1
  · exact maximum_18 x h0 h1
  · exact maximum_19 x h0 h1
  · exact maximum_20 x h0 h1
  · exact maximum_21 x h0 h1
  · exact maximum_22 x h0 h1
  · exact maximum_23 x h0 h1
  · exact maximum_24 x h0 h1
  · exact maximum_25 x h0 h1
  · exact maximum_26 x h0 h1
  · exact maximum_27 x h0 h1
  · exact maximum_28 x h0 h1
  · exact maximum_29 x h0 h1
  · exact maximum_30 x h0 h1
  · exact maximum_31 x h0 h1
  · exact maximum_32 x h0 h1
  · exact maximum_33 x h0 h1
  · exact maximum_34 x h0 h1
  · exact maximum_35 x h0 h1
  · exact maximum_36 x h0 h1
  · exact maximum_37 x h0 h1
  · exact maximum_38 x h0 h1
  · exact maximum_39 x h0 h1
  · exact maximum_40 x h0 h1
  · exact maximum_41 x h0 h1
  · exact maximum_42 x h0 h1
  · exact maximum_43 x h0 h1
  · exact maximum_44 x h0 h1
  · exact maximum_45 x h0 h1
  · exact maximum_46 x h0 h1
  · exact maximum_47 x h0 h1
  · exact maximum_48 x h0 h1
  · exact maximum_49 x h0 h1
  · exact maximum_50 x h0 h1
  · exact maximum_51 x h0 h1
  · exact maximum_52 x h0 h1
  · exact maximum_53 x h0 h1
  · exact maximum_54 x h0 h1
  · exact maximum_55 x h0 h1
  · exact maximum_56 x h0 h1
  · exact maximum_57 x h0 h1
  · exact maximum_58 x h0 h1
  · exact maximum_59 x h0 h1
  · exact maximum_60 x h0 h1
  · exact maximum_61 x h0 h1
  · exact maximum_62 x h0 h1
  · exact maximum_63 x h0 h1
  · exact maximum_64 x h0 h1
  · exact maximum_65 x h0 h1
  · exact maximum_66 x h0 h1
  · exact maximum_67 x h0 h1
  · exact maximum_68 x h0 h1
  · exact maximum_69 x h0 h1
  · exact maximum_70 x h0 h1
  · exact maximum_71 x h0 h1
  · exact maximum_72 x h0 h1
  · exact maximum_73 x h0 h1
  · exact maximum_74 x h0 h1
  · exact maximum_75 x h0 h1
  · exact maximum_76 x h0 h1
  · exact maximum_77 x h0 h1
  · exact maximum_78 x h0 h1
  · exact maximum_79 x h0 h1
  · exact maximum_80 x h0 h1
  · exact maximum_81 x h0 h1
  · exact maximum_82 x h0 h1
  · exact maximum_83 x h0 h1
  · exact maximum_84 x h0 h1
  · exact maximum_85 x h0 h1
  · exact maximum_86 x h0 h1
  · exact maximum_87 x h0 h1
  · exact maximum_88 x h0 h1
  · exact maximum_89 x h0 h1
  · exact maximum_90 x h0 h1
  · exact maximum_91 x h0 h1
  · exact maximum_92 x h0 h1
  · exact maximum_93 x h0 h1
  · exact maximum_94 x h0 h1
  · exact maximum_95 x h0 h1
  · exact maximum_96 x h0 h1
  · exact maximum_97 x h0 h1
  · exact maximum_98 x h0 h1
  · exact maximum_99 x h0 h1
  · exact maximum_100 x h0 h1
  · exact maximum_101 x h0 h1
  · exact maximum_102 x h0 h1
  · exact maximum_103 x h0 h1
  · exact maximum_104 x h0 h1
  · exact maximum_105 x h0 h1
  · exact maximum_106 x h0 h1
  · exact maximum_107 x h0 h1
  · exact maximum_108 x h0 h1
  · exact maximum_109 x h0 h1
  · exact maximum_110 x h0 h1
  · exact maximum_111 x h0 h1
  · exact maximum_112 x h0 h1
  · exact maximum_113 x h0 h1
  · exact maximum_114 x h0 h1
  · exact maximum_115 x h0 h1
  · exact maximum_116 x h0 h1
  · exact maximum_117 x h0 h1
  · exact maximum_118 x h0 h1
  · exact maximum_119 x h0 h1
  · exact maximum_120 x h0 h1
  · exact maximum_121 x h0 h1
  · exact maximum_122 x h0 h1
  · exact maximum_123 x h0 h1
  · exact maximum_124 x h0 h1
  · exact maximum_125 x h0 h1
  · exact maximum_126 x h0 h1
  · exact maximum_127 x h0 h1
  · exact maximum_128 x h0 h1
theorem lowMaximum_all (x : ℝ) (h0 : 0 ≤ x) (h1 : x ≤ 1) (j : Fin 129) :
    lowMaximum j (Data.weight j) (1 - x / 6250) ≤
      (selectedLowPolynomial j (Data.weight j)).eval x := by
  fin_cases j
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · exact lowMaximum_1 x h0 h1
  · exact lowMaximum_2 x h0 h1
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
  · change (0 : ℝ) ≤ zero.eval x
    norm_num [zero, RatPoly.eval_constant]
theorem data_all_valid (j : Fin 129) : checkData j (Data.weight j) = true := by
  fin_cases j
  · exact Data.weight0_valid
  · exact Data.weight1_valid
  · exact Data.weight2_valid
  · exact Data.weight3_valid
  · exact Data.weight4_valid
  · exact Data.weight5_valid
  · exact Data.weight6_valid
  · exact Data.weight7_valid
  · exact Data.weight8_valid
  · exact Data.weight9_valid
  · exact Data.weight10_valid
  · exact Data.weight11_valid
  · exact Data.weight12_valid
  · exact Data.weight13_valid
  · exact Data.weight14_valid
  · exact Data.weight15_valid
  · exact Data.weight16_valid
  · exact Data.weight17_valid
  · exact Data.weight18_valid
  · exact Data.weight19_valid
  · exact Data.weight20_valid
  · exact Data.weight21_valid
  · exact Data.weight22_valid
  · exact Data.weight23_valid
  · exact Data.weight24_valid
  · exact Data.weight25_valid
  · exact Data.weight26_valid
  · exact Data.weight27_valid
  · exact Data.weight28_valid
  · exact Data.weight29_valid
  · exact Data.weight30_valid
  · exact Data.weight31_valid
  · exact Data.weight32_valid
  · exact Data.weight33_valid
  · exact Data.weight34_valid
  · exact Data.weight35_valid
  · exact Data.weight36_valid
  · exact Data.weight37_valid
  · exact Data.weight38_valid
  · exact Data.weight39_valid
  · exact Data.weight40_valid
  · exact Data.weight41_valid
  · exact Data.weight42_valid
  · exact Data.weight43_valid
  · exact Data.weight44_valid
  · exact Data.weight45_valid
  · exact Data.weight46_valid
  · exact Data.weight47_valid
  · exact Data.weight48_valid
  · exact Data.weight49_valid
  · exact Data.weight50_valid
  · exact Data.weight51_valid
  · exact Data.weight52_valid
  · exact Data.weight53_valid
  · exact Data.weight54_valid
  · exact Data.weight55_valid
  · exact Data.weight56_valid
  · exact Data.weight57_valid
  · exact Data.weight58_valid
  · exact Data.weight59_valid
  · exact Data.weight60_valid
  · exact Data.weight61_valid
  · exact Data.weight62_valid
  · exact Data.weight63_valid
  · exact Data.weight64_valid
  · exact Data.weight65_valid
  · exact Data.weight66_valid
  · exact Data.weight67_valid
  · exact Data.weight68_valid
  · exact Data.weight69_valid
  · exact Data.weight70_valid
  · exact Data.weight71_valid
  · exact Data.weight72_valid
  · exact Data.weight73_valid
  · exact Data.weight74_valid
  · exact Data.weight75_valid
  · exact Data.weight76_valid
  · exact Data.weight77_valid
  · exact Data.weight78_valid
  · exact Data.weight79_valid
  · exact Data.weight80_valid
  · exact Data.weight81_valid
  · exact Data.weight82_valid
  · exact Data.weight83_valid
  · exact Data.weight84_valid
  · exact Data.weight85_valid
  · exact Data.weight86_valid
  · exact Data.weight87_valid
  · exact Data.weight88_valid
  · exact Data.weight89_valid
  · exact Data.weight90_valid
  · exact Data.weight91_valid
  · exact Data.weight92_valid
  · exact Data.weight93_valid
  · exact Data.weight94_valid
  · exact Data.weight95_valid
  · exact Data.weight96_valid
  · exact Data.weight97_valid
  · exact Data.weight98_valid
  · exact Data.weight99_valid
  · exact Data.weight100_valid
  · exact Data.weight101_valid
  · exact Data.weight102_valid
  · exact Data.weight103_valid
  · exact Data.weight104_valid
  · exact Data.weight105_valid
  · exact Data.weight106_valid
  · exact Data.weight107_valid
  · exact Data.weight108_valid
  · exact Data.weight109_valid
  · exact Data.weight110_valid
  · exact Data.weight111_valid
  · exact Data.weight112_valid
  · exact Data.weight113_valid
  · exact Data.weight114_valid
  · exact Data.weight115_valid
  · exact Data.weight116_valid
  · exact Data.weight117_valid
  · exact Data.weight118_valid
  · exact Data.weight119_valid
  · exact Data.weight120_valid
  · exact Data.weight121_valid
  · exact Data.weight122_valid
  · exact Data.weight123_valid
  · exact Data.weight124_valid
  · exact Data.weight125_valid
  · exact Data.weight126_valid
  · exact Data.weight127_valid
  · exact Data.weight128_valid
end Spin.Structured.SparsePolynomial
namespace Spin.Imt.Occupation.Sparse
def numericalMatrix (β z : ℝ) : Transfer 5 :=
  matrix 128 count 524287 β (fun j => row j (Spin.Structured.SparsePolynomial.Data.weight j) z)
end Spin.Imt.Occupation.Sparse
