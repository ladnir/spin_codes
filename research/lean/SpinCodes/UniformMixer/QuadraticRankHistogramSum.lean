import SpinCodes.UniformMixer.QuadraticRankHistogramPage0
import SpinCodes.UniformMixer.QuadraticRankHistogramPage1
import SpinCodes.UniformMixer.QuadraticRankHistogramPage2
import SpinCodes.UniformMixer.QuadraticRankHistogramPage3
import SpinCodes.UniformMixer.QuadraticRankHistogramPage4
import SpinCodes.UniformMixer.QuadraticRankHistogramPage5
import SpinCodes.UniformMixer.QuadraticRankHistogramPage6
import SpinCodes.UniformMixer.QuadraticRankHistogramPage7

/-!
# Aggregating bounded quadratic-rank certificates

Every 512-form block is checked independently by the kernel. These
lemmas combine the checked blocks without re-evaluating their ranks.
-/

namespace Spin.UniformMixer

def rankHistogramChunks : Nat → Nat → RankHistogram
  | _, 0 => (0, 0, 0, 0)
  | start, n + 1 =>
      addRankHistogram (fastRankHistogramBlock start 512)
        (rankHistogramChunks (start + 512) n)

theorem rankHistogramChunks_eq_block (start n : Nat) :
    rankHistogramChunks start n =
      fastRankHistogramBlock start (512 * n) := by
  induction n generalizing start with
  | zero => simp [rankHistogramChunks, fastRankHistogramBlock,
      rankHistogramList]
  | succ n ih =>
      simp only [rankHistogramChunks]
      rw [ih (start + 512)]
      have hlen : 512 * (n + 1) = 512 + 512 * n := by omega
      rw [hlen]
      exact (fastRankHistogramBlock_add start 512 (512 * n)).symm

theorem rank_histogram_page_0 :
    fastRankHistogramBlock 0 4096 = (1, 203, 2548, 1344) := by
  have hc : rankHistogramChunks 0 8 = (1, 203, 2548, 1344) := by
    change addRankHistogram (fastRankHistogramBlock 0 512) (addRankHistogram (fastRankHistogramBlock 512 512) (addRankHistogram (fastRankHistogramBlock 1024 512) (addRankHistogram (fastRankHistogramBlock 1536 512) (addRankHistogram (fastRankHistogramBlock 2048 512) (addRankHistogram (fastRankHistogramBlock 2560 512) (addRankHistogram (fastRankHistogramBlock 3072 512) (addRankHistogram (fastRankHistogramBlock 3584 512) ((0, 0, 0, 0))))))))) = _
    rw [first512_fast_histogram, rank_histogram_block_1, rank_histogram_block_2, rank_histogram_block_3, rank_histogram_block_4, rank_histogram_block_5, rank_histogram_block_6, rank_histogram_block_7]
    norm_num [addRankHistogram]
  have h := rankHistogramChunks_eq_block 0 8
  norm_num at h
  exact h.symm.trans hc

theorem rank_histogram_page_1 :
    fastRankHistogramBlock 4096 4096 = (0, 64, 2240, 1792) := by
  have hc : rankHistogramChunks 4096 8 = (0, 64, 2240, 1792) := by
    change addRankHistogram (fastRankHistogramBlock 4096 512) (addRankHistogram (fastRankHistogramBlock 4608 512) (addRankHistogram (fastRankHistogramBlock 5120 512) (addRankHistogram (fastRankHistogramBlock 5632 512) (addRankHistogram (fastRankHistogramBlock 6144 512) (addRankHistogram (fastRankHistogramBlock 6656 512) (addRankHistogram (fastRankHistogramBlock 7168 512) (addRankHistogram (fastRankHistogramBlock 7680 512) ((0, 0, 0, 0))))))))) = _
    rw [rank_histogram_block_8, rank_histogram_block_9, rank_histogram_block_10, rank_histogram_block_11, rank_histogram_block_12, rank_histogram_block_13, rank_histogram_block_14, rank_histogram_block_15]
    norm_num [addRankHistogram]
  have h := rankHistogramChunks_eq_block 4096 8
  norm_num at h
  exact h.symm.trans hc

theorem rank_histogram_page_2 :
    fastRankHistogramBlock 8192 4096 = (0, 64, 2240, 1792) := by
  have hc : rankHistogramChunks 8192 8 = (0, 64, 2240, 1792) := by
    change addRankHistogram (fastRankHistogramBlock 8192 512) (addRankHistogram (fastRankHistogramBlock 8704 512) (addRankHistogram (fastRankHistogramBlock 9216 512) (addRankHistogram (fastRankHistogramBlock 9728 512) (addRankHistogram (fastRankHistogramBlock 10240 512) (addRankHistogram (fastRankHistogramBlock 10752 512) (addRankHistogram (fastRankHistogramBlock 11264 512) (addRankHistogram (fastRankHistogramBlock 11776 512) ((0, 0, 0, 0))))))))) = _
    rw [rank_histogram_block_16, rank_histogram_block_17, rank_histogram_block_18, rank_histogram_block_19, rank_histogram_block_20, rank_histogram_block_21, rank_histogram_block_22, rank_histogram_block_23]
    norm_num [addRankHistogram]
  have h := rankHistogramChunks_eq_block 8192 8
  norm_num at h
  exact h.symm.trans hc

theorem rank_histogram_page_3 :
    fastRankHistogramBlock 12288 4096 = (0, 64, 2240, 1792) := by
  have hc : rankHistogramChunks 12288 8 = (0, 64, 2240, 1792) := by
    change addRankHistogram (fastRankHistogramBlock 12288 512) (addRankHistogram (fastRankHistogramBlock 12800 512) (addRankHistogram (fastRankHistogramBlock 13312 512) (addRankHistogram (fastRankHistogramBlock 13824 512) (addRankHistogram (fastRankHistogramBlock 14336 512) (addRankHistogram (fastRankHistogramBlock 14848 512) (addRankHistogram (fastRankHistogramBlock 15360 512) (addRankHistogram (fastRankHistogramBlock 15872 512) ((0, 0, 0, 0))))))))) = _
    rw [rank_histogram_block_24, rank_histogram_block_25, rank_histogram_block_26, rank_histogram_block_27, rank_histogram_block_28, rank_histogram_block_29, rank_histogram_block_30, rank_histogram_block_31]
    norm_num [addRankHistogram]
  have h := rankHistogramChunks_eq_block 12288 8
  norm_num at h
  exact h.symm.trans hc

theorem rank_histogram_page_4 :
    fastRankHistogramBlock 16384 4096 = (0, 64, 2240, 1792) := by
  have hc : rankHistogramChunks 16384 8 = (0, 64, 2240, 1792) := by
    change addRankHistogram (fastRankHistogramBlock 16384 512) (addRankHistogram (fastRankHistogramBlock 16896 512) (addRankHistogram (fastRankHistogramBlock 17408 512) (addRankHistogram (fastRankHistogramBlock 17920 512) (addRankHistogram (fastRankHistogramBlock 18432 512) (addRankHistogram (fastRankHistogramBlock 18944 512) (addRankHistogram (fastRankHistogramBlock 19456 512) (addRankHistogram (fastRankHistogramBlock 19968 512) ((0, 0, 0, 0))))))))) = _
    rw [rank_histogram_block_32, rank_histogram_block_33, rank_histogram_block_34, rank_histogram_block_35, rank_histogram_block_36, rank_histogram_block_37, rank_histogram_block_38, rank_histogram_block_39]
    norm_num [addRankHistogram]
  have h := rankHistogramChunks_eq_block 16384 8
  norm_num at h
  exact h.symm.trans hc

theorem rank_histogram_page_5 :
    fastRankHistogramBlock 20480 4096 = (0, 64, 2240, 1792) := by
  have hc : rankHistogramChunks 20480 8 = (0, 64, 2240, 1792) := by
    change addRankHistogram (fastRankHistogramBlock 20480 512) (addRankHistogram (fastRankHistogramBlock 20992 512) (addRankHistogram (fastRankHistogramBlock 21504 512) (addRankHistogram (fastRankHistogramBlock 22016 512) (addRankHistogram (fastRankHistogramBlock 22528 512) (addRankHistogram (fastRankHistogramBlock 23040 512) (addRankHistogram (fastRankHistogramBlock 23552 512) (addRankHistogram (fastRankHistogramBlock 24064 512) ((0, 0, 0, 0))))))))) = _
    rw [rank_histogram_block_40, rank_histogram_block_41, rank_histogram_block_42, rank_histogram_block_43, rank_histogram_block_44, rank_histogram_block_45, rank_histogram_block_46, rank_histogram_block_47]
    norm_num [addRankHistogram]
  have h := rankHistogramChunks_eq_block 20480 8
  norm_num at h
  exact h.symm.trans hc

theorem rank_histogram_page_6 :
    fastRankHistogramBlock 24576 4096 = (0, 64, 2240, 1792) := by
  have hc : rankHistogramChunks 24576 8 = (0, 64, 2240, 1792) := by
    change addRankHistogram (fastRankHistogramBlock 24576 512) (addRankHistogram (fastRankHistogramBlock 25088 512) (addRankHistogram (fastRankHistogramBlock 25600 512) (addRankHistogram (fastRankHistogramBlock 26112 512) (addRankHistogram (fastRankHistogramBlock 26624 512) (addRankHistogram (fastRankHistogramBlock 27136 512) (addRankHistogram (fastRankHistogramBlock 27648 512) (addRankHistogram (fastRankHistogramBlock 28160 512) ((0, 0, 0, 0))))))))) = _
    rw [rank_histogram_block_48, rank_histogram_block_49, rank_histogram_block_50, rank_histogram_block_51, rank_histogram_block_52, rank_histogram_block_53, rank_histogram_block_54, rank_histogram_block_55]
    norm_num [addRankHistogram]
  have h := rankHistogramChunks_eq_block 24576 8
  norm_num at h
  exact h.symm.trans hc

theorem rank_histogram_page_7 :
    fastRankHistogramBlock 28672 4096 = (0, 64, 2240, 1792) := by
  have hc : rankHistogramChunks 28672 8 = (0, 64, 2240, 1792) := by
    change addRankHistogram (fastRankHistogramBlock 28672 512) (addRankHistogram (fastRankHistogramBlock 29184 512) (addRankHistogram (fastRankHistogramBlock 29696 512) (addRankHistogram (fastRankHistogramBlock 30208 512) (addRankHistogram (fastRankHistogramBlock 30720 512) (addRankHistogram (fastRankHistogramBlock 31232 512) (addRankHistogram (fastRankHistogramBlock 31744 512) (addRankHistogram (fastRankHistogramBlock 32256 512) ((0, 0, 0, 0))))))))) = _
    rw [rank_histogram_block_56, rank_histogram_block_57, rank_histogram_block_58, rank_histogram_block_59, rank_histogram_block_60, rank_histogram_block_61, rank_histogram_block_62, rank_histogram_block_63]
    norm_num [addRankHistogram]
  have h := rankHistogramChunks_eq_block 28672 8
  norm_num at h
  exact h.symm.trans hc

def rankHistogramPages : Nat → RankHistogram
  | 0 => (0, 0, 0, 0)
  | n + 1 => addRankHistogram (rankHistogramPages n)
      (fastRankHistogramBlock (4096 * n) 4096)

theorem rankHistogramPages_eq_block (n : Nat) :
    rankHistogramPages n = fastRankHistogramBlock 0 (4096 * n) := by
  induction n with
  | zero => simp [rankHistogramPages, fastRankHistogramBlock,
      rankHistogramList]
  | succ n ih =>
      simp only [rankHistogramPages]
      rw [ih]
      have hlen : 4096 * (n + 1) = 4096 * n + 4096 := by omega
      rw [hlen]
      simpa using (fastRankHistogramBlock_add 0 (4096 * n) 4096).symm

theorem rankHistogramPages_8 :
    rankHistogramPages 8 = (1, 651, 18228, 13888) := by
  change addRankHistogram (addRankHistogram (addRankHistogram (addRankHistogram (addRankHistogram (addRankHistogram (addRankHistogram (addRankHistogram ((0, 0, 0, 0)) (fastRankHistogramBlock 0 4096)) (fastRankHistogramBlock 4096 4096)) (fastRankHistogramBlock 8192 4096)) (fastRankHistogramBlock 12288 4096)) (fastRankHistogramBlock 16384 4096)) (fastRankHistogramBlock 20480 4096)) (fastRankHistogramBlock 24576 4096)) (fastRankHistogramBlock 28672 4096) = _
  rw [rank_histogram_page_0, rank_histogram_page_1,
    rank_histogram_page_2, rank_histogram_page_3,
    rank_histogram_page_4, rank_histogram_page_5,
    rank_histogram_page_6, rank_histogram_page_7]
  norm_num [addRankHistogram]

theorem full_fast_rank_histogram :
    fastRankHistogramBlock 0 (2 ^ 15) =
      (1, 651, 18228, 13888) := by
  have h := rankHistogramPages_eq_block 8
  norm_num at h
  exact h.symm.trans rankHistogramPages_8

theorem packedRadicalCount_64 : packedRadicalCount 64 = 1 := by
  rw [packedRadicalCount_eq_fastRange]
  exact (fastRankHistogramBlock_first 0 (2 ^ 15)).symm.trans
    (congrArg Prod.fst full_fast_rank_histogram)

theorem packedRadicalCount_16 : packedRadicalCount 16 = 651 := by
  rw [packedRadicalCount_eq_fastRange]
  exact (fastRankHistogramBlock_second 0 (2 ^ 15)).symm.trans
    (congrArg (fun h : RankHistogram => h.2.1) full_fast_rank_histogram)

theorem packedRadicalCount_4 : packedRadicalCount 4 = 18228 := by
  rw [packedRadicalCount_eq_fastRange]
  exact (fastRankHistogramBlock_third 0 (2 ^ 15)).symm.trans
    (congrArg (fun h : RankHistogram => h.2.2.1) full_fast_rank_histogram)

theorem packedRadicalCount_1 : packedRadicalCount 1 = 13888 := by
  rw [packedRadicalCount_eq_fastRange]
  exact (fastRankHistogramBlock_fourth 0 (2 ^ 15)).symm.trans
    (congrArg (fun h : RankHistogram => h.2.2.2) full_fast_rank_histogram)

theorem actualQuadraticRadicalCount_64 :
    (Finset.univ.filter (fun c : Fin 15 → ZMod 2 =>
      (quadraticRadical c).card = 64)).card = 1 := by
  rw [← packedRadicalCount_eq_actual, packedRadicalCount_64]

theorem actualQuadraticRadicalCount_16 :
    (Finset.univ.filter (fun c : Fin 15 → ZMod 2 =>
      (quadraticRadical c).card = 16)).card = 651 := by
  rw [← packedRadicalCount_eq_actual, packedRadicalCount_16]

theorem actualQuadraticRadicalCount_4 :
    (Finset.univ.filter (fun c : Fin 15 → ZMod 2 =>
      (quadraticRadical c).card = 4)).card = 18228 := by
  rw [← packedRadicalCount_eq_actual, packedRadicalCount_4]

theorem actualQuadraticRadicalCount_1 :
    (Finset.univ.filter (fun c : Fin 15 → ZMod 2 =>
      (quadraticRadical c).card = 1)).card = 13888 := by
  rw [← packedRadicalCount_eq_actual, packedRadicalCount_1]

end Spin.UniformMixer
