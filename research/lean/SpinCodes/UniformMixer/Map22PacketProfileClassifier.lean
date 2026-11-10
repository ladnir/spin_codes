import SpinCodes.UniformMixer.Map22PacketProfileBank
namespace Spin.UniformMixer
open Spin.Structured
set_option maxRecDepth 10000
set_option maxHeartbeats 0
def profileZeroClassifier : List Nat := [
171,0,106,156,80,166,146,116,41,169,153,109,163,83,119,143,15,170,155,107,165,81,117,145,168,42,110,152,84,162,142,120,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,106,0,106,0,43,95,43,95,20,103,20,103,87,51,87,51,4,105,4,105,93,45,93,45,101,22,101,22,53,85,53,85,156,0,0,156,64,136,136,64,33,150,150,33,128,72,72,128,11,154,154,11,134,66,66,134,148,35,35,148,74,126,126,74,80,0,43,64,80,0,43,64,16,77,60,47,16,77,60,47,2,79,63,44,2,79,63,44,75,18,48,59,75,18,48,59,166,0,95,136,0,166,136,95,37,160,132,99,160,37,99,132,13,164,135,96,164,13,96,135,158,39,100,131,39,158,131,100,146,0,43,136,43,136,146,0,29,140,124,55,124,55,29,140,9,144,133,46,133,46,9,144,138,31,58,121,58,121,138,31,116,0,95,64,64,95,0,116,24,113,68,91,91,68,113,24,6,115,65,94,94,65,115,6,111,26,90,69,69,90,26,111,41,0,20,33,16,37,29,24,41,0,20,33,16,37,29,24,1,40,32,21,36,17,25,28,1,40,32,21,36,17,25,28,169,0,103,150,77,160,140,113,0,169,150,103,160,77,113,140,14,167,149,104,159,78,114,139,167,14,104,149,78,159,139,114,153,0,20,150,60,132,124,68,20,150,153,0,124,68,60,132,10,151,147,23,130,62,70,122,147,23,10,151,70,122,130,62,109,0,103,33,47,99,55,91,33,103,0,109,91,55,99,47,5,108,34,102,97,49,89,57,102,34,108,5,57,89,49,97,163,0,87,128,16,160,124,91,16,160,124,91,163,0,87,128,12,161,127,88,157,19,92,123,157,19,92,123,12,161,127,88,83,0,51,72,77,37,55,68,37,77,68,55,0,83,72,51,3,82,71,52,38,76,67,56,76,38,56,67,82,3,52,71,119,0,87,72,60,99,29,113,29,113,60,99,87,72,119,0,7,118,73,86,98,61,112,30,112,30,98,61,73,86,7,118,143,0,51,128,47,132,140,24,24,140,132,47,128,51,0,143,8,141,125,54,129,50,27,137,137,27,50,129,54,125,141,8,15,0,4,11,2,13,9,6,1,14,10,5,12,3,7,8,15,0,4,11,2,13,9,6,1,14,10,5,12,3,7,8,170,0,105,154,79,164,144,115,40,167,151,108,161,82,118,141,0,170,154,105,164,79,115,144,167,40,108,151,82,161,141,118,155,0,4,154,63,135,133,65,32,149,147,34,127,71,73,125,4,154,155,0,133,65,63,135,147,34,32,149,73,125,127,71,107,0,105,11,44,96,46,94,21,104,23,102,88,52,86,54,11,105,0,107,94,46,96,44,102,23,104,21,54,86,52,88,165,0,93,134,2,164,133,94,36,159,130,97,157,38,98,129,2,164,133,94,165,0,93,134,157,38,98,129,36,159,130,97,81,0,45,66,79,13,46,65,17,78,62,49,19,76,61,50,13,79,65,46,0,81,66,45,76,19,50,61,78,17,49,62,117,0,93,66,63,96,9,115,25,114,70,89,92,67,112,27,9,115,63,96,93,66,117,0,112,27,92,67,70,89,25,114,145,0,45,134,44,135,144,6,28,139,122,57,123,56,30,137,6,144,135,44,134,45,0,145,137,30,56,123,57,122,139,28,168,0,101,148,75,158,138,111,1,167,147,102,157,76,112,137,1,167,147,102,157,76,112,137,168,0,101,148,75,158,138,111,42,0,22,35,18,39,31,26,40,14,23,34,19,38,30,27,14,40,34,23,38,19,27,30,0,42,35,22,39,18,26,31,110,0,101,35,48,100,58,90,32,104,10,108,92,56,98,50,10,108,32,104,98,50,92,56,101,35,110,0,58,90,48,100,152,0,22,148,59,131,121,69,21,149,151,5,123,67,61,129,5,151,149,21,129,61,67,123,148,22,0,152,69,121,131,59,84,0,53,74,75,39,58,69,36,78,70,57,12,82,73,54,12,82,73,54,36,78,70,57,75,39,58,69,84,0,53,74,162,0,85,126,18,158,121,90,17,159,122,89,161,3,86,125,3,161,125,86,159,17,89,122,158,18,90,121,0,162,126,85,142,0,53,126,48,131,138,26,25,139,130,49,127,52,7,141,7,141,127,52,130,49,25,139,138,26,48,131,53,126,142,0,120,0,85,74,59,100,31,111,28,114,62,97,88,71,118,8,8,118,71,88,97,62,114,28,111,31,100,59,74,85,0,120
]
def profileZeroChunks : List (List Nat) := [
  [171,0,106,156,80,166,146,116,41,169,153,109,163,83,119,143,15,170,155,107,165,81,117,145,168,42,110,152,84,162,142,120],
  [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
  [106,0,106,0,43,95,43,95,20,103,20,103,87,51,87,51,4,105,4,105,93,45,93,45,101,22,101,22,53,85,53,85],
  [156,0,0,156,64,136,136,64,33,150,150,33,128,72,72,128,11,154,154,11,134,66,66,134,148,35,35,148,74,126,126,74],
  [80,0,43,64,80,0,43,64,16,77,60,47,16,77,60,47,2,79,63,44,2,79,63,44,75,18,48,59,75,18,48,59],
  [166,0,95,136,0,166,136,95,37,160,132,99,160,37,99,132,13,164,135,96,164,13,96,135,158,39,100,131,39,158,131,100],
  [146,0,43,136,43,136,146,0,29,140,124,55,124,55,29,140,9,144,133,46,133,46,9,144,138,31,58,121,58,121,138,31],
  [116,0,95,64,64,95,0,116,24,113,68,91,91,68,113,24,6,115,65,94,94,65,115,6,111,26,90,69,69,90,26,111],
  [41,0,20,33,16,37,29,24,41,0,20,33,16,37,29,24,1,40,32,21,36,17,25,28,1,40,32,21,36,17,25,28],
  [169,0,103,150,77,160,140,113,0,169,150,103,160,77,113,140,14,167,149,104,159,78,114,139,167,14,104,149,78,159,139,114],
  [153,0,20,150,60,132,124,68,20,150,153,0,124,68,60,132,10,151,147,23,130,62,70,122,147,23,10,151,70,122,130,62],
  [109,0,103,33,47,99,55,91,33,103,0,109,91,55,99,47,5,108,34,102,97,49,89,57,102,34,108,5,57,89,49,97],
  [163,0,87,128,16,160,124,91,16,160,124,91,163,0,87,128,12,161,127,88,157,19,92,123,157,19,92,123,12,161,127,88],
  [83,0,51,72,77,37,55,68,37,77,68,55,0,83,72,51,3,82,71,52,38,76,67,56,76,38,56,67,82,3,52,71],
  [119,0,87,72,60,99,29,113,29,113,60,99,87,72,119,0,7,118,73,86,98,61,112,30,112,30,98,61,73,86,7,118],
  [143,0,51,128,47,132,140,24,24,140,132,47,128,51,0,143,8,141,125,54,129,50,27,137,137,27,50,129,54,125,141,8],
  [15,0,4,11,2,13,9,6,1,14,10,5,12,3,7,8,15,0,4,11,2,13,9,6,1,14,10,5,12,3,7,8],
  [170,0,105,154,79,164,144,115,40,167,151,108,161,82,118,141,0,170,154,105,164,79,115,144,167,40,108,151,82,161,141,118],
  [155,0,4,154,63,135,133,65,32,149,147,34,127,71,73,125,4,154,155,0,133,65,63,135,147,34,32,149,73,125,127,71],
  [107,0,105,11,44,96,46,94,21,104,23,102,88,52,86,54,11,105,0,107,94,46,96,44,102,23,104,21,54,86,52,88],
  [165,0,93,134,2,164,133,94,36,159,130,97,157,38,98,129,2,164,133,94,165,0,93,134,157,38,98,129,36,159,130,97],
  [81,0,45,66,79,13,46,65,17,78,62,49,19,76,61,50,13,79,65,46,0,81,66,45,76,19,50,61,78,17,49,62],
  [117,0,93,66,63,96,9,115,25,114,70,89,92,67,112,27,9,115,63,96,93,66,117,0,112,27,92,67,70,89,25,114],
  [145,0,45,134,44,135,144,6,28,139,122,57,123,56,30,137,6,144,135,44,134,45,0,145,137,30,56,123,57,122,139,28],
  [168,0,101,148,75,158,138,111,1,167,147,102,157,76,112,137,1,167,147,102,157,76,112,137,168,0,101,148,75,158,138,111],
  [42,0,22,35,18,39,31,26,40,14,23,34,19,38,30,27,14,40,34,23,38,19,27,30,0,42,35,22,39,18,26,31],
  [110,0,101,35,48,100,58,90,32,104,10,108,92,56,98,50,10,108,32,104,98,50,92,56,101,35,110,0,58,90,48,100],
  [152,0,22,148,59,131,121,69,21,149,151,5,123,67,61,129,5,151,149,21,129,61,67,123,148,22,0,152,69,121,131,59],
  [84,0,53,74,75,39,58,69,36,78,70,57,12,82,73,54,12,82,73,54,36,78,70,57,75,39,58,69,84,0,53,74],
  [162,0,85,126,18,158,121,90,17,159,122,89,161,3,86,125,3,161,125,86,159,17,89,122,158,18,90,121,0,162,126,85],
  [142,0,53,126,48,131,138,26,25,139,130,49,127,52,7,141,7,141,127,52,130,49,25,139,138,26,48,131,53,126,142,0],
  [120,0,85,74,59,100,31,111,28,114,62,97,88,71,118,8,8,118,71,88,97,62,114,28,111,31,100,59,74,85,0,120]]
def profileZeroIndex (p : Fin 1024) : Nat := (profileZeroChunks.getD (p.val/32) []).getD (p.val%32) 0
def profilePackedAffine (p : Nat) : Nat := PackedMap.eval (profileQuad4Rows.take 5) p
def profileClassWitness (p : Fin 1024) : ProfileZeroWitness :=
  profileZeroBank.getD ((profileZeroIndex p)-1) profileZeroWitness0
def profileZeroClass (p : Fin 1024) : Nat :=
  if profileZeroIndex p = 0 then 0 else (profileClassWitness p).kind.val+1
def profileZeroClassifierValid (p : Fin 1024) : Prop :=
  if profileZeroIndex p = 0 then
    ∀ v : Fin 16, ((profilePackedAffine (p.val%32)).testBit v.val ||
      (profilePackedAffine (p.val/32)).testBit v.val) = true
  else (profileZeroIndex p)-1 < profileZeroBank.length ∧
    ∀ v : Fin 16, (profileClassWitness p).zeroMask.testBit v.val =
      !((profilePackedAffine (p.val%32)).testBit v.val ||
        (profilePackedAffine (p.val/32)).testBit v.val)
instance : DecidablePred profileZeroClassifierValid := fun p => by
  unfold profileZeroClassifierValid; infer_instance
theorem profileZeroBlock0_checked : ∀ p : Fin 32,
  profileZeroClassifierValid ⟨0*32+p.val, by omega⟩ := by decide
theorem profileZeroBlock1_checked : ∀ p : Fin 32,
  profileZeroClassifierValid ⟨1*32+p.val, by omega⟩ := by decide
theorem profileZeroBlock2_checked : ∀ p : Fin 32,
  profileZeroClassifierValid ⟨2*32+p.val, by omega⟩ := by decide
theorem profileZeroBlock3_checked : ∀ p : Fin 32,
  profileZeroClassifierValid ⟨3*32+p.val, by omega⟩ := by decide
theorem profileZeroBlock4_checked : ∀ p : Fin 32,
  profileZeroClassifierValid ⟨4*32+p.val, by omega⟩ := by decide
theorem profileZeroBlock5_checked : ∀ p : Fin 32,
  profileZeroClassifierValid ⟨5*32+p.val, by omega⟩ := by decide
theorem profileZeroBlock6_checked : ∀ p : Fin 32,
  profileZeroClassifierValid ⟨6*32+p.val, by omega⟩ := by decide
theorem profileZeroBlock7_checked : ∀ p : Fin 32,
  profileZeroClassifierValid ⟨7*32+p.val, by omega⟩ := by decide
theorem profileZeroBlock8_checked : ∀ p : Fin 32,
  profileZeroClassifierValid ⟨8*32+p.val, by omega⟩ := by decide
theorem profileZeroBlock9_checked : ∀ p : Fin 32,
  profileZeroClassifierValid ⟨9*32+p.val, by omega⟩ := by decide
theorem profileZeroBlock10_checked : ∀ p : Fin 32,
  profileZeroClassifierValid ⟨10*32+p.val, by omega⟩ := by decide
theorem profileZeroBlock11_checked : ∀ p : Fin 32,
  profileZeroClassifierValid ⟨11*32+p.val, by omega⟩ := by decide
theorem profileZeroBlock12_checked : ∀ p : Fin 32,
  profileZeroClassifierValid ⟨12*32+p.val, by omega⟩ := by decide
theorem profileZeroBlock13_checked : ∀ p : Fin 32,
  profileZeroClassifierValid ⟨13*32+p.val, by omega⟩ := by decide
theorem profileZeroBlock14_checked : ∀ p : Fin 32,
  profileZeroClassifierValid ⟨14*32+p.val, by omega⟩ := by decide
theorem profileZeroBlock15_checked : ∀ p : Fin 32,
  profileZeroClassifierValid ⟨15*32+p.val, by omega⟩ := by decide
theorem profileZeroBlock16_checked : ∀ p : Fin 32,
  profileZeroClassifierValid ⟨16*32+p.val, by omega⟩ := by decide
theorem profileZeroBlock17_checked : ∀ p : Fin 32,
  profileZeroClassifierValid ⟨17*32+p.val, by omega⟩ := by decide
theorem profileZeroBlock18_checked : ∀ p : Fin 32,
  profileZeroClassifierValid ⟨18*32+p.val, by omega⟩ := by decide
theorem profileZeroBlock19_checked : ∀ p : Fin 32,
  profileZeroClassifierValid ⟨19*32+p.val, by omega⟩ := by decide
theorem profileZeroBlock20_checked : ∀ p : Fin 32,
  profileZeroClassifierValid ⟨20*32+p.val, by omega⟩ := by decide
theorem profileZeroBlock21_checked : ∀ p : Fin 32,
  profileZeroClassifierValid ⟨21*32+p.val, by omega⟩ := by decide
theorem profileZeroBlock22_checked : ∀ p : Fin 32,
  profileZeroClassifierValid ⟨22*32+p.val, by omega⟩ := by decide
theorem profileZeroBlock23_checked : ∀ p : Fin 32,
  profileZeroClassifierValid ⟨23*32+p.val, by omega⟩ := by decide
theorem profileZeroBlock24_checked : ∀ p : Fin 32,
  profileZeroClassifierValid ⟨24*32+p.val, by omega⟩ := by decide
theorem profileZeroBlock25_checked : ∀ p : Fin 32,
  profileZeroClassifierValid ⟨25*32+p.val, by omega⟩ := by decide
theorem profileZeroBlock26_checked : ∀ p : Fin 32,
  profileZeroClassifierValid ⟨26*32+p.val, by omega⟩ := by decide
theorem profileZeroBlock27_checked : ∀ p : Fin 32,
  profileZeroClassifierValid ⟨27*32+p.val, by omega⟩ := by decide
theorem profileZeroBlock28_checked : ∀ p : Fin 32,
  profileZeroClassifierValid ⟨28*32+p.val, by omega⟩ := by decide
theorem profileZeroBlock29_checked : ∀ p : Fin 32,
  profileZeroClassifierValid ⟨29*32+p.val, by omega⟩ := by decide
theorem profileZeroBlock30_checked : ∀ p : Fin 32,
  profileZeroClassifierValid ⟨30*32+p.val, by omega⟩ := by decide
theorem profileZeroBlock31_checked : ∀ p : Fin 32,
  profileZeroClassifierValid ⟨31*32+p.val, by omega⟩ := by decide
theorem profileZeroClassifier_checked (p : Fin 1024) : profileZeroClassifierValid p := by
  have hb : p.val / 32 < 32 := by omega
  let hc := (⟨p.val / 32, hb⟩ : Fin 32)
  have hr : ∀ b : Fin 32, ∀ r : Fin 32, profileZeroClassifierValid ⟨b.val*32+r.val, by omega⟩ := by
    intro b; fin_cases b
    · exact profileZeroBlock0_checked
    · exact profileZeroBlock1_checked
    · exact profileZeroBlock2_checked
    · exact profileZeroBlock3_checked
    · exact profileZeroBlock4_checked
    · exact profileZeroBlock5_checked
    · exact profileZeroBlock6_checked
    · exact profileZeroBlock7_checked
    · exact profileZeroBlock8_checked
    · exact profileZeroBlock9_checked
    · exact profileZeroBlock10_checked
    · exact profileZeroBlock11_checked
    · exact profileZeroBlock12_checked
    · exact profileZeroBlock13_checked
    · exact profileZeroBlock14_checked
    · exact profileZeroBlock15_checked
    · exact profileZeroBlock16_checked
    · exact profileZeroBlock17_checked
    · exact profileZeroBlock18_checked
    · exact profileZeroBlock19_checked
    · exact profileZeroBlock20_checked
    · exact profileZeroBlock21_checked
    · exact profileZeroBlock22_checked
    · exact profileZeroBlock23_checked
    · exact profileZeroBlock24_checked
    · exact profileZeroBlock25_checked
    · exact profileZeroBlock26_checked
    · exact profileZeroBlock27_checked
    · exact profileZeroBlock28_checked
    · exact profileZeroBlock29_checked
    · exact profileZeroBlock30_checked
    · exact profileZeroBlock31_checked
  have hh := hr hc ⟨p.val%32, Nat.mod_lt _ (by decide)⟩
  have he : (⟨hc.val*32+(p.val%32), by omega⟩ : Fin 1024) = p := by
    apply Fin.ext; dsimp [hc]; omega
  rw [he] at hh; exact hh
def profileZeroKinds : List Nat := [3,0,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,0,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,0,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,0,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,0,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,0,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,0,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,0,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,2,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,2,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,2,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,2,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,2,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,2,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,2]
def profileZeroKindChunks : List (List Nat) := [[3,0,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[2,0,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[2,0,0,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[2,0,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[2,0,1,1,0,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[2,0,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[2,0,1,1,1,1,0,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[2,0,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[2,0,1,1,1,1,1,1,0,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[2,0,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[2,0,1,1,1,1,1,1,1,1,0,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[2,0,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[2,0,1,1,1,1,1,1,1,1,1,1,0,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[2,0,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[2,0,1,1,1,1,1,1,1,1,1,1,1,1,0,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1],[2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,2,1,1,1,1,1,1,1,1,1,1,1,1],[2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1],[2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,2,1,1,1,1,1,1,1,1,1,1],[2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1],[2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,2,1,1,1,1,1,1,1,1],[2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1],[2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,2,1,1,1,1,1,1],[2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1],[2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,2,1,1,1,1],[2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1],[2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,2,1,1],[2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0],[2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,2]]
theorem profileZeroKinds_block0 :
  List.ofFn (fun p : Fin 32 => profileZeroClass ⟨0*32+p.val, by omega⟩) =
  [3,0,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2] := by decide
theorem profileZeroKinds_block1 :
  List.ofFn (fun p : Fin 32 => profileZeroClass ⟨1*32+p.val, by omega⟩) =
  [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] := by decide
theorem profileZeroKinds_block2 :
  List.ofFn (fun p : Fin 32 => profileZeroClass ⟨2*32+p.val, by omega⟩) =
  [2,0,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1] := by decide
theorem profileZeroKinds_block3 :
  List.ofFn (fun p : Fin 32 => profileZeroClass ⟨3*32+p.val, by omega⟩) =
  [2,0,0,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1] := by decide
theorem profileZeroKinds_block4 :
  List.ofFn (fun p : Fin 32 => profileZeroClass ⟨4*32+p.val, by omega⟩) =
  [2,0,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1] := by decide
theorem profileZeroKinds_block5 :
  List.ofFn (fun p : Fin 32 => profileZeroClass ⟨5*32+p.val, by omega⟩) =
  [2,0,1,1,0,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1] := by decide
theorem profileZeroKinds_block6 :
  List.ofFn (fun p : Fin 32 => profileZeroClass ⟨6*32+p.val, by omega⟩) =
  [2,0,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1] := by decide
theorem profileZeroKinds_block7 :
  List.ofFn (fun p : Fin 32 => profileZeroClass ⟨7*32+p.val, by omega⟩) =
  [2,0,1,1,1,1,0,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1] := by decide
theorem profileZeroKinds_block8 :
  List.ofFn (fun p : Fin 32 => profileZeroClass ⟨8*32+p.val, by omega⟩) =
  [2,0,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1] := by decide
theorem profileZeroKinds_block9 :
  List.ofFn (fun p : Fin 32 => profileZeroClass ⟨9*32+p.val, by omega⟩) =
  [2,0,1,1,1,1,1,1,0,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1] := by decide
theorem profileZeroKinds_block10 :
  List.ofFn (fun p : Fin 32 => profileZeroClass ⟨10*32+p.val, by omega⟩) =
  [2,0,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1] := by decide
theorem profileZeroKinds_block11 :
  List.ofFn (fun p : Fin 32 => profileZeroClass ⟨11*32+p.val, by omega⟩) =
  [2,0,1,1,1,1,1,1,1,1,0,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1] := by decide
theorem profileZeroKinds_block12 :
  List.ofFn (fun p : Fin 32 => profileZeroClass ⟨12*32+p.val, by omega⟩) =
  [2,0,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1] := by decide
theorem profileZeroKinds_block13 :
  List.ofFn (fun p : Fin 32 => profileZeroClass ⟨13*32+p.val, by omega⟩) =
  [2,0,1,1,1,1,1,1,1,1,1,1,0,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1] := by decide
theorem profileZeroKinds_block14 :
  List.ofFn (fun p : Fin 32 => profileZeroClass ⟨14*32+p.val, by omega⟩) =
  [2,0,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1] := by decide
theorem profileZeroKinds_block15 :
  List.ofFn (fun p : Fin 32 => profileZeroClass ⟨15*32+p.val, by omega⟩) =
  [2,0,1,1,1,1,1,1,1,1,1,1,1,1,0,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1] := by decide
theorem profileZeroKinds_block16 :
  List.ofFn (fun p : Fin 32 => profileZeroClass ⟨16*32+p.val, by omega⟩) =
  [2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1] := by decide
theorem profileZeroKinds_block17 :
  List.ofFn (fun p : Fin 32 => profileZeroClass ⟨17*32+p.val, by omega⟩) =
  [2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1] := by decide
theorem profileZeroKinds_block18 :
  List.ofFn (fun p : Fin 32 => profileZeroClass ⟨18*32+p.val, by omega⟩) =
  [2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1,1,1] := by decide
theorem profileZeroKinds_block19 :
  List.ofFn (fun p : Fin 32 => profileZeroClass ⟨19*32+p.val, by omega⟩) =
  [2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,2,1,1,1,1,1,1,1,1,1,1,1,1] := by decide
theorem profileZeroKinds_block20 :
  List.ofFn (fun p : Fin 32 => profileZeroClass ⟨20*32+p.val, by omega⟩) =
  [2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1,1,1] := by decide
theorem profileZeroKinds_block21 :
  List.ofFn (fun p : Fin 32 => profileZeroClass ⟨21*32+p.val, by omega⟩) =
  [2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,2,1,1,1,1,1,1,1,1,1,1] := by decide
theorem profileZeroKinds_block22 :
  List.ofFn (fun p : Fin 32 => profileZeroClass ⟨22*32+p.val, by omega⟩) =
  [2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1,1,1] := by decide
theorem profileZeroKinds_block23 :
  List.ofFn (fun p : Fin 32 => profileZeroClass ⟨23*32+p.val, by omega⟩) =
  [2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,2,1,1,1,1,1,1,1,1] := by decide
theorem profileZeroKinds_block24 :
  List.ofFn (fun p : Fin 32 => profileZeroClass ⟨24*32+p.val, by omega⟩) =
  [2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1,1,1] := by decide
theorem profileZeroKinds_block25 :
  List.ofFn (fun p : Fin 32 => profileZeroClass ⟨25*32+p.val, by omega⟩) =
  [2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,2,1,1,1,1,1,1] := by decide
theorem profileZeroKinds_block26 :
  List.ofFn (fun p : Fin 32 => profileZeroClass ⟨26*32+p.val, by omega⟩) =
  [2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1,1,1] := by decide
theorem profileZeroKinds_block27 :
  List.ofFn (fun p : Fin 32 => profileZeroClass ⟨27*32+p.val, by omega⟩) =
  [2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,2,1,1,1,1] := by decide
theorem profileZeroKinds_block28 :
  List.ofFn (fun p : Fin 32 => profileZeroClass ⟨28*32+p.val, by omega⟩) =
  [2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0,1,1] := by decide
theorem profileZeroKinds_block29 :
  List.ofFn (fun p : Fin 32 => profileZeroClass ⟨29*32+p.val, by omega⟩) =
  [2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,2,1,1] := by decide
theorem profileZeroKinds_block30 :
  List.ofFn (fun p : Fin 32 => profileZeroClass ⟨30*32+p.val, by omega⟩) =
  [2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,0] := by decide
theorem profileZeroKinds_block31 :
  List.ofFn (fun p : Fin 32 => profileZeroClass ⟨31*32+p.val, by omega⟩) =
  [2,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,2] := by decide
theorem profileZeroKinds_counts : ∀ k : Fin 4,
  profileZeroKinds.count k.val = ([93,840,90,1] : List Nat).getD k.val 0 := by decide
end Spin.UniformMixer
