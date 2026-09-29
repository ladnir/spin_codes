import SpinCodes.Structured.SparseResiduals
import SpinCodes.Structured.SparseContributionSound
import SpinCodes.Structured.PolyPacked
noncomputable section
namespace Spin.Structured.SparsePolynomial.Residuals
theorem residual0_neg {x : ℝ} (h0 : 0 < x) (h1 : x ≤ 1) :
    packedValue terms0 x < 0 := by
  rw [← checkSumProducts_sound _ _ _ _ row0_checked x]
  change listEval (0 :: SparseData.row0) x / ((19383127018156436864577562173735346694589352901095787354573417703738701064321748192759872788736542462769824605938564968468946278861747241027326030322883069046934946944734453664138271591007394468111442610883779975025061919811395861669029927062903246299723626161766551040916505899358707046015496894601596574106248028513604460927317300826439620144199830354976510065198998329489152406998624613583049110557811427747118871154124659565113923301777225669109521374890654498845056602073493136771413915414780833785426317739441068626016835804080310642356039961818645559741790878022344208502815432471066969810635782778263092041015625000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 : ℕ) : ℝ) < 0
  rw [listEval, Int.cast_zero, zero_add]
  apply div_neg_of_neg_of_pos
  · exact residual_neg_of_tailBound _ SparseData.row0_cert h0 h1
  · positivity
def sourceResidual0 (x : ℝ) : ℝ :=
  [Row0Part0.sourceValue x, Row0Part1.sourceValue x, Row0Part2.sourceValue x, Row0Part3.sourceValue x, Row0Part4.sourceValue x, Row0Part5.sourceValue x, Row0Part6.sourceValue x, Row0Part7.sourceValue x, Row0Part8.sourceValue x, Row0Part9.sourceValue x, Row0Part10.sourceValue x, Row0Part11.sourceValue x, Row0Part12.sourceValue x, Row0Part13.sourceValue x, Row0Part14.sourceValue x, Row0Part15.sourceValue x, Row0Part16.sourceValue x].sum -
    (1 - 96 * (x / 10000)) * (witness 0).eval x
theorem sourceResidual0_neg {x : ℝ} (h0 : 0 < x) (h1 : x ≤ 1) :
    sourceResidual0 x < 0 := by
  have hc : negativeContraction.eval x = -(1 - 96 * (x / 10000)) := by
    norm_num [negativeContraction, RatPoly.eval, listEval]
    ring
  have hw : vector0.eval x = (witness 0).eval x :=
    (RatPoly.checkEq_sound _ _ vector0_checked x).symm
  have he : packedValue terms0 x = sourceResidual0 x := by
    simp only [sourceResidual0, terms0, packedValue, RatPoly.eval_constant,
      Int.cast_one, Nat.cast_one, div_one, mul_one, List.sum_cons, List.sum_nil,
      Row0Part0.expected_eval, Row0Part1.expected_eval, Row0Part2.expected_eval, Row0Part3.expected_eval, Row0Part4.expected_eval, Row0Part5.expected_eval, Row0Part6.expected_eval, Row0Part7.expected_eval, Row0Part8.expected_eval, Row0Part9.expected_eval, Row0Part10.expected_eval, Row0Part11.expected_eval, Row0Part12.expected_eval, Row0Part13.expected_eval, Row0Part14.expected_eval, Row0Part15.expected_eval, Row0Part16.expected_eval, hc, hw]
    ring
  exact he ▸ residual0_neg h0 h1
theorem residual1_neg {x : ℝ} (h0 : 0 < x) (h1 : x ≤ 1) :
    packedValue terms1 x < 0 := by
  rw [← checkSumProducts_sound _ _ _ _ row1_checked x]
  change listEval (0 :: SparseData.row1) x / ((20324643029936367628837552678762367424932336128913614129534466895280104729820112789474990849576631276356406066347424823247356875411177751512987364919786791242820765057627985016420123995268987647005487824264848703531929277504314604253741386740056708585486401578948199491177984256914156882068653640759974526096884920250624284012400811396783098253084192914639139065105974474345758496016175805563232138054046560038507423201565110778829767012237752628760847266154593150471956381462610992396948574942138398005719615697316679093456977586427707651497842246924024449160688628131401556086631189287916576802217605290934443473815917968750000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 : ℕ) : ℝ) < 0
  rw [listEval, Int.cast_zero, zero_add]
  apply div_neg_of_neg_of_pos
  · exact residual_neg_of_tailBound _ SparseData.row1_cert h0 h1
  · positivity
def sourceResidual1 (x : ℝ) : ℝ :=
  [Row1Part0.sourceValue x, Row1Part1.sourceValue x, Row1Part2.sourceValue x, Row1Part3.sourceValue x, Row1Part4.sourceValue x, Row1Part5.sourceValue x, Row1Part6.sourceValue x, Row1Part7.sourceValue x, Row1Part8.sourceValue x, Row1Part9.sourceValue x, Row1Part10.sourceValue x, Row1Part11.sourceValue x, Row1Part12.sourceValue x, Row1Part13.sourceValue x, Row1Part14.sourceValue x, Row1Part15.sourceValue x, Row1Part16.sourceValue x].sum -
    (1 - 96 * (x / 10000)) * (witness 1).eval x
theorem sourceResidual1_neg {x : ℝ} (h0 : 0 < x) (h1 : x ≤ 1) :
    sourceResidual1 x < 0 := by
  have hc : negativeContraction.eval x = -(1 - 96 * (x / 10000)) := by
    norm_num [negativeContraction, RatPoly.eval, listEval]
    ring
  have hw : vector1.eval x = (witness 1).eval x :=
    (RatPoly.checkEq_sound _ _ vector1_checked x).symm
  have he : packedValue terms1 x = sourceResidual1 x := by
    simp only [sourceResidual1, terms1, packedValue, RatPoly.eval_constant,
      Int.cast_one, Nat.cast_one, div_one, mul_one, List.sum_cons, List.sum_nil,
      Row1Part0.expected_eval, Row1Part1.expected_eval, Row1Part2.expected_eval, Row1Part3.expected_eval, Row1Part4.expected_eval, Row1Part5.expected_eval, Row1Part6.expected_eval, Row1Part7.expected_eval, Row1Part8.expected_eval, Row1Part9.expected_eval, Row1Part10.expected_eval, Row1Part11.expected_eval, Row1Part12.expected_eval, Row1Part13.expected_eval, Row1Part14.expected_eval, Row1Part15.expected_eval, Row1Part16.expected_eval, hc, hw]
    ring
  exact he ▸ residual1_neg h0 h1
theorem residual2_neg {x : ℝ} (h0 : 0 < x) (h1 : x ≤ 1) :
    packedValue terms2 x < 0 := by
  rw [← checkSumProducts_sound _ _ _ _ row2_checked x]
  change listEval (0 :: SparseData.row2) x / ((52498552946325637585287398569243195058600224220983865296587527990508510517125351335213901364456438586828596869375398318447922809187072132158046363587809281780206036143853085297413180279779795092215175050076104201222973323793644622787414001949566478276311375278423199285712733335609267226383332354083014200908253749007362525604031295837890742787716470298512896205168732067235094195209782105769828612593602264579464674129642681141717288192610115040089268488477314107669063333317924193361318169075543482048773767346168982098399373105742768863818926523804755152182058726463410219371768361930688517880128074466483667492866516113281250000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 : ℕ) : ℝ) < 0
  rw [listEval, Int.cast_zero, zero_add]
  apply div_neg_of_neg_of_pos
  · exact residual_neg_of_tailBound _ SparseData.row2_cert h0 h1
  · positivity
def sourceResidual2 (x : ℝ) : ℝ :=
  [Row2Part0.sourceValue x, Row2Part1.sourceValue x, Row2Part2.sourceValue x, Row2Part3.sourceValue x, Row2Part4.sourceValue x, Row2Part5.sourceValue x, Row2Part6.sourceValue x, Row2Part7.sourceValue x, Row2Part8.sourceValue x, Row2Part9.sourceValue x, Row2Part10.sourceValue x, Row2Part11.sourceValue x, Row2Part12.sourceValue x, Row2Part13.sourceValue x, Row2Part14.sourceValue x, Row2Part15.sourceValue x, Row2Part16.sourceValue x].sum -
    (1 - 96 * (x / 10000)) * (witness 2).eval x
theorem sourceResidual2_neg {x : ℝ} (h0 : 0 < x) (h1 : x ≤ 1) :
    sourceResidual2 x < 0 := by
  have hc : negativeContraction.eval x = -(1 - 96 * (x / 10000)) := by
    norm_num [negativeContraction, RatPoly.eval, listEval]
    ring
  have hw : vector2.eval x = (witness 2).eval x :=
    (RatPoly.checkEq_sound _ _ vector2_checked x).symm
  have he : packedValue terms2 x = sourceResidual2 x := by
    simp only [sourceResidual2, terms2, packedValue, RatPoly.eval_constant,
      Int.cast_one, Nat.cast_one, div_one, mul_one, List.sum_cons, List.sum_nil,
      Row2Part0.expected_eval, Row2Part1.expected_eval, Row2Part2.expected_eval, Row2Part3.expected_eval, Row2Part4.expected_eval, Row2Part5.expected_eval, Row2Part6.expected_eval, Row2Part7.expected_eval, Row2Part8.expected_eval, Row2Part9.expected_eval, Row2Part10.expected_eval, Row2Part11.expected_eval, Row2Part12.expected_eval, Row2Part13.expected_eval, Row2Part14.expected_eval, Row2Part15.expected_eval, Row2Part16.expected_eval, hc, hw]
    ring
  exact he ▸ residual2_neg h0 h1
theorem residual3_neg {x : ℝ} (h0 : 0 < x) (h1 : x ≤ 1) :
    packedValue terms3 x < 0 := by
  rw [← checkSumProducts_sound _ _ _ _ row3_checked x]
  change listEval (0 :: SparseData.row3) x / ((140097764405351382065577250614708998660058592936601542194881080309165761902650037457851111926131719387924707015332799306644030942209248241179021906392090352036763533542229700718183914699389131850808827572657602113445588509837240567121039378799210892279757766083689939092689845482909283388099229545758504408385827755287553189697478792958025896258509341760607585575775482051665313313039499827747359127606542938345431668128388308598473584015354828870048520205603610586203195337421777570592166527076159977453425311001603868991198946503246188841774626608047300528064626713709750926105148787761608963897685953270411118865013122558593750000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 : ℕ) : ℝ) < 0
  rw [listEval, Int.cast_zero, zero_add]
  apply div_neg_of_neg_of_pos
  · exact residual_neg_of_tailBound _ SparseData.row3_cert h0 h1
  · positivity
def sourceResidual3 (x : ℝ) : ℝ :=
  [Row3Part0.sourceValue x, Row3Part1.sourceValue x, Row3Part2.sourceValue x, Row3Part3.sourceValue x, Row3Part4.sourceValue x, Row3Part5.sourceValue x, Row3Part6.sourceValue x, Row3Part7.sourceValue x, Row3Part8.sourceValue x, Row3Part9.sourceValue x, Row3Part10.sourceValue x, Row3Part11.sourceValue x, Row3Part12.sourceValue x, Row3Part13.sourceValue x, Row3Part14.sourceValue x, Row3Part15.sourceValue x, Row3Part16.sourceValue x].sum -
    (1 - 96 * (x / 10000)) * (witness 3).eval x
theorem sourceResidual3_neg {x : ℝ} (h0 : 0 < x) (h1 : x ≤ 1) :
    sourceResidual3 x < 0 := by
  have hc : negativeContraction.eval x = -(1 - 96 * (x / 10000)) := by
    norm_num [negativeContraction, RatPoly.eval, listEval]
    ring
  have hw : vector3.eval x = (witness 3).eval x :=
    (RatPoly.checkEq_sound _ _ vector3_checked x).symm
  have he : packedValue terms3 x = sourceResidual3 x := by
    simp only [sourceResidual3, terms3, packedValue, RatPoly.eval_constant,
      Int.cast_one, Nat.cast_one, div_one, mul_one, List.sum_cons, List.sum_nil,
      Row3Part0.expected_eval, Row3Part1.expected_eval, Row3Part2.expected_eval, Row3Part3.expected_eval, Row3Part4.expected_eval, Row3Part5.expected_eval, Row3Part6.expected_eval, Row3Part7.expected_eval, Row3Part8.expected_eval, Row3Part9.expected_eval, Row3Part10.expected_eval, Row3Part11.expected_eval, Row3Part12.expected_eval, Row3Part13.expected_eval, Row3Part14.expected_eval, Row3Part15.expected_eval, Row3Part16.expected_eval, hc, hw]
    ring
  exact he ▸ residual3_neg h0 h1
theorem residual4_neg {x : ℝ} (h0 : 0 < x) (h1 : x ≤ 1) :
    packedValue terms4 x < 0 := by
  rw [← checkSumProducts_sound _ _ _ _ row4_checked x]
  change listEval (0 :: SparseData.row4) x / ((1192873624069995352504104804269242106536703739742068926876507396550884626697872239727076687952502066240633828439996710301210622374757433414048741434507206564832393521997244068598713497406332153990399085893928235258992461226005728438256335729160668283590782395070048776336727074022548781565491350829843664911152272854429389852971816021688596819571764366353085710870134747873826911889685374204311657414530046655220039175123057916720297855715245939534602886897879226594349591984422101754769308811929044717353689964891213212674083471525028589774059859314217918945689976273660088728280471130497111809098953472130233421921730041503906250000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 : ℕ) : ℝ) < 0
  rw [listEval, Int.cast_zero, zero_add]
  apply div_neg_of_neg_of_pos
  · exact residual_neg_of_tailBound _ SparseData.row4_cert h0 h1
  · positivity
def sourceResidual4 (x : ℝ) : ℝ :=
  [Row4Part0.sourceValue x, Row4Part1.sourceValue x, Row4Part2.sourceValue x, Row4Part3.sourceValue x, Row4Part4.sourceValue x, Row4Part5.sourceValue x, Row4Part6.sourceValue x, Row4Part7.sourceValue x, Row4Part8.sourceValue x, Row4Part9.sourceValue x, Row4Part10.sourceValue x, Row4Part11.sourceValue x, Row4Part12.sourceValue x, Row4Part13.sourceValue x, Row4Part14.sourceValue x, Row4Part15.sourceValue x, Row4Part16.sourceValue x].sum -
    (1 - 96 * (x / 10000)) * (witness 4).eval x
theorem sourceResidual4_neg {x : ℝ} (h0 : 0 < x) (h1 : x ≤ 1) :
    sourceResidual4 x < 0 := by
  have hc : negativeContraction.eval x = -(1 - 96 * (x / 10000)) := by
    norm_num [negativeContraction, RatPoly.eval, listEval]
    ring
  have hw : vector4.eval x = (witness 4).eval x :=
    (RatPoly.checkEq_sound _ _ vector4_checked x).symm
  have he : packedValue terms4 x = sourceResidual4 x := by
    simp only [sourceResidual4, terms4, packedValue, RatPoly.eval_constant,
      Int.cast_one, Nat.cast_one, div_one, mul_one, List.sum_cons, List.sum_nil,
      Row4Part0.expected_eval, Row4Part1.expected_eval, Row4Part2.expected_eval, Row4Part3.expected_eval, Row4Part4.expected_eval, Row4Part5.expected_eval, Row4Part6.expected_eval, Row4Part7.expected_eval, Row4Part8.expected_eval, Row4Part9.expected_eval, Row4Part10.expected_eval, Row4Part11.expected_eval, Row4Part12.expected_eval, Row4Part13.expected_eval, Row4Part14.expected_eval, Row4Part15.expected_eval, Row4Part16.expected_eval, hc, hw]
    ring
  exact he ▸ residual4_neg h0 h1
theorem residual5_neg {x : ℝ} (h0 : 0 < x) (h1 : x ≤ 1) :
    packedValue terms5 x < 0 := by
  rw [← checkSumProducts_sound _ _ _ _ row5_checked x]
  change listEval (0 :: SparseData.row5) x / ((139894517975052018389288875087921374985809269575312406053585735640212960855351836329956362017635953075161142954669325058411557373455136463663892032742892484124335325891653420868019713459436441974338772694414953626410269217062197421078501964931810325193902902067900457097778065640340141819278543009350904663124858906085046946857354784844058065275978499831461194185124422306921855728079338069691726806226002472745046593896372657490685286345232451343760911732942064654698475773607151460668197041326738593473368114844630702200264376727381911765259648185578060283573019827428436910544282475868729798129663777217501774430274963378906250000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 : ℕ) : ℝ) < 0
  rw [listEval, Int.cast_zero, zero_add]
  apply div_neg_of_neg_of_pos
  · exact residual_neg_of_tailBound _ SparseData.row5_cert h0 h1
  · positivity
def sourceResidual5 (x : ℝ) : ℝ :=
  [Row5Part0.sourceValue x, Row5Part1.sourceValue x, Row5Part2.sourceValue x, Row5Part3.sourceValue x, Row5Part4.sourceValue x, Row5Part5.sourceValue x, Row5Part6.sourceValue x, Row5Part7.sourceValue x, Row5Part8.sourceValue x, Row5Part9.sourceValue x, Row5Part10.sourceValue x, Row5Part11.sourceValue x, Row5Part12.sourceValue x, Row5Part13.sourceValue x, Row5Part14.sourceValue x, Row5Part15.sourceValue x, Row5Part16.sourceValue x].sum -
    (1 - 96 * (x / 10000)) * (witness 5).eval x
theorem sourceResidual5_neg {x : ℝ} (h0 : 0 < x) (h1 : x ≤ 1) :
    sourceResidual5 x < 0 := by
  have hc : negativeContraction.eval x = -(1 - 96 * (x / 10000)) := by
    norm_num [negativeContraction, RatPoly.eval, listEval]
    ring
  have hw : vector5.eval x = (witness 5).eval x :=
    (RatPoly.checkEq_sound _ _ vector5_checked x).symm
  have he : packedValue terms5 x = sourceResidual5 x := by
    simp only [sourceResidual5, terms5, packedValue, RatPoly.eval_constant,
      Int.cast_one, Nat.cast_one, div_one, mul_one, List.sum_cons, List.sum_nil,
      Row5Part0.expected_eval, Row5Part1.expected_eval, Row5Part2.expected_eval, Row5Part3.expected_eval, Row5Part4.expected_eval, Row5Part5.expected_eval, Row5Part6.expected_eval, Row5Part7.expected_eval, Row5Part8.expected_eval, Row5Part9.expected_eval, Row5Part10.expected_eval, Row5Part11.expected_eval, Row5Part12.expected_eval, Row5Part13.expected_eval, Row5Part14.expected_eval, Row5Part15.expected_eval, Row5Part16.expected_eval, hc, hw]
    ring
  exact he ▸ residual5_neg h0 h1
theorem residual6_neg {x : ℝ} (h0 : 0 < x) (h1 : x ≤ 1) :
    packedValue terms6 x < 0 := by
  rw [← checkSumProducts_sound _ _ _ _ row6_checked x]
  change listEval (0 :: SparseData.row6) x / ((426817503628663720205588606254009715923579058707185896720223804800882199326222368578974807841109256803484527393295921288194494383634732781772734663315522616099236066210187685344822603900648740587115244309561822774170514827590606689328569121541190880295214433157912189314737669395197294523441726455959465048034583325263109964260417039332445063314768051207421920367225463961260928416339691916827874899134977760808655887232867326355425107256992805203977792589246456159911084010714830840335920073784906358120111929643650260962596529314981860681454687185404513432374461190759432677819254975046248112846569711109623312950134277343750000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 : ℕ) : ℝ) < 0
  rw [listEval, Int.cast_zero, zero_add]
  apply div_neg_of_neg_of_pos
  · exact residual_neg_of_tailBound _ SparseData.row6_cert h0 h1
  · positivity
def sourceResidual6 (x : ℝ) : ℝ :=
  [Row6Part0.sourceValue x, Row6Part1.sourceValue x, Row6Part2.sourceValue x, Row6Part3.sourceValue x, Row6Part4.sourceValue x, Row6Part5.sourceValue x, Row6Part6.sourceValue x, Row6Part7.sourceValue x, Row6Part8.sourceValue x, Row6Part9.sourceValue x, Row6Part10.sourceValue x, Row6Part11.sourceValue x, Row6Part12.sourceValue x, Row6Part13.sourceValue x, Row6Part14.sourceValue x, Row6Part15.sourceValue x, Row6Part16.sourceValue x].sum -
    (1 - 96 * (x / 10000)) * (witness 6).eval x
theorem sourceResidual6_neg {x : ℝ} (h0 : 0 < x) (h1 : x ≤ 1) :
    sourceResidual6 x < 0 := by
  have hc : negativeContraction.eval x = -(1 - 96 * (x / 10000)) := by
    norm_num [negativeContraction, RatPoly.eval, listEval]
    ring
  have hw : vector6.eval x = (witness 6).eval x :=
    (RatPoly.checkEq_sound _ _ vector6_checked x).symm
  have he : packedValue terms6 x = sourceResidual6 x := by
    simp only [sourceResidual6, terms6, packedValue, RatPoly.eval_constant,
      Int.cast_one, Nat.cast_one, div_one, mul_one, List.sum_cons, List.sum_nil,
      Row6Part0.expected_eval, Row6Part1.expected_eval, Row6Part2.expected_eval, Row6Part3.expected_eval, Row6Part4.expected_eval, Row6Part5.expected_eval, Row6Part6.expected_eval, Row6Part7.expected_eval, Row6Part8.expected_eval, Row6Part9.expected_eval, Row6Part10.expected_eval, Row6Part11.expected_eval, Row6Part12.expected_eval, Row6Part13.expected_eval, Row6Part14.expected_eval, Row6Part15.expected_eval, Row6Part16.expected_eval, hc, hw]
    ring
  exact he ▸ residual6_neg h0 h1
end Spin.Structured.SparsePolynomial.Residuals
