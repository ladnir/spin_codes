import SpinCodes.Structured.ConcreteFixedNumericDefs

namespace Spin.Structured.ConcreteFixedNumeric
open DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace Segment00
def lo : QInput := ⟨13, 125⟩
def hi : QInput := ⟨18296026121, 120000000000⟩
def s : QInput := ⟨833, 500⟩
def c : QInput := ⟨-43311, 250000⟩
def u : QInput := ⟨96838234476759233, 100000000000000000⟩
def v : QInput := ⟨39357831448094154198102851, 39321525000000000000000000⟩
theorem check_lo : rowCheck lo s c u v=true := by decide +kernel
theorem check_hi : rowCheck hi s c u v=true := by decide +kernel
end Segment00
namespace Segment01
def lo : QInput := ⟨18296026121, 120000000000⟩
def hi : QInput := ⟨791599208623, 5000000000000⟩
def s : QInput := ⟨83, 50⟩
def c : QInput := ⟨-3446583973879, 20000000000000⟩
def u : QInput := ⟨96838234476759233, 100000000000000000⟩
def v : QInput := ⟨39357831448094154198102851, 39321525000000000000000000⟩
theorem check_lo : rowCheck lo s c u v=true := by decide +kernel
theorem check_hi : rowCheck hi s c u v=true := by decide +kernel
end Segment01
namespace Segment02
def lo : QInput := ⟨791599208623, 5000000000000⟩
def hi : QInput := ⟨1609032935677, 10000000000000⟩
def s : QInput := ⟨33, 20⟩
def c : QInput := ⟨-5335812508647, 31250000000000⟩
def u : QInput := ⟨96838234476759233, 100000000000000000⟩
def v : QInput := ⟨39357831448094154198102851, 39321525000000000000000000⟩
theorem check_lo : rowCheck lo s c u v=true := by decide +kernel
theorem check_hi : rowCheck hi s c u v=true := by decide +kernel
end Segment02
namespace Segment03
def lo : QInput := ⟨1609032935677, 10000000000000⟩
def hi : QInput := ⟨2469788513329, 15000000000000⟩
def s : QInput := ⟨163, 100⟩
def c : QInput := ⟨-3350558688107, 20000000000000⟩
def u : QInput := ⟨96838234476759233, 100000000000000000⟩
def v : QInput := ⟨39357831448094154198102851, 39321525000000000000000000⟩
theorem check_lo : rowCheck lo s c u v=true := by decide +kernel
theorem check_hi : rowCheck hi s c u v=true := by decide +kernel
end Segment03
namespace Segment04
def lo : QInput := ⟨2469788513329, 15000000000000⟩
def hi : QInput := ⟨4271577070219, 25000000000000⟩
def s : QInput := ⟨8, 5⟩
def c : QInput := ⟨-40647089344673, 250000000000000⟩
def u : QInput := ⟨96838234476759233, 100000000000000000⟩
def v : QInput := ⟨39357831448094154198102851, 39321525000000000000000000⟩
theorem check_lo : rowCheck lo s c u v=true := by decide +kernel
theorem check_hi : rowCheck hi s c u v=true := by decide +kernel
end Segment04
namespace Segment05
def lo : QInput := ⟨4271577070219, 25000000000000⟩
def hi : QInput := ⟨13956016256281, 75000000000000⟩
def s : QInput := ⟨31, 20⟩
def c : QInput := ⟨-77022601619127, 500000000000000⟩
def u : QInput := ⟨96838234476759233, 100000000000000000⟩
def v : QInput := ⟨39357831448094154198102851, 39321525000000000000000000⟩
theorem check_lo : rowCheck lo s c u v=true := by decide +kernel
theorem check_hi : rowCheck hi s c u v=true := by decide +kernel
end Segment05
namespace Segment06
def lo : QInput := ⟨13956016256281, 75000000000000⟩
def hi : QInput := ⟨428613422794653, 2000000000000000⟩
def s : QInput := ⟨7, 5⟩
def c : QInput := ⟨-31533292681423, 250000000000000⟩
def u : QInput := ⟨96838234476759233, 100000000000000000⟩
def v : QInput := ⟨39357831448094154198102851, 39321525000000000000000000⟩
theorem check_lo : rowCheck lo s c u v=true := by decide +kernel
theorem check_hi : rowCheck hi s c u v=true := by decide +kernel
end Segment06
namespace Segment07
def lo : QInput := ⟨428613422794653, 2000000000000000⟩
def hi : QInput := ⟨499776708257287, 2000000000000000⟩
def s : QInput := ⟨6, 5⟩
def c : QInput := ⟨-832718284462267, 10000000000000000⟩
def u : QInput := ⟨96838234476759233, 100000000000000000⟩
def v : QInput := ⟨39357831448094154198102851, 39321525000000000000000000⟩
theorem check_lo : rowCheck lo s c u v=true := by decide +kernel
theorem check_hi : rowCheck hi s c u v=true := by decide +kernel
end Segment07
namespace Segment08
def lo : QInput := ⟨499776708257287, 2000000000000000⟩
def hi : QInput := ⟨144603079418107, 500000000000000⟩
def s : QInput := ⟨1, 1⟩
def c : QInput := ⟨-16647078810249, 500000000000000⟩
def u : QInput := ⟨96838234476759233, 100000000000000000⟩
def v : QInput := ⟨39357831448094154198102851, 39321525000000000000000000⟩
theorem check_lo : rowCheck lo s c u v=true := by decide +kernel
theorem check_hi : rowCheck hi s c u v=true := by decide +kernel
end Segment08
namespace Segment09
def lo : QInput := ⟨144603079418107, 500000000000000⟩
def hi : QInput := ⟨663870278625033, 2000000000000000⟩
def s : QInput := ⟨4, 5⟩
def c : QInput := ⟨30683842683431, 1250000000000000⟩
def u : QInput := ⟨96838234476759233, 100000000000000000⟩
def v : QInput := ⟨39357831448094154198102851, 39321525000000000000000000⟩
theorem check_lo : rowCheck lo s c u v=true := by decide +kernel
theorem check_hi : rowCheck hi s c u v=true := by decide +kernel
end Segment09
namespace Segment10
def lo : QInput := ⟨663870278625033, 2000000000000000⟩
def hi : QInput := ⟨755265762055919, 2000000000000000⟩
def s : QInput := ⟨3, 5⟩
def c : QInput := ⟨909341020092481, 10000000000000000⟩
def u : QInput := ⟨96838234476759233, 100000000000000000⟩
def v : QInput := ⟨39357831448094154198102851, 39321525000000000000000000⟩
theorem check_lo : rowCheck lo s c u v=true := by decide +kernel
theorem check_hi : rowCheck hi s c u v=true := by decide +kernel
end Segment10
namespace Segment11
def lo : QInput := ⟨755265762055919, 2000000000000000⟩
def hi : QInput := ⟨10640568152347, 25000000000000⟩
def s : QInput := ⟨2, 5⟩
def c : QInput := ⟨4161516955371, 25000000000000⟩
def u : QInput := ⟨96838234476759233, 100000000000000000⟩
def v : QInput := ⟨39357831448094154198102851, 39321525000000000000000000⟩
theorem check_lo : rowCheck lo s c u v=true := by decide +kernel
theorem check_hi : rowCheck hi s c u v=true := by decide +kernel
end Segment11
namespace Segment12
def lo : QInput := ⟨10640568152347, 25000000000000⟩
def hi : QInput := ⟨46255924070167, 100000000000000⟩
def s : QInput := ⟨1, 5⟩
def c : QInput := ⟨15724076464601, 62500000000000⟩
def u : QInput := ⟨96838234476759233, 100000000000000000⟩
def v : QInput := ⟨39357831448094154198102851, 39321525000000000000000000⟩
theorem check_lo : rowCheck lo s c u v=true := by decide +kernel
theorem check_hi : rowCheck hi s c u v=true := by decide +kernel
end Segment12
namespace Segment13
def lo : QInput := ⟨46255924070167, 100000000000000⟩
def hi : QInput := ⟨12032207560909, 25000000000000⟩
def s : QInput := ⟨1, 10⟩
def c : QInput := ⟨297841147503783, 1000000000000000⟩
def u : QInput := ⟨96838234476759233, 100000000000000000⟩
def v : QInput := ⟨39357831448094154198102851, 39321525000000000000000000⟩
theorem check_lo : rowCheck lo s c u v=true := by decide +kernel
theorem check_hi : rowCheck hi s c u v=true := by decide +kernel
end Segment13
namespace Segment14
def lo : QInput := ⟨12032207560909, 25000000000000⟩
def hi : QInput := ⟨246680276543717, 500000000000000⟩
def s : QInput := ⟨1, 20⟩
def c : QInput := ⟨109535552428011023053662565657799459079, 340282366920938463463374607431768211456⟩
def u : QInput := ⟨96838234476759233, 100000000000000000⟩
def v : QInput := ⟨39357831448094154198102851, 39321525000000000000000000⟩
theorem check_lo : rowCheck lo s c u v=true := by decide +kernel
theorem check_hi : rowCheck hi s c u v=true := by decide +kernel
end Segment14
namespace Segment15
def lo : QInput := ⟨246680276543717, 500000000000000⟩
def hi : QInput := ⟨253319723456283, 500000000000000⟩
def s : QInput := ⟨0, 1⟩
def c : QInput := ⟨3465735902799727, 10000000000000000⟩
def u : QInput := ⟨49900171601753829, 50000000000000000⟩
def v : QInput := ⟨6658451237321793314805221, 6553587500000000000000000⟩
theorem check_lo : rowCheck lo s c u v=true := by decide +kernel
theorem check_hi : rowCheck hi s c u v=true := by decide +kernel
end Segment15
namespace Segment16
def lo : QInput := ⟨253319723456283, 500000000000000⟩
def hi : QInput := ⟨12967792439091, 25000000000000⟩
def s : QInput := ⟨-1, 20⟩
def c : QInput := ⟨632748353870289731134156480146939348259, 1701411834604692317316873037158841057280⟩
def u : QInput := ⟨10499206513492543, 10000000000000000⟩
def v : QInput := ⟨4098984794519733140991421, 3932152500000000000000000⟩
theorem check_lo : rowCheck lo s c u v=true := by decide +kernel
theorem check_hi : rowCheck hi s c u v=true := by decide +kernel
end Segment16
namespace Segment17
def lo : QInput := ⟨12967792439091, 25000000000000⟩
def hi : QInput := ⟨53744075929833, 100000000000000⟩
def s : QInput := ⟨-1, 10⟩
def c : QInput := ⟨397841147503783, 1000000000000000⟩
def u : QInput := ⟨11167521538448069, 10000000000000000⟩
def v : QInput := ⟨4232750283614970196036943, 3932152500000000000000000⟩
theorem check_lo : rowCheck lo s c u v=true := by decide +kernel
theorem check_hi : rowCheck hi s c u v=true := by decide +kernel
end Segment17
namespace Segment18
def lo : QInput := ⟨53744075929833, 100000000000000⟩
def hi : QInput := ⟨14359431847653, 25000000000000⟩
def s : QInput := ⟨-1, 5⟩
def c : QInput := ⟨28224076464601, 62500000000000⟩
def u : QInput := ⟨12492915353143481, 10000000000000000⟩
def v : QInput := ⟨4498032291719354693380907, 3932152500000000000000000⟩
theorem check_lo : rowCheck lo s c u v=true := by decide +kernel
theorem check_hi : rowCheck hi s c u v=true := by decide +kernel
end Segment18
namespace Segment19
def lo : QInput := ⟨14359431847653, 25000000000000⟩
def hi : QInput := ⟨1244734237944081, 2000000000000000⟩
def s : QInput := ⟨-2, 5⟩
def c : QInput := ⟨14161516955371, 25000000000000⟩
def u : QInput := ⟨1486892840550039, 1000000000000000⟩
def v : QInput := ⟨165786641888809208943511, 131071750000000000000000⟩
theorem check_lo : rowCheck lo s c u v=true := by decide +kernel
theorem check_hi : rowCheck hi s c u v=true := by decide +kernel
end Segment19
namespace Segment20
def lo : QInput := ⟨1244734237944081, 2000000000000000⟩
def hi : QInput := ⟨1336129721374967, 2000000000000000⟩
def s : QInput := ⟨-3, 5⟩
def c : QInput := ⟨6909341020092481, 10000000000000000⟩
def u : QInput := ⟨907492662852861, 500000000000000⟩
def v : QInput := ⟨93838122449314232558589, 65535875000000000000000⟩
theorem check_lo : rowCheck lo s c u v=true := by decide +kernel
theorem check_hi : rowCheck hi s c u v=true := by decide +kernel
end Segment20
namespace Segment21
def lo : QInput := ⟨1336129721374967, 2000000000000000⟩
def hi : QInput := ⟨355396920581893, 500000000000000⟩
def s : QInput := ⟨-4, 5⟩
def c : QInput := ⟨1030683842683431, 1250000000000000⟩
def u : QInput := ⟨4430894293586297, 2000000000000000⟩
def v : QInput := ⟨1286365017064497977886059, 786430500000000000000000⟩
theorem check_lo : rowCheck lo s c u v=true := by decide +kernel
theorem check_hi : rowCheck hi s c u v=true := by decide +kernel
end Segment21
namespace Segment22
def lo : QInput := ⟨355396920581893, 500000000000000⟩
def hi : QInput := ⟨1500223291742713, 2000000000000000⟩
def s : QInput := ⟨-1, 1⟩
def c : QInput := ⟨483352921189751, 500000000000000⟩
def u : QInput := ⟨6761285149803551, 2500000000000000⟩
def v : QInput := ⟨1852677220254577125135197, 983038125000000000000000⟩
theorem check_lo : rowCheck lo s c u v=true := by decide +kernel
theorem check_hi : rowCheck hi s c u v=true := by decide +kernel
end Segment22
namespace Segment23
def lo : QInput := ⟨1500223291742713, 2000000000000000⟩
def hi : QInput := ⟨1571386577205347, 2000000000000000⟩
def s : QInput := ⟨-6, 5⟩
def c : QInput := ⟨11167281715537733, 10000000000000000⟩
def u : QInput := ⟨8254819521857547, 2500000000000000⟩
def v : QInput := ⟨717204374560242563086603, 327679375000000000000000⟩
theorem check_lo : rowCheck lo s c u v=true := by decide +kernel
theorem check_hi : rowCheck hi s c u v=true := by decide +kernel
end Segment23
namespace Segment24
def lo : QInput := ⟨1571386577205347, 2000000000000000⟩
def hi : QInput := ⟨61043983743719, 75000000000000⟩
def s : QInput := ⟨-7, 5⟩
def c : QInput := ⟨318466707318577, 250000000000000⟩
def u : QInput := ⟨7974382895176137, 2000000000000000⟩
def v : QInput := ⟨665202040243017991693513, 262143500000000000000000⟩
theorem check_lo : rowCheck lo s c u v=true := by decide +kernel
theorem check_hi : rowCheck hi s c u v=true := by decide +kernel
end Segment24
namespace Segment25
def lo : QInput := ⟨61043983743719, 75000000000000⟩
def hi : QInput := ⟨20728422929781, 25000000000000⟩
def s : QInput := ⟨-31, 20⟩
def c : QInput := ⟨697977398380873, 500000000000000⟩
def u : QInput := ⟨918788069754107, 200000000000000⟩
def v : QInput := ⟨223849176744953983046129, 78643050000000000000000⟩
theorem check_lo : rowCheck lo s c u v=true := by decide +kernel
theorem check_hi : rowCheck hi s c u v=true := by decide +kernel
end Segment25
namespace Segment26
def lo : QInput := ⟨20728422929781, 25000000000000⟩
def hi : QInput := ⟨12530211486671, 15000000000000⟩
def s : QInput := ⟨-8, 5⟩
def c : QInput := ⟨359352910655327, 250000000000000⟩
def u : QInput := ⟨49510688747396809, 10000000000000000⟩
def v : QInput := ⟨11907263535066708858469723, 3932152500000000000000000⟩
theorem check_lo : rowCheck lo s c u v=true := by decide +kernel
theorem check_hi : rowCheck hi s c u v=true := by decide +kernel
end Segment26
namespace Segment27
def lo : QInput := ⟨12530211486671, 15000000000000⟩
def hi : QInput := ⟨8390967064323, 10000000000000⟩
def s : QInput := ⟨-163, 100⟩
def c : QInput := ⟨29249441311893, 20000000000000⟩
def u : QInput := ⟨51330699179736179, 10000000000000000⟩
def v : QInput := ⟨12271544714674350803841113, 3932152500000000000000000⟩
theorem check_lo : rowCheck lo s c u v=true := by decide +kernel
theorem check_hi : rowCheck hi s c u v=true := by decide +kernel
end Segment27
namespace Segment28
def lo : QInput := ⟨8390967064323, 10000000000000⟩
def hi : QInput := ⟨4208400791377, 5000000000000⟩
def s : QInput := ⟨-33, 20⟩
def c : QInput := ⟨46226687491353, 31250000000000⟩
def u : QInput := ⟨52546977145860421, 10000000000000000⟩
def v : QInput := ⟨12514986820476470457979087, 3932152500000000000000000⟩
theorem check_lo : rowCheck lo s c u v=true := by decide +kernel
theorem check_hi : rowCheck hi s c u v=true := by decide +kernel
end Segment28
namespace Segment29
def lo : QInput := ⟨4208400791377, 5000000000000⟩
def hi : QInput := ⟨101703973879, 120000000000⟩
def s : QInput := ⟨-83, 50⟩
def c : QInput := ⟨29753416026121, 20000000000000⟩
def u : QInput := ⟨10848858740933581, 2000000000000000⟩
def v : QInput := ⟨2570942082128062141845607, 786430500000000000000000⟩
theorem check_lo : rowCheck lo s c u v=true := by decide +kernel
theorem check_hi : rowCheck hi s c u v=true := by decide +kernel
end Segment29
namespace Segment30
def lo : QInput := ⟨101703973879, 120000000000⟩
def hi : QInput := ⟨112, 125⟩
def s : QInput := ⟨-833, 500⟩
def c : QInput := ⟨373189, 250000⟩
def u : QInput := ⟨16961761337730057, 2500000000000000⟩
def v : QInput := ⟨1298112223420622763683593, 327679375000000000000000⟩
theorem check_lo : rowCheck lo s c u v=true := by decide +kernel
theorem check_hi : rowCheck hi s c u v=true := by decide +kernel
end Segment30
end Spin.Structured.ConcreteFixedNumeric
