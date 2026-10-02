#pragma once
#include "../kernel/K16NativeOuter.h"
namespace spin::research::k16codesign::outervariants {
// Exact native-field map, replacing two-input VBMI permutations by qword
// unpacks plus single-input VBMI. No changes to setup or the coefficient law.
void fieldUnaryPack(const Block*,Block*,const rs::Plan&,const NativeOuterTables&);
void fieldUnaryUnpack(const Block*,Block*,const rs::Plan&,const NativeOuterTables&);
void fieldUnaryBoth(const Block*,Block*,const rs::Plan&,const NativeOuterTables&);
}
