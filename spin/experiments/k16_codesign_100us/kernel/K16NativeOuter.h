#pragma once
#include "K16CodeSign.h"
namespace spin::research::k16codesign {
struct NativeOuterTables {
    std::vector<std::uint64_t> gl;
    std::vector<std::array<std::uint8_t,3>> field;
};
// Exact same binary map as the Plan's scalar matrices; only the compact
// coefficient representation absorbs the interleaved input-coordinate layout.
void prepareNativeGl16(const rs::Plan&,NativeOuterTables&);
// New transitive family Q^T M_c^T, with fixed symbol-coordinate permutation Q.
// Also updates literal Plan matrices and old compact tables for scalar checks.
void customizeNativeField16(rs::Plan&,std::uint64_t seed,NativeOuterTables&);
void outerNativeGl16(const Block* scratch,Block*,const rs::Plan&,const NativeOuterTables&);
void outerNativeField16(const Block* scratch,Block*,const rs::Plan&,const NativeOuterTables&);
}
