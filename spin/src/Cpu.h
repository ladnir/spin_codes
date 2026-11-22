#pragma once
#include <bit>
// Architecture selection for the private kernels. SPIN_PORTABLE=1 (CMake
// option SPIN_PORTABLE) compiles only the portable C++ kernel: no ISA headers,
// intrinsics or feature detection. Otherwise x86-64 adds the AVX2/AVX-512
// kernels and AArch64 the NEON kernel; every other target is portable.
#if !defined(SPIN_PORTABLE) || !SPIN_PORTABLE
#if defined(__x86_64__) || defined(_M_X64)
#define SPIN_X86 1
#endif
#if defined(__aarch64__) || defined(_M_ARM64)
#define SPIN_AARCH64 1
#endif
#endif
#if defined(SPIN_PORTABLE) && SPIN_PORTABLE
#define SPIN_PORTABLE_BUILD 1
#else
#define SPIN_PORTABLE_BUILD 0
#endif
#ifndef SPIN_X86
#define SPIN_X86 0
#endif
#ifndef SPIN_AARCH64
#define SPIN_AARCH64 0
#endif
// Tables are written as 64-bit words and read as bytes, and records are byte
// views of 64-bit lanes. Every supported target (x86, ARM, RISC-V, wasm) is
// little-endian.
static_assert(std::endian::native == std::endian::little, "SPIN requires a little-endian target");

namespace spin::detail {
bool cpu_avx2() noexcept;
// AVX2 with VEX-encoded (256-bit) GFNI.
bool cpu_avx2gfni() noexcept;
// The Turbo AVX-512 kernels' ISA: AVX2, AVX-512 F/VL/BW/DQ, VBMI and GFNI.
// False when the build excludes the AVX-512 kernels.
bool cpu_packet512() noexcept;
}
