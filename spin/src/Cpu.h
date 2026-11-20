#pragma once
namespace spin::detail {
bool cpu_avx2() noexcept;
// AVX2 with VEX-encoded (256-bit) GFNI.
bool cpu_avx2gfni() noexcept;
// The Turbo AVX-512 kernels' ISA: AVX2, AVX-512 F/VL/BW/DQ, VBMI and GFNI.
// False when the build excludes the AVX-512 kernels.
bool cpu_packet512() noexcept;
}
