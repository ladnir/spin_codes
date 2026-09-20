#pragma once
namespace spin::detail {
bool cpu_avx2() noexcept;
bool cpu_avx512f() noexcept;
bool cpu_avx512vl() noexcept;
bool cpu_wide512() noexcept;
}
