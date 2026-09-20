#pragma once
#include <immintrin.h>
#include "Spin.h"
namespace spin::detail::kernel {
void bchForward2(const block*,const block*,block*,block*);
}
