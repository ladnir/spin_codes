// AVX2 + GFNI Turbo (PacketRs8T64) kernels: Avx2.cpp built with GFNI byte
// maps and transposes (transposeAvx2Gfni, forwardAvx2Gfni, stage entry points).
// Compiled with -mavx2 -mgfni; selected only when the CPU has both.
#define SPIN_RS8_GFNI 1
#include "Avx2.cpp"
