// Exact store-policy control: same paired outer, cached instead of NT.
#include "../../src/packet/PacketPlan.h"
#include "../../src/packet/PacketRsField.h"
#define outerFast outerFastCached
#define _mm512_stream_si512 _mm512_store_si512
#define _mm_sfence() ((void)0)
#include "../../src/packet/PacketOuterFast.cpp"
#undef _mm_sfence
#undef _mm512_stream_si512
#undef outerFast
