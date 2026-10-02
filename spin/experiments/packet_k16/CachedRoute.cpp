// Exact store-policy control: same generated schedule, cached instead of NT.
// Parse intrinsics and shared helpers before rebinding the two store operations.
#include "../../src/packet/PacketPlan.h"
#include "../../src/packet/PacketRsInner.h"
#define reverseRoute reverseRouteCached
#define transposeFast transposeFastCachedRoute
#define _mm512_stream_si512 _mm512_store_si512
#define _mm_sfence() ((void)0)
#include "../../src/packet/PacketFast.cpp"
#undef _mm_sfence
#undef _mm512_stream_si512
#undef transposeFast
#undef reverseRoute
