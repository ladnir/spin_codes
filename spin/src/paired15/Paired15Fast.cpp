#include "Paired15FastCommon.h"
#include "Paired15TransposeOuter.h"

namespace spin::detail::paired15 {
namespace {
using namespace fast;
struct ByteRoute {
    std::byte* scratch;
    const std::uint32_t* offsets;
    SPIN_FORCEINLINE void operator()(std::size_t packet, __m512i value) {
        _mm512_store_si512(reinterpret_cast<__m512i*>(scratch+offsets[packet]), value);
    }
};

// Exact selected mode-52 loop: peeled boundaries, wide state, raw feedback,
// and cached routing stores. Input and scratch remain disjoint under aliasing.
static SPIN_NOINLINE void reverseRoute(const Block* __restrict input,
                                      Block* __restrict scratch, const Plan& p) {
    const auto epochs = p.n/64;
    ByteRoute emit{reinterpret_cast<std::byte*>(scratch), p.routeBytes.data()};
    alignas(64) __m512i packets[16], high[16];
    PackedState before;
    WordPackets16 state, feedback;
    const auto first = epochs-1;
    loadPackets(input+64*first, packets, std::make_index_sequence<16>{});
    for(unsigned h=16; h-->0;) emit(16*first+h, packets[h]);
    highFromPackets(packets, high);
    finishWide(high, state);
    for(std::size_t epoch=first; --epoch;) {
        packWordGroup<0>(state, before);
        packWordGroup<1>(state, before);
        streamHigh<false,true>(state, input+64*epoch, nullptr, 16*epoch, emit, high);
        finishWide(high, feedback);
        updateWide<false>(before, p.updates.data()+4*epoch, feedback, state);
    }
    streamHigh<false,true>(state, input, nullptr, 0, emit, high);
}

template<bool Stream>
static SPIN_FORCEINLINE void storeOutput(Block* output, __m512i value) {
    if constexpr(Stream) _mm512_stream_si512(reinterpret_cast<__m512i*>(output), value);
    else _mm512_storeu_si512(output, value);
}
template<bool Stream> struct DirectOutput {
    Block* output;
    SPIN_FORCEINLINE void operator()(std::size_t packet, __m512i value) {
        storeOutput<Stream>(output+4*packet, value);
    }
};
struct EpochOutput {
    __m512i* output;
    SPIN_FORCEINLINE void operator()(std::size_t packet, __m512i value) {
        output[packet&15] = value;
    }
};

template<unsigned... Packet>
static SPIN_FORCEINLINE void gatherPackets(const Block* input, const std::uint32_t* route,
                                          __m512i* packets, std::integer_sequence<unsigned,Packet...>) {
    ((packets[Packet] = _mm512_loadu_si512(input+route[Packet])), ...);
}
template<class Emit>
static SPIN_FORCEINLINE void firstForward(const Block* input, const Plan& p,
                                         WordPackets16& state, Emit& emit) {
    alignas(64) __m512i packets[16], high[16];
    gatherPackets(input, p.route.data(), packets, std::make_integer_sequence<unsigned,16>{});
    // The initial state is zero, so the first epoch passes through unchanged.
    for(unsigned packet=0; packet<16; ++packet) emit(packet, packets[packet]);
    highFromPackets(packets, high);
    finishWide(high, state);
}
template<bool Last, class Emit>
static SPIN_FORCEINLINE void forwardStep(const Block* input, const Plan& p, std::size_t epoch,
                                        WordPackets16& state, Emit& emit) {
    alignas(64) __m512i high[16];
    PackedState before;
    if constexpr(!Last) {
        packWordGroup<0>(state, before);
        packWordGroup<1>(state, before);
    }
    streamHigh<true,true>(state, input, p.route.data()+16*epoch, 16*epoch, emit, high);
    if constexpr(!Last) {
        WordPackets16 feedback;
        finishWide(high, feedback);
        updateWide<false>(before, p.forwardUpdates.data()+4*epoch, feedback, state);
    }
}
template<bool Stream>
static SPIN_NOINLINE void forwardNarrow(const Block* input, Block* output, const Plan& p) {
    WordPackets16 state;
    DirectOutput<Stream> emit{output};
    firstForward(input, p, state, emit);
    for(std::size_t epoch=1; epoch+1<p.n/64; ++epoch)
        forwardStep<false>(input, p, epoch, state, emit);
    forwardStep<true>(input, p, p.n/64-1, state, emit);
}

template<bool Stream,unsigned L>
static SPIN_FORCEINLINE void storeEpoch(Block* output, const __m512i (&v)[L][16]) {
    const auto low = _mm512_setr_epi64(0,1,8,9,2,3,10,11);
    const auto high = _mm512_setr_epi64(4,5,12,13,6,7,14,15);
    for(unsigned packet=0; packet<16; ++packet) {
        const auto a = _mm512_permutex2var_epi64(v[0][packet], low, v[1][packet]);
        const auto b = _mm512_permutex2var_epi64(v[0][packet], high, v[1][packet]);
        if constexpr(L==2) {
            storeOutput<Stream>(output+8*packet, a);
            storeOutput<Stream>(output+8*packet+4, b);
        } else {
            const auto c = _mm512_permutex2var_epi64(v[2][packet], low, v[3][packet]);
            const auto d = _mm512_permutex2var_epi64(v[2][packet], high, v[3][packet]);
            storeOutput<Stream>(output+16*packet, _mm512_shuffle_i32x4(a,c,0x44));
            storeOutput<Stream>(output+16*packet+4, _mm512_shuffle_i32x4(a,c,0xee));
            storeOutput<Stream>(output+16*packet+8, _mm512_shuffle_i32x4(b,d,0x44));
            storeOutput<Stream>(output+16*packet+12, _mm512_shuffle_i32x4(b,d,0xee));
        }
    }
}
template<unsigned Lane,unsigned L>
static SPIN_FORCEINLINE void firstLane(const Block* input, const Plan& p,
                                      WordPackets16 (&states)[L], __m512i (&v)[L][16]) {
    EpochOutput emit{v[Lane]};
    firstForward(input+Lane*p.scratchBlocks(), p, states[Lane], emit);
}
template<bool Last,unsigned Lane,unsigned L>
static SPIN_FORCEINLINE void stepLane(const Block* input, const Plan& p, std::size_t epoch,
                                     WordPackets16 (&states)[L], __m512i (&v)[L][16]) {
    EpochOutput emit{v[Lane]};
    forwardStep<Last>(input+Lane*p.scratchBlocks(), p, epoch, states[Lane], emit);
}
template<bool Stream,unsigned L,unsigned... Lane>
static SPIN_NOINLINE void forwardWide(const Block* input, Block* output, const Plan& p,
                                      std::integer_sequence<unsigned,Lane...>) {
    WordPackets16 states[L];
    alignas(64) __m512i values[L][16];
    (firstLane<Lane>(input,p,states,values), ...);
    storeEpoch<Stream>(output, values);
    for(std::size_t epoch=1; epoch+1<p.n/64; ++epoch) {
        (stepLane<false,Lane>(input,p,epoch,states,values), ...);
        storeEpoch<Stream>(output+64*L*epoch, values);
    }
    (stepLane<true,Lane>(input,p,p.n/64-1,states,values), ...);
    storeEpoch<Stream>(output+L*(p.n-64), values);
}
template<bool Stream>
static SPIN_FORCEINLINE void forwardInner(const Block* scratch, Block* output, const Plan& p, unsigned lanes) {
    if(lanes==1) forwardNarrow<Stream>(scratch, output, p);
    else if(lanes==2) forwardWide<Stream,2>(scratch, output, p, std::make_integer_sequence<unsigned,2>{});
    else forwardWide<Stream,4>(scratch, output, p, std::make_integer_sequence<unsigned,4>{});
    if constexpr(Stream) _mm_sfence();
}
}
void transposeFast(const Block* input, Block* output, Block* scratch, const Plan& p) {
    reverseRoute(input, scratch, p);
    outerfast::run(scratch, output, p);
}
void forwardFast(const Block* input, Block* output, Block* scratch, const Plan& p, unsigned lanes, bool stream) {
    forwardOuterFast(input, scratch, p, lanes);
    if(stream && (reinterpret_cast<std::uintptr_t>(output)&63U)==0)
        forwardInner<true>(scratch, output, p, lanes);
    else forwardInner<false>(scratch, output, p, lanes);
}
void forwardFastCached(const Block* input, Block* output, Block* scratch, const Plan& p, unsigned lanes) {
    forwardOuterFast(input, scratch, p, lanes);
    forwardInner<false>(scratch, output, p, lanes);
}
}
