#pragma once
#include "LiteralEncoder.h"
#include "SparseByteField128MulYSimd.h"

// Isolated gamma30 construction. The measured literal encoder is unchanged.
// Each record is one sparse-tower GF(2^128) element, not eight subfield
// symbols. The fixed byte maps act on all 16 coefficients of that element.
// For each packet: y=x+A*q; q_next=D2*H30*D1*(q+Y*C*x), C=A^T.
namespace spin::research::char2::field128_lift::gamma30 {
namespace literal=spin::research::char2::field128_lift::encoder;
namespace field=sparse_byte;
namespace old=spin::research::char2::hdgroup;
namespace sr=spin::research::char2::hdsumrank;
using Element=field::Element;
using Block=spin::research::char2::Block;
using Kind=hdbw::detail::Kind;
inline constexpr unsigned stateSize=30,packetWidth=64;

#if defined(_MSC_VER)
#define GAMMA30_INLINE __forceinline
#define GAMMA30_NOINLINE __declspec(noinline)
#else
#define GAMMA30_INLINE inline __attribute__((always_inline))
#define GAMMA30_NOINLINE inline __attribute__((noinline))
#endif

struct alignas(64) State {__m512i values[8]{};};
struct alignas(64) Packet {__m512i values[16]{};};
struct alignas(64) Refresh {std::array<Element,32> before{},after{};};
static_assert(sizeof(Block)==sizeof(Element));
static_assert(sizeof(State)==512&&sizeof(Packet)==1024&&sizeof(Refresh)==1024);

// State coordinates are the 30 Lagrange values, not FFT coefficients.
// Coordinates 30/31 are inert storage padding; their diagonal keys are one.
struct Plan {
    literal::Plan outerPlan;
    std::vector<Refresh> refresh;
    explicit Plan(std::size_t k,std::uint64_t seed=113):outerPlan(k,seed,false){
        // Reuse the frozen outer and route setup byte-for-byte. Its unused
        // state16 diagonals are discarded, not used as gamma30 randomness.
        decltype(outerPlan.before){}.swap(outerPlan.before);
        decltype(outerPlan.after){}.swap(outerPlan.after);
        decltype(outerPlan.refresh){}.swap(outerPlan.refresh);
        // Deterministic research setup ONLY. MT19937 is neither ideal
        // independent setup nor a cryptographic seed-expansion guarantee.
        std::mt19937_64 random(seed^0x97a75b63c012e48dULL);
        const auto sample=[&](){Element x;do{x={random(),random()};}while(!literal::nonzero(x));return x;};
        refresh.resize(outerPlan.n/64);
        for(auto& d:refresh){
            for(unsigned j=0;j!=30;++j){d.before[j]=sample();d.after[j]=sample();}
            d.before[30]=d.before[31]=d.after[30]=d.after[31]=field::embedByte(1);
        }
    }
    std::size_t scratchSize() const noexcept{return outerPlan.scratchSize();}
    std::size_t hotKeyBytes() const noexcept{
        return outerPlan.outerKeys.size()*sizeof(Element)+refresh.size()*sizeof(Refresh);
    }
};

namespace detail {
// Compile-time coefficients of the two missing evaluations. Translation by
// 128 leaves all differences unchanged, so the same tails serve nodes 0..29
// and nodes 128..157. This is not zero-padding the missing evaluations.
constexpr std::uint8_t completionCoefficient(unsigned target,unsigned j) noexcept {
    if(j>=30)return 0;
    std::uint8_t numerator=1,denominator=1;
    for(unsigned k=0;k!=30;++k)if(k!=j){
        numerator=sfMultiply(numerator,static_cast<std::uint8_t>(target^k));
        denominator=sfMultiply(denominator,static_cast<std::uint8_t>(j^k));
    }
    return sfMultiply(numerator,hdbw::gf8Inverse(denominator));
}
inline constexpr auto completion=[] {
    std::array<std::array<std::uint8_t,32>,2> result{};
    for(unsigned i=0;i!=2;++i)for(unsigned j=0;j!=30;++j)result[i][j]=completionCoefficient(30+i,j);
    return result;
}();

template<unsigned Tail,unsigned Vector>
GAMMA30_INLINE __m512i tailProduct(__m512i x) noexcept {
    return sr::byteProduct4<true,completion[Tail][4*Vector],completion[Tail][4*Vector+1],
        completion[Tail][4*Vector+2],completion[Tail][4*Vector+3]>(x);
}
GAMMA30_INLINE void clearPadding(__m512i* x) noexcept {
    x[7]=_mm512_maskz_mov_epi64(0x0f,x[7]);
}
GAMMA30_INLINE void complete(__m512i* x) noexcept {
    clearPadding(x);
    auto a=_mm512_setzero_si512(),b=a;
    c8Unroll(std::make_index_sequence<8>{},[&](auto i){
        a=_mm512_xor_si512(a,tailProduct<0,i>(x[i]));
        b=_mm512_xor_si512(b,tailProduct<1,i>(x[i]));
    });
    a=hdbw::detail::repeatRecordSum(a);b=hdbw::detail::repeatRecordSum(b);
    x[7]=_mm512_mask_mov_epi64(x[7],0x30,a);
    x[7]=_mm512_mask_mov_epi64(x[7],0xc0,b);
}
GAMMA30_INLINE void completeTranspose(__m512i* x) noexcept {
    const auto a=sr::repeat<2>(x[7]),b=sr::repeat<3>(x[7]);
    c8Unroll(std::make_index_sequence<8>{},[&](auto i){
        x[i]=_mm512_ternarylogic_epi64(x[i],tailProduct<0,i>(a),tailProduct<1,i>(b),0x96);
    });
    clearPadding(x);
}

// Eight four-record vectors, in natural coordinate order. Add just the
// size-32 layer to the existing fixed-16 transform and its exact adjoints.
template<unsigned Offset,Kind K>
GAMMA30_INLINE void fft32(__m512i* x) noexcept {
    static_assert(Offset%32==0);
    if constexpr(K==Kind::Forward||K==Kind::InverseTranspose){
        c8Unroll(std::make_index_sequence<4>{},[&](auto i){
            hdbw::detail::pair<hdbw::detail::normalized(4,Offset),K>(x[i],x[i+4]);
        });
    }
    hdbw::detail::transform16<Offset,K>(x[0],x[1],x[2],x[3]);
    hdbw::detail::transform16<Offset+16,K>(x[4],x[5],x[6],x[7]);
    if constexpr(K==Kind::Inverse||K==Kind::ForwardTranspose){
        c8Unroll(std::make_index_sequence<4>{},[&](auto i){
            hdbw::detail::pair<hdbw::detail::normalized(4,Offset),K>(x[i],x[i+4]);
        });
    }
}

template<unsigned Offset>
GAMMA30_INLINE void emitCoset(const State& coefficients,Packet& output) noexcept {
    static_assert(Offset==0||Offset==32);
    State values=coefficients;
    fft32<Offset,Kind::Forward>(values.values);
    c8Unroll(std::make_index_sequence<8>{},[&](auto i){
        output.values[Offset/4+i]=_mm512_xor_si512(output.values[Offset/4+i],values.values[i]);
    });
}
} // namespace detail

// Add A*q to a natural-order packet. A is exactly interpolation from
// nodes 128..157 to nodes 0..63, via completion to nodes 128..159.
GAMMA30_INLINE void emitAdd(const State& state,Packet& output) noexcept {
    State coefficients=state;
    detail::complete(coefficients.values);
    detail::fft32<128,Kind::Inverse>(coefficients.values);
    detail::emitCoset<0>(coefficients,output);
    detail::emitCoset<32>(coefficients,output);
}

// Consume the packet storage after it has been emitted. The sequence is
// precisely the transpose of emitAdd's linear A map, including completion.
GAMMA30_INLINE State feedback(Packet& input) noexcept {
    detail::fft32<0,Kind::ForwardTranspose>(input.values);
    detail::fft32<32,Kind::ForwardTranspose>(input.values+8);
    State result;
    c8Unroll(std::make_index_sequence<8>{},[&](auto i){
        result.values[i]=_mm512_xor_si512(input.values[i],input.values[i+8]);
    });
    detail::fft32<128,Kind::InverseTranspose>(result.values);
    detail::completeTranspose(result.values);
    return result;
}

// H30: values at 0..29 -> values at 30..59. The missing two values are also
// the first two outputs; nodes 32..59 come from the other size-32 coset.
GAMMA30_INLINE void refreshMap(State& state) noexcept {
    detail::complete(state.values);
    const auto tail=state.values[7];
    detail::fft32<0,Kind::Inverse>(state.values);
    detail::fft32<32,Kind::Forward>(state.values);
    // Descend so each previous vector remains available for the splice.
    c8Unroll(std::make_index_sequence<7>{},[&](auto index){
        constexpr unsigned i=7-index;
        state.values[i]=_mm512_shuffle_i32x4(state.values[i-1],state.values[i],0x4e);
    });
    state.values[0]=_mm512_shuffle_i32x4(tail,state.values[0],0x4e);
    detail::clearPadding(state.values);
}

GAMMA30_INLINE void finish(State& state,State syndrome,const Refresh& diagonal) noexcept {
    c8Unroll(std::make_index_sequence<8>{},[&](auto i){
        syndrome.values[i]=field::multiply4(
            _mm512_xor_si512(state.values[i],field::multiplyY4(syndrome.values[i])),
            diagonal.before.data()+4*i);
    });
    detail::clearPadding(syndrome.values);
    refreshMap(syndrome);
    c8Unroll(std::make_index_sequence<8>{},[&](auto i){
        state.values[i]=field::multiply4(syndrome.values[i],diagonal.after.data()+4*i);
    });
    detail::clearPadding(state.values);
}

GAMMA30_INLINE Packet loadPacket(const Block* input,const std::uint16_t* route) noexcept {
    Packet result;
    // Preserve the old paired route: lane l of innerInput4<J> is row16*l+J.
    // Four fixed transposes recover natural row order for the full-RS map.
    c8Unroll(std::make_index_sequence<4>{},[&](auto group){
        constexpr unsigned j=4*group;
        auto a=old::innerInput4<j>(input,route),b=old::innerInput4<j+1>(input,route);
        auto c=old::innerInput4<j+2>(input,route),d=old::innerInput4<j+3>(input,route);
        sr::transposeRecords4(a,b,c,d);
        result.values[group]=a;result.values[4+group]=b;
        result.values[8+group]=c;result.values[12+group]=d;
    });
    return result;
}
template<bool Stream>
GAMMA30_INLINE void storePacket(const Packet& packet,Block* output) noexcept {
    c8Unroll(std::make_index_sequence<4>{},[&](auto group){
        auto a=packet.values[group],b=packet.values[4+group];
        auto c=packet.values[8+group],d=packet.values[12+group];
        sr::transposeRecords4(a,b,c,d);
        if constexpr(Stream){
            _mm512_stream_si512(reinterpret_cast<__m512i*>(output+16*group),a);
            _mm512_stream_si512(reinterpret_cast<__m512i*>(output+16*group+4),b);
            _mm512_stream_si512(reinterpret_cast<__m512i*>(output+16*group+8),c);
            _mm512_stream_si512(reinterpret_cast<__m512i*>(output+16*group+12),d);
        }else{
            _mm512_storeu_si512(output+16*group,a);_mm512_storeu_si512(output+16*group+4,b);
            _mm512_storeu_si512(output+16*group+8,c);_mm512_storeu_si512(output+16*group+12,d);
        }
    });
}
template<bool Stream>
GAMMA30_INLINE void innerPacket(const Block* input,const std::uint16_t* route,Block* output,
                               State& state,const Refresh& diagonal) noexcept {
    Packet y=loadPacket(input,route);
    emitAdd(state,y);
    storePacket<Stream>(y,output);
    // C*y=C*x because C*A=0. Gamma multiplies only feedback, not old state.
    finish(state,feedback(y),diagonal);
}
template<bool Stream>
GAMMA30_NOINLINE void inner(const Plan& plan,const Block* planes,Block* output){
    const auto& p=plan.outerPlan;
    if constexpr(Stream)if(reinterpret_cast<std::uintptr_t>(output)&63U){inner<false>(plan,planes,output);return;}
    State state{};
    for(unsigned region=0;region!=256;++region){
        const auto begin=region*p.rows;const auto* source=planes+region*p.stride;
        for(std::size_t j=0;j!=p.rows;j+=64){
            if(region+1<256){
                const auto* next=reinterpret_cast<const char*>(planes+(region+1)*p.stride+j);
                for(unsigned h=0;h!=16;++h)_mm_prefetch(next+64*h,_MM_HINT_T1);
            }
            innerPacket<Stream>(source,p.packetRoute.data()+(begin+j)/2,output+begin+j,
                                state,plan.refresh[(begin+j)/64]);
        }
    }
    if constexpr(Stream)_mm_sfence();
}
template<bool Stream=true>
inline void forward(const Plan& plan,const Block* input,Block* output,Block* scratch){
    // Frozen outer/route, including completing outer before output mutation.
    // Exact/partial input-output alias is allowed; scratch must be disjoint.
    // Buffers hold k, 2*k, and scratchSize() Block records respectively.
    // Normal Block alignment suffices; NT stores fall back for addresses
    // that are not 64-byte aligned.
    // State starts at zero and persists across every region; no flush is added.
    literal::outer<true>(plan.outerPlan,input,scratch);
    inner<Stream>(plan,scratch,output);
}

#undef GAMMA30_INLINE
#undef GAMMA30_NOINLINE
} // namespace spin::research::char2::field128_lift::gamma30
