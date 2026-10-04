#include <algorithm>
#ifndef SPIN_WORKSPACE_ROUTING_OPT
#define SPIN_WORKSPACE_ROUTING_OPT 1
#endif
#include "WorkspaceRouting.h"
// Included by separate AVX2/AVX-512 translation units, with no runtime dispatch
// inside the XOR circuit or routing loops.
#include "Spin.h"
#include "generated/WideCircuit.h"
#include <stdexcept>
#include <type_traits>
#include "../../Cpu.h"
// The explicit-register helpers have their own ISA target. AVX2 callers can
// select them after one runtime capability check outside the encoding loops.
#if HC_WIDE_BITS == 256 && defined(__GNUC__) && SPIN_BCH_AVX512
#include "generated/BchRegisterSchedule.h"
#include "feedbackRegister.h"
#define SPIN_WIDE_REGISTERS 1
#endif

#if HC_WIDE_BITS == 256
#define HC_NS wide256
#define HC_NAME(x) x##256
#else
#define HC_NS wide512
#define HC_NAME(x) x##512
#endif

namespace spin::detail::kernel::HC_NS {
struct Ops {
#if HC_WIDE_BITS == 256
    using Vec=__m256i;
    static SPIN_FORCEINLINE Vec zero() { return _mm256_setzero_si256(); }
    static SPIN_FORCEINLINE Vec vx(Vec a,Vec b) { return _mm256_xor_si256(a,b); }
    static SPIN_FORCEINLINE Vec load(const Vec* p) { return _mm256_loadu_si256(p); }
    static SPIN_FORCEINLINE void store(Vec* p,Vec v) { _mm256_storeu_si256(p,v); }
    static SPIN_FORCEINLINE void stream(Vec* p,Vec v) { _mm256_stream_si256(p,v); }
#else
    using Vec=__m512i;
    static SPIN_FORCEINLINE Vec zero() { return _mm512_setzero_si512(); }
    static SPIN_FORCEINLINE Vec vx(Vec a,Vec b) { return _mm512_xor_si512(a,b); }
    static SPIN_FORCEINLINE Vec load(const Vec* p) { return _mm512_loadu_si512(p); }
    static SPIN_FORCEINLINE void store(Vec* p,Vec v) { _mm512_storeu_si512(p,v); }
    static SPIN_FORCEINLINE void stream(Vec* p,Vec v) { _mm512_stream_si512(p,v); }
#endif
};
using Vec=Ops::Vec;
using Circuit=WideCircuit<Ops>;

static_assert(sizeof(Vec)==HC_WIDE_BITS/8 && alignof(Vec)==HC_WIDE_BITS/8);
#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
#endif
struct Workspace {
    std::vector<Vec> buckets;
    std::vector<u32> direct_slots;
    bool registers=false;
    explicit Workspace(const Spin& code):buckets(code.codeBlocks()),direct_slots(code.codeBlocks()) {
        const auto view=code.wideForwardView();
        auto packed=[](const u8* p) {u32 x;std::memcpy(&x,p,4);return x&0xffffff;};
        // Compose the tile permutation and final gather once, including a
        // partial final tile and either packed-24 or full-32 setup indices.
        for(std::size_t i=0;i<view.n;++i) {
            const auto slot=view.slots?packed(view.slots+3*i):view.slots32[i];
            const auto offset=view.offsets?packed(view.offsets+3*slot):view.offsets32[slot];
            direct_slots[i]=static_cast<u32>((slot/view.tile)*view.tile+offset);
        }
#if defined(SPIN_WIDE_REGISTERS)
        registers=code.bchBackend()!=BchBackend::Avx2 && spin::detail::cpu_avx512vl();
#endif
        if(workspace_routing::eligible(true,view.n/2)) {
            workspace_routing::adviseOwned(buckets.data(),buckets.size()*sizeof(Vec));
        }
    }
};
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
#if defined(SPIN_WIDE_REGISTERS)
template<class Map> struct RegisterFeedback:Map {
    static SPIN_FORCEINLINE void feedback(const Vec* in,Vec* out) {
        if constexpr(std::is_same_v<Map,WideMap128S19<Ops>>) feedbackRegister(in,out);
        else Map::feedback(in,out);
    }
};
template<class Map,bool Stream> __attribute__((noinline)) static void registerInner(
    Spin::WideView view,Workspace& w,Vec* out) {
    Circuit::innerForward<RegisterFeedback<Map>,Stream>(view.n,view.fieldRows,[&](std::size_t i) {
        return w.buckets[w.direct_slots[i]];
    },out);
}
#endif
template<class Map,bool Stream> static SPIN_FORCEINLINE void inner(
    Spin::WideView view,Workspace& w,Vec* out) {
    Circuit::innerForward<Map,Stream>(view.n,view.fieldRows,[&](std::size_t i) {
        return w.buckets[w.direct_slots[i]];
    },out);
}
template<class Map,bool Registers> static void encode(Spin::WideView view,const Vec* in,Vec* out,Workspace& w,bool stream) {
    auto* values=w.buckets.data();
    for(std::size_t base=0;base<view.n;base+=view.tile) {
        const auto active=std::min(view.tile,view.n-base);
        for(std::size_t j=0;j<active;j+=256) {
#if defined(SPIN_WIDE_REGISTERS)
            // Input and workspace output are disjoint. The register schedule
            // may reuse already-written BCH output as dependency backing.
            if constexpr(Registers) bchForwardRegister(in+(base+j)/2,values+base+j);
            else
#endif
            Circuit::bchForward(in+(base+j)/2,values+base+j);
        }
    }
#if defined(SPIN_WIDE_REGISTERS)
    if constexpr(Registers) {
        if(stream) registerInner<Map,true>(view,w,out);
        else registerInner<Map,false>(view,w,out);
        return;
    }
#endif
    if(stream) inner<Map,true>(view,w,out);
    else inner<Map,false>(view,w,out);
}
template<class Map> static void dispatch(Spin::WideView v,const Vec* in,Vec* out,Workspace& w,bool stream) {
#if defined(SPIN_WIDE_REGISTERS)
    if(w.registers) {encode<Map,true>(v,in,out,w,stream);return;}
#endif
    encode<Map,false>(v,in,out,w,stream);
}

}

extern "C" void* HC_NAME(spin_internal_wide_workspace)(const void* code) noexcept {
    try { return new spin::detail::kernel::HC_NS::Workspace(*static_cast<const spin::detail::kernel::Spin*>(code)); }
    catch(...) { return nullptr; }
}
extern "C" void HC_NAME(spin_internal_wide_destroy)(void* work) noexcept {
    delete static_cast<spin::detail::kernel::HC_NS::Workspace*>(work);
}
extern "C" std::size_t HC_NAME(spin_internal_wide_bytes)(const void* work) noexcept {
    const auto& w=*static_cast<const spin::detail::kernel::HC_NS::Workspace*>(work);
    return w.buckets.capacity()*(HC_WIDE_BITS/8)+w.direct_slots.capacity()*sizeof(std::uint32_t);
}
extern "C" int HC_NAME(spin_internal_wide_encode)(const void* code,void* work,const void* in,
                                            std::size_t ni,void* out,std::size_t no,bool stream) noexcept {
    try {
        const auto v=static_cast<const spin::detail::kernel::Spin*>(code)->wideForwardView();
        auto& w=*static_cast<spin::detail::kernel::HC_NS::Workspace*>(work);
        constexpr auto lanes=HC_WIDE_BITS/128;
        if(!in || !out || ni!=v.n/2*lanes || no!=v.n*lanes || w.buckets.size()!=v.n || w.direct_slots.size()!=v.n)
            return 1;
        const auto a=reinterpret_cast<std::uintptr_t>(in), b=reinterpret_cast<std::uintptr_t>(out);
        if(a<=b ? b-a<ni*16 : a-b<no*16) return 1;
        using namespace spin::detail::kernel;
        using Ops=HC_NS::Ops;
        const auto* input=static_cast<const HC_NS::Vec*>(in);auto* output=static_cast<HC_NS::Vec*>(out);
        // Streaming stores require vector alignment. Unaligned caller buffers
        // use the cached storeu circuit, selected once before encoding.
        stream=stream && (b&(HC_WIDE_BITS/8-1))==0;
        if(v.configuration==Configuration::T64S12R2) HC_NS::dispatch<WideMap64S12R2<Ops>>(v,input,output,w,stream);
        else HC_NS::dispatch<WideMap128S19<Ops>>(v,input,output,w,stream);
        return 0;
    } catch(...) { return 1; }
}
