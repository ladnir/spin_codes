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
#else
    using Vec=__m512i;
    static SPIN_FORCEINLINE Vec zero() { return _mm512_setzero_si512(); }
    static SPIN_FORCEINLINE Vec vx(Vec a,Vec b) { return _mm512_xor_si512(a,b); }
    static SPIN_FORCEINLINE Vec load(const Vec* p) { return _mm512_loadu_si512(p); }
    static SPIN_FORCEINLINE void store(Vec* p,Vec v) { _mm512_storeu_si512(p,v); }
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
    std::vector<Vec> buckets,tile;
    explicit Workspace(Spin::WideView view):buckets(view.n),tile(view.tile) {
        if(workspace_routing::eligible(true,view.n/2)) {
            workspace_routing::adviseOwned(buckets.data(),buckets.size()*sizeof(Vec));
            workspace_routing::adviseOwned(tile.data(),tile.size()*sizeof(Vec));
        }
    }
};
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
static SPIN_FORCEINLINE u32 unpack(const u8* p) {
    u32 x; std::memcpy(&x,p,4); return x&0xffffff;
}
template<class Map,bool Packed> static void encode(Spin::WideView view,const Vec* in,Vec* out,Workspace& w) {
    auto* values=w.buckets.data(); auto* tile=w.tile.data();
    for(std::size_t base=0;base<view.n;base+=view.tile) {
        // One vector is now an entire 256/512-row record. Apply the identical
        // BCH DAG to one outer group at a time, without half-vector packing.
        for(std::size_t j=0;j<view.tile;j+=256)
            Circuit::bchForward(in+(base+j)/2,tile+j);
        for(std::size_t j=0;j<view.tile;++j)
            values[base+j]=tile[Packed?unpack(view.offsets+3*(base+j)):view.offsets32[base+j]];
    }
    Circuit::innerForward<Map>(view.n,view.fieldRows,[&](std::size_t i) {
        // Each routing lookup loads one complete 32/64-byte row record.
        return values[Packed?unpack(view.slots+3*i):view.slots32[i]];
    },out);
}
template<class Map,bool Packed> static void encodeTail(Spin::WideView view,const Vec* in,Vec* out,Workspace& w) {
    auto* values=w.buckets.data(); auto* tile=w.tile.data();
    for(std::size_t base=0;base<view.n;base+=view.tile) {
        const auto activeSize=std::min(view.tile,view.n-base);
        // One vector is now an entire 256/512-row record. Apply the identical
        // BCH DAG to one outer group at a time, without half-vector packing.
        for(std::size_t j=0;j<activeSize;j+=256)
            Circuit::bchForward(in+(base+j)/2,tile+j);
        for(std::size_t j=0;j<activeSize;++j)
            values[base+j]=tile[Packed?unpack(view.offsets+3*(base+j)):view.offsets32[base+j]];
    }
    Circuit::innerForward<Map>(view.n,view.fieldRows,[&](std::size_t i) {
        // Each routing lookup loads one complete 32/64-byte row record.
        return values[Packed?unpack(view.slots+3*i):view.slots32[i]];
    },out);
}
template<class Map> static void dispatch(Spin::WideView v,const Vec* in,Vec* out,Workspace& w) {
    if(v.n%v.tile) {
        if(v.slots) encodeTail<Map,true>(v,in,out,w);else encodeTail<Map,false>(v,in,out,w);
    } else {
        if(v.slots) encode<Map,true>(v,in,out,w);else encode<Map,false>(v,in,out,w);
    }
}

}

extern "C" void* HC_NAME(spin_internal_wide_workspace)(const void* code) noexcept {
    try { return new spin::detail::kernel::HC_NS::Workspace(static_cast<const spin::detail::kernel::Spin*>(code)->wideForwardView()); }
    catch(...) { return nullptr; }
}
extern "C" void HC_NAME(spin_internal_wide_destroy)(void* work) noexcept {
    delete static_cast<spin::detail::kernel::HC_NS::Workspace*>(work);
}
extern "C" std::size_t HC_NAME(spin_internal_wide_bytes)(const void* work) noexcept {
    const auto& w=*static_cast<const spin::detail::kernel::HC_NS::Workspace*>(work);
    return (w.buckets.capacity()+w.tile.capacity())*(HC_WIDE_BITS/8);
}
extern "C" int HC_NAME(spin_internal_wide_encode)(const void* code,void* work,const void* in,
                                            std::size_t ni,void* out,std::size_t no) noexcept {
    try {
        const auto v=static_cast<const spin::detail::kernel::Spin*>(code)->wideForwardView();
        auto& w=*static_cast<spin::detail::kernel::HC_NS::Workspace*>(work);
        constexpr auto lanes=HC_WIDE_BITS/128;
        if(!in || !out || ni!=v.n/2*lanes || no!=v.n*lanes || w.buckets.size()!=v.n || w.tile.size()!=v.tile)
            return 1;
        const auto a=reinterpret_cast<std::uintptr_t>(in), b=reinterpret_cast<std::uintptr_t>(out);
        if(a<=b ? b-a<ni*16 : a-b<no*16) return 1;
        using namespace spin::detail::kernel;
        using Ops=HC_NS::Ops;
        const auto* input=static_cast<const HC_NS::Vec*>(in);auto* output=static_cast<HC_NS::Vec*>(out);
        if(v.configuration==Configuration::T64S12R2) HC_NS::dispatch<WideMap64S12R2<Ops>>(v,input,output,w);
        else HC_NS::dispatch<WideMap128S19<Ops>>(v,input,output,w);
        return 0;
    } catch(...) { return 1; }
}
