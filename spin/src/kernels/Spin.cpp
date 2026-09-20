#ifndef SPIN_WORKSPACE_ROUTING_OPT
#define SPIN_WORKSPACE_ROUTING_OPT 1
#endif
#include "WorkspaceRouting.h"
#include "Spin.h"
#include "Inner.h"
#include "K16Inner.h"
#include "K16R2Inner.h"
#include "generated/BchCircuit.h"
#include "generated/BchForward.h"
#include <algorithm>
#include <array>
#include <bit>
#include <cstring>
#include <numeric>
#include <stdexcept>

namespace spin::detail::kernel {
bool bchAvx512Available() noexcept {
#if SPIN_BCH_AVX512 && !SPIN_TEST_NO_AVX512
    return cpu_avx512f() && cpu_avx512vl();
#else
    return false;
#endif
}
static u32 packFour(u32 x) {return (x&~1023U)|((x&255U)<<2)|((x>>8)&3U);}

u64 splitmix(u64& s) noexcept {
    u64 v=(s+=0x9e3779b97f4a7c15ULL);
    v=(v^(v>>30))*0xbf58476d1ce4e5b9ULL;
    v=(v^(v>>27))*0x94d049bb133111ebULL;
    return v^(v>>31);
}
static u64 below(u64& s,u64 bound) {
    const u64 threshold=(0-bound)%bound;
    for(;;) { const auto v=splitmix(s); if(v>=threshold) return v%bound; }
}
template<class T> static void shuffle(T* p,std::size_t n,u64& seed) {
    std::iota(p,p+n,T{0});
    for(std::size_t i=n;i>1;--i) std::swap(p[i-1],p[below(seed,i)]);
}
static void pack(u8* p,u32 v) { p[0]=u8(v); p[1]=u8(v>>8); p[2]=u8(v>>16); }
static SPIN_FORCEINLINE u32 unpack(const u8* p) { u32 v; std::memcpy(&v,p,4); return v&0xffffff; }

Spin::Spin(Configuration c,unsigned exponent,u64 routeSeed,u64 coefficientSeed,unsigned tileRows,BchBackend backend)
    :Spin(c,MessageLength{length_geometry::fromExponent(exponent)},routeSeed,coefficientSeed,tileRows,backend) {}
Spin::Spin(Configuration c,MessageLength length,u64 routeSeed,u64 coefficientSeed,unsigned tileRows,BchBackend backend)
    :mConfig(c),mK(0),mTileRows(tileRows) {
    if((tileRows && (tileRows<2 || !std::has_single_bit(tileRows))))
        throw std::invalid_argument("tile rows must be a power of two");
    if(static_cast<unsigned>(c)>static_cast<unsigned>(Configuration::T64S12R2))
        throw std::invalid_argument("unknown configuration");
    mK=length.value;
    length_geometry::check(mK,step());
    const std::size_t n=codeBlocks(),rows=mK/128;
    const auto selectedTile=tileRows?tileRows:(mK<=(1U<<18)?256U:2048U);
    mTileRows=static_cast<unsigned>(std::min<std::size_t>(selectedTile,std::bit_floor(rows)));
    if(backend!=BchBackend::Auto && backend!=BchBackend::Avx2 && backend!=BchBackend::Avx512)
        throw std::invalid_argument("unknown BCH backend");
    const bool four=(c==Configuration::T128S19 || (c==Configuration::T64S12 || c==Configuration::T64S12R2)) && mTileRows>=4 && bchAvx512Available();
    if(backend==BchBackend::Avx512 && !four) throw std::invalid_argument("unsupported AVX-512 BCH configuration");
    mBchBackend=backend!=BchBackend::Avx2 && four?BchBackend::Avx512:BchBackend::Avx2;
    if(mBchBackend==BchBackend::Avx512) {if(packed24Available()) mBchOffsets24.resize(3*n+4);mBchOffsets32.resize(n);}
    mPartialTile=(rows%mTileRows)!=0;
    if(rows%step()) throw std::invalid_argument("incomplete region epochs");
    std::vector<u8> regionCoordinate(n);
    std::array<u8,256> coordinateRegion;
    for(std::size_t row=0;row<rows;++row) {
        shuffle(coordinateRegion.data(),256,routeSeed);
        for(unsigned j=0;j<256;++j) regionCoordinate[coordinateRegion[j]*rows+row]=u8(j);
    }
    mRoute.resize(n);
    std::vector<u32> positions(rows);
    for(unsigned region=0;region<256;++region) {
        shuffle(positions.data(),rows,routeSeed);
        for(std::size_t row=0;row<rows;++row)
            mRoute[region*rows+positions[row]]=u32(row*256+regionCoordinate[region*rows+row]);
    }
    if(mK==(1U<<16) && mBchBackend==BchBackend::Avx512) {
        mSmallRoute32.resize(n);
        for(std::size_t i=0;i<n;++i) mSmallRoute32[i]=packFour(mRoute[i]);
    }
    if(mK==(1U<<18) && mBchBackend==BchBackend::Avx512) {
        mK18Route32.resize(n);
        for(std::size_t i=0;i<n;++i) mK18Route32[i]=packFour(mRoute[i]);
    }
    if(mK<=(1U<<18) && mBchBackend==BchBackend::Avx512) {
        mForwardDirect.resize(n);
        for(std::size_t i=0;i<n;++i) mForwardDirect[i]=packFour(mRoute[i]);
    }
    // Standalone transpose's measured range-direct schedule. Reuse the
    // forward direct table where present, without extending forward dispatch.
    if(mK<=458752 && mBchBackend==BchBackend::Avx512 && mForwardDirect.empty()) {
        mRangeRoute32.resize(n);
        for(std::size_t i=0;i<n;++i) mRangeRoute32[i]=packFour(mRoute[i]);
    }
    if(packed24Available()) {mSlots24.resize(3*n+4); mOffsets24.resize(3*n+4);}
    mSlots32.resize(n); mOffsets32.resize(n);
    const auto tile=tileBlocks();
    std::vector<u32> counts((n+tile-1)/tile,u32(tile));
    if(n%tile) counts.back()=u32(n%tile);
    for(std::size_t inner=n;inner-->0;) {
        const u32 outer=mRoute[inner],bucket=outer/u32(tile),local=outer&u32(tile-1);
        const u32 slot=bucket*u32(tile)+--counts[bucket];
        mSlots32[inner]=slot; mOffsets32[slot]=local;
        if(packed24Available()) {pack(mSlots24.data()+3*inner,slot); pack(mOffsets24.data()+3*slot,local);}
        if(mBchBackend==BchBackend::Avx512) {
            mBchOffsets32[slot]=packFour(local);if(packed24Available()) pack(mBchOffsets24.data()+3*slot,packFour(local));
        }
    }
    switch(c) {
        case Configuration::T64S16: setupInner<Map64S16>(coefficientSeed); break;
        case Configuration::T64S20: setupInner<Map64S20>(coefficientSeed); break;
        case Configuration::T64S12R2: setupInner<Map64S12R2>(coefficientSeed); break;
        case Configuration::T64S12: setupInner<Map64S12>(coefficientSeed); break;
        case Configuration::T128S19: setupInner<Map128S19>(coefficientSeed); break;
        case Configuration::T256S14: setupInner<Map256S14>(coefficientSeed); break;
    }
}
template<class Map> void Spin::setupInner(u64 seed) {
    const auto epochs=codeBlocks()/Map::T;
    constexpr unsigned rounds=isTwoRoundMap<Map>?2:1;
    if constexpr(isImtMap<Map>) {
        mCoefficients.resize(2*rounds*epochs);mFieldRows.resize(2*rounds*epochs);mForwardFieldRows.resize(2*rounds*epochs);
        for(std::size_t e=0;e<rounds*epochs;++e) {
            u32 u;do u=u32(splitmix(seed))&((1U<<Map::S)-1);while(!u);
            u32 v=u32(splitmix(seed))&((1U<<Map::S)-1);
            if(std::popcount(u&v)&1) v^=u&-u;
            mCoefficients[2*e]=mFieldRows[2*e]=mForwardFieldRows[2*e]=u;
            mCoefficients[2*e+1]=mFieldRows[2*e+1]=mForwardFieldRows[2*e+1]=v;
        }
        return;
    }
    mCoefficients.resize(epochs); mFieldRows.resize(epochs*Map::S);
    mForwardFieldRows.resize(epochs*Map::S);
    for(std::size_t e=0;e<epochs;++e) {
        u32 a;
        do a=u32(splitmix(seed))&((1U<<Map::S)-1); while(!a);
        mCoefficients[e]=a;
        for(unsigned j=0;j<Map::S;++j) mFieldRows[e*Map::S+j]=fieldMultiply<Map>(1U<<j,a);
        for(unsigned j=0;j<Map::S;++j) {
            u32 mask=0;
            for(unsigned k=0;k<Map::S;++k) mask|=((mFieldRows[e*Map::S+k]>>j)&1U)<<k;
            mForwardFieldRows[e*Map::S+j]=mask;
        }
    }
}
Spin::Workspace::Workspace(const Spin& code):buckets(code.codeBlocks()),tile(code.tileBlocks()) {
    if(workspace_routing::eligible(true,code.messageBlocks())) {
        workspace_routing::adviseOwned(buckets.data(),buckets.size()*sizeof(block));
        workspace_routing::adviseOwned(tile.data(),tile.size()*sizeof(block));
    }
}
Spin::WideView Spin::wideForwardView() const {
    if((mConfig!=Configuration::T128S19 && mConfig!=Configuration::T64S12R2) || (mSlots24.empty() && mSlots32.empty()))
        throw std::invalid_argument("wide forward requires S19 or K16 R2 routing");
    return {codeBlocks(),tileBlocks(),mSlots24.empty()?nullptr:mSlots24.data(),mOffsets24.data(),mForwardFieldRows.data(),mConfig,mSlots32.data(),mOffsets32.data()};
}
Spin::WideView Spin::wideView() const {
    if(mConfig!=Configuration::T128S19 || mSlots24.empty())
        throw std::invalid_argument("wide encoding requires the S19 packed schedule");
    return {codeBlocks(),tileBlocks(),mSlots24.data(),mOffsets24.data(),mForwardFieldRows.data()};
}
std::size_t Spin::Workspace::bytes() const noexcept { return (buckets.capacity()+tile.capacity())*sizeof(block); }
unsigned Spin::step() const noexcept {
    switch(mConfig) { case Configuration::T64S16: case Configuration::T64S20:case Configuration::T64S12:case Configuration::T64S12R2:return 64;
        case Configuration::T128S19:return 128; default:return 256; }
}
unsigned Spin::state() const noexcept {
    switch(mConfig) { case Configuration::T64S12:case Configuration::T64S12R2:return 12; case Configuration::T64S16:return 16; case Configuration::T64S20:return 20;
        case Configuration::T128S19:return 19; default:return 14; }
}
const char* Spin::name() const noexcept {
    if(mConfig==Configuration::T64S12R2) return "imt_t64_s12_subspace_v1_r2";
    if(mConfig==Configuration::T64S12) return Map64S12::name;
    switch(mConfig) { case Configuration::T64S16:return "t64_s16"; case Configuration::T64S20:return "t64_s20";
        case Configuration::T128S19:return "imt_t128_s19_weight5"; default:return "t256_s14"; }
}
std::size_t Spin::setupBytes() const noexcept {
    return 4*mRangeRoute32.capacity()+4*mForwardDirect.capacity()+mK18Route32.capacity()*sizeof(u32)+mSmallRoute32.capacity()*sizeof(u32)+mBchOffsets24.capacity()+4*mBchOffsets32.capacity()+mSlots24.capacity()+mOffsets24.capacity()+4*(mSlots32.capacity()+mOffsets32.capacity()
        +mRoute.capacity()+mCoefficients.capacity()+mFieldRows.capacity()+mForwardFieldRows.capacity());
}
u64 Spin::routeHash() const noexcept {
    u64 h=0xcbf29ce484222325ULL;
    for(u32 v:mRoute) { h^=v;h*=0x100000001b3ULL; }
    return h;
}
void Spin::encode(const block* in,std::size_t ni,block* out,std::size_t no,Workspace& w,Layout layout) const {
    if(layout==Layout::Auto) layout=preferredLayout();
    if(layout==Layout::Packed24 && !packed24Available()) throw std::invalid_argument("Packed24 index range exceeded; use Auto or Indices32");

    if(!in || !out || ni!=codeBlocks() || no!=messageBlocks() ||
       w.buckets.size()!=codeBlocks() || w.tile.size()!=tileBlocks())
        throw std::invalid_argument("input, output, or workspace geometry mismatch");
    const auto a=reinterpret_cast<std::uintptr_t>(in),b=reinterpret_cast<std::uintptr_t>(out);
    if(a<=b ? b-a<ni*sizeof(block) : a-b<no*sizeof(block))
        throw std::invalid_argument("input and output must not overlap");
    if(layout!=Layout::Packed24 && layout!=Layout::Indices32) throw std::invalid_argument("unknown layout");
    if(mCompacted && layout!=mRetainedLayout) throw std::invalid_argument("routing layout was discarded");
    encodeUnchecked(in,out,w,layout);
}
void Spin::compact(Layout layout) {
    if(layout==Layout::Auto) layout=preferredLayout();
    if(layout==Layout::Packed24 && !packed24Available()) throw std::invalid_argument("Packed24 index range exceeded; use Auto or Indices32");

    if(layout!=Layout::Packed24 && layout!=Layout::Indices32) throw std::invalid_argument("unknown layout");
    if(mCompacted && layout!=mRetainedLayout) throw std::invalid_argument("routing layout was discarded");
    if(layout==Layout::Packed24) {std::vector<u32>().swap(mSlots32);std::vector<u32>().swap(mOffsets32);}
    else {std::vector<u8>().swap(mSlots24);std::vector<u8>().swap(mOffsets24);}
    if(layout==Layout::Packed24) std::vector<u32>().swap(mBchOffsets32);
    else std::vector<u8>().swap(mBchOffsets24);
    // The packed single-row path uses a direct bit permutation, not SIMD buckets.
    if(mConfig!=Configuration::T128S19 && (mConfig!=Configuration::T64S12 && mConfig!=Configuration::T64S12R2)) std::vector<u32>().swap(mRoute);
    std::vector<u32>().swap(mCoefficients);
    if(!mSmallRoute32.empty()) {
        /* Retained for four-row forward routing. */
        /* Retained for four-row forward routing. */
    }
    if(!mK18Route32.empty()) {
        /* Retained for four-row forward routing. */
        /* Retained for four-row forward routing. */
    }
    mRetainedLayout=layout;mCompacted=true;
}
template<class Map,bool Packed> void Spin::run(const block* in,block* out,Workspace& w) const {
    block* values=w.buckets.data();
    innerReverse<Map>(in,codeBlocks(),mFieldRows.data(),[&](std::size_t i,block v) {
        if constexpr(Packed) values[unpack(mSlots24.data()+3*i)]=v;
        else values[mSlots32[i]]=v;
    });
    block* tile=w.tile.data();
    const auto tileSize=tileBlocks(),n=codeBlocks();
    for(std::size_t base=0;base<n;base+=tileSize) {
        for(std::size_t j=0;j<tileSize;++j) {
            if(j+32<tileSize) {
                const auto future=Packed?unpack(mOffsets24.data()+3*(base+j+32)):mOffsets32[base+j+32];
                _mm_prefetch(reinterpret_cast<const char*>(tile+future),_MM_HINT_T0);
            }
            const auto offset=Packed?unpack(mOffsets24.data()+3*(base+j)):mOffsets32[base+j];
            tile[offset]=values[base+j];
        }
        for(std::size_t j=0;j<tileSize;j+=512)
            bchTranspose2(tile+j,tile+j+256,out+base/2+j/2,out+base/2+j/2+128);
    }
}
template<class Map,bool Packed> void Spin::runTail(const block* in,block* out,Workspace& w) const {
    block* values=w.buckets.data();
    innerReverse<Map>(in,codeBlocks(),mFieldRows.data(),[&](std::size_t i,block v) {
        if constexpr(Packed) values[unpack(mSlots24.data()+3*i)]=v;
        else values[mSlots32[i]]=v;
    });
    block* tile=w.tile.data();
    const auto tileSize=tileBlocks(),n=codeBlocks();
    for(std::size_t base=0;base<n;base+=tileSize) {
        const auto activeSize=std::min(tileSize,n-base);
        for(std::size_t j=0;j<activeSize;++j) {
            if(j+32<activeSize) {
                const auto future=Packed?unpack(mOffsets24.data()+3*(base+j+32)):mOffsets32[base+j+32];
                _mm_prefetch(reinterpret_cast<const char*>(tile+future),_MM_HINT_T0);
            }
            const auto offset=Packed?unpack(mOffsets24.data()+3*(base+j)):mOffsets32[base+j];
            tile[offset]=values[base+j];
        }
        for(std::size_t j=0;j<activeSize;j+=512)
            bchTranspose2(tile+j,tile+j+256,out+base/2+j/2,out+base/2+j/2+128);
    }
}
void Spin::encodeUnchecked(const block* in,block* out,Workspace& w,Layout layout) const {
#if SPIN_BCH_AVX512
    if(mK<=458752 && mBchBackend==BchBackend::Avx512 && mK!=(1U<<16) &&
       (mConfig!=Configuration::T128S19 || mK!=(1U<<18))) {
        switch(mConfig) {
        case Configuration::T128S19:runRange<Map128S19>(in,out,w);return;
        case Configuration::T64S12:runRange<Map64S12>(in,out,w);return;
        case Configuration::T64S12R2:runRange<Map64S12R2>(in,out,w);return;
        default:break;
        }
    }
#endif

    if(layout==Layout::Auto) layout=preferredLayout();
    if(mPartialTile) {
#if SPIN_BCH_AVX512
        if(mBchBackend==BchBackend::Avx512) {
            if(mConfig==Configuration::T128S19) {if(layout==Layout::Packed24) runFourTail<true>(in,out,w);else runFourTail<false>(in,out,w);return;}
            if(mConfig==Configuration::T64S12) {if(layout==Layout::Packed24) runK16FourTail<true>(in,out,w);else runK16FourTail<false>(in,out,w);return;}
            if(mConfig==Configuration::T64S12R2) {if(layout==Layout::Packed24) runK16R2FourTail<true>(in,out,w);else runK16R2FourTail<false>(in,out,w);return;}
        }
#endif
        switch(mConfig) {
            case Configuration::T64S16:if(layout==Layout::Packed24) runTail<Map64S16,true>(in,out,w);else runTail<Map64S16,false>(in,out,w);return;
            case Configuration::T64S20:if(layout==Layout::Packed24) runTail<Map64S20,true>(in,out,w);else runTail<Map64S20,false>(in,out,w);return;
            case Configuration::T128S19:if(layout==Layout::Packed24) runTail<Map128S19,true>(in,out,w);else runTail<Map128S19,false>(in,out,w);return;
            case Configuration::T256S14:if(layout==Layout::Packed24) runTail<Map256S14,true>(in,out,w);else runTail<Map256S14,false>(in,out,w);return;
            case Configuration::T64S12:if(layout==Layout::Packed24) runTail<Map64S12,true>(in,out,w);else runTail<Map64S12,false>(in,out,w);return;
            case Configuration::T64S12R2:if(layout==Layout::Packed24) runTail<Map64S12R2,true>(in,out,w);else runTail<Map64S12R2,false>(in,out,w);return;
        }
    }

#if SPIN_BCH_AVX512
    if(mBchBackend==BchBackend::Avx512) {
        if(mConfig==Configuration::T128S19 && mK==(1U<<18)) {runK18(in,out,w);return;}
        if(mConfig==Configuration::T64S12R2 && mK==(1U<<16)) {runSmallK16R2(in,out,w);return;}
        if(mConfig==Configuration::T64S12 && mK==(1U<<16)) {runSmallK16(in,out,w);return;}
        if(mK==(1U<<16)) {runSmall(in,out,w);return;}
        if(mConfig==Configuration::T64S12) {if(layout==Layout::Packed24) runK16Four<true>(in,out,w);else runK16Four<false>(in,out,w);return;}
        if(mConfig==Configuration::T64S12R2) {if(layout==Layout::Packed24) runK16R2Four<true>(in,out,w);else runK16R2Four<false>(in,out,w);return;}
        if(layout==Layout::Packed24) runFour<true>(in,out,w); else runFour<false>(in,out,w);
        return;
    }
#endif
#define RUN_CASE(Enum,Map) case Configuration::Enum: if(layout==Layout::Packed24) run<Map,true>(in,out,w); else run<Map,false>(in,out,w); break
    switch(mConfig) {
        RUN_CASE(T64S16,Map64S16); RUN_CASE(T64S20,Map64S20);
        RUN_CASE(T64S12R2,Map64S12R2); RUN_CASE(T64S12,Map64S12); RUN_CASE(T128S19,Map128S19); RUN_CASE(T256S14,Map256S14);
    }
#undef RUN_CASE
}
void Spin::encodeInplace(block* buffer,std::size_t count,Workspace& w,Layout layout) const {
    if(layout==Layout::Auto) layout=preferredLayout();
    if(layout==Layout::Packed24 && !packed24Available()) throw std::invalid_argument("Packed24 index range exceeded; use Auto or Indices32");

    if(!buffer || count!=codeBlocks() || w.buckets.size()!=codeBlocks() || w.tile.size()!=tileBlocks())
        throw std::invalid_argument("buffer or workspace geometry mismatch");
    if(layout!=Layout::Packed24 && layout!=Layout::Indices32) throw std::invalid_argument("unknown layout");
    if(mCompacted && layout!=mRetainedLayout) throw std::invalid_argument("routing layout was discarded");
    encodeUnchecked(buffer,buffer,w,layout);
}
void Spin::validateSetup() const {
    if(!mCompacted && !mForwardDirect.empty()) {
        for(std::size_t i=0;i<codeBlocks();++i)
            if(mForwardDirect[i]!=packFour(mRoute[i])) throw std::runtime_error("forward direct route mismatch");
    }
    if(!mCompacted && !mK18Route32.empty()) {
        for(std::size_t i=0;i<codeBlocks();++i)
            if(mK18Route32[i]!=packFour(mRoute[i])) throw std::runtime_error("K18 direct route mismatch");
    }
    if(!mCompacted && (mConfig==Configuration::T64S12 || mConfig==Configuration::T64S12R2)) {
        const auto count=2*(mConfig==Configuration::T64S12R2?2:1)*codeBlocks()/64;
        if(mCoefficients.size()!=count || mFieldRows.size()!=count || mForwardFieldRows.size()!=count)
            throw std::runtime_error("invalid IMT schedule length");
        for(std::size_t i=0;i<count;i+=2) {
            const auto u=mCoefficients[i],v=mCoefficients[i+1];
            if(!u || (u>>12) || (v>>12) || (std::popcount(u&v)&1) ||
               mFieldRows[i]!=u || mFieldRows[i+1]!=v || mForwardFieldRows[i]!=u || mForwardFieldRows[i+1]!=v)
                throw std::runtime_error("invalid IMT schedule sample");
        }
    }
    if(!mCompacted && !mSmallRoute32.empty()) {
        for(std::size_t i=0;i<codeBlocks();++i)
            if(mSmallRoute32[i]!=packFour(mRoute[i]))
                throw std::runtime_error("small direct route mismatch");
    }

    if(mCompacted) throw std::logic_error("oracle schedules were discarded");
    std::vector<bool> seen(codeBlocks()),slots(codeBlocks());
    for(std::size_t i=0;i<codeBlocks();++i) {
        const auto outer=mRoute[i],slot=mSlots32[i];
        if(mBchBackend==BchBackend::Avx512 &&
           (slot>=mBchOffsets32.size() || mBchOffsets32[slot]!=packFour(mOffsets32[slot]) ||
            (packed24Available() && unpack(mBchOffsets24.data()+3*slot)!=mBchOffsets32[slot])))
            throw std::runtime_error("packed BCH offset mismatch");
        if(outer>=seen.size() || seen[outer] || slot>=slots.size() || slots[slot]) throw std::runtime_error("invalid route");
        seen[outer]=true; slots[slot]=true;
        if((packed24Available() && (unpack(mSlots24.data()+3*i)!=slot || unpack(mOffsets24.data()+3*slot)!=mOffsets32[slot])) ||
           (slot/tileBlocks())*tileBlocks()+mOffsets32[slot]!=outer)
            throw std::runtime_error("packed route mismatch");
    }
    const auto rows=mK/128;
    std::vector<bool> rowRegion(codeBlocks());
    for(std::size_t i=0;i<codeBlocks();++i) {
        const auto key=(mRoute[i]/256)*256+i/rows;
        if(rowRegion[key]) throw std::runtime_error("route has duplicate row/region");
        rowRegion[key]=true;
    }
}

// Independent polynomial multiplication for the reference (not the shift kernel).
template<class Map> static u32 referenceMultiply(u32 a,u32 b) {
    u64 product=0;
    for(unsigned j=0;j<Map::S;++j) if((b>>j)&1) product^=u64(a)<<j;
    const u64 modulus=(u64{1}<<Map::S)|Map::modulusLow;
    for(int j=2*Map::S-2;j>=int(Map::S);--j) if((product>>j)&1) product^=modulus<<(j-Map::S);
    return u32(product);
}
template<class Map> void Spin::oracle(const block* in,block* out) const {
    const auto n=codeBlocks(),epochs=n/Map::T;
    std::vector<block> routed(n);
    std::array<block,Map::S> state{},syndrome{},next{};
    for(auto& v:state) v=block(0,0);
    for(std::size_t e=epochs;e-->0;) {
        for(auto& v:syndrome) v=block(0,0);
        for(unsigned p=0;p<Map::T;++p) {
            block value=in[e*Map::T+p];
            for(unsigned j=0;j<Map::S;++j) {
                u32 column=Map::columns[p];
                if constexpr(isImtMap<Map>) column=Map::feedbackColumns[p];
                if((column>>j)&1) value^=state[j];
            }
            routed[mRoute[e*Map::T+p]]=value;
            for(unsigned j=0;j<Map::S;++j) if((Map::columns[p]>>j)&1) {
                if constexpr(isImtMap<Map>) syndrome[j]^=in[e*Map::T+p];
                else syndrome[j]^=value;
            }
        }
        for(unsigned j=0;j<Map::S;++j) {
            block v=syndrome[j];
            u32 mask;
            if constexpr(isTwoRoundMap<Map>) mask=imtR2ReferenceRow(j,mCoefficients.data()+4*e);
            else if constexpr(isImtMap<Map>) mask=(1U<<j)^(((mCoefficients[2*e+1]>>j)&1)?mCoefficients[2*e]:0);
            else mask=referenceMultiply<Map>(1U<<j,mCoefficients[e]);
            for(unsigned b=0;b<Map::S;++b) if((mask>>b)&1) v^=state[b];
            next[j]=v;
        }
        state=next;
    }
    for(std::size_t r=0;r<mK/128;++r) for(unsigned j=0;j<128;++j) {
        block v(0,0);
        for(unsigned c=0;c<256;++c) if((BchRows[j][c/64]>>(c%64))&1) v^=routed[r*256+c];
        out[r*128+j]=v;
    }
}
void Spin::reference(const block* in,block* out) const {
    if(mCompacted) throw std::logic_error("oracle schedules were discarded");
    switch(mConfig) {
        case Configuration::T64S16:oracle<Map64S16>(in,out);break;
        case Configuration::T64S20:oracle<Map64S20>(in,out);break;
        case Configuration::T64S12R2:oracle<Map64S12R2>(in,out);break;
        case Configuration::T64S12:oracle<Map64S12>(in,out);break;
        case Configuration::T128S19:oracle<Map128S19>(in,out);break;
        case Configuration::T256S14:oracle<Map256S14>(in,out);break;
    }
}

void Spin::forward(const block* in,std::size_t ni,block* out,std::size_t no,Workspace& w,Layout layout) const {
    if(layout==Layout::Auto) layout=preferredLayout();
    if(layout==Layout::Packed24 && !packed24Available()) throw std::invalid_argument("Packed24 index range exceeded; use Auto or Indices32");

    if(!in || !out || ni!=messageBlocks() || no!=codeBlocks() ||
       w.buckets.size()!=codeBlocks() || w.tile.size()!=tileBlocks())
        throw std::invalid_argument("forward input, output, or workspace geometry mismatch");
    const auto a=reinterpret_cast<std::uintptr_t>(in),b=reinterpret_cast<std::uintptr_t>(out);
    if(a<=b ? b-a<ni*sizeof(block) : a-b<no*sizeof(block))
        throw std::invalid_argument("input and output must not overlap");
    if(layout!=Layout::Packed24 && layout!=Layout::Indices32) throw std::invalid_argument("unknown layout");
    if(mCompacted && layout!=mRetainedLayout) throw std::invalid_argument("routing layout was discarded");
    forwardUnchecked(in,out,w,layout);
}
template<class Map,bool Packed> void Spin::runForward(const block* in,block* out,Workspace& w) const {
    block* values=w.buckets.data(); block* tile=w.tile.data();
    const auto tileSize=tileBlocks(),n=codeBlocks();
    for(std::size_t base=0;base<n;base+=tileSize) {
        for(std::size_t j=0;j<tileSize;j+=512)
            bchForward2(in+base/2+j/2,in+base/2+j/2+128,tile+j,tile+j+256);
        for(std::size_t j=0;j<tileSize;++j) {
            const auto offset=Packed?unpack(mOffsets24.data()+3*(base+j)):mOffsets32[base+j];
            values[base+j]=tile[offset];
        }
    }
    innerForward<Map>(n,mForwardFieldRows.data(),[&](std::size_t i) {
        if constexpr(Packed) return values[unpack(mSlots24.data()+3*i)];
        else return values[mSlots32[i]];
    },out);
}
template<class Map,bool Packed> void Spin::runForwardTail(const block* in,block* out,Workspace& w) const {
    block* values=w.buckets.data(); block* tile=w.tile.data();
    const auto tileSize=tileBlocks(),n=codeBlocks();
    for(std::size_t base=0;base<n;base+=tileSize) {
        const auto activeSize=std::min(tileSize,n-base);
        for(std::size_t j=0;j<activeSize;j+=512)
            bchForward2(in+base/2+j/2,in+base/2+j/2+128,tile+j,tile+j+256);
        for(std::size_t j=0;j<activeSize;++j) {
            const auto offset=Packed?unpack(mOffsets24.data()+3*(base+j)):mOffsets32[base+j];
            values[base+j]=tile[offset];
        }
    }
    innerForward<Map>(n,mForwardFieldRows.data(),[&](std::size_t i) {
        if constexpr(Packed) return values[unpack(mSlots24.data()+3*i)];
        else return values[mSlots32[i]];
    },out);
}
void Spin::forwardUnchecked(const block* in,block* out,Workspace& w,Layout layout) const {
    if(layout==Layout::Auto) layout=preferredLayout();
    if(mPartialTile && mForwardDirect.empty()) {
#if SPIN_BCH_AVX512
        if(mBchBackend==BchBackend::Avx512) {
            if(mConfig==Configuration::T128S19) {if(layout==Layout::Packed24) runForwardFourTail<Map128S19,true>(in,out,w);else runForwardFourTail<Map128S19,false>(in,out,w);return;}
            if(mConfig==Configuration::T64S12) {if(layout==Layout::Packed24) runForwardFourTail<Map64S12,true>(in,out,w);else runForwardFourTail<Map64S12,false>(in,out,w);return;}
            if(mConfig==Configuration::T64S12R2) {if(layout==Layout::Packed24) runForwardFourTail<Map64S12R2,true>(in,out,w);else runForwardFourTail<Map64S12R2,false>(in,out,w);return;}
        }
#endif
        switch(mConfig) {
            case Configuration::T64S16:if(layout==Layout::Packed24) runForwardTail<Map64S16,true>(in,out,w);else runForwardTail<Map64S16,false>(in,out,w);return;
            case Configuration::T64S20:if(layout==Layout::Packed24) runForwardTail<Map64S20,true>(in,out,w);else runForwardTail<Map64S20,false>(in,out,w);return;
            case Configuration::T128S19:if(layout==Layout::Packed24) runForwardTail<Map128S19,true>(in,out,w);else runForwardTail<Map128S19,false>(in,out,w);return;
            case Configuration::T256S14:if(layout==Layout::Packed24) runForwardTail<Map256S14,true>(in,out,w);else runForwardTail<Map256S14,false>(in,out,w);return;
            case Configuration::T64S12:if(layout==Layout::Packed24) runForwardTail<Map64S12,true>(in,out,w);else runForwardTail<Map64S12,false>(in,out,w);return;
            case Configuration::T64S12R2:if(layout==Layout::Packed24) runForwardTail<Map64S12R2,true>(in,out,w);else runForwardTail<Map64S12R2,false>(in,out,w);return;
        }
    }

#if SPIN_BCH_AVX512
    if(mBchBackend==BchBackend::Avx512) {
#define FAST_FORWARD(Enum,Map) case Configuration::Enum: if(!mForwardDirect.empty()) runForwardDirect<Map>(in,out,w); else if(layout==Layout::Packed24) runForwardFour<Map,true>(in,out,w); else runForwardFour<Map,false>(in,out,w); return
        switch(mConfig) {
            FAST_FORWARD(T128S19,Map128S19);FAST_FORWARD(T64S12,Map64S12);FAST_FORWARD(T64S12R2,Map64S12R2);
            default:break;
        }
#undef FAST_FORWARD
    }
#endif
#define FORWARD_CASE(Enum,Map) case Configuration::Enum: if(layout==Layout::Packed24) runForward<Map,true>(in,out,w); else runForward<Map,false>(in,out,w); break
    switch(mConfig) {
        FORWARD_CASE(T64S16,Map64S16); FORWARD_CASE(T64S20,Map64S20);
        FORWARD_CASE(T64S12R2,Map64S12R2); FORWARD_CASE(T64S12,Map64S12); FORWARD_CASE(T128S19,Map128S19); FORWARD_CASE(T256S14,Map256S14);
    }
#undef FORWARD_CASE
}
template<class Map> void Spin::forwardOracle(const block* in,block* out) const {
    const auto n=codeBlocks(),epochs=n/Map::T;
    std::vector<block> outer(n,block(0,0));
    for(std::size_t r=0;r<mK/128;++r) for(unsigned j=0;j<128;++j)
        for(unsigned c=0;c<256;++c) if((BchRows[j][c/64]>>(c%64))&1)
            outer[r*256+c]^=in[r*128+j];
    std::array<block,Map::S> state{},syndrome{},next{};
    for(auto& v:state) v=block(0,0);
    for(std::size_t e=0;e<epochs;++e) {
        for(auto& v:syndrome) v=block(0,0);
        for(unsigned p=0;p<Map::T;++p) {
            const block x=outer[mRoute[e*Map::T+p]];
            block y=x;
            for(unsigned j=0;j<Map::S;++j) if((Map::columns[p]>>j)&1) {
                y^=state[j];
                if constexpr(!isImtMap<Map>) syndrome[j]^=x;
            }
            if constexpr(isImtMap<Map>) for(unsigned j=0;j<Map::S;++j)
                if((Map::feedbackColumns[p]>>j)&1) syndrome[j]^=x;
            out[e*Map::T+p]=y;
        }
        next=syndrome;
        for(unsigned j=0;j<Map::S;++j) {
            u32 mask;
            if constexpr(isTwoRoundMap<Map>) mask=imtR2ReferenceRow(j,mCoefficients.data()+4*e);
            else if constexpr(isImtMap<Map>) mask=(1U<<j)^(((mCoefficients[2*e+1]>>j)&1)?mCoefficients[2*e]:0);
            else mask=referenceMultiply<Map>(1U<<j,mCoefficients[e]);
            for(unsigned k=0;k<Map::S;++k) if((mask>>k)&1) next[k]^=state[j];
        }
        state=next;
    }
}
void Spin::forwardReference(const block* in,block* out) const {
    if(mCompacted) throw std::logic_error("oracle schedules were discarded");
    switch(mConfig) {
        case Configuration::T64S16:forwardOracle<Map64S16>(in,out);break;
        case Configuration::T64S20:forwardOracle<Map64S20>(in,out);break;
        case Configuration::T64S12R2:forwardOracle<Map64S12R2>(in,out);break;
        case Configuration::T64S12:forwardOracle<Map64S12>(in,out);break;
        case Configuration::T128S19:forwardOracle<Map128S19>(in,out);break;
        case Configuration::T256S14:forwardOracle<Map256S14>(in,out);break;
    }
}

void Spin::forwardBits(const u64* in,std::size_t ni,u64* out,std::size_t no,
                       u64* scratch,std::size_t ns) const {
    if((mConfig!=Configuration::T128S19 && (mConfig!=Configuration::T64S12 && mConfig!=Configuration::T64S12R2)) || !in || !out || !scratch ||
       ni!=mK/64 || no!=codeBlocks()/64 || ns!=no)
        throw std::invalid_argument("packed bit encoder geometry mismatch");
    const auto overlaps=[](const void* p,std::size_t pn,const void* q,std::size_t qn) {
        const auto a=reinterpret_cast<std::uintptr_t>(p),b=reinterpret_cast<std::uintptr_t>(q);
        return a<=b ? b-a<pn : a-b<qn;
    };
    if(overlaps(in,ni*8,out,no*8) || overlaps(in,ni*8,scratch,ns*8) || overlaps(out,no*8,scratch,ns*8))
        throw std::invalid_argument("packed bit buffers overlap");
    if(mConfig==Configuration::T64S12R2) forwardBitsMap<Map64S12R2>(in,out,scratch);
    else if(mConfig==Configuration::T64S12) forwardBitsMap<Map64S12>(in,out,scratch);
    else forwardBitsMap<Map128S19>(in,out,scratch);
}
template<class Map> void Spin::forwardBitsMap(const u64* in,u64* out,u64* scratch) const {
    constexpr unsigned S=Map::S,T=Map::T,W=T/64,G=(S+3)/4;
    // 16 KiB of four-bit BCH tables. This is the same generator as the SIMD DAG.
    static const auto bch=[] {
        std::array<std::array<std::array<u64,4>,16>,32> t{};
        for(unsigned g=0;g<32;++g) for(unsigned b=0;b<4;++b)
            for(unsigned x=0;x<(1U<<b);++x) for(unsigned w=0;w<4;++w)
                t[g][x+(1U<<b)][w]=t[g][x][w]^BchRows[4*g+b][w];
        return t;
    }();
    static const auto masks=[] {
        std::array<std::array<u64,W>,S> t{};
        for(unsigned j=0;j<S;++j) for(unsigned p=0;p<T;++p)
            t[j][p/64]|=u64((Map::columns[p]>>j)&1)<<(p%64);
        return t;
    }();
    static const auto feedback=[] {
        std::array<std::array<std::array<u64,W>,16>,G> t{};
        for(unsigned g=0;g<G;++g) for(unsigned b=0;b<4 && 4*g+b<S;++b)
            for(unsigned x=0;x<(1U<<b);++x) for(unsigned w=0;w<W;++w)
                t[g][x+(1U<<b)][w]=t[g][x][w]^masks[4*g+b][w];
        return t;
    }();
    static const auto feedbackMasks=[] {
        std::array<std::array<u64,W>,S> t{};
        for(unsigned j=0;j<S;++j) for(unsigned p=0;p<T;++p)
            t[j][p/64]|=u64((Map::feedbackColumns[p]>>j)&1)<<(p%64);
        return t;
    }();
    for(std::size_t r=0;r<mK/128;++r) {
        u64 a=0,b=0,c=0,d=0;
        for(unsigned g=0;g<32;++g) {
            const auto& t=bch[g][(in[2*r+g/16]>>(4*(g%16)))&15];
            a^=t[0]; b^=t[1]; c^=t[2]; d^=t[3];
        }
        scratch[4*r]=a; scratch[4*r+1]=b; scratch[4*r+2]=c; scratch[4*r+3]=d;
    }
    u32 state=0;
    for(std::size_t base=0;base<codeBlocks();base+=T) {
        u64 x0=0,x1=0;
        for(unsigned p=0;p<64;++p) {
            const auto a=mRoute[base+p];
            x0|=((scratch[a/64]>>(a%64))&1)<<p;
            if constexpr(W==2) {const auto b=mRoute[base+64+p];x1|=((scratch[b/64]>>(b%64))&1)<<p;}
        }
        u64 y0=x0,y1=x1;
        for(unsigned g=0;g<G;++g) {
            const auto& t=feedback[g][(state>>(4*g))&15]; y0^=t[0]; if constexpr(W==2) y1^=t[1];
        }
        out[base/64]=y0; if constexpr(W==2) out[base/64+1]=y1;
        if(base+T==codeBlocks()) break;
        u32 mixed=state;
        if constexpr(isTwoRoundMap<Map>) {
            const auto* pair=mForwardFieldRows.data()+4*(base/T);
            mixed^=(0U-u32(std::popcount(mixed&pair[1])&1))&pair[0];
            mixed^=(0U-u32(std::popcount(mixed&pair[3])&1))&pair[2];
        }
        u32 next=0;
        for(unsigned j=0;j<S;++j) {
            const auto parity=std::popcount(x0&feedbackMasks[j][0])^([&] {if constexpr(W==2) return std::popcount(x1&feedbackMasks[j][1]); else return 0;}())^
                ([&] {if constexpr(isTwoRoundMap<Map>) return int((mixed>>j)&1);
                else return int((state>>j)&1)^int(((mForwardFieldRows[2*(base/T)]>>j)&1) &
                (std::popcount(state&mForwardFieldRows[2*(base/T)+1])&1));}());
            next|=(u32(parity)&1)<<j;
        }
        state=next;
    }
}
}
