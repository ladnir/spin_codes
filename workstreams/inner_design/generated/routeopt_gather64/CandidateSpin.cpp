#include "Spin.h"
#include "Inner.h"
#include "AsymmetricMap.h"
#include "../../asymmetric/AsymmetricInner.h"
#include "generated/BchCircuit.h"
#include "QuarterCircuit.h"
#include <algorithm>
#include <array>
#include <bit>
#include <cstring>
#include <numeric>
#include <stdexcept>

namespace bare_spin {
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
static OC_FORCEINLINE u32 unpack(const u8* p) { u32 v; std::memcpy(&v,p,4); return v&0xffffff; }

Spin::Spin(Configuration c,unsigned exponent,u64 routeSeed,u64 coefficientSeed,unsigned tileRows,Outer outer)
    :mConfig(c),mOuter(outer),mK(0),mTileRows(tileRows) {
    if(exponent<16 || exponent>20 || (tileRows && (tileRows<2 || !std::has_single_bit(tileRows))))
        throw std::invalid_argument("supported message exponents: 16..20; tile rows must be a power of two >=2");
    if(static_cast<unsigned>(c)>static_cast<unsigned>(Configuration::T256S14))
        throw std::invalid_argument("unknown configuration");
    if(outer!=Outer::Bch256x128 && outer!=Outer::Bch128x32)
        throw std::invalid_argument("unknown outer");
    if(outer==Outer::Bch128x32 && c!=Configuration::T128S19)
        throw std::invalid_argument("quarter-rate implementation selects t128_s19");
    mK=std::size_t{1}<<exponent;
    const auto b=outerLength();
    const std::size_t n=codeBlocks(),rows=mK/outerDimension();
    const auto defaultTile=outer==Outer::Bch128x32?
        (exponent<=16?256U:exponent<=18?2048U:4096U):(exponent<=18?256U:2048U);
    const auto selectedTile=tileRows?tileRows:defaultTile;
    mTileRows=static_cast<unsigned>(std::min<std::size_t>(selectedTile,rows));
    if(rows%step()) throw std::invalid_argument("incomplete region epochs");
    std::vector<u8> regionCoordinate(n);
    std::array<u8,256> coordinateRegion;
    for(std::size_t row=0;row<rows;++row) {
        shuffle(coordinateRegion.data(),b,routeSeed);
        for(unsigned j=0;j<b;++j) regionCoordinate[coordinateRegion[j]*rows+row]=u8(j);
    }
    mRoute.resize(n);
    std::vector<u32> positions(rows);
    for(unsigned region=0;region<b;++region) {
        shuffle(positions.data(),rows,routeSeed);
        for(std::size_t row=0;row<rows;++row)
            mRoute[region*rows+positions[row]]=u32(row*b+regionCoordinate[region*rows+row]);
    }
    mSlots24.resize(3*n+4); mOffsets24.resize(3*n+4);
    mSlots32.resize(n); mOffsets32.resize(n);
    const auto tile=tileBlocks();
    std::vector<u32> counts(n/tile,u32(tile));
    for(std::size_t inner=n;inner-->0;) {
        const u32 outer=mRoute[inner],bucket=outer/u32(tile),local=outer&u32(tile-1);
        const u32 slot=bucket*u32(tile)+--counts[bucket];
        mSlots32[inner]=slot;
        const auto offsetIndex=mOuter==Outer::Bch128x32?outer:slot;
        const auto offsetValue=mOuter==Outer::Bch128x32?slot-bucket*u32(tile):local;
        mOffsets32[offsetIndex]=offsetValue;
        pack(mSlots24.data()+3*inner,slot);pack(mOffsets24.data()+3*offsetIndex,offsetValue);
    }
    switch(c) {
        case Configuration::T64S16: setupInner<Map64S16>(coefficientSeed); break;
        case Configuration::T64S20: setupInner<Map64S20>(coefficientSeed); break;
        case Configuration::T128S19:
            if(mOuter==Outer::Bch128x32) setupInner<AsymmetricMap>(coefficientSeed);
            else setupInner<Map128S19>(coefficientSeed);
            break;
        case Configuration::T256S14: setupInner<Map256S14>(coefficientSeed); break;
    }
}
template<class Map> void Spin::setupInner(u64 seed) {
    const auto epochs=codeBlocks()/Map::T;
    if constexpr(Map::T==128 && Map::S==19) {
        if(mOuter==Outer::Bch128x32) {
            mCoefficients.resize(2*epochs);mFieldRows.resize(2*epochs);
            for(std::size_t e=0;e<epochs;++e) {
                u32 u;
                do u=u32(splitmix(seed))&((1U<<Map::S)-1); while(!u);
                u32 v=u32(splitmix(seed))&((1U<<Map::S)-1);
                if(std::popcount(u&v)&1) v^=u&-u;
                mCoefficients[2*e]=u;mCoefficients[2*e+1]=v;
                u32 grouped=0;
                if constexpr(SPIN_GROUPED_A) {
                    for(unsigned j=0;j<Map::S;++j) grouped|=((u>>Map::groupOrder[j])&1)<<j;
                } else grouped=u;
                mFieldRows[2*e]=grouped;mFieldRows[2*e+1]=v;
            }
            return;
        }
    }

    mCoefficients.resize(epochs); mFieldRows.resize(epochs*Map::S);
    for(std::size_t e=0;e<epochs;++e) {
        u32 a;
        do a=u32(splitmix(seed))&((1U<<Map::S)-1); while(!a);
        mCoefficients[e]=a;
        for(unsigned j=0;j<Map::S;++j) mFieldRows[e*Map::S+j]=fieldMultiply<Map>(1U<<j,a);
    }
}
Spin::Workspace::Workspace(const Spin& code):buckets(code.codeBlocks()),tile(code.tileBlocks()) {}
std::size_t Spin::Workspace::bytes() const noexcept { return (buckets.capacity()+tile.capacity())*sizeof(block); }
unsigned Spin::step() const noexcept {
    switch(mConfig) { case Configuration::T64S16: case Configuration::T64S20:return 64;
        case Configuration::T128S19:return 128; default:return 256; }
}
unsigned Spin::state() const noexcept {
    switch(mConfig) { case Configuration::T64S16:return 16; case Configuration::T64S20:return 20;
        case Configuration::T128S19:return 19; default:return 14; }
}
const char* Spin::name() const noexcept {
    if(mOuter==Outer::Bch128x32) return "routeopt_gather64";
    switch(mConfig) { case Configuration::T64S16:return "t64_s16"; case Configuration::T64S20:return "t64_s20";
        case Configuration::T128S19:return "t128_s19"; default:return "t256_s14"; }
}
std::size_t Spin::setupBytes() const noexcept {
    return mSlots24.capacity()+mOffsets24.capacity()+4*(mSlots32.capacity()+mOffsets32.capacity()
        +mRoute.capacity()+mCoefficients.capacity()+mFieldRows.capacity());
}
u64 Spin::routeHash() const noexcept {
    u64 h=0xcbf29ce484222325ULL;
    for(u32 v:mRoute) { h^=v;h*=0x100000001b3ULL; }
    return h;
}
void Spin::encode(const block* in,std::size_t ni,block* out,std::size_t no,Workspace& w,Layout layout) const {
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
    if(layout!=Layout::Packed24 && layout!=Layout::Indices32) throw std::invalid_argument("unknown layout");
    if(mCompacted && layout!=mRetainedLayout) throw std::invalid_argument("routing layout was discarded");
    if(layout==Layout::Packed24) {std::vector<u32>().swap(mSlots32);std::vector<u32>().swap(mOffsets32);}
    else {std::vector<u8>().swap(mSlots24);std::vector<u8>().swap(mOffsets24);}
    std::vector<u32>().swap(mRoute);std::vector<u32>().swap(mCoefficients);
    mRetainedLayout=layout;mCompacted=true;
}

void Spin::encodeInplace(block* buffer,std::size_t count,Workspace& w,Layout layout) const {
    if(!buffer || count!=codeBlocks() || w.buckets.size()!=codeBlocks() || w.tile.size()!=tileBlocks())
        throw std::invalid_argument("buffer or workspace geometry mismatch");
    if(layout!=Layout::Packed24 && layout!=Layout::Indices32) throw std::invalid_argument("unknown layout");
    if(mCompacted && layout!=mRetainedLayout) throw std::invalid_argument("routing layout was discarded");
    encodeUnchecked(buffer,buffer,w,layout);
}
template<class Map,bool Packed,bool Quarter> void Spin::run(const block* in,block* out,Workspace& w) const {
    block* values=w.buckets.data();
    asymmetricReverse<Map,Quarter,false>(in,codeBlocks(),mFieldRows.data(),[&](std::size_t i,block v) {
        if constexpr(Packed) values[unpack(mSlots24.data()+3*i)]=v;
        else values[mSlots32[i]]=v;
    });
    block* tile=w.tile.data();
    const auto tileSize=tileBlocks(),n=codeBlocks();
    for(std::size_t base=0;base<n;base+=tileSize) {
        if constexpr(Quarter) {
            for(std::size_t j=0;j<tileSize;j+=256) {
                for(unsigned p=0;p<256;p+=8) {
                    auto point=[&]<unsigned K>() {
                        const auto pos=base+j+p+K;
                        if(j+p+K+64<tileSize) {
                            const auto off=Packed?unpack(mOffsets24.data()+3*(pos+64)):mOffsets32[pos+64];
                            _mm_prefetch(reinterpret_cast<const char*>(values+base+off),_MM_HINT_T0);
                        }
                        const auto off=Packed?unpack(mOffsets24.data()+3*pos):mOffsets32[pos];
                        tile[p+K]=values[base+off];
                    };
                    point.template operator()<0>();
                    point.template operator()<1>();
                    point.template operator()<2>();
                    point.template operator()<3>();
                    point.template operator()<4>();
                    point.template operator()<5>();
                    point.template operator()<6>();
                    point.template operator()<7>();
                }
                quarterTranspose2(tile,tile+128,out+base/4+j/4,out+base/4+j/4+32);
            }
            continue;
        }
        for(std::size_t j=0;j<tileSize;++j) {
            if(j+32<tileSize) {
                const auto future=Packed?unpack(mOffsets24.data()+3*(base+j+32)):mOffsets32[base+j+32];
                _mm_prefetch(reinterpret_cast<const char*>(tile+future),_MM_HINT_T0);
            }
            const auto offset=Packed?unpack(mOffsets24.data()+3*(base+j)):mOffsets32[base+j];
            tile[offset]=values[base+j];
        }
        if constexpr(Quarter) {
            for(std::size_t j=0;j<tileSize;j+=256)
                quarterTranspose2(tile+j,tile+j+128,out+base/4+j/4,out+base/4+j/4+32);
        } else {
            for(std::size_t j=0;j<tileSize;j+=512)
                bchTranspose2(tile+j,tile+j+256,out+base/2+j/2,out+base/2+j/2+128);
        }
    }
}
void Spin::encodeUnchecked(const block* in,block* out,Workspace& w,Layout layout) const {
    if(mOuter==Outer::Bch128x32) {
        if(layout==Layout::Packed24) run<AsymmetricMap,true,true>(in,out,w);
        else run<AsymmetricMap,false,true>(in,out,w);
        return;
    }
#define RUN_CASE(Enum,Map) case Configuration::Enum: if(layout==Layout::Packed24) run<Map,true,false>(in,out,w); else run<Map,false,false>(in,out,w); break
    switch(mConfig) {
        RUN_CASE(T64S16,Map64S16); RUN_CASE(T64S20,Map64S20);
        RUN_CASE(T128S19,Map128S19); RUN_CASE(T256S14,Map256S14);
    }
#undef RUN_CASE
}
void Spin::validateSetup() const {
    if(!mCompacted && mOuter==Outer::Bch128x32) {
        for(std::size_t e=0;e<codeBlocks()/step();++e) {
            const auto u=mCoefficients[2*e],v=mCoefficients[2*e+1];
            if(!u || (u>>19) || (v>>19) || (std::popcount(u&v)&1))
                throw std::runtime_error("invalid transvection sample");
        }
    }
    if(mCompacted) throw std::logic_error("oracle schedules were discarded");
    std::vector<bool> seen(codeBlocks()),slots(codeBlocks());
    for(std::size_t i=0;i<codeBlocks();++i) {
        const auto outer=mRoute[i],slot=mSlots32[i];
        if(outer>=seen.size() || seen[outer] || slot>=slots.size() || slots[slot]) throw std::runtime_error("invalid route");
        seen[outer]=true; slots[slot]=true;
        const auto index=mOuter==Outer::Bch128x32?outer:slot;
        const auto target=mOuter==Outer::Bch128x32?slot:outer;
        if(unpack(mSlots24.data()+3*i)!=slot || unpack(mOffsets24.data()+3*index)!=mOffsets32[index] ||
           (slot/tileBlocks())*tileBlocks()+mOffsets32[index]!=target)
            throw std::runtime_error("packed route mismatch");
    }
    const auto rows=mK/outerDimension();
    std::vector<bool> rowRegion(codeBlocks());
    for(std::size_t i=0;i<codeBlocks();++i) {
        const auto key=(mRoute[i]/outerLength())*outerLength()+i/rows;
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
template<class Map,bool Quarter> void Spin::oracle(const block* in,block* out) const {
    const auto n=codeBlocks(),epochs=n/Map::T;
    std::vector<block> routed(n);
    std::array<block,Map::S> state{},syndrome{},next{};
    for(auto& v:state) v=block(0,0);
    for(std::size_t e=epochs;e-->0;) {
        for(auto& v:syndrome) v=block(0,0);
        for(unsigned p=0;p<Map::T;++p) {
            block value=in[e*Map::T+p];
            for(unsigned j=0;j<Map::S;++j) if((Map::columns[p]>>j)&1) value^=state[j];
            routed[mRoute[e*Map::T+p]]=value;
            if constexpr(Quarter) {
                for(unsigned j=0;j<Map::S;++j) if((Map::feedbackColumns[p]>>j)&1) syndrome[j]^=in[e*Map::T+p];
            } else {
                for(unsigned j=0;j<Map::S;++j) if((Map::columns[p]>>j)&1) syndrome[j]^=value;
            }
        }
        for(unsigned j=0;j<Map::S;++j) {
            block v=syndrome[j];
            u32 mask;
            if constexpr(Quarter) {
                const u32 u=mCoefficients[2*e],v=mCoefficients[2*e+1];
                mask=1U<<j;
                if((v>>j)&1) mask^=u;
            } else mask=referenceMultiply<Map>(1U<<j,mCoefficients[e]);
            for(unsigned b=0;b<Map::S;++b) if((mask>>b)&1) v^=state[b];
            next[j]=v;
        }
        state=next;
    }
    constexpr unsigned B=Quarter?128:256,D=Quarter?32:128;
    for(std::size_t r=0;r<mK/D;++r) for(unsigned j=0;j<D;++j) {
        block v(0,0);
        for(unsigned c=0;c<B;++c) {
            const auto word=Quarter?QuarterRows[j][c/64]:BchRows[j][c/64];
            if((word>>(c%64))&1) v^=routed[r*B+c];
        }
        out[r*D+j]=v;
    }
}
void Spin::reference(const block* in,block* out) const {
    if(mCompacted) throw std::logic_error("oracle schedules were discarded");
    if(mOuter==Outer::Bch128x32) {oracle<AsymmetricMap,true>(in,out);return;}
    switch(mConfig) {
        case Configuration::T64S16:oracle<Map64S16,false>(in,out);break;
        case Configuration::T64S20:oracle<Map64S20,false>(in,out);break;
        case Configuration::T128S19:oracle<Map128S19,false>(in,out);break;
        case Configuration::T256S14:oracle<Map256S14,false>(in,out);break;
    }
}
}
