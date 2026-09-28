#include "Spin.h"
#include "Inner.h"
#include "SetupRandom.h"
#include "WorkspaceRouting.h"
#include "generated/BchCircuit.h"
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
namespace k=spin::detail::kernel;
namespace spin::detail::kernel {
// Research-only map: same fixed circuits, two independent mixing updates.
struct LocalMap128S19R2 : Map128S19 {
    template<class Emit> static SPIN_FORCEINLINE void emitShared(
        const block* in,__m128i* raw,const __m128i* state,std::size_t base,Emit& emit) {
        imtReversePoints<Map128S19>(in,raw,state,base,emit,std::make_index_sequence<128>{});
    }
};
template<> inline constexpr bool isTwoRoundMap<LocalMap128S19R2> = true;
void bchTranspose4(const block*,block*);
void bchTranspose4Mapped(const block*,block*,const unsigned*);
void bchTranspose4Stride16(const block*,block*);
void bchTranspose4Stride20(const block*,block*);
void bchTranspose4Greedy(const block*,block*);
void bchTranspose4Release(const block*,block*);
void bchTranspose4Reverse(const block*,block*);
void bchTranspose4Restrict(const block*,block*);
void bchTranspose4Share6(const block*,block*);
void bchTranspose4Share8(const block*,block*);
void bchTranspose4Gfni1(const block*,block*);
void bchTranspose4Gfni2(const block*,block*);
void bchTranspose4GfniStatic(const block*,block*);
void bchTranspose4GfniTernary(const block*,block*);
void bchTranspose4GfniUnroll(const block*,block*);
void bchTranspose4GfniBlend(const block*,block*);
void bchTranspose4GfniIndependent(const block*,block*,const std::uint16_t*);
void bchTranspose4GfniMapped(const block*,block*,const unsigned*);
void bchTranspose4GfniMapped1(const block*,block*,const unsigned*);
void bchTranspose4GfniPacked(const block*,block*,const unsigned*);
void bchTranspose4GfniPrepared(const block*,block*,const unsigned*);
}
using Clock=std::chrono::steady_clock;

struct Sampler {
    k::setup::Words words;
    std::vector<k::setup::Divisor> divisors;
    Sampler(unsigned n,std::uint64_t seed):words(seed),divisors(std::max(n,256U)+1) {
        for(unsigned i=2;i<divisors.size();++i)divisors[i]=k::setup::Divisor(i);
    }
    template<class T> void shuffle(T& v) {
        std::iota(v.begin(),v.end(),0);
        for(unsigned i=unsigned(v.size());i>1;--i)std::swap(v[i-1],v[divisors[i].sample(words,i)]);
    }
};
// Inverse route: inner position -> outer coordinate. Fixed contiguous row groups.
static std::vector<unsigned> route(unsigned rows,unsigned group,bool shared,std::uint64_t seed,unsigned columns=1,bool globalColumns=false,bool xorLanes=false,bool preserveBlocks=false) {
    if(!group || rows%group || !columns || 256%columns)throw std::invalid_argument("invalid rectangular geometry");
    Sampler random(rows,seed);
    std::vector<unsigned> coordinates(rows*256),out(rows*256),groups(rows/group),lanes(group);
    std::array<unsigned,256> row;
    for(unsigned r=0;r<rows;++r) {
        if(globalColumns ? r==0 : (!shared || r%group==0)) {
            if(preserveBlocks) {
                std::vector<unsigned> blocks(256/columns),within(columns);
                random.shuffle(blocks);
                for(unsigned b=0;b<blocks.size();++b) {
                    random.shuffle(within);
                    for(unsigned c=0;c<columns;++c)row[b*columns+c]=blocks[b]*columns+within[c];
                }
            } else random.shuffle(row);
        }
        for(unsigned c=0;c<256;++c)coordinates[row[c]*rows+r]=c;
    }
    for(unsigned region=0;region<256;region+=columns) {
        random.shuffle(groups);
        for(unsigned b=0;b<rows/group;++b) {
            for(unsigned c=0;c<columns;++c) {
                if(xorLanes){const unsigned shift=random.words()&(group-1);for(unsigned l=0;l<group;++l)lanes[l]=l^shift;}
                else random.shuffle(lanes);
                for(unsigned l=0;l<group;++l) {
                    const auto r=b*group+l;
                    out[region*rows+groups[b]*group*columns+c*group+lanes[l]]=r*256+coordinates[(region+c)*rows+r];
                }
            }
        }
    }
    return out;
}
static void checkRoute(const std::vector<unsigned>& r,unsigned rows,unsigned columns=1) {
    std::vector<bool> seen(r.size());
    for(unsigned region=0;region<256;region+=columns) {
        std::vector<unsigned> rowSeen(rows);
        for(unsigned p=0;p<rows*columns;++p) {
            const auto x=r[region*rows+p];
            if(x>=r.size() || seen[x])throw std::runtime_error("route bijection");
            seen[x]=true;++rowSeen[x/256];
        }
        for(auto count:rowSeen)if(count!=columns)throw std::runtime_error("macroregion occupancy");
    }
}
// Each shared G-row group contributes a 16-bit chunk to each of G/16 bands
// within a region. Band group permutations are independent; lane assignment
// is a fresh uniform permutation of G rows in each region and group.
static std::vector<unsigned> splitRoute(unsigned rows,unsigned group,std::uint64_t seed) {
    if((group!=32 && group!=64 && group!=128) || rows%group)throw std::invalid_argument("split geometry");
    const unsigned count=rows/group,bands=group/16,bandSize=rows/bands;
    Sampler random(rows,seed);
    std::vector<unsigned> out(rows*256),columns(count*256),positions(count*bands),permutation(count),lanes(group);
    std::array<unsigned,256> coordinate;
    for(unsigned b=0;b<count;++b){random.shuffle(coordinate);for(unsigned c=0;c<256;++c)columns[b*256+c]=coordinate[c];}
    for(unsigned region=0;region<256;++region) {
        for(unsigned band=0;band<bands;++band) {
            random.shuffle(permutation);
            for(unsigned b=0;b<count;++b)positions[band*count+b]=permutation[b];
        }
        for(unsigned b=0;b<count;++b) {
            random.shuffle(lanes);
            for(unsigned band=0;band<bands;++band)for(unsigned l=0;l<16;++l)
                out[region*rows+band*bandSize+16*positions[band*count+b]+l]=(b*group+lanes[band*16+l])*256+columns[b*256+region];
        }
    }
    return out;
}
static unsigned packed(unsigned x) {return (x&~1023U)|((x&255U)<<2)|((x>>8)&3U);}
// Keep the existing unrolled IMT recurrence. Only the destination schedule changes.
static SPIN_NOINLINE void direct(const k::block* in,k::block* out,k::block* scratch,
                               std::size_t n,const unsigned* offsets,const unsigned* masks) {
    k::innerReverse<k::Map128S19>(in,n,masks,[&](std::size_t i,k::block v){scratch[offsets[i]]=v;});
    for(std::size_t i=0;i<n;i+=1024)k::bchTranspose4(scratch+i,out+i/2);
}
// Cost diagnostic only: the route below makes packed destinations sequential.
static SPIN_NOINLINE void sequentialControl(const k::block* in,k::block* out,k::block* scratch,
                               std::size_t n,const unsigned* masks,double* phases=nullptr) {
    auto begin=Clock::time_point{};if(phases)begin=Clock::now();
    k::innerReverse<k::Map128S19>(in,n,masks,[&](std::size_t i,k::block v){scratch[i]=v;});
    auto middle=Clock::time_point{};if(phases)middle=Clock::now();
    for(std::size_t i=0;i<n;i+=1024)k::bchTranspose4(scratch+i,out+i/2);
    if(phases){phases[0]+=std::chrono::duration<double,std::milli>(middle-begin).count();phases[1]+=std::chrono::duration<double,std::milli>(Clock::now()-middle).count();}
}
// Still the invalid-distance sequential route. Keep only one input tile hot,
// then buffer compressed output: in-place writes would destroy unread input.
template<unsigned TileBlocks> static SPIN_NOINLINE void fusedControl(
    const k::block* in,k::block* out,k::block* compressed,k::block* tile,
    std::size_t n,const unsigned* masks) {
    static_assert(std::has_single_bit(TileBlocks) && TileBlocks>=1024);
    k::innerReverse<k::Map128S19>(in,n,masks,[&](std::size_t i,k::block v) {
        tile[i&(TileBlocks-1)]=v;
        if((i&(TileBlocks-1))==0)
            for(unsigned j=0;j<TileBlocks;j+=1024)k::bchTranspose4(tile+j,compressed+(i+j)/2);
    });
    std::memcpy(out,compressed,(n/2)*sizeof(k::block));
}
static void runFused(unsigned tileBlocks,const k::block* in,k::block* out,k::block* compressed,
    k::block* tile,std::size_t n,const unsigned* masks) {
    switch(tileBlocks) {
        case 1024:fusedControl<1024>(in,out,compressed,tile,n,masks);break;
        case 4096:fusedControl<4096>(in,out,compressed,tile,n,masks);break;
        case 16384:fusedControl<16384>(in,out,compressed,tile,n,masks);break;
        default:throw std::invalid_argument("fused diagnostic tile size");
    }
}
template<bool Streaming> static SPIN_NOINLINE void lineDirect(
    const k::block* in,k::block* out,k::block* scratch,std::size_t n,
    const unsigned* offsets,const unsigned* masks) {
    alignas(64) __m128i pending[4];
    k::innerReverse<k::Map128S19>(in,n,masks,[&](std::size_t i,k::block v) {
        const auto at=offsets[i];pending[at&3]=v.mData;
        if((i&3)==0) {
            auto* destination=reinterpret_cast<__m256i*>(scratch+(at&~3U));
            const auto lo=_mm256_load_si256(reinterpret_cast<const __m256i*>(pending));
            const auto hi=_mm256_load_si256(reinterpret_cast<const __m256i*>(pending)+1);
            if constexpr(Streaming) {_mm256_stream_si256(destination,lo);_mm256_stream_si256(destination+1,hi);}
            else {_mm256_store_si256(destination,lo);_mm256_store_si256(destination+1,hi);}
        }
    });
    if constexpr(Streaming)_mm_sfence();
    for(std::size_t i=0;i<n;i+=1024)k::bchTranspose4(scratch+i,out+i/2);
}
static k::block* alignedScratch(std::vector<k::block>& x) {
    return reinterpret_cast<k::block*>((reinterpret_cast<std::uintptr_t>(x.data())+63)&~std::uintptr_t{63});
}
struct MappedLayout {
    // One cache-line pad rotates cache-set indices between four-row blocks.
    // Without it, a fixed shuffled column lands at one offset in every 16 KiB block.
    static constexpr unsigned stride=1028;
    std::vector<unsigned> offsets,permutation;
    MappedLayout(const std::vector<unsigned>& r,unsigned rows,unsigned columns):offsets(r.size()),permutation(r.size()/4,~0U) {
        for(unsigned i=0;i<r.size();++i) {
            const unsigned row=r[i]/256,c=r[i]%256;
            const unsigned shuffled=(i/(rows*columns))*columns+(i%(4*columns))/4;
            offsets[i]=(row/4)*stride+4*shuffled+(row&3);
            auto& p=permutation[(row/4)*256+c];
            if(p!=~0U && p!=4*shuffled)throw std::runtime_error("mapped BCH needs shared four-row permutation");
            p=4*shuffled;
        }
        for(auto p:permutation)if(p==~0U)throw std::runtime_error("incomplete mapped permutation");
    }
};
struct RouteBundle {unsigned base;std::uint64_t lanes;};
static SPIN_NOINLINE void repackGfni(const k::block* in,k::block* out,std::size_t n,const unsigned* permutation) {
    alignas(64) k::block tile[1024];
    for(std::size_t group=0;group<n/1024;++group) {
        const auto* source=in+group*MappedLayout::stride;
        const auto* p=permutation+group*256;
        for(unsigned i=0;i<256;i+=4) {
            const auto x0=_mm512_loadu_si512(source+p[i]);
            const auto x1=_mm512_loadu_si512(source+p[i+1]);
            const auto x2=_mm512_loadu_si512(source+p[i+2]);
            const auto x3=_mm512_loadu_si512(source+p[i+3]);
            _mm512_store_si512(tile+4*i,x0);_mm512_store_si512(tile+4*(i+1),x1);
            _mm512_store_si512(tile+4*(i+2),x2);_mm512_store_si512(tile+4*(i+3),x3);
        }
        k::bchTranspose4GfniBlend(tile,out+group*512);
    }
}
template<bool Copy,bool Prefetch,bool Prepared,bool Manual=false,bool Packed=false,unsigned Tile=0,unsigned Ahead=1,bool L1=false> static SPIN_NOINLINE void warmGfni(
    const k::block* in,k::block* out,std::size_t n,const unsigned* permutation) {
    alignas(64) k::block tile[Copy?1024:1];
    for(std::size_t group=0;group<n/1024;++group) {
        if constexpr(Prefetch)if(group+Ahead<n/1024) {
            const auto* next=in+(group+Ahead)*MappedLayout::stride;
            for(unsigned i=0;i<1024;i+=16) {
                _mm_prefetch(reinterpret_cast<const char*>(next+i),L1?_MM_HINT_T0:_MM_HINT_T1);
                _mm_prefetch(reinterpret_cast<const char*>(next+i+4),L1?_MM_HINT_T0:_MM_HINT_T1);
                _mm_prefetch(reinterpret_cast<const char*>(next+i+8),L1?_MM_HINT_T0:_MM_HINT_T1);
                _mm_prefetch(reinterpret_cast<const char*>(next+i+12),L1?_MM_HINT_T0:_MM_HINT_T1);
            }
        }
        const auto* source=in+group*MappedLayout::stride;
        if constexpr(Copy) {
            if constexpr(Manual) {
                for(unsigned i=0;i<1024;i+=16) {
                    const auto a=_mm512_loadu_si512(source+i);
                    const auto b=_mm512_loadu_si512(source+i+4);
                    const auto c=_mm512_loadu_si512(source+i+8);
                    const auto d=_mm512_loadu_si512(source+i+12);
                    _mm512_store_si512(tile+i,a);_mm512_store_si512(tile+i+4,b);
                    _mm512_store_si512(tile+i+8,c);_mm512_store_si512(tile+i+12,d);
                }
            } else std::memcpy(tile,source,sizeof(tile));
            source=tile;
        }
        if constexpr(Tile==2)k::bchTranspose4GfniBlend(source,out+group*512);
        else if constexpr(Tile==1)k::bchTranspose4GfniMapped1(source,out+group*512,permutation+group*256);
        else if constexpr(Packed)k::bchTranspose4GfniPacked(source,out+group*512,permutation+group*256);
        else if constexpr(Prepared)k::bchTranspose4GfniPrepared(source,out+group*512,permutation+group*256);
        else k::bchTranspose4GfniMapped(source,out+group*512,permutation+group*256);
    }
}
struct CanonicalBundles {
    std::vector<unsigned> bases;
    std::vector<unsigned char> lanes;
    std::vector<std::uint64_t> shuffle;
    CanonicalBundles()=default;
    // A single complete cache line per destination, already in BCH order.
    // This changes storage only, not the sampled route or the lane permutation.
    CanonicalBundles(const std::vector<unsigned>& route):bases(route.size()/4),shuffle(route.size()/4) {
        for(unsigned i=0;i<route.size();i+=4) {
            const unsigned row=route[i]/256,column=route[i]%256;
            bases[i/4]=(row/4)*MappedLayout::stride+4*column;
            unsigned seen=0;
            for(unsigned source=0;source<4;++source) {
                const unsigned r=route[i+source]/256,c=route[i+source]%256,lane=r%4;
                if(r/4!=row/4 || c!=column || (seen&(1U<<lane)))throw std::runtime_error("physical column is not a row permutation");
                seen|=1U<<lane;
                for(unsigned half=0;half<2;++half)
                    shuffle[i/4]|=std::uint64_t(2*source+half)<<(8*(2*lane+half));
            }
            if(seen!=15)throw std::runtime_error("incomplete physical column");
        }
    }
    CanonicalBundles(const MappedLayout& mapped,unsigned width):bases(mapped.offsets.size()/width),lanes(mapped.offsets.size()),shuffle(mapped.offsets.size()/4) {
        for(unsigned i=0;i<mapped.offsets.size();i+=width) {
            const unsigned base=mapped.offsets[i]&~3U;
            if(base%4)throw std::runtime_error("unaligned mapped bundle");
            bases[i/width]=base;unsigned seen=0;
            for(unsigned j=0;j<width;++j) {
                const auto at=mapped.offsets[i+j];
                if(at<base || at>=base+width)throw std::runtime_error("mapped bundle not contiguous");
                lanes[i+j]=static_cast<unsigned char>(at-base);seen|=1U<<(at-base);
            }
            if(seen!=((1U<<width)-1))throw std::runtime_error("mapped bundle lane bijection");
            for(unsigned column=0;column<width;column+=4)for(unsigned input=0;input<4;++input) {
                const unsigned destination=lanes[i+column+input];
                if(destination/4!=column/4)throw std::runtime_error("mapped column changed during lane shuffle");
                for(unsigned half=0;half<2;++half)
                    shuffle[(i+column)/4]|=std::uint64_t(2*input+half)<<(8*(2*(destination%4)+half));
            }
        }
    }
    CanonicalBundles(const std::vector<unsigned>& route,unsigned width):bases(route.size()/width),lanes(route.size()) {
        for(unsigned i=0;i<route.size();i+=width) {
            const unsigned base=packed(route[i])&~(width-1);
            bases[i/width]=base;
            unsigned seen=0;
            for(unsigned j=0;j<width;++j) {
                const unsigned destination=packed(route[i+j]);
                if((destination&~(width-1))!=base)throw std::runtime_error("canonical bundle not contiguous");
                lanes[i+j]=static_cast<unsigned char>(destination&(width-1));
                seen|=1U<<lanes[i+j];
            }
            if(seen!=((1U<<width)-1))throw std::runtime_error("canonical bundle lane bijection");
        }
    }
};
// Keep the long-range route cache-line sized. Row zero chooses the physical
// column; the other rows are restored locally immediately before BCH.
struct IndependentPackets {
    CanonicalBundles bundles;
    std::vector<std::uint16_t> offsets;
    explicit IndependentPackets(const std::vector<unsigned>& route):offsets(route.size(),65535) {
        bundles.bases.resize(route.size()/4);bundles.shuffle.resize(route.size()/4);
        for(std::size_t i=0;i<route.size();i+=4) {
            const unsigned group=route[i]/1024;
            unsigned anchor=~0U,seen=0;
            for(unsigned l=0;l<4;++l)if((route[i+l]/256)%4==0)anchor=route[i+l]%256;
            if(anchor==~0U)throw std::runtime_error("packet missing row zero");
            bundles.bases[i/4]=group*MappedLayout::stride+4*anchor;
            for(unsigned l=0;l<4;++l) {
                const unsigned row=route[i+l]/256,lane=row%4,column=route[i+l]%256;
                if(row/4!=group || (seen&(1U<<lane)))throw std::runtime_error("packet row membership");
                seen|=1U<<lane;
                auto& offset=offsets[group*1024+4*column+lane];
                if(offset!=65535)throw std::runtime_error("duplicate packet coordinate");
                offset=4*anchor+lane;
                for(unsigned half=0;half<2;++half)
                    bundles.shuffle[i/4]|=std::uint64_t(2*l+half)<<(8*(2*lane+half));
            }
            if(seen!=15)throw std::runtime_error("incomplete packet");
        }
        for(auto offset:offsets)if(offset>=1024)throw std::runtime_error("missing packet coordinate");
    }
};
template<bool Copy=false,bool Prefetch=false> static SPIN_NOINLINE void independentRepack(const k::block* in,k::block* out,std::size_t n,const std::uint16_t* offsets) {
    alignas(64) k::block tile[1024],local[Copy?1024:1];
    for(std::size_t group=0;group<n/1024;++group) {
        const auto* source=in+group*MappedLayout::stride;
        if constexpr(Copy) {
            for(unsigned i=0;i<1024;i+=16) {
                const auto a=_mm512_load_si512(source+i),b=_mm512_load_si512(source+i+4);
                const auto c=_mm512_load_si512(source+i+8),d=_mm512_load_si512(source+i+12);
                _mm512_store_si512(local+i,a);_mm512_store_si512(local+i+4,b);
                _mm512_store_si512(local+i+8,c);_mm512_store_si512(local+i+12,d);
            }
            source=local;
        }
        if constexpr(Prefetch)for(unsigned i=0;i<1024;i+=4)_mm_prefetch(reinterpret_cast<const char*>(source+i),_MM_HINT_T0);
        const auto* p=offsets+group*1024;
#if defined(__GNUC__)
#pragma GCC unroll 4
#endif
        for(unsigned c=0;c<256;++c) {
            auto x=_mm512_castsi128_si512(_mm_load_si128(reinterpret_cast<const __m128i*>(source+4*c)));
            x=_mm512_inserti32x4(x,_mm_load_si128(reinterpret_cast<const __m128i*>(source+p[4*c+1])),1);
            x=_mm512_inserti32x4(x,_mm_load_si128(reinterpret_cast<const __m128i*>(source+p[4*c+2])),2);
            x=_mm512_inserti32x4(x,_mm_load_si128(reinterpret_cast<const __m128i*>(source+p[4*c+3])),3);
            _mm512_store_si512(tile+4*c,x);
        }
        k::bchTranspose4GfniBlend(tile,out+group*512);
    }
}
static SPIN_NOINLINE void independentFusedCopy(const k::block* in,k::block* out,std::size_t n,const std::uint16_t* offsets) {
    alignas(64) k::block local[1024];
    for(std::size_t group=0;group<n/1024;++group) {
        const auto* source=in+group*MappedLayout::stride;
        for(unsigned i=0;i<1024;i+=16) {
            const auto a=_mm512_load_si512(source+i),b=_mm512_load_si512(source+i+4);
            const auto c=_mm512_load_si512(source+i+8),d=_mm512_load_si512(source+i+12);
            _mm512_store_si512(local+i,a);_mm512_store_si512(local+i+4,b);
            _mm512_store_si512(local+i+8,c);_mm512_store_si512(local+i+12,d);
        }
        k::bchTranspose4GfniIndependent(local,out+group*512,offsets+group*1024);
    }
}
template<unsigned Width,bool Streaming,bool VectorRoute=false,class Map=k::Map128S19> static SPIN_NOINLINE void canonicalDirect(
    const k::block* in,k::block* out,k::block* scratch,std::size_t n,
    const unsigned* masks,const CanonicalBundles& route,double* phases=nullptr,unsigned schedule=0,const unsigned* permutation=nullptr,const std::uint16_t* rowOffsets=nullptr) {
    static_assert(Width==4 || Width==8 || Width==16);
    alignas(64) __m128i pending[Width];
    __m128i v1{},v2{},v3{};
    auto start=Clock::time_point{};if(phases)start=Clock::now();
    k::innerReverse<Map>(in,n,masks,[&](std::size_t i,k::block value) {
        if constexpr(VectorRoute) {
            // Reverse traversal produces lanes 3,2,1,0. Keep one whole column
            // in registers, then permute its four 128-bit elements together.
            switch(i&3) {
                case 3:v3=value.mData;break;
                case 2:v2=value.mData;break;
                case 1:v1=value.mData;break;
                case 0: {
                    auto z=_mm512_castsi128_si512(value.mData);
                    z=_mm512_inserti32x4(z,v1,1);z=_mm512_inserti32x4(z,v2,2);z=_mm512_inserti32x4(z,v3,3);
                    const auto control=_mm512_cvtepu8_epi64(_mm_loadl_epi64(reinterpret_cast<const __m128i*>(route.shuffle.data()+i/4)));
                    z=_mm512_permutexvar_epi64(control,z);
                    auto* destination=reinterpret_cast<__m512i*>(scratch+route.bases[i/Width]+(i&(Width-1)));
                    if constexpr(Streaming)_mm512_stream_si512(destination,z);
                    else _mm512_store_si512(destination,z);
                    break;
                }
            }
        } else {
        pending[route.lanes[i]]=value.mData;
        if((i&(Width-1))==0) {
            auto* destination=reinterpret_cast<__m512i*>(scratch+route.bases[i/Width]);
            const auto* source=reinterpret_cast<const __m512i*>(pending);
            if constexpr(Streaming) {
                _mm512_stream_si512(destination,_mm512_load_si512(source));
                _mm512_stream_si512(destination+1,_mm512_load_si512(source+1));
                if constexpr(Width==16) {
                    _mm512_stream_si512(destination+2,_mm512_load_si512(source+2));
                    _mm512_stream_si512(destination+3,_mm512_load_si512(source+3));
                }
            } else {
                _mm512_store_si512(destination,_mm512_load_si512(source));
                _mm512_store_si512(destination+1,_mm512_load_si512(source+1));
                if constexpr(Width==16) {
                    _mm512_store_si512(destination+2,_mm512_load_si512(source+2));
                    _mm512_store_si512(destination+3,_mm512_load_si512(source+3));
                }
            }
        }
        }
    });
    if constexpr(Streaming)_mm_sfence();
    auto middle=Clock::time_point{};if(phases)middle=Clock::now();
    switch(schedule) {
        case 0:for(std::size_t i=0;i<n;i+=1024)k::bchTranspose4(scratch+i,out+i/2);break;
        case 1:for(std::size_t i=0;i<n;i+=1024)k::bchTranspose4Greedy(scratch+i,out+i/2);break;
        case 2:for(std::size_t i=0;i<n;i+=1024)k::bchTranspose4Release(scratch+i,out+i/2);break;
        case 3:for(std::size_t i=0;i<n;i+=1024)k::bchTranspose4Reverse(scratch+i,out+i/2);break;
        case 4:for(std::size_t i=0;i<n;i+=1024)k::bchTranspose4Restrict(scratch+i,out+i/2);break;
        case 5:for(std::size_t i=0;i<n;i+=1024)k::bchTranspose4Share6(scratch+i,out+i/2);break;
        case 6:for(std::size_t i=0;i<n;i+=1024)k::bchTranspose4Share8(scratch+i,out+i/2);break;
        case 7:for(std::size_t i=0;i<n;i+=1024)k::bchTranspose4Gfni1(scratch+i,out+i/2);break;
        case 8:for(std::size_t i=0;i<n;i+=1024)k::bchTranspose4Gfni2(scratch+i,out+i/2);break;
        case 9:for(std::size_t i=0;i<n;i+=1024)k::bchTranspose4GfniStatic(scratch+i,out+i/2);break;
        case 10:for(std::size_t i=0;i<n;i+=1024)k::bchTranspose4GfniTernary(scratch+i,out+i/2);break;
        case 11:for(std::size_t i=0;i<n;i+=1024)k::bchTranspose4GfniUnroll(scratch+i,out+i/2);break;
        case 12:for(std::size_t i=0;i<n;i+=1024)k::bchTranspose4GfniBlend(scratch+i,out+i/2);break;
        case 13:for(std::size_t g=0;g<n/1024;++g)k::bchTranspose4GfniMapped(scratch+g*MappedLayout::stride,out+g*512,permutation+g*256);break;
        case 14:repackGfni(scratch,out,n,permutation);break;
        case 15:for(std::size_t g=0;g<n/1024;++g)k::bchTranspose4GfniPrepared(scratch+g*MappedLayout::stride,out+g*512,permutation+g*256);break;
        case 16:warmGfni<false,true,false>(scratch,out,n,permutation);break;
        case 17:warmGfni<true,false,true>(scratch,out,n,permutation);break;
        case 18:warmGfni<false,true,true>(scratch,out,n,permutation);break;
        case 19:warmGfni<true,true,true>(scratch,out,n,permutation);break;
        case 20:warmGfni<true,false,false>(scratch,out,n,permutation);break;
        case 21:warmGfni<true,false,false,true>(scratch,out,n,permutation);break;
        case 22:warmGfni<true,true,false>(scratch,out,n,permutation);break;
        case 23:warmGfni<true,false,false,false,true>(scratch,out,n,permutation);break;
        case 24:warmGfni<true,false,false,true,true>(scratch,out,n,permutation);break;
        case 25:warmGfni<true,false,false,true,true>(scratch,out,n,permutation);break;
        case 26:warmGfni<true,false,false,true>(scratch,out,n,permutation);break;
        case 27:warmGfni<true,true,false,true,false,0,1>(scratch,out,n,permutation);break;
        case 28:warmGfni<true,true,false,true,false,0,2>(scratch,out,n,permutation);break;
        case 29:warmGfni<true,true,false,true,false,0,4>(scratch,out,n,permutation);break;
        case 30:warmGfni<true,false,false,true,false,1>(scratch,out,n,permutation);break;
        case 31:warmGfni<true,true,false,true,false,1,2>(scratch,out,n,permutation);break;
        case 32:warmGfni<false,false,false,false,false,2>(scratch,out,n,permutation);break;
        case 33:warmGfni<true,false,false,true,false,2>(scratch,out,n,permutation);break;
        case 34:warmGfni<false,true,false,false,false,2,0,true>(scratch,out,n,permutation);break;
        case 35:warmGfni<false,true,false,false,false,0,0,true>(scratch,out,n,permutation);break;
        case 36:independentRepack<>(scratch,out,n,rowOffsets);break;
        case 37:for(std::size_t g=0;g<n/1024;++g)k::bchTranspose4GfniIndependent(scratch+g*MappedLayout::stride,out+g*512,rowOffsets+g*1024);break;
        case 38:independentRepack<true>(scratch,out,n,rowOffsets);break;
        case 39:independentRepack<false,true>(scratch,out,n,rowOffsets);break;
        case 40:independentFusedCopy(scratch,out,n,rowOffsets);break;
        default:throw std::invalid_argument("BCH schedule");
    }
    if(phases){phases[0]+=std::chrono::duration<double,std::milli>(middle-start).count();phases[1]+=std::chrono::duration<double,std::milli>(Clock::now()-middle).count();}
}
static void runCanonical(unsigned columns,bool streaming,const k::block* in,k::block* out,k::block* scratch,
    std::size_t n,const unsigned* masks,const CanonicalBundles& route,double* phases=nullptr,unsigned schedule=0,const unsigned* permutation=nullptr) {
    if(schedule>=25) {
        if(!streaming || route.shuffle.size()!=n/4)throw std::invalid_argument("vector route requires mapped shuffle controls");
        if(schedule>=32 && schedule<=34)canonicalDirect<4,true,true>(in,out,scratch,n,masks,route,phases,schedule,permutation);
        else if(columns==2)canonicalDirect<8,true,true>(in,out,scratch,n,masks,route,phases,schedule,permutation);
        else if(columns==4)canonicalDirect<16,true,true>(in,out,scratch,n,masks,route,phases,schedule,permutation);
        else throw std::invalid_argument("vector route column count");
        return;
    }
    if(columns==2) {
        if(streaming)canonicalDirect<8,true>(in,out,scratch,n,masks,route,phases,schedule,permutation);
        else canonicalDirect<8,false>(in,out,scratch,n,masks,route,phases,schedule,permutation);
    } else if(columns==4) {
        if(streaming)canonicalDirect<16,true>(in,out,scratch,n,masks,route,phases,schedule,permutation);
        else canonicalDirect<16,false>(in,out,scratch,n,masks,route,phases,schedule,permutation);
    } else throw std::invalid_argument("canonical bundle columns");
}
struct WideLayout {
    unsigned group,stride;
    std::vector<unsigned> bases;
    std::vector<unsigned char> lanes;
    WideLayout(const std::vector<unsigned>& route,unsigned g,unsigned columnPad=0):group(g),stride(256*(g+columnPad)+4),bases(route.size()/g),lanes(route.size()) {
        for(unsigned i=0;i<route.size();i+=g) {
            const unsigned row=route[i]/256,column=route[i]&255U;
            bases[i/g]=(row/g)*stride+column*(g+columnPad);
            std::uint64_t seen=0;
            for(unsigned j=0;j<g;++j) {
                const unsigned r=route[i+j]/256;
                if(r/g!=row/g || (route[i+j]&255U)!=column)throw std::runtime_error("wide group not contiguous");
                lanes[i+j]=static_cast<unsigned char>(r&(g-1));seen|=std::uint64_t{1}<<(r&(g-1));
            }
            if(seen!=((std::uint64_t{1}<<g)-1))throw std::runtime_error("wide group lane bijection");
        }
    }
};
struct SplitLayout {
    std::vector<unsigned> pendingBase,destination;
    std::vector<unsigned char> lanes;
    SplitLayout()=default;
    SplitLayout(const std::vector<unsigned>& route,unsigned g):pendingBase(route.size()/16),destination(route.size()/16),lanes(route.size()) {
        const unsigned rows=route.size()/256,bandSize=rows/(g/16),stride=256*g+4;
        std::vector<std::uint64_t> seenLo(rows/g),seenHi(rows/g);
        std::vector<unsigned> groupColumn(rows/g,~0U);
        for(unsigned i=0;i<route.size();i+=16) {
            if(i%rows==0){std::fill(seenLo.begin(),seenLo.end(),0);std::fill(seenHi.begin(),seenHi.end(),0);std::fill(groupColumn.begin(),groupColumn.end(),~0U);}
            const unsigned row=route[i]/256,b=row/g,column=route[i]&255;
            if(groupColumn[b]!=~0U && groupColumn[b]!=column)throw std::runtime_error("split column disagreement");groupColumn[b]=column;
            pendingBase[i/16]=b*g;destination[i/16]=b*stride+column*g;
            for(unsigned j=0;j<16;++j) {
                const auto r=route[i+j]/256,lane=r%g;
                if(r/g!=b || (route[i+j]&255)!=column)throw std::runtime_error("split chunk geometry");
                lanes[i+j]=static_cast<unsigned char>(lane);
                auto& seen=lane<64?seenLo[b]:seenHi[b];const auto bit=std::uint64_t{1}<<(lane%64);
                if(seen&bit)throw std::runtime_error("split duplicate lane");seen|=bit;
            }
        }
        // checkRoute already checks every row exactly once per region.
        if(bandSize%128)throw std::invalid_argument("split proof tests require epoch-aligned bands");
    }
};
template<unsigned G> static SPIN_NOINLINE void splitDirect(const k::block* in,k::block* out,
    k::block* scratch,k::block* tile,k::block* pending,std::size_t n,const unsigned* masks,const SplitLayout& layout) {
    constexpr unsigned stride=256*G+4;
    const auto rows=n/256,regionMask=rows-1,bandSize=rows/(G/16);
    unsigned pendingBase=0;
    k::innerReverse<k::Map128S19>(in,n,masks,[&](std::size_t i,k::block v) {
        if((i&15)==15)pendingBase=layout.pendingBase[i/16];
        pending[pendingBase+layout.lanes[i]]=v;
        if((i&15)==0 && (i&regionMask)<bandSize) {
            auto* dest=scratch+layout.destination[i/16];
            for(unsigned j=0;j<G;j+=4)_mm512_stream_si512(reinterpret_cast<__m512i*>(dest+j),_mm512_load_si512(pending+pendingBase+j));
        }
    });
    _mm_sfence();
    for(std::size_t b=0;b<n/(256*G);++b)for(unsigned row=0;row<G;row+=4) {
        for(unsigned c=0;c<256;++c)_mm512_storeu_si512(tile+4*c,_mm512_load_si512(scratch+b*stride+G*c+row));
        k::bchTranspose4(tile,out+b*G*128+row*128);
    }
}
static void runSplit(unsigned group,const k::block* in,k::block* out,k::block* scratch,k::block* tile,
    k::block* pending,std::size_t n,const unsigned* masks,const SplitLayout& layout) {
    switch(group) {
        case 32:splitDirect<32>(in,out,scratch,tile,pending,n,masks,layout);break;
        case 64:splitDirect<64>(in,out,scratch,tile,pending,n,masks,layout);break;
        case 128:splitDirect<128>(in,out,scratch,tile,pending,n,masks,layout);break;
        default:throw std::invalid_argument("split group");
    }
}
template<unsigned G,bool Streaming,bool XorLanes=false,unsigned ColumnPad=0,bool StridedBch=false> static SPIN_NOINLINE void wideDirect(
    const k::block* in,k::block* out,k::block* scratch,k::block* tile,std::size_t n,
    const unsigned* masks,const unsigned* bases,const unsigned char* lanes,double* phases=nullptr) {
    constexpr unsigned stride=256*(G+ColumnPad)+4;
    alignas(64) k::block pending[G];
    auto begin=Clock::time_point{};if(phases)begin=Clock::now();
    k::innerReverse<k::Map128S19>(in,n,masks,[&](std::size_t i,k::block v) {
        if constexpr(XorLanes)pending[i&(G-1)]=v;
        else pending[lanes[i]]=v;
        if((i&(G-1))==0) {
            auto* dest=scratch+bases[i/G];
            const unsigned shift=XorLanes?lanes[i]:0;
            const auto indices=_mm512_xor_si512(_mm512_setr_epi64(0,1,2,3,4,5,6,7),_mm512_set1_epi64(2*(shift&3)));
            for(unsigned j=0;j<G;j+=4) {
                auto value=_mm512_load_si512(pending+j);
                if constexpr(XorLanes)value=_mm512_permutexvar_epi64(indices,value);
                auto* at=dest+(j^(shift&~3U));
                if constexpr(Streaming)_mm512_stream_si512(reinterpret_cast<__m512i*>(at),value);
                else _mm512_store_si512(at,value);
            }
        }
    });
    if constexpr(Streaming)_mm_sfence();
    auto middle=Clock::time_point{};if(phases)middle=Clock::now();
    for(std::size_t group=0;group<n/(256*G);++group) {
        const auto* src=scratch+group*stride;
        for(unsigned row=0;row<G;row+=4) {
            auto* result=out+group*G*128+row*128;
            if constexpr(StridedBch) {
                static_assert(G==16);
                if constexpr(ColumnPad==0)k::bchTranspose4Stride16(src+row,result);
                else {static_assert(ColumnPad==4);k::bchTranspose4Stride20(src+row,result);}
            } else {
                for(unsigned c=0;c<256;++c)_mm512_storeu_si512(tile+4*c,_mm512_load_si512(src+(G+ColumnPad)*c+row));
                k::bchTranspose4(tile,result);
            }
        }
    }
    if(phases){phases[0]+=std::chrono::duration<double,std::milli>(middle-begin).count();phases[1]+=std::chrono::duration<double,std::milli>(Clock::now()-middle).count();}
}
template<bool Streaming> static void runWide(unsigned group,const k::block* in,k::block* out,
    k::block* scratch,k::block* tile,std::size_t n,const unsigned* masks,const WideLayout& layout,double* phases=nullptr,bool xorLanes=false,unsigned stridedBch=0) {
    if(stridedBch) {
        if(group!=16 || !Streaming || xorLanes)throw std::invalid_argument("strided BCH geometry");
        if(stridedBch==16)wideDirect<16,true,false,0,true>(in,out,scratch,tile,n,masks,layout.bases.data(),layout.lanes.data(),phases);
        else if(stridedBch==20)wideDirect<16,true,false,4,true>(in,out,scratch,tile,n,masks,layout.bases.data(),layout.lanes.data(),phases);
        else throw std::invalid_argument("strided BCH stride");
        return;
    }
    if(xorLanes) {
        if(group!=16 || !Streaming)throw std::invalid_argument("xor wide geometry");
        wideDirect<16,true,true>(in,out,scratch,tile,n,masks,layout.bases.data(),layout.lanes.data(),phases);return;
    }
    switch(group) {
        case 8:wideDirect<8,Streaming>(in,out,scratch,tile,n,masks,layout.bases.data(),layout.lanes.data(),phases);break;
        case 16:wideDirect<16,Streaming>(in,out,scratch,tile,n,masks,layout.bases.data(),layout.lanes.data(),phases);break;
        case 32:wideDirect<32,Streaming>(in,out,scratch,tile,n,masks,layout.bases.data(),layout.lanes.data(),phases);break;
        default:throw std::invalid_argument("wide group must be 8,16,32");
    }
}
// Save the exact reverse state before each epoch's emission. Padding to twenty
// blocks aligns every state to a cache line without a 512-byte power-of-two stride.
static SPIN_NOINLINE void reverseStates(const k::block* in,k::block* saved,std::size_t n,const unsigned* masks) {
    using Map=k::Map128S19;
    alignas(32) __m128i state[Map::S]{},syndrome[Map::S],raw[Map::T];
    const auto epochs=n/Map::T;
    for(std::size_t epoch=epochs;epoch-->0;) {
        std::memcpy(saved+20*epoch,state,sizeof(state));
        if(epoch==0)break;
        std::memcpy(raw,in+epoch*Map::T,sizeof(raw));
        k::zeta<Map::T>(raw);Map::finish(raw,syndrome);
        if(epoch+1==epochs)std::memcpy(state,syndrome,sizeof(state));
        else {
            auto u=masks[2*epoch];auto dot=_mm_setzero_si128();
            while(u){const auto j=std::countr_zero(u);u&=u-1;dot=_mm_xor_si128(dot,state[j]);}
            auto v=masks[2*epoch+1];
            while(v){const auto j=std::countr_zero(v);v&=v-1;state[j]=_mm_xor_si128(state[j],dot);}
            for(unsigned j=0;j<Map::S;++j)state[j]=_mm_xor_si128(state[j],syndrome[j]);
        }
    }
}
template<unsigned Start,std::size_t... J> static SPIN_FORCEINLINE void gatherPoints(
    const k::block* in,k::block* tile,const __m128i* state,const unsigned* offsets,std::index_sequence<J...>) {
    ((tile[offsets[J]&1023U]=k::block(_mm_xor_si128(in[J].mData,
        k::imtSparseSum<k::Map128S19::feedbackColumns[Start+J]>(state)))),...);
}
template<unsigned Width,unsigned Start=0> static SPIN_FORCEINLINE void gatherDispatch(
    unsigned start,const k::block* in,k::block* tile,const __m128i* state,const unsigned* offsets) {
    if(start==Start)gatherPoints<Start>(in,tile,state,offsets,std::make_index_sequence<Width>{});
    else if constexpr(Start+Width<128)gatherDispatch<Width,Start+Width>(start,in,tile,state,offsets);
    else throw std::logic_error("gather epoch alignment");
}
static std::vector<unsigned> gatherSources(const std::vector<unsigned>& offsets,unsigned columns) {
    const unsigned width=4*columns;
    std::vector<unsigned> sources(offsets.size()/width,~0U),counts(offsets.size()/1024);
    for(unsigned i=0;i<offsets.size();i+=width) {
        const unsigned group=offsets[i]/1024;
        for(unsigned j=1;j<width;++j)if(offsets[i+j]/1024!=group)throw std::runtime_error("gather group split");
        auto& count=counts[group];if(count>=1024/width)throw std::runtime_error("gather group overflow");
        sources[group*(1024/width)+count++]=i;
    }
    for(auto c:counts)if(c!=1024/width)throw std::runtime_error("gather group incomplete");
    return sources;
}
template<unsigned Columns> static SPIN_NOINLINE void stateGather(
    const k::block* in,k::block* out,k::block* compressed,k::block* tile,k::block* states,
    std::size_t n,const unsigned* masks,const unsigned* offsets,const unsigned* sources) {
    constexpr unsigned width=4*Columns;
    reverseStates(in,states,n,masks);
    for(unsigned group=0;group<n/1024;++group) {
        for(unsigned b=0;b<1024/width;++b) {
            const unsigned source=sources[group*(1024/width)+b];
            const auto* state=reinterpret_cast<const __m128i*>(states+20*(source/128));
            gatherDispatch<width>(source&127,in+source,tile,state,offsets+source);
        }
        k::bchTranspose4(tile,compressed+group*512);
    }
    std::memcpy(out,compressed,(n/2)*sizeof(k::block));
}
static void runStateGather(unsigned columns,const k::block* in,k::block* out,k::block* compressed,
    k::block* tile,k::block* states,std::size_t n,const unsigned* masks,const unsigned* offsets,const unsigned* sources) {
    switch(columns) {
        case 1:stateGather<1>(in,out,compressed,tile,states,n,masks,offsets,sources);break;
        case 4:stateGather<4>(in,out,compressed,tile,states,n,masks,offsets,sources);break;
        case 8:stateGather<8>(in,out,compressed,tile,states,n,masks,offsets,sources);break;
        default:throw std::invalid_argument("state-gather columns must be 1,4,8");
    }
}
static std::vector<RouteBundle> bundles(const MappedLayout& mapped,unsigned columns) {
    const unsigned width=4*columns;
    if(width>32 || !std::has_single_bit(width))throw std::invalid_argument("bundle width");
    std::vector<RouteBundle> result(mapped.offsets.size()/width);
    for(unsigned i=0;i<mapped.offsets.size();i+=width) {
        auto& b=result[i/width];b.base=mapped.offsets[i]&~3U;b.lanes=0;
        for(unsigned j=0;j<width;++j) {
            if((mapped.offsets[i+j]&~3U)!=b.base+(j&~3U))throw std::runtime_error("noncontiguous route bundle");
            b.lanes|=std::uint64_t(mapped.offsets[i+j]&3)<<(2*j);
        }
    }
    return result;
}
static SPIN_NOINLINE void repackBch(const k::block* scratch,k::block* out,std::size_t n,const unsigned* permutation) {
    alignas(64) k::block tile[1024];
    for(std::size_t group=0;group<n/1024;++group) {
        const auto* src=scratch+group*MappedLayout::stride;
        const auto* perm=permutation+group*256;
        for(unsigned c=0;c<256;++c) {
            const auto* at=reinterpret_cast<const __m256i*>(src+perm[c]);
            auto* to=reinterpret_cast<__m256i*>(tile+4*c);
            _mm256_store_si256(to,_mm256_load_si256(at));
            _mm256_store_si256(to+1,_mm256_load_si256(at+1));
        }
        k::bchTranspose4(tile,out+group*512);
    }
}
template<unsigned Columns> static SPIN_NOINLINE void bundledDirect(
    const k::block* in,k::block* out,k::block* scratch,std::size_t n,
    const RouteBundle* route,const unsigned* masks,const unsigned* permutation) {
    constexpr unsigned width=4*Columns;
    unsigned destination=0;std::uint64_t lanes=0;
    k::innerReverse<k::Map128S19>(in,n,masks,[&](std::size_t i,k::block v) {
        const unsigned local=unsigned(i)&(width-1);
        if(local==width-1){destination=route[i/width].base;lanes=route[i/width].lanes;}
        scratch[destination+(local&~3U)+unsigned((lanes>>(2*local))&3)]=v;
    });
    repackBch(scratch,out,n,permutation);
}
static void runBundled(unsigned columns,const k::block* in,k::block* out,k::block* scratch,
                       std::size_t n,const RouteBundle* route,const unsigned* masks,const unsigned* permutation) {
    switch(columns) {
        case 1:bundledDirect<1>(in,out,scratch,n,route,masks,permutation);break;
        case 2:bundledDirect<2>(in,out,scratch,n,route,masks,permutation);break;
        case 4:bundledDirect<4>(in,out,scratch,n,route,masks,permutation);break;
        case 8:bundledDirect<8>(in,out,scratch,n,route,masks,permutation);break;
        default:throw std::invalid_argument("bundled columns");
    }
}
// A shared global coordinate order permits region-local writes into a
// column-major matrix. Repack a fixed 32-row tile before the original BCH circuit.
// Scratch is a separate allocation; in and out may alias for in-place encoding.
static SPIN_NOINLINE void columnDirect(const k::block* in,k::block* out,k::block* __restrict scratch,
    k::block* tile,std::size_t n,const unsigned* route,const unsigned* masks,const unsigned* order,double* phases=nullptr) {
    auto start=Clock::time_point{};if(phases)start=Clock::now();
    const auto rows=n/256;
    const auto regionMask=~(rows-1); // Experimental lengths are powers of two.
    k::innerReverse<k::Map128S19>(in,n,masks,[&](std::size_t i,k::block v) {
        scratch[(i&regionMask)+(route[i]>>8)]=v;
    });
    auto middle=Clock::time_point{};if(phases)middle=Clock::now();
    constexpr unsigned tileRows=32;
    for(std::size_t first=0;first<rows;first+=tileRows) {
        for(unsigned c=0;c<256;++c) {
            const auto* from=scratch+c*rows+first;
            const unsigned column=order[c]*4;
            for(unsigned group=0;group<tileRows/4;++group) {
                const auto* src=reinterpret_cast<const __m256i*>(from+4*group);
                auto* dst=reinterpret_cast<__m256i*>(tile+group*1024+column);
                _mm256_storeu_si256(dst,_mm256_load_si256(src));
                _mm256_storeu_si256(dst+1,_mm256_load_si256(src+1));
            }
        }
        for(unsigned group=0;group<tileRows/4;++group)
            k::bchTranspose4(tile+group*1024,out+first*128+group*512);
    }
    if(phases){phases[0]+=std::chrono::duration<double,std::milli>(middle-start).count();phases[1]+=std::chrono::duration<double,std::milli>(Clock::now()-middle).count();}
}
template<bool Streaming> static SPIN_NOINLINE void mappedDirect(
    const k::block* in,k::block* out,k::block* scratch,std::size_t n,
    const unsigned* offsets,const unsigned* masks,const unsigned* permutation) {
    alignas(64) __m128i pending[4];
    k::innerReverse<k::Map128S19>(in,n,masks,[&](std::size_t i,k::block v) {
        const auto at=offsets[i];
        if constexpr(!Streaming)scratch[at]=v;
        else {
            pending[at&3]=v.mData;
            if((i&3)==0) {
                auto* destination=reinterpret_cast<__m256i*>(scratch+(at&~3U));
                _mm256_stream_si256(destination,_mm256_load_si256(reinterpret_cast<const __m256i*>(pending)));
                _mm256_stream_si256(destination+1,_mm256_load_si256(reinterpret_cast<const __m256i*>(pending)+1));
            }
        }
    });
    if constexpr(Streaming)_mm_sfence();
    for(std::size_t group=0;group<n/1024;++group)
        k::bchTranspose4Mapped(scratch+group*MappedLayout::stride,out+group*512,permutation+group*256);
}
static void fill(std::vector<k::block>& x) {
    k::setup::Words random(913);
    for(auto& v:x){auto a=random(),b=random();v=k::block(a,b);}
}
static double median(std::vector<double> x){std::sort(x.begin(),x.end());return x[x.size()/2];}
static std::uint64_t hash(const std::vector<k::block>& x,std::size_t n) {
    std::uint64_t h=0;
    for(std::size_t i=0;i<n;++i){std::uint64_t a[2];std::memcpy(a,&x[i],16);for(auto v:a)h=(h^v)*0x100000001b3ULL;}
    return h;
}
// Independent dense transpose reference. Use explicit binary matrix rows,
// not the optimized sparse mixing or zeta circuits.
static void r2Reference(const k::block* in,k::block* out,std::size_t n,const unsigned* masks) {
    using Map=k::Map128S19;
    std::array<k::block,19> state{},next{};
    for(std::size_t epoch=n/128;epoch-->0;) {
        const auto base=128*epoch;
        for(unsigned p=0;p<128;++p) {
            auto value=in[base+p].mData;
            for(unsigned j=0;j<19;++j)if((Map::feedbackColumns[p]>>j)&1)
                value=_mm_xor_si128(value,state[j].mData);
            out[base+p]=k::block(value);
        }
        for(unsigned j=0;j<19;++j) {
            auto value=_mm_setzero_si128();
            const auto row=k::imtR2ReferenceRow(j,masks+4*epoch);
            for(unsigned i=0;i<19;++i)if((row>>i)&1)value=_mm_xor_si128(value,state[i].mData);
            for(unsigned p=0;p<128;++p)if((Map::columns[p]>>j)&1)value=_mm_xor_si128(value,in[base+p].mData);
            next[j]=k::block(value);
        }
        state=next;
    }
}
static SPIN_NOINLINE void independentScatter(const k::block* in,k::block* out,k::block* scratch,
    std::size_t n,const unsigned* masks,const unsigned* offsets,double* phases) {
    auto start=Clock::time_point{};if(phases)start=Clock::now();
    k::innerReverse<k::LocalMap128S19R2>(in,n,masks,[&](std::size_t i,k::block v){scratch[offsets[i]]=v;});
    auto middle=Clock::time_point{};if(phases)middle=Clock::now();
    for(std::size_t g=0;g<n/1024;++g)k::bchTranspose4GfniBlend(scratch+g*MappedLayout::stride,out+g*512);
    if(phases){phases[0]+=std::chrono::duration<double,std::milli>(middle-start).count();phases[1]+=std::chrono::duration<double,std::milli>(Clock::now()-middle).count();}
}
static void r2Experiment(const std::vector<unsigned>& route,unsigned exponent,
                         std::uint64_t seed,unsigned calls,bool profile,const std::string& mode,bool shared) {
    const auto n=route.size();
    k::setup::Words random(2);
    std::vector<unsigned> masks(4*(n/128));
    for(std::size_t i=0;i<masks.size();i+=2) {
        unsigned u;do u=unsigned(random())&((1U<<19)-1);while(!u);
        unsigned v=unsigned(random())&((1U<<19)-1);
        if(std::popcount(u&v)&1)v^=u&-u;
        masks[i]=u;masks[i+1]=v;
    }
    const unsigned schedule=mode=="r2-fused-copy"?40:mode=="r2-repack-copy"?38:mode=="r2-repack-prefetch"?39:mode.starts_with("r2-repack")?36:mode.starts_with("r2-fused")?37:32;
    const bool scatter=mode.starts_with("r2-scatter");
    const IndependentPackets packets(route);
    const auto bundles=(schedule==32 && !scatter)?CanonicalBundles(route):packets.bundles;
    std::vector<unsigned> scattered(n);
    if(scatter)for(std::size_t i=0;i<n;++i) {
        const auto row=route[i]/256;
        scattered[i]=(row/4)*MappedLayout::stride+4*(route[i]%256)+row%4;
    }
    std::vector<k::block> input(n),expected(n),inner(n),forward(n),actual(n),packedReference(n);
    fill(input);expected=input;
    r2Reference(input.data(),inner.data(),n,masks.data());
    k::innerReverse<k::LocalMap128S19R2>(input.data(),n,masks.data(),
        [&](std::size_t i,k::block v){actual[i]=v;});
    if(std::memcmp(actual.data(),inner.data(),n*sizeof(k::block)))throw std::runtime_error("r2 dense inner mismatch");
    k::innerForward<k::LocalMap128S19R2>(n,masks.data(),[&](std::size_t i){return input[i];},forward.data());
    auto lhs=_mm_setzero_si128(),rhs=lhs;
    for(std::size_t i=0;i<n;++i) {
        lhs=_mm_xor_si128(lhs,_mm_and_si128(input[i].mData,forward[i].mData));
        rhs=_mm_xor_si128(rhs,_mm_and_si128(input[i].mData,inner[i].mData));
    }
    if(_mm_movemask_epi8(_mm_cmpeq_epi8(lhs,rhs))!=65535)throw std::runtime_error("r2 inner adjoint mismatch");
    for(std::size_t i=0;i<n;++i)packedReference[packed(route[i])]=inner[i];
    for(std::size_t i=0;i<n;i+=1024)k::bchTranspose4(packedReference.data()+i,expected.data()+i/2);
    const auto scratchBlocks=(n/1024)*MappedLayout::stride;
    std::vector<k::block> scratch(scratchBlocks+4);auto* aligned=alignedScratch(scratch);
    k::workspace_routing::adviseOwned(aligned,scratchBlocks*sizeof(k::block));
    double phases[2]{};
    auto run=[&] {
        if(scatter)independentScatter(input.data(),input.data(),aligned,n,masks.data(),scattered.data(),profile?phases:nullptr);
        else canonicalDirect<4,true,true,k::LocalMap128S19R2>(input.data(),input.data(),aligned,n,
            masks.data(),bundles,profile?phases:nullptr,schedule,nullptr,packets.offsets.data());
    };
    run();
    if(std::memcmp(input.data(),expected.data(),n*sizeof(k::block)))throw std::runtime_error("r2 full encoder or suffix mismatch");
    std::cout<<"r2 checks: dense inner, adjoint, materialized route+BCH, unchanged suffix\n";
    fill(input);phases[0]=phases[1]=0;
    for(unsigned i=0;i<3;++i)run();
    std::vector<double> samples;
    for(unsigned i=0;i<calls;++i){auto begin=Clock::now();run();samples.push_back(std::chrono::duration<double,std::milli>(Clock::now()-begin).count());}
    std::cout<<"exponent,group,columns,shared,mode,seed,calls,median_ms,checksum\n"
        <<exponent<<",4,1,"<<shared<<','<<mode<<','<<seed<<','<<calls<<','<<std::fixed<<std::setprecision(6)<<median(samples)<<','<<std::hex<<hash(input,n)<<'\n';
    if(profile)std::cout<<std::dec<<"phase_means_including_warmups,inner_route,"<<phases[0]/(calls+3)<<",transpose_bch,"<<phases[1]/(calls+3)<<'\n';
}
// Exact diagnostic: rank of all feedback equations B*x_epoch=0, restricted to
// messages in the first row group. A rank deficit yields a zero-state subcode.
static unsigned feedbackRank(const std::vector<unsigned>& r,unsigned group,std::vector<std::uint64_t>* witness=nullptr,unsigned prefixRegions=0) {
    const unsigned variables=128*group,words=(variables+63)/64;
    std::vector<std::vector<std::uint64_t>> basis(variables);
    unsigned rank=0;
    const auto end=prefixRegions?(r.size()/256)*prefixRegions:r.size()-128;
    for(unsigned base=0;base<end;base+=128) {
        std::vector<std::vector<std::uint64_t>> equations(19,std::vector<std::uint64_t>(words));
        bool active=false;
        for(unsigned p=0;p<128;++p) {
            const auto x=r[base+p],row=x/256,c=x%256;
            if(row>=group)continue;
            active=true;
            auto syndrome=k::Map128S19::feedbackColumns[p];
            for(unsigned j=0;j<128;++j)if((k::BchRows[j][c/64]>>(c%64))&1) {
                const unsigned v=row*128+j;
                for(unsigned s=0;s<19;++s)if((syndrome>>s)&1)equations[s][v/64]^=std::uint64_t{1}<<(v%64);
            }
        }
        if(!active)continue;
        for(auto& eq:equations) {
            bool inserted=false;
            for(unsigned w=words;w-->0 && !inserted;) {
                while(eq[w]) {
                    const unsigned pivot=64*w+63-std::countl_zero(eq[w]);
                    if(basis[pivot].empty()){basis[pivot]=eq;++rank;inserted=true;break;}
                    for(unsigned j=0;j<=w;++j)eq[j]^=basis[pivot][j];
                }
            }
            if(rank==variables)return rank;
        }
    }
    if(witness && rank<variables) {
        witness->assign(words,0);
        unsigned free=0;while(!basis[free].empty())++free;
        (*witness)[free/64]|=std::uint64_t{1}<<(free%64);
        for(unsigned pivot=0;pivot<variables;++pivot)if(!basis[pivot].empty()) {
            unsigned parity=0;
            for(unsigned w=0;w<words;++w)parity^=std::popcount(basis[pivot][w]&(*witness)[w])&1;
            if(parity)(*witness)[pivot/64]^=std::uint64_t{1}<<(pivot%64);
        }
        for(const auto& eq:basis)if(!eq.empty()) {
            unsigned parity=0;for(unsigned w=0;w<words;++w)parity^=std::popcount(eq[w]&(*witness)[w])&1;
            if(parity)throw std::runtime_error("nullspace back-substitution");
        }
    }
    return rank;
}
static unsigned verifyWitness(const std::vector<unsigned>& r,unsigned group,std::uint64_t seed,
                              const std::vector<std::uint64_t>& witness,unsigned prefixRegions=0) {
    const auto n=r.size();
    std::vector<std::uint64_t> message(n/128),output(n/64),scratch(n/64),outer(group*4);
    std::copy(witness.begin(),witness.end(),message.begin());
    if(std::all_of(witness.begin(),witness.end(),[](auto x){return x==0;}))throw std::runtime_error("zero witness");
    for(unsigned row=0;row<group;++row)for(unsigned j=0;j<128;++j)
        if((witness[2*row+j/64]>>(j%64))&1)
            for(unsigned w=0;w<4;++w)outer[4*row+w]^=k::BchRows[j][w];
    k::Spin code(k::Configuration::T128S19,k::MessageLength{n/2},seed,2,0,k::BchBackend::Avx512,false,r);
    code.forwardBits(message.data(),message.size(),output.data(),output.size(),scratch.data(),scratch.size());
    unsigned weight=0;
    const auto end=prefixRegions?(n/256)*prefixRegions:n;
    for(unsigned i=0;i<n;++i) {
        const unsigned x=r[i];
        const bool expected=x<256*group && ((outer[x/64]>>(x%64))&1);
        const bool actual=(output[i/64]>>(i%64))&1;
        if(i<end && expected!=actual)throw std::runtime_error("witness inner is not identity on constrained prefix");
        weight+=actual;
    }
    const auto support=prefixRegions?(n-end)+prefixRegions*group:256*group;
    if(!weight || weight>support)throw std::runtime_error("invalid witness weight");
    return weight;
}
// Exact local feedback distribution when a group occupies one aligned window
// in an IMT epoch and its active lanes are a uniformly sampled subset.
static void feedbackCensus(unsigned group) {
    if(!std::has_single_bit(group) || group>16)throw std::invalid_argument("census group must be 1,2,4,8,16");
    std::vector<std::uint32_t> keys;
    keys.reserve((128/group)*(1U<<group));
    for(unsigned base=0;base<128;base+=group) {
        unsigned syndrome=0,previous=0;
        for(unsigned j=0;j<(1U<<group);++j) {
            const unsigned mask=j^(j>>1);
            if(j)syndrome^=k::Map128S19::feedbackColumns[base+std::countr_zero(mask^previous)];
            keys.push_back((std::popcount(mask)<<19)|syndrome);previous=mask;
        }
    }
    std::sort(keys.begin(),keys.end());
    std::vector<unsigned> totals(group+1),zeros(group+1),maximum(group+1);
    for(std::size_t i=0;i<keys.size();) {
        std::size_t end=i+1;while(end<keys.size() && keys[end]==keys[i])++end;
        const unsigned a=keys[i]>>19,count=unsigned(end-i);
        totals[a]+=count;maximum[a]=std::max(maximum[a],count);
        if(!(keys[i]&((1U<<19)-1)))zeros[a]=count;
        i=end;
    }
    std::cout<<"group,active,zero_count,max_syndrome_count,total\n";
    for(unsigned a=0;a<=group;++a)std::cout<<group<<','<<a<<','<<zeros[a]<<','<<maximum[a]<<','<<totals[a]<<'\n';
}
// Enumerate two nonempty supports in distinct aligned windows. These are
// local feedback cancellations, not complete low-weight SPIN codewords.
static void pairCensus(unsigned group) {
    if(!std::has_single_bit(group) || group>8)throw std::invalid_argument("pair census group must be 1,2,4,8");
    struct Entry {unsigned syndrome,window,mask;};
    std::vector<Entry> entries;
    for(unsigned window=0;window<128/group;++window)
        for(unsigned mask=1;mask<(1U<<group);++mask) {
            unsigned syndrome=0;
            for(unsigned j=0;j<group;++j)if(mask&(1U<<j))syndrome^=k::Map128S19::feedbackColumns[window*group+j];
            entries.push_back({syndrome,window,mask});
        }
    std::sort(entries.begin(),entries.end(),[](const Entry& a,const Entry& b){return a.syndrome<b.syndrome;});
    std::array<std::array<unsigned,9>,9> count{};
    unsigned total=0;
    for(std::size_t begin=0;begin<entries.size();) {
        std::size_t end=begin+1;while(end<entries.size() && entries[end].syndrome==entries[begin].syndrome)++end;
        for(auto i=begin;i<end;++i)for(auto j=i+1;j<end;++j)if(entries[i].window!=entries[j].window) {
            unsigned a=std::popcount(entries[i].mask),b=std::popcount(entries[j].mask);
            if(a>b)std::swap(a,b);++count[a][b];++total;
        }
        begin=end;
    }
    std::cout<<"group,active_a,active_b,unordered_zero_feedback_pairs\n";
    for(unsigned a=1;a<=group;++a)for(unsigned b=a;b<=group;++b)if(count[a][b])
        std::cout<<group<<','<<a<<','<<b<<','<<count[a][b]<<'\n';
    std::cout<<"total,"<<total<<'\n';
    if(group==4) {
        std::array<unsigned,13> triples{};
        for(std::size_t i=0;i<entries.size();++i)for(auto j=i+1;j<entries.size();++j) {
            if(entries[i].window==entries[j].window)continue;
            const auto target=entries[i].syndrome^entries[j].syndrome;
            auto it=std::lower_bound(entries.begin()+j+1,entries.end(),target,
                [](const Entry& a,unsigned b){return a.syndrome<b;});
            for(;it!=entries.end() && it->syndrome==target;++it)
                if(it->window!=entries[i].window && it->window!=entries[j].window)
                    ++triples[std::popcount(entries[i].mask)+std::popcount(entries[j].mask)+std::popcount(it->mask)];
        }
        std::cout<<"group,total_weight,unordered_zero_feedback_triples\n";
        for(unsigned a=0;a<triples.size();++a)if(triples[a])std::cout<<group<<','<<a<<','<<triples[a]<<'\n';
        std::cout<<"triple_total,"<<std::accumulate(triples.begin(),triples.end(),0U)<<'\n';
    }
}
// Exact distribution for each support orbit under the 16 XOR translations.
// A support is fixed before choosing the window and the fresh XOR offset.
static void xorCensus() {
    std::array<std::vector<unsigned>,8> syndrome;
    for(unsigned window=0;window<8;++window) {
        auto& table=syndrome[window];table.resize(65536);
        for(unsigned mask=1;mask<65536;++mask)table[mask]=table[mask&(mask-1)]^k::Map128S19::feedbackColumns[window*16+std::countr_zero(mask)];
    }
    std::vector<bool> seen(65536);
    std::array<unsigned,17> orbitCount{},maxZero{},maxAtom{},badSupports{};
    for(unsigned mask=0;mask<65536;++mask)if(!seen[mask]) {
        std::array<unsigned,16> translations{};
        unsigned orbitSize=0;
        for(unsigned shift=0;shift<16;++shift) {
            unsigned translated=0;
            for(unsigned j=0;j<16;++j)if(mask&(1U<<j))translated|=1U<<(j^shift);
            translations[shift]=translated;
            if(!seen[translated]){seen[translated]=true;++orbitSize;}
        }
        std::array<unsigned,128> distribution{};
        for(unsigned window=0;window<8;++window)for(unsigned shift=0;shift<16;++shift)
            distribution[16*window+shift]=syndrome[window][translations[shift]];
        std::sort(distribution.begin(),distribution.end());
        const auto a=std::popcount(mask);++orbitCount[a];
        const auto zeros=unsigned(std::count(distribution.begin(),distribution.end(),0U));
        maxZero[a]=std::max(maxZero[a],zeros);if(zeros)badSupports[a]+=orbitSize;
        for(unsigned i=0;i<128;) {
            unsigned end=i+1;while(end<128 && distribution[end]==distribution[i])++end;
            maxAtom[a]=std::max(maxAtom[a],end-i);i=end;
        }
    }
    if(std::accumulate(orbitCount.begin(),orbitCount.end(),0U)!=4336)throw std::runtime_error("XOR orbit coverage");
    std::cout<<"active,orbits,worst_zero_count,worst_syndrome_count,random_choices,supports_with_zero_probability\n";
    for(unsigned a=0;a<=16;++a)std::cout<<a<<','<<orbitCount[a]<<','<<maxZero[a]<<','<<maxAtom[a]<<",128,"<<badSupports[a]<<'\n';
}
int main(int argc,char** argv) {try {
    if(argc<6)throw std::invalid_argument("usage: locality exponent group shared mode seed [calls] [columns] [prefix_macroregions]");
    const unsigned exponent=std::stoul(argv[1]),group=std::stoul(argv[2]);
    const bool shared=std::stoul(argv[3])!=0;
    const std::string requestedMode=argv[4];const bool xorLanes=requestedMode.starts_with("xor-"),split=requestedMode.starts_with("split-"),canonical=requestedMode.starts_with("block-");
    const std::string mode=xorLanes?requestedMode.substr(4):split?requestedMode.substr(6):canonical?requestedMode.substr(6):requestedMode;const auto seed=std::stoull(argv[5]);
    const unsigned calls=argc>6?std::stoul(argv[6]):31;
    const unsigned columns=argc>7?std::stoul(argv[7]):1;
    const bool gfniCheck=canonical && mode=="gfni-check";
    const bool randomGfni=mode.starts_with("random-gfni");
    const bool profile=mode=="profile" || mode.ends_with("-profile");
    if(randomGfni && (group!=4 || !shared || (columns!=1 && columns!=2 && columns!=4)))throw std::invalid_argument("random GFNI geometry");
    unsigned bchSchedule=mode=="greedy"?1:mode=="release"?2:mode=="reverse"?3:mode=="restrict"?4:mode=="share6"?5:mode=="share8"?6:mode=="gfni1"?7:mode=="gfni2"?8:mode=="gfni-static"?9:mode=="gfni-ternary"?10:mode=="gfni-unroll"?11:(mode=="gfni-blend" || mode=="gfni-profile")?12:0;
    if(randomGfni) {
        const auto selection=profile?mode.substr(0,mode.size()-8):mode;
        if(selection=="random-gfni")bchSchedule=13;
        else if(selection=="random-gfni-repack")bchSchedule=14;
        else if(selection=="random-gfni-prepared")bchSchedule=15;
        else if(selection=="random-gfni-prefetch")bchSchedule=16;
        else if(selection=="random-gfni-copy")bchSchedule=17;
        else if(selection=="random-gfni-prepared-prefetch")bchSchedule=18;
        else if(selection=="random-gfni-copy-prefetch")bchSchedule=19;
        else if(selection=="random-gfni-copy-mapped")bchSchedule=20;
        else if(selection=="random-gfni-manual-copy")bchSchedule=21;
        else if(selection=="random-gfni-copy-mapped-prefetch")bchSchedule=22;
        else if(selection=="random-gfni-copy-packed")bchSchedule=23;
        else if(selection=="random-gfni-manual-packed")bchSchedule=24;
        else if(selection=="random-gfni-vector-packed")bchSchedule=25;
        else if(selection=="random-gfni-vector-mapped")bchSchedule=26;
        else if(selection=="random-gfni-vector-prefetch1")bchSchedule=27;
        else if(selection=="random-gfni-vector-prefetch2")bchSchedule=28;
        else if(selection=="random-gfni-vector-prefetch4")bchSchedule=29;
        else if(selection=="random-gfni-vector-tile1")bchSchedule=30;
        else if(selection=="random-gfni-vector-tile1-prefetch2")bchSchedule=31;
        else if(selection=="random-gfni-physical")bchSchedule=32;
        else if(selection=="random-gfni-physical-copy")bchSchedule=33;
        else if(selection=="random-gfni-physical-prefetch")bchSchedule=34;
        else if(selection=="random-gfni-vector-current-prefetch")bchSchedule=35;
        else throw std::invalid_argument("random GFNI mode");
        if(columns==1 && (bchSchedule<32 || bchSchedule>34))throw std::invalid_argument("single columns require physical cache-line routing");
    }
#if defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
    if((bchSchedule>=7 || gfniCheck || mode.starts_with("r2-")) && (!__builtin_cpu_supports("gfni") || !__builtin_cpu_supports("avx512bw")))throw std::runtime_error("GFNI experiment requires GFNI/AVX512BW");
#endif
    if(gfniCheck) {
        std::array<unsigned,4> order{0,1,2,3};
        do {
            std::uint64_t control=0;
            for(unsigned source=0;source<4;++source)for(unsigned half=0;half<2;++half)
                control|=std::uint64_t(2*source+half)<<(8*(2*order[source]+half));
            for(unsigned source=0;source<4;++source)for(unsigned bit=0;bit<128;++bit) {
                alignas(64) k::block input[4]{},expected[4]{},actual[4];
                input[source]=k::block(bit>=64?std::uint64_t{1}<<(bit-64):0,bit<64?std::uint64_t{1}<<bit:0);
                expected[order[source]]=input[source];
                const auto indices=_mm512_cvtepu8_epi64(_mm_loadl_epi64(reinterpret_cast<const __m128i*>(&control)));
                _mm512_store_si512(actual,_mm512_permutexvar_epi64(indices,_mm512_load_si512(input)));
                if(std::memcmp(actual,expected,sizeof(actual)))throw std::runtime_error("vector route permutation mismatch");
            }
        } while(std::next_permutation(order.begin(),order.end()));
        std::cout<<"All 24 column permutations match on all 512 input basis vectors\n";
        alignas(64) k::block input[1024],shuffled[1024],expected[512],actual[512];
        std::array<unsigned,256> permutation;Sampler random(256,seed);random.shuffle(permutation);
        for(auto& offset:permutation)offset*=4;
        std::memset(input,0,sizeof(input));
        std::memset(shuffled,0,sizeof(shuffled));
        for(unsigned coordinate=0;coordinate<1024;++coordinate)for(unsigned bit=0;bit<128;++bit) {
            input[coordinate]=k::block(bit>=64?std::uint64_t{1}<<(bit-64):0,bit<64?std::uint64_t{1}<<bit:0);
            k::bchTranspose4(input,expected);
            // Explicit calls: dispatch belongs in this correctness test only.
#define CHECK_GFNI(NAME) k::bchTranspose4Gfni##NAME(input,actual); if(std::memcmp(actual,expected,sizeof(actual)))throw std::runtime_error("GFNI basis mismatch: " #NAME)
            CHECK_GFNI(1);CHECK_GFNI(2);CHECK_GFNI(Static);CHECK_GFNI(Ternary);CHECK_GFNI(Unroll);CHECK_GFNI(Blend);
#undef CHECK_GFNI
            const auto at=permutation[coordinate/4]+coordinate%4;shuffled[at]=input[coordinate];
            k::bchTranspose4GfniMapped(shuffled,actual,permutation.data());
            if(std::memcmp(actual,expected,sizeof(actual)))throw std::runtime_error("mapped GFNI basis mismatch");
            k::bchTranspose4GfniMapped1(shuffled,actual,permutation.data());
            if(std::memcmp(actual,expected,sizeof(actual)))throw std::runtime_error("mapped tile1 GFNI basis mismatch");
            k::bchTranspose4GfniPrepared(shuffled,actual,permutation.data());
            if(std::memcmp(actual,expected,sizeof(actual)))throw std::runtime_error("prepared GFNI basis mismatch");
            k::bchTranspose4GfniPacked(shuffled,actual,permutation.data());
            if(std::memcmp(actual,expected,sizeof(actual)))throw std::runtime_error("packed GFNI basis mismatch");
            shuffled[at]=k::block(0,0);
            input[coordinate]=k::block(0,0);
        }
        std::cout<<"All ten GFNI implementations match BCH on all 131072 input basis vectors\n";
        return 0;
    }
    if(canonical && (group!=4 || !shared || (columns!=2 && columns!=4) || (mode!="stream" && mode!="cached" && mode!="profile" && mode!="tiled" && mode!="verify" && mode!="rank" && mode!="prefix-rank" && !bchSchedule)))throw std::invalid_argument("canonical block mode/geometry");
    if(split && (!shared || columns!=1 || (group!=32 && group!=64 && group!=128) || (mode!="stream" && mode!="tiled" && mode!="rank" && mode!="verify")))throw std::invalid_argument("split mode/geometry");
    if(xorLanes && (group!=16 || !shared || columns!=1 || (mode!="wide-stream" && mode!="wide-profile" && mode!="tiled" && mode!="rank" && mode!="prefix-rank" && mode!="verify" && mode!="census")))throw std::invalid_argument("xor experiment geometry/mode");
    if(exponent<14 || exponent>22 || !calls)throw std::invalid_argument("experimental size/call limit");
    if(!k::bchAvx512Available())throw std::runtime_error("requires AVX512 BCH backend");
    const unsigned n=2U<<exponent,rows=n/256;
    if(split && (rows/(group/16))%128)throw std::invalid_argument("split requires epoch-aligned bands");
    if(mode=="census"){if(xorLanes)xorCensus();else feedbackCensus(group);return 0;}
    if(mode=="pair-census"){pairCensus(group);return 0;}
    const unsigned fused=mode=="fused-1k"?1024:mode=="fused-4k"?4096:mode=="fused-16k"?16384:0;
    const bool sequential=mode=="sequential-control" || mode=="sequential-profile" || fused;
    const bool localControl=mode=="local-control" || sequential;
    const unsigned stridedBch=mode=="wide-stride16"?16:mode=="wide-stride20"?20:0;
    const bool wideMode=mode=="wide" || mode=="wide-stream" || mode=="wide-profile" || stridedBch;
    if(stridedBch && group!=16)throw std::invalid_argument("strided BCH group");
    if(wideMode && (!shared || columns!=1 || (group!=8 && group!=16 && group!=32)))throw std::invalid_argument("wide geometry");
    const bool stateGatherMode=mode=="state-gather";
    if(stateGatherMode && (group!=4 || !shared || (columns!=1 && columns!=4 && columns!=8)))throw std::invalid_argument("state-gather needs g=4 shared c=1,4,8");
    const bool columnMode=mode=="column" || mode=="verify-column" || mode=="column-profile";
    if(columnMode && (group!=4 || !shared || columns!=1))throw std::invalid_argument("column experiment uses shared g=4,c=1");
    auto r=split?splitRoute(rows,group,seed):route(rows,group,shared,seed,columns,columnMode,xorLanes,canonical);checkRoute(r,rows,columns);
    if(mode.starts_with("r2-")) {
        const auto selection=profile?mode.substr(0,mode.size()-8):mode;
        if(selection!="r2-physical" && selection!="r2-repack" && selection!="r2-repack-copy" && selection!="r2-repack-prefetch" && selection!="r2-fused-copy" && selection!="r2-fused" && selection!="r2-scatter")throw std::invalid_argument("r2 mode");
        if(group!=4 || (selection=="r2-physical" && !shared) || columns!=1 || split || canonical || xorLanes)throw std::invalid_argument("r2 geometry");
        r2Experiment(r,exponent,seed,calls,profile,selection,shared);return 0;
    }
    if(localControl)std::iota(r.begin(),r.end(),0U); // Invalid-distance cost diagnostic only.
    if(sequential)for(unsigned i=0;i<n;++i)r[i]=(i&~1023U)|((i&3U)<<8)|((i>>2)&255U);
    if(mode=="rank" || mode=="prefix-rank") {
        const unsigned prefix=mode=="prefix-rank"?(argc>8?std::stoul(argv[8]):((128*group-1)/19))*columns:0;
        if(mode=="prefix-rank" && (!prefix || prefix>=256 || group*columns>128 || 128%(group*columns)))throw std::invalid_argument("prefix rank geometry");
        std::vector<std::uint64_t> witness;
        const auto rank=feedbackRank(r,group,&witness,prefix);
        std::cout<<"exponent,group,columns,shared,seed,variables,rank,nullity\n"<<exponent<<','<<group<<','<<columns<<','<<shared<<','<<seed<<','<<128*group<<','<<rank<<','<<128*group-rank<<'\n';
        if(prefix)std::cout<<"zero_state_prefix_regions,"<<prefix<<",support_bound,"<<((256-prefix)*(n/256)+prefix*group)<<'\n';
        if(!witness.empty())std::cout<<(prefix?"verified_late_activation_witness_weight,":"verified_zero_state_witness_weight,")<<verifyWitness(r,group,seed,witness,prefix)<<",length,"<<n<<'\n';
        return 0;
    }
    if(mode!="tiled" && mode!="direct" && mode!="verify" && mode!="line" && mode!="stream" && mode!="mapped" && mode!="mapped-stream" && mode!="plain-bch" && mode!="mapped-bch" && mode!="bundled" && mode!="repack-bch" && !localControl && !columnMode && !stateGatherMode && !wideMode && !canonical && !randomGfni)throw std::invalid_argument("mode");
    const bool lines=group==4 && shared;
    if((mode=="line" || (mode=="stream" && !split) || mode=="mapped" || mode=="mapped-stream" || mode=="mapped-bch" || mode=="bundled" || mode=="repack-bch") && !lines)throw std::invalid_argument("line kernels need g=4, shared shuffle");
    k::Spin code(k::Configuration::T128S19,k::MessageLength{n/2},seed,2,0,k::BchBackend::Avx512,false,r);
    // The library's validateSetup also requires its original 256-region law.
    // Rectangular routes deliberately replace that law; checkRoute validates
    // their bijection and macroregion occupancy before the same linear kernels.
    if(columns==1 && !localControl)code.validateSetup();
    if(group==1 && columns==1 && !localControl) {
        k::Spin original(k::Configuration::T128S19,k::MessageLength{n/2},seed,2,0,k::BchBackend::Avx512);
        if(code.routeHash()!=original.routeHash())throw std::runtime_error("g=1 differs from original seed-to-route");
    }
    std::vector<unsigned> offsets(n);std::transform(r.begin(),r.end(),offsets.begin(),packed);
    const auto gather=stateGatherMode?gatherSources(offsets,columns):std::vector<unsigned>{};
    const WideLayout wide=wideMode?WideLayout(r,group,stridedBch==20?4:0):WideLayout(route(rows,8,true,seed),8);
    if(wideMode && xorLanes)for(unsigned i=0;i<n;i+=16)for(unsigned j=0;j<16;++j)
        if(wide.lanes[i+j]!=(j^wide.lanes[i]))throw std::runtime_error("non-XOR lane route");
    // Reserve the same additional buffers for controls as for mapped candidates.
    MappedLayout mapped=lines?MappedLayout(r,rows,columns):MappedLayout(route(rows,4,true,seed,columns),rows,columns);
    const auto canonicalRoute=(bchSchedule>=32 && bchSchedule<=34)?CanonicalBundles(r):randomGfni?CanonicalBundles(mapped,4*columns):canonical?CanonicalBundles(r,4*columns):CanonicalBundles{};
    const auto compactRoute=columns<=8?bundles(mapped,columns):std::vector<RouteBundle>{};
    std::array<unsigned,256> columnOrder{};
    if(columnMode)for(unsigned c=0;c<256;++c) {
        columnOrder[c]=r[c*rows]%256;
        for(unsigned j=0;j<rows;++j)if(r[c*rows+j]%256!=columnOrder[c])throw std::runtime_error("non-global column order");
    }
    if(lines)for(unsigned i=0;i<n;i+=4) {
        unsigned lanes=0;
        for(unsigned j=0;j<4;++j){if((offsets[i+j]&~3U)!=(offsets[i]&~3U))throw std::runtime_error("not a cache-line group");lanes|=1U<<(offsets[i+j]&3);}
        if(lanes!=15)throw std::runtime_error("incomplete cache-line group");
    }
    if(mode=="verify" || mode=="verify-column") {
        std::vector<k::block> input(n),expected(n),actual(n),message(n/2),forward(n);
        fill(input);fill(message);actual=input;
        k::Spin::Workspace w(code);
        code.reference(input.data(),expected.data());
        code.encodeInplace(actual.data(),n,w);
        if(std::memcmp(actual.data(),expected.data(),n/2*sizeof(k::block)))throw std::runtime_error("dense transpose mismatch");
        code.forward(message.data(),n/2,forward.data(),n,w);
        __m128i lhs=_mm_setzero_si128(),rhs=lhs;
        for(unsigned i=0;i<n;++i)lhs=_mm_xor_si128(lhs,_mm_and_si128(input[i].mData,forward[i].mData));
        for(unsigned i=0;i<n/2;++i)rhs=_mm_xor_si128(rhs,_mm_and_si128(message[i].mData,actual[i].mData));
        if(_mm_movemask_epi8(_mm_cmpeq_epi8(lhs,rhs))!=65535)throw std::runtime_error("adjoint mismatch");
        std::vector<k::block> scratch(n);
        direct(input.data(),input.data(),scratch.data(),n,offsets.data(),code.wideView().fieldRows);
        if(std::memcmp(input.data(),actual.data(),n*sizeof(k::block)))throw std::runtime_error("direct/suffix mismatch");
        if(lines) {
            scratch.resize((n/1024)*MappedLayout::stride+4);auto* aligned=alignedScratch(scratch);
            fill(input);lineDirect<false>(input.data(),input.data(),aligned,n,offsets.data(),code.wideView().fieldRows);
            if(std::memcmp(input.data(),actual.data(),n*sizeof(k::block)))throw std::runtime_error("line mismatch");
            fill(input);lineDirect<true>(input.data(),input.data(),aligned,n,offsets.data(),code.wideView().fieldRows);
            if(std::memcmp(input.data(),actual.data(),n*sizeof(k::block)))throw std::runtime_error("stream mismatch");
            if(canonical)for(bool streaming:{false,true}) {
                fill(input);runCanonical(columns,streaming,input.data(),input.data(),aligned,n,code.wideView().fieldRows,canonicalRoute);
                if(std::memcmp(input.data(),actual.data(),n*sizeof(k::block)))throw std::runtime_error("canonical bundle mismatch");
            }
            fill(input);mappedDirect<false>(input.data(),input.data(),aligned,n,mapped.offsets.data(),code.wideView().fieldRows,mapped.permutation.data());
            if(std::memcmp(input.data(),actual.data(),n*sizeof(k::block)))throw std::runtime_error("mapped mismatch");
            fill(input);mappedDirect<true>(input.data(),input.data(),aligned,n,mapped.offsets.data(),code.wideView().fieldRows,mapped.permutation.data());
            if(std::memcmp(input.data(),actual.data(),n*sizeof(k::block)))throw std::runtime_error("mapped streaming mismatch");
            if(columns<=8) {
                fill(input);runBundled(columns,input.data(),input.data(),aligned,n,compactRoute.data(),code.wideView().fieldRows,mapped.permutation.data());
                if(std::memcmp(input.data(),actual.data(),n*sizeof(k::block)))throw std::runtime_error("bundled mismatch");
            }
            if(columnMode) {
                fill(input);columnDirect(input.data(),input.data(),aligned,w.tile.data(),n,r.data(),code.wideView().fieldRows,columnOrder.data());
                if(std::memcmp(input.data(),actual.data(),n*sizeof(k::block)))throw std::runtime_error("column-major mismatch");
            }
        }
        std::cout<<"verified dense transpose, adjoint, direct route and suffix: "<<exponent<<','<<group<<','<<columns<<','<<shared<<','<<seed<<'\n';return 0;
    }
    code.compact();k::Spin::Workspace workspace(code);
    const auto splitLayout=split?SplitLayout(r,group):SplitLayout{};
    std::vector<k::block> pending(rows+4);auto* pendingBuffer=alignedScratch(pending);
    std::vector<k::block> states(20*(n/128)+4);auto* stateBuffer=alignedScratch(states);
    const auto scratchBlocks=std::max((n/1024)*MappedLayout::stride,(n/4096)*(256*20+4));
    std::vector<k::block> input(n),expected(n),scratch(scratchBlocks+4);fill(input);expected=input;
    auto* aligned=alignedScratch(scratch);
    k::workspace_routing::adviseOwned(aligned,scratchBlocks*sizeof(k::block));
    const auto* masks=code.wideView().fieldRows;
    code.encodeInplace(expected.data(),n,workspace);
    if(canonical || randomGfni)runCanonical(columns,mode!="cached",input.data(),input.data(),aligned,n,masks,canonicalRoute,nullptr,bchSchedule,mapped.permutation.data());
    else if(split && mode=="stream")runSplit(group,input.data(),input.data(),aligned,workspace.tile.data(),pendingBuffer,n,masks,splitLayout);
    else if(mode=="wide")runWide<false>(group,input.data(),input.data(),aligned,workspace.tile.data(),n,masks,wide);
    else if(wideMode)runWide<true>(group,input.data(),input.data(),aligned,workspace.tile.data(),n,masks,wide,nullptr,xorLanes,stridedBch);
    else if(stateGatherMode)runStateGather(columns,input.data(),input.data(),aligned,workspace.tile.data(),stateBuffer,n,masks,offsets.data(),gather.data());
    else if(fused)runFused(fused,input.data(),input.data(),aligned,workspace.tile.data(),n,masks);
    else if(sequential)sequentialControl(input.data(),input.data(),aligned,n,masks);
    else if(columnMode)columnDirect(input.data(),input.data(),aligned,workspace.tile.data(),n,r.data(),masks,columnOrder.data());
    else if(mode=="bundled" || mode=="repack-bch")runBundled(columns,input.data(),input.data(),aligned,n,compactRoute.data(),masks,mapped.permutation.data());
    else if(mode=="mapped" || mode=="mapped-bch")mappedDirect<false>(input.data(),input.data(),aligned,n,mapped.offsets.data(),masks,mapped.permutation.data());
    else if(mode=="mapped-stream")mappedDirect<true>(input.data(),input.data(),aligned,n,mapped.offsets.data(),masks,mapped.permutation.data());
    else if(mode=="stream")lineDirect<true>(input.data(),input.data(),aligned,n,offsets.data(),masks);
    else if(mode=="line")lineDirect<false>(input.data(),input.data(),aligned,n,offsets.data(),masks);
    else direct(input.data(),input.data(),aligned,n,offsets.data(),masks);
    if(std::memcmp(input.data(),expected.data(),n*sizeof(k::block)))throw std::runtime_error("candidate differs from materialized reference");
    fill(input);
    double phases[2]{};
    auto run=[&] {
        if((canonical && mode!="tiled") || randomGfni)runCanonical(columns,mode!="cached",input.data(),input.data(),aligned,n,masks,canonicalRoute,profile?phases:nullptr,bchSchedule,mapped.permutation.data());
        else if(split && mode=="stream")runSplit(group,input.data(),input.data(),aligned,workspace.tile.data(),pendingBuffer,n,masks,splitLayout);
        else if(mode=="wide")runWide<false>(group,input.data(),input.data(),aligned,workspace.tile.data(),n,masks,wide);
        else if(wideMode)runWide<true>(group,input.data(),input.data(),aligned,workspace.tile.data(),n,masks,wide,mode=="wide-profile"?phases:nullptr,xorLanes,stridedBch);
        else if(stateGatherMode)runStateGather(columns,input.data(),input.data(),aligned,workspace.tile.data(),stateBuffer,n,masks,offsets.data(),gather.data());
        else if(fused)runFused(fused,input.data(),input.data(),aligned,workspace.tile.data(),n,masks);
        else if(sequential)sequentialControl(input.data(),input.data(),aligned,n,masks,mode=="sequential-profile"?phases:nullptr);
        else if(columnMode)columnDirect(input.data(),input.data(),aligned,workspace.tile.data(),n,r.data(),masks,columnOrder.data(),mode=="column-profile"?phases:nullptr);
        else if(mode=="plain-bch")for(std::size_t i=0;i<n;i+=1024)k::bchTranspose4(aligned+i,input.data()+i/2);
        else if(mode=="repack-bch")repackBch(aligned,input.data(),n,mapped.permutation.data());
        else if(mode=="bundled")runBundled(columns,input.data(),input.data(),aligned,n,compactRoute.data(),masks,mapped.permutation.data());
        else if(mode=="mapped-bch")for(std::size_t b=0;b<n/1024;++b)k::bchTranspose4Mapped(aligned+b*MappedLayout::stride,input.data()+b*512,mapped.permutation.data()+b*256);
        else if(mode=="direct" || localControl)direct(input.data(),input.data(),aligned,n,offsets.data(),masks);
        else if(mode=="mapped")mappedDirect<false>(input.data(),input.data(),aligned,n,mapped.offsets.data(),masks,mapped.permutation.data());
        else if(mode=="mapped-stream")mappedDirect<true>(input.data(),input.data(),aligned,n,mapped.offsets.data(),masks,mapped.permutation.data());
        else if(mode=="line")lineDirect<false>(input.data(),input.data(),aligned,n,offsets.data(),masks);
        else if(mode=="stream")lineDirect<true>(input.data(),input.data(),aligned,n,offsets.data(),masks);
        else code.encodeInplace(input.data(),n,workspace);
    };
    for(unsigned j=0;j<3;++j)run();
    std::vector<double> samples;
    for(unsigned j=0;j<calls;++j){auto a=Clock::now();run();auto b=Clock::now();samples.push_back(std::chrono::duration<double,std::milli>(b-a).count());}
    std::cout<<"exponent,group,columns,shared,mode,seed,calls,median_ms,checksum\n";
    std::cout<<exponent<<','<<group<<','<<columns<<','<<shared<<','<<requestedMode<<','<<seed<<','<<calls<<','<<std::fixed<<std::setprecision(6)<<median(samples)<<','<<std::hex<<hash(input,n)<<'\n';
    if(mode=="column-profile" || mode=="sequential-profile" || mode=="wide-profile" || ((canonical || randomGfni) && profile))std::cout<<std::dec<<"phase_means_including_warmups,inner_route,"<<phases[0]/(calls+3)<<",transpose_bch,"<<phases[1]/(calls+3)<<'\n';
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
