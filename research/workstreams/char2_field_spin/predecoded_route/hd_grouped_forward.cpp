#include "HdGroupedOuter.h"
#include "HdGroupedInner.h"
#include "HdGroupedWideInner.h"
#include "HdSumrankInner.h"
#include "BenchPayload.h"
#include "Punctured8Oracle.h"
#include <atomic>
#include <iomanip>
#include <iostream>
#ifndef HDGROUP_PREFETCH
#define HDGROUP_PREFETCH 1
#endif
#ifndef HDGROUP_OUTER_TILE
#define HDGROUP_OUTER_TILE 64
#endif
static_assert(!HDGROUP_WIDE_INNER||HDGROUP_SUBFIELD,"Wide striped inner requires GF256 refresh");
static_assert(!HDGROUP_XOR_REFRESH||(HDGROUP_SUBFIELD&&!HDGROUP_WIDE_INNER),
              "The XOR-refresh research candidate is GF256/h4 only");

namespace c2=spin::research::char2;
namespace sr=c2::hdsumrank;
using Block=c2::Block;
using Clock=std::chrono::steady_clock;
static void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}

// This is only a representation change of a sampled 64-record route.
static void decodePacketRoute(const std::uint16_t* packed,std::uint16_t* decoded){
    auto* masks=reinterpret_cast<std::uint8_t*>(decoded+32);
    for(unsigned j=0;j!=16;++j){
        const unsigned a=packed[2*j],b=packed[2*j+1];
        decoded[2*j]=std::uint16_t((a&0x7fffU)<<HDGROUP_ROUTE_ADDRESS_SHIFT);
        decoded[2*j+1]=std::uint16_t((b&0x7fffU)<<HDGROUP_ROUTE_ADDRESS_SHIFT);
        masks[j]=std::uint8_t(((0U-(a>>15))&0x0fU)|((0U-(b>>15))&0xf0U));
    }
}

// Isolated 35% research candidate; no production/default profile changes.
// MT is reproducible benchmark setup, not the ideal independent randomness
// in the probability theorem. Setup and all allocations precede encoding.
template<unsigned Group=4>
struct Plan {
#if HDGROUP_LONG_RS
    static_assert(Group==1,"Long RS has one constituent per group");
#elif HDGROUP_TOEPLITZ_COMPACT7
    static_assert(Group==4,"Compact7 supports only four-symbol Toeplitz columns");
#elif HDGROUP_XOR_MIXER==5
    static_assert(Group==16,"The GF16-coset mixer has exactly sixteen constituents");
#elif HDGROUP_XOR_MIXER>=3
    static_assert(Group==8,"The two-sum mixer has exactly eight constituents");
#else
    static_assert(Group==2||Group==4||Group==8||Group==16);
#endif
    std::size_t k,n,rows,stride;
    std::vector<std::uint16_t> route,diagonal,outerBefore;
    std::vector<std::uint16_t> stripedRoute,packetRoute;
    std::vector<std::uint16_t> decodedRoute;
    std::array<std::uint8_t,256> pairParity{};
    std::vector<c2::hdgroup::Coefficient> prepared,preparedBefore;
    std::vector<std::array<std::uint16_t,2*Group-1>> toeplitz;
    std::size_t toeplitzRejected=0;
    using Refresh=std::conditional_t<HDGROUP_SUBFIELD,sr::ByteDiagonalRefresh,sr::DiagonalRefresh>;
    std::vector<Refresh> refresh;
    std::vector<sr::ByteDiagonalRefresh> stripedRefresh;
    std::vector<std::array<std::uint16_t,16>> before,after;
    static bool nonsingular(const std::array<std::uint16_t,2*Group-1>& t){
        std::uint16_t a[Group][Group];
        for(unsigned i=0;i!=Group;++i)for(unsigned j=0;j!=Group;++j)a[i][j]=t[Group-1+i-j];
        const auto mul=[](std::uint16_t x,std::uint16_t y)->std::uint16_t{
            if constexpr(HDGROUP_SUBFIELD)return c2::p8Oracle().mul[x][y];
            else return c2::medium32::quadraticMultiply(x,y);
        };
        // Fraction-free elimination. Only setup uses this scalar matrix loop.
        for(unsigned c=0;c!=Group;++c){
            auto p=c;while(p!=Group&&!a[p][c])++p;if(p==Group)return false;
            for(unsigned j=c;j!=Group;++j)std::swap(a[c][j],a[p][j]);
            for(unsigned i=c+1;i!=Group;++i){
                const auto f=a[i][c];if(!f)continue;
                for(unsigned j=c+1;j!=Group;++j)a[i][j]=mul(a[c][c],a[i][j])^mul(f,a[c][j]);
                a[i][c]=0;
            }
        }
        return true;
    }
    explicit Plan(std::size_t size,std::uint64_t seed):k(size),n(2*size),rows(size/128),stride(rows+4){
        if(size<4096||size%4096||rows%64||rows>65536)throw std::invalid_argument("common-refresh research dimensions");
        if constexpr(HDGROUP_PAIR_SWAPS){
            if(rows>32768)throw std::invalid_argument("pair-swap route bases require at most 32768 rows");
        }
        if constexpr(HDGROUP_ROUTE_ADDRESS_SHIFT){
            if(rows>(65536U>>HDGROUP_ROUTE_ADDRESS_SHIFT))
                throw std::invalid_argument("pre-scaled research route exceeds 16-bit indices");
        }
        std::mt19937_64 rr(seed^0x91af823017ccb91aULL),dr(seed^0xe685639417f293b5ULL),tr(seed^0xb2d174839ae61713ULL);
        std::uniform_int_distribution<unsigned> nz(1,HDGROUP_SUBFIELD?255:65535);
#if HDGROUP_LONG_RS
        // New outer ensemble: deterministic RS[256,128], with one independent
        // uniform nonzero GF256 scalar on every output coordinate/group.
        // No D1 gauge or column mixer. Routing and inner setup streams stay
        // independent and unchanged. The packed hot keys are group-major.
        route.resize(n);diagonal.resize(n);prepared.resize(n);
        for(auto& d:diagonal)d=std::uint16_t(nz(dr));
        for(std::size_t g=0;g!=rows;g+=4)for(unsigned c=0;c!=256;++c)for(unsigned lane=0;lane!=4;++lane)
            prepared[g*256+4*c+lane]=c2::hdgroup::prepareCoefficient(diagonal[c*rows+g+lane]);
#else
        route.resize(n);preparedBefore.resize(rows*c2::hdgroup::beforeRegions<Group>);
        if constexpr(HDGROUP_TOEPLITZ){
            static_assert(Group==2||Group==4);
            std::uniform_int_distribution<unsigned> field(0,HDGROUP_SUBFIELD?255:65535);
            toeplitz.resize(rows*(256/Group));
            for(auto& t:toeplitz){
                do{for(auto& a:t)a=std::uint16_t(field(dr));if(nonsingular(t))break;++toeplitzRejected;}while(true);
            }
            constexpr auto products=c2::hdgroup::toeplitzProducts<Group>;
            for(std::size_t g=0;g!=rows;g+=4)for(unsigned c=0;c!=256/Group;++c)for(unsigned lane=0;lane!=4;++lane){
                const auto& t=toeplitz[c*rows+g+lane];std::array<std::uint16_t,products> coeff;
                const auto pair=[&](unsigned i,std::uint16_t a,std::uint16_t b,std::uint16_t d){coeff[i]=a;coeff[i+1]=a^b;coeff[i+2]=a^d;};
                if constexpr(Group==2)pair(0,t[1],t[0],t[2]);
                else{pair(0,t[3],t[2],t[4]);pair(3,t[5]^t[3],t[4]^t[2],t[6]^t[4]);pair(6,t[1]^t[3],t[0]^t[2],t[2]^t[4]);}
#if HDGROUP_TOEPLITZ_COMPACT7
                // Keep the sampler and its rejection sequence unchanged.
                // The two omitted coefficients are reconstructed after expansion.
                constexpr unsigned stored=c2::hdgroup::toeplitzStoredCoefficients<Group>;
                constexpr std::array<unsigned,stored> selected{0,1,2,3,5,6,7};
                for(unsigned i=0;i!=stored;++i)preparedBefore[g*c2::hdgroup::beforeRegions<Group>+4*(c*stored+i)+lane]=c2::hdgroup::prepareCoefficient(coeff[selected[i]]);
#else
                for(unsigned i=0;i!=products;++i)preparedBefore[g*c2::hdgroup::beforeRegions<Group>+4*(c*products+i)+lane]=c2::hdgroup::prepareCoefficient(coeff[i]);
#endif
            }
        }else{
        if constexpr(HDGROUP_XOR_MIXER!=1){
            diagonal.resize(n);prepared.resize(n);
            for(auto& d:diagonal)d=std::uint16_t(nz(dr));
        }
        outerBefore.resize(n);
        for(auto& d:outerBefore)d=std::uint16_t(nz(dr));
        if constexpr(HDGROUP_NORMALIZE&&!(HDGROUP_FIELD32&&Group==2)){
            // Divide D1 by its first entry and multiply D2 by that entry.
            // This preserves each sampled matrix, not just its distribution.
            std::array<std::uint16_t,65536> inverses{};
            auto inverse=[](std::uint16_t a){
                std::uint16_t r=1;for(unsigned e=65534;e;e>>=1){if(e&1)r=c2::medium32::quadraticMultiply(r,a);a=c2::medium32::quadraticMultiply(a,a);}return r;
            };
            for(unsigned c=0;c!=256;c+=Group)for(std::size_t g=0;g!=rows;++g){
                const auto a=outerBefore[c*rows+g];
                if(!inverses[a])inverses[a]=inverse(a);
                for(unsigned r=0;r!=Group;++r){
                    outerBefore[(c+r)*rows+g]=c2::medium32::quadraticMultiply(outerBefore[(c+r)*rows+g],inverses[a]);
                    diagonal[(c+r)*rows+g]=c2::medium32::quadraticMultiply(diagonal[(c+r)*rows+g],a);
                }
            }
        }
        // Consume all coefficients of four groups contiguously. Region-major
        // coefficient planes cause pathological set conflicts in this kernel.
        for(std::size_t g=0;g!=rows;g+=4)for(unsigned c=0;c!=256;++c)for(unsigned lane=0;lane!=4;++lane){
            if constexpr(HDGROUP_NORMALIZE&&!(HDGROUP_FIELD32&&Group==2)){
                if(c%Group)preparedBefore[g*c2::hdgroup::beforeRegions<Group>+4*((c/Group)*(Group-1)+c%Group-1)+lane]=c2::hdgroup::prepareCoefficient(outerBefore[c*rows+g+lane]);
            }else preparedBefore[g*256+4*c+lane]=c2::hdgroup::prepareCoefficient(outerBefore[c*rows+g+lane]);
            if constexpr(HDGROUP_XOR_MIXER!=1)prepared[g*256+4*c+lane]=c2::hdgroup::prepareCoefficient(diagonal[c*rows+g+lane]);
        }
        if constexpr(HDGROUP_FIELD32&&Group==2){
            for(std::size_t g=0;g!=rows;++g)for(unsigned c=0;c!=128;++c){
                std::uint32_t scalar;do{scalar=std::uint32_t(dr());}while(!scalar);
                const auto a=std::uint16_t(scalar),b=std::uint16_t(scalar>>16);
                outerBefore[(2*c)*rows+g]=a;outerBefore[(2*c+1)*rows+g]=b;
                const auto i=(g/4)*1024+8*c+g%4;
                preparedBefore[i]=c2::hdgroup::prepareCoefficient(a);
                preparedBefore[i+4]=c2::hdgroup::prepareCoefficient(b);
                prepared[i]=c2::hdgroup::prepareCoefficient(c2::medium32::quadraticMultiply(0x2000,b));
                prepared[i+4]=c2::hdgroup::prepareCoefficient(a^b);
            }
        }
        }
#endif
        if constexpr(HDGROUP_PACKET_ROUTE){
            // Independent uniform packet permutations in each region. Keep
            // the fully expanded logical route for the independent oracle;
            // the hot kernel consumes only the compact packet-base table.
            constexpr unsigned width=HDGROUP_ROUTE_WIDTH;
            const auto packets=rows/width;packetRoute.resize(n/width);
            if constexpr(HDGROUP_SHIFTED_PAIRS==1){
                // Separate setup stream leaves packet permutations unchanged
                // when toggling the independently sampled regional parity.
                std::mt19937_64 parityRandom(seed^0x643aa56472f93731ULL);
                std::uniform_int_distribution<unsigned> parity(0,1);
                for(auto& p:pairParity)p=std::uint8_t(parity(parityRandom));
            }else if constexpr(HDGROUP_SHIFTED_PAIRS==2){
                // Deterministic balanced schedule; consumes no parity RNG.
                for(unsigned c=0;c!=256;++c)pairParity[c]=std::uint8_t(c&1);
            }
#if HDGROUP_PAIR_SWAPS
            // A separate reproducible research stream preserves the packet
            // permutations and all coefficient streams when swaps toggle.
            // The mathematical ensemble uses independent fair setup bits.
            std::mt19937_64 swapRandom(seed^0x28622e73287794afULL);
            std::uniform_int_distribution<unsigned> swapBit(0,1);
#endif
            for(unsigned c=0;c!=256;++c){
                auto* p=packetRoute.data()+c*packets;
                for(std::size_t j=0;j!=packets;++j)p[j]=std::uint16_t(width*j+pairParity[c]);
                for(auto j=packets-1;j;--j){std::uniform_int_distribution<std::size_t> d(0,j);std::swap(p[j],p[d(rr)]);}
#if HDGROUP_PAIR_SWAPS
                for(std::size_t j=0;j!=packets;++j)p[j]|=std::uint16_t(swapBit(swapRandom)<<15);
#endif
                for(std::size_t first=0;first!=rows;first+=64)for(unsigned j=0;j!=16;++j)for(unsigned lane=0;lane!=4;++lane){
#if HDGROUP_PAIR_SWAPS
                    const unsigned entry=p[first/width+(4/width)*j+lane/width];
                    const auto record=std::size_t(entry&0x7fffU)+((lane%width)^(entry>>15));
#else
                    const auto record=std::size_t(p[first/width+(4/width)*j+lane/width])+lane%width;
#endif
                    route[c*rows+first+16*lane+j]=std::uint16_t(HDGROUP_SHIFTED_PAIRS?record%rows:record);
                }
            }
        }else{
            for(unsigned c=0;c!=256;++c){
                auto* p=route.data()+c*rows;std::iota(p,p+rows,std::uint16_t(0));
                for(auto j=rows-1;j;--j){std::uniform_int_distribution<std::size_t> d(0,j);std::swap(p[j],p[d(rr)]);}
            }
        }
        if constexpr(HDGROUP_DECODED_ROUTE){
            static_assert(!HDGROUP_DECODED_ROUTE||!HDGROUP_WIDE_INNER);
            decodedRoute.resize((n/64)*40);
            for(std::size_t first=0;first!=n;first+=64)
                decodePacketRoute(packetRoute.data()+first/2,decodedRoute.data()+(first/64)*40);
        }
        if constexpr((HDGROUP_INNER>=2||HDGROUP_WIDE_INNER)&&!HDGROUP_PACKET_ROUTE){
            stripedRoute.resize(n);
            for(std::size_t first=0;first!=n;first+=64)for(unsigned j=0;j!=16;++j)for(unsigned lane=0;lane!=4;++lane)
                stripedRoute[first+4*j+lane]=route[first+16*lane+j];
        }
        refresh.resize(n/64);before.resize(n/64);after.resize(n/64);
        if constexpr(HDGROUP_INNER>=4||HDGROUP_WIDE_INNER){
            static_assert(HDGROUP_INNER<4||HDGROUP_SUBFIELD,"Persistent striped state is a byte-field research option");
            stripedRefresh.resize(n/64);
        }
        for(std::size_t i=0;i!=n/64;++i){
            for(unsigned j=0;j!=16;++j){before[i][j]=std::uint16_t(nz(tr));after[i][j]=std::uint16_t(nz(tr));}
            refresh[i].prepare(before[i].data(),after[i].data());
            if constexpr(HDGROUP_INNER>=4||HDGROUP_WIDE_INNER){
                std::array<std::uint16_t,16> a,b;
                for(unsigned j=0;j!=16;++j){a[j]=before[i][4*(j%4)+j/4];b[j]=after[i][4*(j%4)+j/4];}
                stripedRefresh[i].prepare(a.data(),b.data());
            }
        }
    }
    std::size_t scratch()const{return 256*stride;}
    std::size_t oddPairRegions()const{return std::count(pairParity.begin(),pairParity.end(),std::uint8_t(1));}
    std::size_t swappedPairs()const{
        if constexpr(HDGROUP_PAIR_SWAPS)return std::count_if(packetRoute.begin(),packetRoute.end(),[](std::uint16_t p){return (p&0x8000U)!=0;});
        else return 0;
    }
};
// The common-refresh instantiation never reads a mixer. An empty type makes
// that compile-time property explicit without allocating unused coefficients.
struct NoMixer {};

// Rank-two controls retain the established kernel and use this explicit
// output conversion. Rank-four native output is stored directly in registers.
template<bool Stream,unsigned Rank,class Refresh>
static inline void innerNativeControl(const Block* input,const std::uint16_t* route,
                                     Block* output,sr::State& state,const Refresh& refresh){
    alignas(64) std::array<Block,64> natural;
    sr::mdsStriped<true,std::uint16_t,true,Rank,Rank,2,false,false>(input,route,natural.data(),state,NoMixer{},refresh);
    c2::c8Unroll(std::make_index_sequence<4>{},[&](auto g){
        constexpr unsigned j=4*g;
        auto a=_mm512_load_si512(natural.data()+j),b=_mm512_load_si512(natural.data()+16+j);
        auto c=_mm512_load_si512(natural.data()+32+j),d=_mm512_load_si512(natural.data()+48+j);
        sr::transposeRecords4(a,b,c,d);
        if constexpr(Stream){
            _mm512_stream_si512(reinterpret_cast<__m512i*>(output+4*j),a);
            _mm512_stream_si512(reinterpret_cast<__m512i*>(output+4*j+4),b);
            _mm512_stream_si512(reinterpret_cast<__m512i*>(output+4*j+8),c);
            _mm512_stream_si512(reinterpret_cast<__m512i*>(output+4*j+12),d);
        }else{
            _mm512_storeu_si512(output+4*j,a);_mm512_storeu_si512(output+4*j+4,b);
            _mm512_storeu_si512(output+4*j+8,c);_mm512_storeu_si512(output+4*j+12,d);
        }
    });
}

template<bool Stream,unsigned Rank=4,unsigned Group>
HD35_NOINLINE static void inner(const Plan<Group>& p,const Block* planes,Block* output){
    static_assert(Rank==2||Rank==4||(HDGROUP_WIDE_INNER&&Rank==16));
    sr::State state{};
    for(unsigned region=0;region!=256;++region){
        const auto begin=region*p.rows;
        const auto* source=planes+region*p.stride;
        for(std::size_t j=0;j!=p.rows;j+=64){
            if(HDGROUP_PREFETCH&&region+1<256){
                const auto* next=reinterpret_cast<const char*>(planes+(region+1)*p.stride+j);
                for(unsigned h=0;h!=16;++h)_mm_prefetch(next+64*h,_MM_HINT_T1);
            }
            if constexpr(Rank==16){
                const auto* route=HDGROUP_PACKET_ROUTE?p.packetRoute.data()+(begin+j)/HDGROUP_ROUTE_WIDTH:
                    p.stripedRoute.data()+begin+j;
                c2::hdgroup::innerPacket16<Stream>(source,route,output+begin+j,state,NoMixer{},p.stripedRefresh[(begin+j)/64]);
            }else if constexpr(HDGROUP_INNER&&Rank==4){
                const auto* route=HDGROUP_DECODED_ROUTE?p.decodedRoute.data()+((begin+j)/64)*40:
                    HDGROUP_PACKET_ROUTE?p.packetRoute.data()+(begin+j)/HDGROUP_ROUTE_WIDTH:
                    (HDGROUP_INNER>=2?p.stripedRoute.data():p.route.data())+begin+j;
                if constexpr(HDGROUP_INNER>=4)c2::hdgroup::innerPacket4<Stream>(source,route,output+begin+j,state,NoMixer{},p.stripedRefresh[(begin+j)/64]);
                else c2::hdgroup::innerPacket4<Stream>(source,route,output+begin+j,state,NoMixer{},p.refresh[(begin+j)/64]);
            }else if constexpr(HDGROUP_NATIVE_OUTPUT){
                innerNativeControl<Stream,Rank>(source,p.route.data()+begin+j,output+begin+j,state,p.refresh[(begin+j)/64]);
            }else sr::mdsStriped<true,std::uint16_t,true,Rank,Rank,2,false,Stream>(source,
                    p.route.data()+begin+j,output+begin+j,state,NoMixer{},p.refresh[(begin+j)/64]);
        }
    }
    if constexpr(Stream)_mm_sfence();
}

#if HDGROUP_COSET_PIPELINE && !HDGROUP_LONG_RS
// This range kernel deliberately leaves the original full inner untouched.
// State belongs to the caller; neither half restarts the feedback chain.
template<unsigned First,unsigned Last,bool Stream,unsigned Rank=4,unsigned Group>
HD35_NOINLINE static void innerRange(const Plan<Group>& p,const Block* planes,Block* output,sr::State& state){
    static_assert(First<Last&&Last<=256);
    static_assert(Rank==2||Rank==4||(HDGROUP_WIDE_INNER&&Rank==16));
    for(unsigned region=First;region!=Last;++region){
        const auto begin=region*p.rows;
        const auto* source=planes+region*p.stride;
        for(std::size_t j=0;j!=p.rows;j+=64){
            if(HDGROUP_PREFETCH&&region+1<Last){
                const auto* next=reinterpret_cast<const char*>(planes+(region+1)*p.stride+j);
                for(unsigned h=0;h!=16;++h)_mm_prefetch(next+64*h,_MM_HINT_T1);
            }
            if constexpr(Rank==16){
                const auto* route=HDGROUP_PACKET_ROUTE?p.packetRoute.data()+(begin+j)/HDGROUP_ROUTE_WIDTH:
                    p.stripedRoute.data()+begin+j;
                c2::hdgroup::innerPacket16<Stream>(source,route,output+begin+j,state,NoMixer{},p.stripedRefresh[(begin+j)/64]);
            }else if constexpr(HDGROUP_INNER&&Rank==4){
                const auto* route=HDGROUP_DECODED_ROUTE?p.decodedRoute.data()+((begin+j)/64)*40:
                    HDGROUP_PACKET_ROUTE?p.packetRoute.data()+(begin+j)/HDGROUP_ROUTE_WIDTH:
                    (HDGROUP_INNER>=2?p.stripedRoute.data():p.route.data())+begin+j;
                if constexpr(HDGROUP_INNER>=4)c2::hdgroup::innerPacket4<Stream>(source,route,output+begin+j,state,NoMixer{},p.stripedRefresh[(begin+j)/64]);
                else c2::hdgroup::innerPacket4<Stream>(source,route,output+begin+j,state,NoMixer{},p.refresh[(begin+j)/64]);
            }else if constexpr(HDGROUP_NATIVE_OUTPUT){
                innerNativeControl<Stream,Rank>(source,p.route.data()+begin+j,output+begin+j,state,p.refresh[(begin+j)/64]);
            }else sr::mdsStriped<true,std::uint16_t,true,Rank,Rank,2,false,Stream>(source,
                    p.route.data()+begin+j,output+begin+j,state,NoMixer{},p.refresh[(begin+j)/64]);
        }
    }
}

// Compare unsigned distances, not potentially wrapping end addresses or byte
// lengths. Dividing the distance is exactly distance < count*sizeof(Block).
static inline bool blockRangesOverlap(std::uintptr_t a,std::size_t aBlocks,
                                     std::uintptr_t b,std::size_t bBlocks){
    if(!aBlocks||!bBlocks)return false;
    return a<=b?(b-a)/sizeof(Block)<aBlocks:(a-b)/sizeof(Block)<bBlocks;
}

template<unsigned Half,unsigned Group>
HD35_NOINLINE static void outerCoset(const Plan<Group>& p,const Block* input,Block* scratch){
    c2::hdgroup::forwardCoset<Group,Half,HDGROUP_COSET_PIPELINE_STREAM,HDGROUP_OUTER_TILE>(
        p.k,p.stride,input,scratch,p.preparedBefore.data(),p.prepared.data());
    if constexpr(HDGROUP_SHIFTED_PAIRS){
        // forwardCoset has fenced any NT stores before we read the wrap value.
        for(unsigned region=128*Half;region!=128*(Half+1);++region){
            auto* plane=scratch+region*p.stride;
            _mm_storeu_si128(reinterpret_cast<__m128i*>(plane+p.rows),
                _mm_loadu_si128(reinterpret_cast<const __m128i*>(plane)));
        }
    }
}

template<bool Stream,unsigned Rank=4,unsigned Group>
HD35_NOINLINE static void forwardPipelined(const Plan<Group>& p,const Block* input,Block* output,Block* scratch){
    sr::State state{};
    outerCoset<0>(p,input,scratch);
    innerRange<0,128,Stream,Rank>(p,scratch,output,state);
    outerCoset<1>(p,input,scratch);
    innerRange<128,256,Stream,Rank>(p,scratch,output,state);
    if constexpr(Stream)_mm_sfence();
}
#endif

template<unsigned Group>
HD35_NOINLINE static void outer(const Plan<Group>& p,const Block* input,Block* scratch){
    c2::hdgroup::forward<Group,true,HDGROUP_OUTER_TILE>(p.k,p.stride,input,scratch,p.preparedBefore.data(),p.prepared.data());
    if constexpr(HDGROUP_SHIFTED_PAIRS){
        // The outer store fence has completed. Fill exactly one record per
        // region so an unaligned YMM load at rows-1 reads the circular pair.
        for(unsigned region=0;region!=256;++region){
            auto* plane=scratch+region*p.stride;
            _mm_storeu_si128(reinterpret_cast<__m128i*>(plane+p.rows),
                _mm_loadu_si128(reinterpret_cast<const __m128i*>(plane)));
        }
    }
}
template<bool Stream,unsigned Rank=4,unsigned Group>
HD35_NOINLINE static void forward(const Plan<Group>& p,const Block* input,Block* output,Block* scratch){
#if HDGROUP_COSET_PIPELINE && !HDGROUP_LONG_RS
    // Half-one output must not clobber input still needed by half two. Long
    // RS is unsupported by this experiment and keeps the complete old path.
    if(!blockRangesOverlap(reinterpret_cast<std::uintptr_t>(input),p.k,
                           reinterpret_cast<std::uintptr_t>(output),p.n)){
        if constexpr(Stream){
            if(reinterpret_cast<std::uintptr_t>(output)&63U)forwardPipelined<false,Rank>(p,input,output,scratch);
            else forwardPipelined<true,Rank>(p,input,output,scratch);
        }else forwardPipelined<false,Rank>(p,input,output,scratch);
        return;
    }
#endif
    outer(p,input,scratch);
    if constexpr(Stream){
        if(reinterpret_cast<std::uintptr_t>(output)&63U)inner<false,Rank>(p,scratch,output);
        else inner<true,Rank>(p,scratch,output);
    }else inner<false,Rank>(p,scratch,output);
}

static std::uint16_t oracleMul(std::uint16_t a,std::uint16_t b){
    const auto& f=c2::p8Oracle().mul;
    const auto h=f[a>>8][b>>8];
    return std::uint16_t((f[a&255][b&255]^f[0x20][h])|((f[a&255][b>>8]^f[a>>8][b&255]^h)<<8));
}
static Block scaled(Block x,std::uint16_t a){
    std::array<std::uint16_t,8> v;std::memcpy(v.data(),&x,16);
    for(auto& c:v)c=oracleMul(c,a);std::memcpy(&x,v.data(),16);return x;
}
static void add(Block& a,Block b){a.lo^=b.lo;a.hi^=b.hi;}
static std::uint8_t byteInverse(std::uint8_t a){
    const auto& f=c2::p8Oracle().mul;
    for(unsigned b=1;b!=256;++b)if(f[a][b]==1)return std::uint8_t(b);
    throw std::runtime_error("zero inverse");
}
static std::uint8_t basis(unsigned degree,unsigned point){
    const auto& f=c2::p8Oracle().mul;std::uint8_t result=1;
    for(unsigned bit=0;bit!=7;++bit)if(degree>>bit&1){
        std::uint8_t value=1,normalizer=1;
        for(unsigned j=0;j!=(1U<<bit);++j){value=f[value][point^j];normalizer=f[normalizer][(1U<<bit)^j];}
        result=f[result][f[value][byteInverse(normalizer)]];
    }
    return result;
}

#if HDGROUP_XOR_MIXER==5
// Independent scalar labels, not generated from the SIMD bit-sum circuit.
static constexpr std::array<std::uint8_t,16> twoSum16Labels{
    0x06,0x07,0x0a,0x0b,0x56,0x57,0x5a,0x5b,
    0xb6,0xb7,0xba,0xbb,0xe6,0xe7,0xea,0xeb};
#endif

template<unsigned Rank=4,unsigned Group>
static void reference(const Plan<Group>& p,const Block* input,Block* output){
    constexpr unsigned length=256/Group,dimension=length/2;
    std::array<std::array<std::uint8_t,dimension>,length> coefficients;
    for(unsigned c=0;c!=length;++c)for(unsigned d=0;d!=dimension;++d)coefficients[c][d]=basis(d,c);
    std::array<std::array<std::uint8_t,Group>,Group> mix;
    const auto& table=c2::p8Oracle().mul;
    for(unsigned i=0;i!=Group;++i)for(unsigned j=0;j!=Group;++j){
        std::uint8_t a=1,b=1;
        for(unsigned l=0;l!=Group;++l)if(l!=j){a=table[a][(Group+i)^l];b=table[b][j^l];}
        mix[i][j]=table[a][byteInverse(b)];
    }
    std::vector<Block> planes(p.n);
    for(std::size_t group=0;group!=p.rows;++group)for(unsigned c=0;c!=length;++c){
        std::array<Block,Group> column{};
        for(unsigned r=0;r!=Group;++r){
            for(unsigned d=0;d!=dimension;++d)add(column[r],scaled(input[(group/4)*512+4*(r*dimension+d)+group%4],coefficients[c][d]));
#if !HDGROUP_LONG_RS
            if constexpr(!HDGROUP_TOEPLITZ&&!(HDGROUP_FIELD32&&Group==2))column[r]=scaled(column[r],p.outerBefore[(c*Group+r)*p.rows+group]);
#endif
        }
#if HDGROUP_LONG_RS
        planes[c*p.rows+group]=scaled(column[0],p.diagonal[c*p.rows+group]);
#else
        if constexpr(HDGROUP_TOEPLITZ){
            const auto& t=p.toeplitz[c*p.rows+group];
            for(unsigned r=0;r!=Group;++r){
                Block y{};for(unsigned a=0;a!=Group;++a)add(y,scaled(column[a],t[Group-1+r-a]));
                planes[(c*Group+r)*p.rows+group]=y;
            }
        }else if constexpr(HDGROUP_FIELD32&&Group==2){
            std::array<std::uint16_t,8> low,high;
            std::memcpy(low.data(),&column[0],16);std::memcpy(high.data(),&column[1],16);
            const auto scalar=std::uint32_t(p.outerBefore[(2*c)*p.rows+group])|(std::uint32_t(p.outerBefore[(2*c+1)*p.rows+group])<<16);
            for(unsigned i=0;i!=8;++i){
                const auto z=c2::medium32::towerMultiplyOracle(std::uint32_t(low[i])|(std::uint32_t(high[i])<<16),scalar);
                low[i]=std::uint16_t(z);high[i]=std::uint16_t(z>>16);
            }
            std::memcpy(&planes[(2*c)*p.rows+group],low.data(),16);
            std::memcpy(&planes[(2*c+1)*p.rows+group],high.data(),16);
        }else{
        for(unsigned r=0;r!=Group;++r){
            Block y{};for(unsigned a=0;a!=Group;++a)add(y,scaled(column[a],mix[r][a]));
            if constexpr(HDGROUP_XOR_MIXER){
#if HDGROUP_XOR_MIXER==5
                y={};for(unsigned a=0;a!=Group;++a)
                    add(y,scaled(column[a],std::uint8_t((a==r)^1^
                        table[twoSum16Labels[r]][twoSum16Labels[a]])));
                y=scaled(y,p.diagonal[(c*Group+r)*p.rows+group]);
#elif HDGROUP_XOR_MIXER>=3
                // Dense scalar matrix oracle, independent of twoSumColumn.
                y={};for(unsigned a=0;a!=Group;++a)
                    add(y,scaled(column[a],std::uint8_t((a==r)^1^
                        table[r+(HDGROUP_XOR_MIXER==4?8:0)][a+(HDGROUP_XOR_MIXER==4?8:0)])));
                y=scaled(y,p.diagonal[(c*Group+r)*p.rows+group]);
#else
                // Scalar H=I+11^T oracle, independent of the shared-sum SIMD
                // circuit. The pre-diagonal was applied to column above.
                y={};for(unsigned a=0;a!=Group;++a)if(a!=r)add(y,column[a]);
                if constexpr(HDGROUP_XOR_MIXER==2)y=scaled(y,p.diagonal[(c*Group+r)*p.rows+group]);
#endif
                planes[(c*Group+r)*p.rows+group]=y;
            }else planes[(c*Group+r)*p.rows+group]=scaled(y,p.diagonal[(c*Group+r)*p.rows+group]);
        }
        }
#endif
    }
    // Independent dense interpolation: Lagrange weights from nodes0..15
    // to nodes16..31, with polynomial-byte arithmetic, not SIMD FFT calls.
    std::array<std::array<std::uint8_t,16>,16> h;
    const auto& f=c2::p8Oracle().mul;
    for(unsigned i=0;i!=16;++i)for(unsigned j=0;j!=16;++j){
        if constexpr(HDGROUP_XOR_REFRESH&&Rank==4)h[i][j]=std::uint8_t(i!=j);
        else{
            std::uint8_t a=1,b=1;
            for(unsigned k=0;k!=16;++k)if(k!=j){a=f[a][(16+i)^k];b=f[b][j^k];}
            h[i][j]=f[a][byteInverse(b)];
        }
    }
    // Scalar basis evaluations are independent of the SIMD FFT schedule.
    // Cache them once so K20 checks do not rebuild the same dense A per step.
    std::array<std::array<std::uint8_t,Rank>,4*Rank> emission;
    for(unsigned i=0;i!=4*Rank;++i)for(unsigned j=0;j!=Rank;++j)emission[i][j]=basis(j,i);
    std::array<Block,16> state{};
    for(std::size_t first=0;first!=p.n;first+=64){
        std::array<Block,64> x;
        std::array<Block,16> feedback{},next{};
        const auto region=first/p.rows,step=first/64;
        for(unsigned i=0;i!=64;++i){
            x[i]=planes[region*p.rows+p.route[first+i]];
            Block y=x[i];
            for(unsigned j=0;j!=Rank;++j){
                const auto a=emission[i%(4*Rank)][j];
                add(y,scaled(state[Rank*(i/(4*Rank))+j],a));
                add(feedback[Rank*(i/(4*Rank))+j],scaled(x[i],a));
            }
            // Explicit scalar permutation, independent of SIMD load/store
            // helpers. Native layout does not change the logical A/C map.
            output[first+(HDGROUP_NATIVE_OUTPUT?4*(i%16)+i/16:i)]=y;
        }
        for(unsigned j=0;j!=16;++j){add(feedback[j],state[j]);feedback[j]=scaled(feedback[j],p.before[step][j]);}
        for(unsigned i=0;i!=16;++i){
            for(unsigned j=0;j!=16;++j)add(next[i],scaled(feedback[j],h[i][j]));
            next[i]=scaled(next[i],p.after[step][i]);
        }
        state=next;
    }
}
static void produce(Block* p,std::size_t n){
    std::uint64_t x=913,y=113;
    for(std::size_t i=0;i!=n;++i){x^=x<<13;x^=x>>7;x^=x<<17;y+=0x9e3779b97f4a7c15ULL;_mm_stream_si128(reinterpret_cast<__m128i*>(p+i),_mm_set_epi64x(y,x));}_mm_sfence();
}
static std::uint64_t hash(const Block* p,std::size_t n){std::uint64_t h=0;for(std::size_t i=0;i!=n;++i){h=(h^p[i].lo)*0x100000001b3ULL;h=(h^p[i].hi)*0x100000001b3ULL;}return h;}

#if HDGROUP_TOEPLITZ_COMPACT7
static_assert(c2::hdgroup::toeplitzProducts<4> == 9);
static_assert(c2::hdgroup::toeplitzStoredCoefficients<4> == 7);
static_assert(c2::hdgroup::beforeRegions<4> == 448);
static void checkCompactToeplitzKeys(const Plan<4>& p){
    require(p.preparedBefore.size()==p.rows*448,"Compact7 prepared key size");
    require(p.prepared.empty()&&p.outerBefore.empty(),"Compact7 has no diagonal keys");
    for(std::size_t g=0;g!=p.rows;g+=4)for(unsigned c=0;c!=64;++c)for(unsigned lane=0;lane!=4;++lane){
        const auto& t=p.toeplitz[c*p.rows+g+lane];
        const auto* keys=p.preparedBefore.data()+g*448+28*c+lane;
        const auto keyXor=std::uint8_t(keys[4]^keys[8]);
        const std::array<std::uint8_t,9> reconstructed{
            keys[0],keys[4],keys[8],keys[12],std::uint8_t(keys[12]^keyXor),keys[16],
            keys[20],keys[24],std::uint8_t(keys[20]^keyXor)};
        const std::array<std::uint8_t,9> expected{
            std::uint8_t(t[3]),std::uint8_t(t[3]^t[2]),std::uint8_t(t[3]^t[4]),
            std::uint8_t(t[5]^t[3]),std::uint8_t(t[5]^t[3]^t[4]^t[2]),std::uint8_t(t[5]^t[3]^t[6]^t[4]),
            std::uint8_t(t[1]^t[3]),std::uint8_t(t[1]^t[3]^t[0]^t[2]),std::uint8_t(t[1]^t[3]^t[2]^t[4])};
        require(reconstructed==expected,"Compact7 sampled coefficient identities");
    }
}
static void checkCompactToeplitzColumn(){
    const auto& f=c2::p8Oracle().mul;
    alignas(64) std::uint8_t input[4][64],expected[4][64],actual[4][64];
    alignas(64) __m512i values[4];
    std::array<std::array<std::uint8_t,7>,4> matrices{};
    std::array<c2::hdgroup::Coefficient,28> keys{};
    std::mt19937_64 rng(0x67e3c18bd9205fa4ULL);
    const auto prepare=[&]{
        for(unsigned lane=0;lane!=4;++lane){
            const auto& t=matrices[lane];
            const std::array<std::uint8_t,7> compact{
                t[3],std::uint8_t(t[3]^t[2]),std::uint8_t(t[3]^t[4]),
                std::uint8_t(t[5]^t[3]),std::uint8_t(t[5]^t[3]^t[6]^t[4]),
                std::uint8_t(t[1]^t[3]),std::uint8_t(t[1]^t[3]^t[0]^t[2])};
            for(unsigned i=0;i!=7;++i)keys[4*i+lane]=compact[i];
        }
    };
    const auto checkColumn=[&]{
        std::memset(expected,0,sizeof(expected));
        for(unsigned i=0;i!=4;++i)for(unsigned j=0;j!=4;++j)for(unsigned b=0;b!=64;++b)
            expected[i][b]^=f[matrices[b/16][3+i-j]][input[j][b]];
        std::memcpy(values,input,sizeof(input));c2::hdgroup::toeplitzColumn<4>(values,keys.data());
        std::memcpy(actual,values,sizeof(actual));
        require(std::memcmp(actual,expected,sizeof(actual))==0,"Compact7 dense scalar matrix mismatch");
    };
    // The map is bilinear in the seven matrix bytes and four input bytes.
    // Exhaust both GF2 bases, with independent coefficient bases per lane.
    for(unsigned coefficient=0;coefficient!=56;++coefficient){
        for(unsigned lane=0;lane!=4;++lane){
            matrices[lane].fill(0);const auto bit=(coefficient+11*lane)%56;
            matrices[lane][bit/8]=std::uint8_t(1U<<(bit%8));
        }
        prepare();
        for(unsigned direction=0;direction!=32;++direction){
            std::memset(input,0,sizeof(input));
            for(unsigned b=0;b!=64;++b)input[direction/8][b]=std::uint8_t(1U<<((direction+b)%8));
            checkColumn();
        }
    }
    // Every individual SIMD bit under two independent dense key samples.
    // Then random full columns also exercise cancellations and mixed bytes.
    for(unsigned sample=0;sample!=2;++sample){
        for(auto& t:matrices)for(auto& byte:t)byte=std::uint8_t(rng());
        prepare();
        for(unsigned direction=0;direction!=2048;++direction){
            std::memset(input,0,sizeof(input));input[direction/512][direction/8%64]=std::uint8_t(1U<<(direction%8));
            checkColumn();
        }
    }
    for(unsigned sample=0;sample!=16;++sample){
        for(auto& t:matrices)for(auto& byte:t)byte=std::uint8_t(rng());
        for(auto& row:input)for(auto& byte:row)byte=std::uint8_t(rng());
        prepare();checkColumn();
    }
    std::cout<<"check,toeplitz_compact7,coefficient_input_basis_pairs,1792,simd_basis_directions,4096,random_columns,16,pass\n"<<std::flush;
}
#endif

#if HDGROUP_COSET_PIPELINE && !HDGROUP_LONG_RS
static void checkCosetRanges(){
    require(!blockRangesOverlap(64,0,64,4),"coset empty input range");
    require(!blockRangesOverlap(64,4,64,0),"coset empty output range");
    require(blockRangesOverlap(64,4,64,8),"coset equal ranges");
    require(!blockRangesOverlap(64,4,128,8),"coset adjacent forward ranges");
    require(!blockRangesOverlap(128,8,64,4),"coset adjacent reverse ranges");
    require(blockRangesOverlap(64,4,127,8),"coset partial forward ranges");
    require(blockRangesOverlap(127,8,64,4),"coset partial reverse ranges");
    const auto top=~std::uintptr_t(0);
    require(blockRangesOverlap(top-31,2,top-15,1),"coset overlap near address limit");
    require(!blockRangesOverlap(top-31,1,top-15,1),"coset adjacency near address limit");
    require(!blockRangesOverlap(top-15,1,0,1),"coset distant ranges do not wrap");
    std::cout<<"check,coset-range-overlap,unsigned_distance,pass\n"<<std::flush;
}
#endif

#if HDGROUP_LONG_RS
static void checkLongOuter(){
    // Exhaust the 128 polynomial-basis inputs in all four interleaved groups.
    // The dense scalar polynomial oracle shares no SIMD butterfly helper.
    c2::BenchPayload input(512),actual(1024);
    std::array<c2::hdgroup::Coefficient,1024> diagonal;
    for(unsigned c=0;c!=256;++c)for(unsigned lane=0;lane!=4;++lane)
        diagonal[4*c+lane]=std::uint8_t(1+(37*c+17*lane)%255);
    const std::array<Block,4> seed{{{0x0123456789abcdefULL,0xfedcba9876543210ULL},
        {0x1020304050607080ULL,0x9080706050403020ULL},
        {0xa5a55a5aff00ff00ULL,0x55aaff0000ffaa55ULL},
        {0xf1e2d3c4b5a69788ULL,0x7766554433221100ULL}}};
    for(unsigned d=0;d!=128;++d){
        std::fill_n(input.p,512,Block{});
        for(unsigned lane=0;lane!=4;++lane)input.p[4*d+lane]=seed[lane];
        c2::hdgroup::forward<1,false>(512,4,input.p,actual.p,nullptr,diagonal.data());
        for(unsigned c=0;c!=256;++c)for(unsigned lane=0;lane!=4;++lane)
            require(actual.p[4*c+lane]==scaled(scaled(seed[lane],basis(d,c)),diagonal[4*c+lane]),"long RS polynomial basis oracle");
    }
    std::cout<<"check,long-rs256,polynomial_basis,128,interleaved_groups,4,outputs,256,pass\n"<<std::flush;
}
#endif

#if HDGROUP_XOR_MIXER==5
static void checkTwoSum16Column(){
    const auto& f=c2::p8Oracle().mul;
    alignas(64) std::uint8_t input[16][64],expected[16][64],actual[16][64];
    alignas(64) __m512i values[16];
    std::mt19937_64 rng(0x67e3c18bd9205fa4ULL);
    for(unsigned test=0;test!=8192+16;++test){
        std::memset(input,0,sizeof(input));std::memset(expected,0,sizeof(expected));
        if(test<8192)input[test/512][test/8%64]=std::uint8_t(1U<<(test%8));
        else for(auto& row:input)for(auto& byte:row)byte=std::uint8_t(rng());
        for(unsigned i=0;i!=16;++i)for(unsigned j=0;j!=16;++j){
            const auto h=std::uint8_t((i==j)^1^f[twoSum16Labels[i]][twoSum16Labels[j]]);
            for(unsigned b=0;b!=64;++b)expected[i][b]^=f[h][input[j][b]];
        }
        std::memcpy(values,input,sizeof(input));c2::hdgroup::twoSum16Column(values);
        std::memcpy(actual,values,sizeof(actual));
        require(std::memcmp(actual,expected,sizeof(actual))==0,"GF16-coset dense matrix mismatch");
        c2::hdgroup::twoSum16Column(values);std::memcpy(actual,values,sizeof(actual));
        require(std::memcmp(actual,input,sizeof(actual))==0,"GF16-coset involution mismatch");
    }
    std::cout<<"two_sum_gf16_basis_directions,8192,random_columns,16,ok\n";
}
#elif HDGROUP_XOR_MIXER>=3
static void checkTwoSumColumn(){
    const auto& f=c2::p8Oracle().mul;
    alignas(64) std::uint8_t input[8][64],expected[8][64],actual[8][64];
    alignas(64) __m512i values[8];
    std::mt19937_64 rng(0x67e3c18bd9205fa4ULL);
    // Every bit of every SIMD byte, then independent random columns.
    for(unsigned test=0;test!=4096+16;++test){
        std::memset(input,0,sizeof(input));std::memset(expected,0,sizeof(expected));
        if(test<4096)input[test/512][test/8%64]=std::uint8_t(1U<<(test%8));
        else for(auto& row:input)for(auto& byte:row)byte=std::uint8_t(rng());
        for(unsigned i=0;i!=8;++i)for(unsigned j=0;j!=8;++j){
            const auto h=std::uint8_t((i==j)^1^
                f[i+(HDGROUP_XOR_MIXER==4?8:0)][j+(HDGROUP_XOR_MIXER==4?8:0)]);
            for(unsigned b=0;b!=64;++b)expected[i][b]^=f[h][input[j][b]];
        }
        std::memcpy(values,input,sizeof(input));c2::hdgroup::twoSumColumn(values);
        std::memcpy(actual,values,sizeof(actual));
        require(std::memcmp(actual,expected,sizeof(actual))==0,"two-sum dense matrix mismatch");
        c2::hdgroup::twoSumColumn(values);std::memcpy(actual,values,sizeof(actual));
        require(std::memcmp(actual,input,sizeof(actual))==0,"two-sum involution mismatch");
    }
    std::cout<<"two_sum_basis_directions,4096,random_columns,16,ok\n";
}
#endif

static void checkPairLoads(){
#if HDGROUP_PAIR_SWAPS
    // Check all four independent swap masks without any FFT/state algebra.
    // Odd bases include the wrapping pair; distinct records expose lane loss.
    alignas(64) std::array<Block,65> input{};
    alignas(64) std::array<Block,4> actual{};
    for(unsigned i=0;i!=64;++i)input[i]={0x12340000U+i,0x98760000U+3*i};
    input[64]=input[0];
    for(unsigned parity=0;parity!=1+bool(HDGROUP_SHIFTED_PAIRS);++parity){
        for(unsigned wrapPart=0;wrapPart!=2;++wrapPart){
            const std::array<unsigned,2> base{wrapPart?2U+parity:62U+parity,wrapPart?62U+parity:4U+parity};
            for(unsigned swaps=0;swaps!=4;++swaps){
                const std::array<std::uint16_t,32> route{
                    std::uint16_t(base[0]|((swaps&1U)<<15)),
                    std::uint16_t(base[1]|((swaps&2U)<<14))};
                std::array<std::uint16_t,40> decoded{};
                decodePacketRoute(route.data(),decoded.data());
                _mm512_store_si512(actual.data(),c2::hdgroup::innerInput4<0>(input.data(),HDGROUP_DECODED_ROUTE?decoded.data():route.data()));
                for(unsigned part=0;part!=2;++part)for(unsigned lane=0;lane!=2;++lane)
                    require(actual[2*part+lane]==input[(base[part]+(lane^((swaps>>part)&1U)))%64],"independent pair swap load");
            }
        }
    }
    std::cout<<"check,pair-swaps,all_masks,4,wrapping,"<<bool(HDGROUP_SHIFTED_PAIRS)<<",pass\n"<<std::flush;
#endif
}

template<unsigned Rank=4>
static void checkPackets(){
#if HDGROUP_INNER || HDGROUP_WIDE_INNER
    static_assert(Rank==4||(HDGROUP_WIDE_INNER&&Rank==16));
    // Exhaust every GF2 basis direction of the 64-record input and 16-record
    // incoming state. Both kernels are binary-linear for each fixed setup.
    alignas(64) std::array<Block,64+(HDGROUP_SHIFTED_PAIRS?1:0)> input{};
    alignas(64) std::array<Block,64> actual{},expected{};
    alignas(64) std::array<Block,16> initial{};
    std::array<std::uint16_t,64> route{},striped{};
    std::array<std::uint16_t,64/HDGROUP_ROUTE_WIDTH> packets{};
    std::array<std::uint16_t,40> decoded{};
    std::array<std::array<std::uint8_t,4>,16> scalarEmission{};
    if constexpr(HDGROUP_XOR_REFRESH)for(unsigned i=0;i!=16;++i)for(unsigned j=0;j!=4;++j)
        scalarEmission[i][j]=basis(j,i);
    std::mt19937_64 rng(391);
    for(unsigned sample=0;sample!=3;++sample){
        if constexpr(HDGROUP_PACKET_ROUTE){
            constexpr unsigned width=HDGROUP_ROUTE_WIDTH;
            // Keep one even route and exercise wrapping odd pairs in both
            // shuffled routes when the shifted-pair experiment is enabled.
            const unsigned parity=HDGROUP_SHIFTED_PAIRS&&sample!=0;
            for(unsigned j=0;j!=64/width;++j)packets[j]=std::uint16_t(width*j+parity);
            if(sample)std::shuffle(packets.begin(),packets.end(),rng);
#if HDGROUP_PAIR_SWAPS
            for(unsigned j=0;j!=64/width;++j){
                // Sample0 cycles through all four ZMM swap masks. Sample1
                // swaps every pair, including the odd wrapping pair;
                // sample2 exercises that same wrapping pair without a swap.
                const unsigned swap=sample==0?((j/2)>>(j%2))&1U:
                    sample==1?1U:(packets[j]==63?0U:unsigned(rng()&1U));
                packets[j]|=std::uint16_t(swap<<15);
            }
#endif
            for(unsigned j=0;j!=16;++j)for(unsigned lane=0;lane!=4;++lane){
#if HDGROUP_PAIR_SWAPS
                const unsigned entry=packets[(4/width)*j+lane/width];
                const unsigned record=(entry&0x7fffU)+((lane%width)^(entry>>15));
#else
                const unsigned record=packets[(4/width)*j+lane/width]+lane%width;
#endif
                route[16*lane+j]=std::uint16_t(HDGROUP_SHIFTED_PAIRS?record%64:record);
            }
        }else{
            std::iota(route.begin(),route.end(),std::uint16_t(0));
            if(sample)std::shuffle(route.begin(),route.end(),rng);
        }
        if constexpr(HDGROUP_DECODED_ROUTE)decodePacketRoute(packets.data(),decoded.data());
        for(unsigned j=0;j!=16;++j)for(unsigned lane=0;lane!=4;++lane)striped[4*j+lane]=route[16*lane+j];
        std::array<std::uint16_t,16> a,b,ta,tb;
        for(unsigned j=0;j!=16;++j){a[j]=std::uint16_t(1+rng()%255);b[j]=std::uint16_t(1+rng()%255);}
        for(unsigned j=0;j!=16;++j){ta[j]=a[4*(j%4)+j/4];tb[j]=b[4*(j%4)+j/4];}
        sr::ByteDiagonalRefresh refresh,stripedRefresh;refresh.prepare(a.data(),b.data());stripedRefresh.prepare(ta.data(),tb.data());
        for(unsigned direction=0;direction!=80*128;++direction){
            input.fill({});initial.fill({});
            auto* bytes=reinterpret_cast<std::uint8_t*>(direction<64*128?input.data():initial.data());
            const auto bit=direction<64*128?direction:direction-64*128;
            bytes[bit/8]=std::uint8_t(1U<<(bit%8));
            if constexpr(HDGROUP_SHIFTED_PAIRS)input[64]=input[0];
            sr::State oldState,newState;
            for(unsigned j=0;j!=4;++j)oldState.values[j]=newState.values[j]=_mm512_load_si512(initial.data()+4*j);
            sr::mdsStriped<true,std::uint16_t,true,Rank,Rank,2,false,false>(input.data(),route.data(),expected.data(),oldState,NoMixer{},refresh);
            if constexpr(HDGROUP_XOR_REFRESH&&Rank==4){
                // The emission is unchanged. Independently compute the new
                // recurrence with scalar A^T and a dense off-diagonal-one H.
                std::array<Block,16> pre{},next{};
                for(unsigned j=0;j!=16;++j){
                    pre[j]=initial[j];
                    for(unsigned i=0;i!=16;++i)
                        add(pre[j],scaled(input[route[16*(j/4)+i]],scalarEmission[i][j%4]));
                    pre[j]=scaled(pre[j],a[j]);
                }
                for(unsigned i=0;i!=16;++i){
                    for(unsigned j=0;j!=16;++j)if(i!=j)add(next[i],pre[j]);
                    next[i]=scaled(next[i],b[i]);
                }
                for(unsigned j=0;j!=4;++j)oldState.values[j]=_mm512_loadu_si512(next.data()+4*j);
            }
            if constexpr(Rank==16||HDGROUP_INNER>=4)sr::transposeRecords4(newState.values[0],newState.values[1],newState.values[2],newState.values[3]);
            const auto* r=HDGROUP_DECODED_ROUTE?decoded.data():HDGROUP_PACKET_ROUTE?packets.data():((Rank==16||HDGROUP_INNER>=2)?striped.data():route.data());
            if constexpr(Rank==16)c2::hdgroup::innerPacket16<false>(input.data(),r,actual.data(),newState,NoMixer{},stripedRefresh);
            else c2::hdgroup::innerPacket4<false>(input.data(),r,actual.data(),newState,NoMixer{},HDGROUP_INNER>=4?stripedRefresh:refresh);
            if constexpr(Rank==16||HDGROUP_INNER>=4)sr::transposeRecords4(newState.values[0],newState.values[1],newState.values[2],newState.values[3]);
            if constexpr(HDGROUP_NATIVE_OUTPUT){
                for(unsigned i=0;i!=64;++i)require(actual[4*(i%16)+i/16]==expected[i],"packet basis native output");
            }else require(actual==expected,"packet basis output");
            for(unsigned j=0;j!=4;++j)require(_mm512_cmpeq_epi8_mask(oldState.values[j],newState.values[j])==~__mmask64(0),"packet basis next state");
        }
    }
    std::cout<<"check,packet-basis,variant,"<<HDGROUP_INNER<<",rank,"<<Rank
        <<",xor_refresh,"<<(HDGROUP_XOR_REFRESH&&Rank==4)
        <<",route_mode,"<<(HDGROUP_PACKET_ROUTE?(HDGROUP_ROUTE_WIDTH==2?"packet2":"packet4"):"record")
        <<",packet_width,"<<(HDGROUP_PACKET_ROUTE?HDGROUP_ROUTE_WIDTH:1)
        <<",shifted_pairs,"<<HDGROUP_SHIFTED_PAIRS<<",odd_pair_routes,"<<(HDGROUP_SHIFTED_PAIRS?2:0)
#if HDGROUP_PAIR_SWAPS
        <<",pair_swaps,"<<HDGROUP_PAIR_SWAPS
#endif
        <<",output_layout,"<<(HDGROUP_NATIVE_OUTPUT?"native_striped4":"natural")
        <<",directions,30720,pass\n"<<std::flush;
#endif
}

template<unsigned Group>
static void checkRoutes(const Plan<Group>& p){
    if constexpr(HDGROUP_SHIFTED_PAIRS==2)require(p.oddPairRegions()==128,"balanced alternating pair parity");
    std::vector<bool> seen(p.rows);
    for(unsigned c=0;c!=256;++c){
        std::fill(seen.begin(),seen.end(),false);
        const auto* route=p.route.data()+c*p.rows;
        for(std::size_t i=0;i!=p.rows;++i){
            require(route[i]<p.rows,"route range");
            require(!seen[route[i]],"route permutation");seen[route[i]]=true;
        }
        if constexpr(HDGROUP_PACKET_ROUTE){
            constexpr unsigned width=HDGROUP_ROUTE_WIDTH;
            require(p.pairParity[c]<2,"pair parity range");
            if constexpr(!HDGROUP_SHIFTED_PAIRS)require(!p.pairParity[c],"unshifted packet parity");
            if constexpr(HDGROUP_SHIFTED_PAIRS==2)require(p.pairParity[c]==(c&1),"alternating pair parity");
            const auto* packets=p.packetRoute.data()+c*p.rows/width;
            for(std::size_t first=0;first!=p.rows;first+=64)for(unsigned j=0;j!=16;++j)for(unsigned part=0;part!=4/width;++part){
#if HDGROUP_PAIR_SWAPS
                const unsigned entry=packets[first/width+(4/width)*j+part];
                const auto base=entry&0x7fffU;
#else
                const auto base=packets[first/width+(4/width)*j+part];
#endif
                require(base<p.rows,"packet base range");require(base%width==p.pairParity[c],"packet base parity");
                for(unsigned lane=0;lane!=width;++lane){
#if HDGROUP_PAIR_SWAPS
                    const auto record=std::size_t(base)+(lane^(entry>>15));
#else
                    const auto record=std::size_t(base)+lane;
#endif
                    require(route[first+16*(part*width+lane)+j]==(HDGROUP_SHIFTED_PAIRS?record%p.rows:record),"packet lane identity");
                }
            }
        }else if constexpr(HDGROUP_INNER>=2||HDGROUP_WIDE_INNER){
            const auto* striped=p.stripedRoute.data()+c*p.rows;
            for(std::size_t first=0;first!=p.rows;first+=64)for(unsigned j=0;j!=16;++j)for(unsigned lane=0;lane!=4;++lane)
                require(striped[first+4*j+lane]==route[first+16*lane+j],"striped route layout");
        }
    }
}

template<unsigned Group=4,unsigned Rank=4>
static void check(std::size_t k){
    Plan<Group> p(k,113);c2::BenchPayload input(k+4),actual(p.n+8),expected(p.n),scratch(p.scratch()+8);
#if HDGROUP_TOEPLITZ_COMPACT7
    checkCompactToeplitzKeys(p);
#endif
#if HDGROUP_LONG_RS
    require(p.preparedBefore.empty()&&p.outerBefore.empty()&&p.toeplitz.empty(),"long RS has no column-mixer keys");
    require(p.diagonal.size()==p.n&&p.prepared.size()==p.n,"long RS diagonal dimensions");
    for(std::size_t g=0;g!=p.rows;g+=4)for(unsigned c=0;c!=256;++c)for(unsigned lane=0;lane!=4;++lane){
        const auto scalar=p.diagonal[c*p.rows+g+lane];
        require(scalar>=1&&scalar<=255,"long RS nonzero GF256 diagonal");
        require(p.prepared[g*256+4*c+lane]==scalar,"long RS packed diagonal identity");
    }
#endif
    checkRoutes(p);
    produce(input.p,k);reference<Rank>(p,input.p,expected.p);
#if HDGROUP_COSET_PIPELINE
    // Independent schedule control: all outer work, then the original full
    // inner. Both this control and every pipelined case use the dense oracle.
    outer(p,input.p,scratch.p);inner<false,Rank>(p,scratch.p,actual.p);
    require(std::equal(expected.p,expected.p+p.n,actual.p),"unfused full-path oracle");
#endif
    const Block guard{0xa1369,0x73495};
    for(unsigned offset=0;offset!=4;++offset){
        std::fill_n(actual.p,p.n+8,guard);
        std::fill_n(scratch.p,p.scratch()+8,guard);
        forward<true,Rank>(p,input.p,actual.p+offset,scratch.p+offset);
        require(std::equal(expected.p,expected.p+p.n,actual.p+offset),"full independent dense oracle");
        for(unsigned j=0;j!=offset;++j)require(actual.p[j]==guard,"prefix guard");
        for(auto j=p.n+offset;j!=p.n+8;++j)require(actual.p[j]==guard,"suffix guard");
        for(unsigned j=0;j!=offset;++j)require(scratch.p[j]==guard,"scratch prefix guard");
        for(auto j=p.scratch()+offset;j!=p.scratch()+8;++j)require(scratch.p[j]==guard,"scratch suffix guard");
        for(unsigned region=0;region!=256;++region){
            const auto* plane=scratch.p+offset+region*p.stride;
            if constexpr(HDGROUP_SHIFTED_PAIRS)require(plane[p.rows]==plane[0],"wrapped pair padding");
            for(auto j=p.rows+(HDGROUP_SHIFTED_PAIRS?1:0);j!=p.stride;++j)
                require(plane[j]==Block{},"unused region padding");
        }
    }
    forward<false,Rank>(p,input.p,actual.p,scratch.p);
    require(std::equal(expected.p,expected.p+p.n,actual.p),"cached output oracle");
    std::copy_n(input.p,k,actual.p);forward<true,Rank>(p,actual.p,actual.p,scratch.p);
    require(std::equal(expected.p,expected.p+p.n,actual.p),"full alias");
#if HDGROUP_COSET_PIPELINE
    // Both overlap directions must fall back before first-half output can
    // destroy source records required by the second outer coset.
    std::copy_n(input.p,k,actual.p+1);forward<true,Rank>(p,actual.p+1,actual.p,scratch.p);
    require(std::equal(expected.p,expected.p+p.n,actual.p),"partial alias output before input");
    std::copy_n(input.p,k,actual.p);forward<true,Rank>(p,actual.p,actual.p+1,scratch.p);
    require(std::equal(expected.p,expected.p+p.n,actual.p+1),"partial alias output after input");
#endif
    // Four interleaved outer rows must not impose input alignment.
    std::move_backward(input.p,input.p+k,input.p+k+1);
    forward<true,Rank>(p,input.p+1,actual.p,scratch.p);
    require(std::equal(expected.p,expected.p+p.n,actual.p),"unaligned input oracle");
    std::cout<<"check,grouped,rows,"<<Group<<",rank,"<<Rank<<",K,"<<k
        <<",xor_refresh,"<<(HDGROUP_XOR_REFRESH&&Rank==4)
#if HDGROUP_COSET_PIPELINE
        <<",coset_pipeline,"<<!HDGROUP_LONG_RS<<",coset_pipeline_stream,"<<HDGROUP_COSET_PIPELINE_STREAM
#endif
#if HDGROUP_LONG_RS
        <<",outer_construction,rs256-diagonal-gf256,long_rs,1"
#endif
#if HDGROUP_TOEPLITZ_COMPACT7
        <<",toeplitz_compact7,1,toeplitz_keys_per_column,7,toeplitz_gfni_per_column,9"
#endif
        <<",route_mode,"<<(HDGROUP_PACKET_ROUTE?(HDGROUP_ROUTE_WIDTH==2?"packet2":"packet4"):"record")
        <<",packet_width,"<<(HDGROUP_PACKET_ROUTE?HDGROUP_ROUTE_WIDTH:1)
        <<",shifted_pairs,"<<HDGROUP_SHIFTED_PAIRS<<",odd_pair_regions,"<<p.oddPairRegions()
#if HDGROUP_PAIR_SWAPS
        <<",pair_swaps,"<<HDGROUP_PAIR_SWAPS<<",swapped_pairs,"<<p.swappedPairs()
#endif
        <<",output_layout,"<<(HDGROUP_NATIVE_OUTPUT?"native_striped4":"natural")
        <<",pass,checksum,"<<hash(expected.p,p.n)<<'\n'<<std::flush;
}
template<bool Stream,unsigned Group=4,unsigned Rank=4>
static void bench(unsigned calls,unsigned logK,const std::string& mode){
    Plan<Group> p(std::size_t(1)<<logK,113);c2::BenchPayload input(p.k),out(p.n),scratch(p.scratch());
    produce(input.p,p.k);outer(p,input.p,scratch.p);
    auto run=[&]{
        if(mode=="outer")outer(p,input.p,scratch.p);
        else if(mode=="inner")inner<Stream,Rank>(p,scratch.p,out.p);
        else forward<Stream,Rank>(p,input.p,out.p,scratch.p);
    };
    for(unsigned j=0;j!=3;++j)run();
    std::vector<double> times(calls);std::uint64_t consume=0;
    const auto* result=mode=="outer"?scratch.p:out.p;
    for(unsigned j=0;j!=calls;++j){
        produce(input.p,p.k);const auto start=Clock::now();run();
        times[j]=std::chrono::duration<double,std::micro>(Clock::now()-start).count();
        std::atomic_signal_fence(std::memory_order_seq_cst);consume^=result[(j*104729U)%p.n].lo;
    }
    std::sort(times.begin(),times.end());
    std::cout<<"profile,grouped,rows,"<<Group<<",rank,"<<Rank<<",K,"<<p.k<<",N,"<<p.n<<",stream,"<<Stream<<",calls,"<<calls
#if HDGROUP_LONG_RS
        <<",outer_construction,rs256-diagonal-gf256,long_rs,1"
        <<",field32,0,normalized,0"
#else
        <<",field32,"<<(Group==2?HDGROUP_FIELD32:0)<<",normalized,"<<(!HDGROUP_TOEPLITZ&&HDGROUP_NORMALIZE)
#endif
        <<",toeplitz,"<<HDGROUP_TOEPLITZ<<",setup_rejections,"<<p.toeplitzRejected
#if HDGROUP_TOEPLITZ_COMPACT7
        <<",toeplitz_compact7,1,toeplitz_keys_per_column,7,toeplitz_gfni_per_column,9"
        <<",toeplitz_key_expansions_per_column,7,toeplitz_key_xors_per_column,3"
        <<",outer_hot_key_bytes,"<<p.preparedBefore.size()*sizeof(c2::hdgroup::Coefficient)
#endif
        <<",xor_mixer,"<<HDGROUP_XOR_MIXER
        <<",byte_keys,"<<HDGROUP_BYTE_KEYS
        <<",decoded_route,"<<HDGROUP_DECODED_ROUTE
        <<",route_address_shift,"<<HDGROUP_ROUTE_ADDRESS_SHIFT
        <<",inner_variant,"<<HDGROUP_INNER<<",inner_byte_keys,"<<HD_SUMRANK_BYTE_KEYS
        <<",xor_refresh,"<<(HDGROUP_XOR_REFRESH&&Rank==4)
#if HDGROUP_COSET_PIPELINE
        <<",coset_pipeline,"<<(!HDGROUP_LONG_RS&&mode!="outer"&&mode!="inner")
        <<",coset_pipeline_stream,"<<HDGROUP_COSET_PIPELINE_STREAM
#endif
        <<",wide_inner,"<<(HDGROUP_WIDE_INNER&&Rank==16)
        <<",route_mode,"<<(HDGROUP_PACKET_ROUTE?(HDGROUP_ROUTE_WIDTH==2?"packet2":"packet4"):"record")
        <<",packet_width,"<<(HDGROUP_PACKET_ROUTE?HDGROUP_ROUTE_WIDTH:1)
        <<",shifted_pairs,"<<HDGROUP_SHIFTED_PAIRS<<",odd_pair_regions,"<<p.oddPairRegions()
#if HDGROUP_PAIR_SWAPS
        <<",pair_swaps,"<<HDGROUP_PAIR_SWAPS<<",swapped_pairs,"<<p.swappedPairs()
#endif
        <<",output_layout,"<<(HDGROUP_NATIVE_OUTPUT?"native_striped4":"natural")
        <<",fused_refresh,"<<HDGROUP_FUSE_REFRESH<<",prefetch,"<<HDGROUP_PREFETCH
        <<",outer_tile,"<<HDGROUP_OUTER_TILE
        <<",field_order,"<<(HDGROUP_SUBFIELD?256:65536)
        <<",fused_last,"<<(HDGROUP_FUSE_LAST&&Group<=8)<<",checksum,"<<hash(result,p.n)<<",consume,"<<consume<<"\nname,median_us,p10_us,p90_us\n"
        <<std::fixed<<std::setprecision(3)<<mode<<','<<times[calls/2]<<','<<times[calls/10]<<','<<times[9*calls/10]<<'\n';
}
template<unsigned Rank>
static void dispatch(unsigned group,bool nt,unsigned calls,unsigned logK,const std::string& mode){
    switch(group){
#define GROUP_CASE(G) case G: if(nt)bench<true,G,Rank>(calls,logK,mode);else bench<false,G,Rank>(calls,logK,mode);break
#if HDGROUP_LONG_RS
    GROUP_CASE(1);
#elif HDGROUP_TOEPLITZ_COMPACT7
    GROUP_CASE(4);
#elif HDGROUP_XOR_MIXER==5
    GROUP_CASE(16);
#elif HDGROUP_XOR_MIXER>=3
    GROUP_CASE(8);
#else
    GROUP_CASE(2);GROUP_CASE(4);
#if !HDGROUP_TOEPLITZ
    GROUP_CASE(8);GROUP_CASE(16);
#endif
#endif
#undef GROUP_CASE
    default:throw std::invalid_argument("invalid group width");
    }
}
int main(int argc,char** argv){try{
    if(!c2::gfni512Available())return 77;
    if(argc==2&&std::string(argv[1])=="--check"){
#if HDGROUP_COSET_PIPELINE && !HDGROUP_LONG_RS
        checkCosetRanges();
#endif
#if HDGROUP_LONG_RS
        checkLongOuter();
#endif
#if HDGROUP_TOEPLITZ_COMPACT7
        checkCompactToeplitzColumn();
#endif
#if HDGROUP_XOR_MIXER==5
        checkTwoSum16Column();
#elif HDGROUP_XOR_MIXER>=3
        checkTwoSumColumn();
#endif
        checkPairLoads();
#if HDGROUP_INNER
        checkPackets<4>();
#endif
#if HDGROUP_WIDE_INNER
        checkPackets<16>();
#endif
#if HDGROUP_LONG_RS
        for(auto k:{8192U,24576U}){check<1,2>(k);check<1,4>(k);}
#elif HDGROUP_XOR_MIXER==5
        for(auto k:{8192U,24576U}){check<16,2>(k);check<16,4>(k);}
#elif HDGROUP_XOR_MIXER>=3
        for(auto k:{8192U,24576U}){check<8,2>(k);check<8,4>(k);}
#elif HDGROUP_TOEPLITZ_COMPACT7
        for(auto k:{8192U,24576U}){check<4,2>(k);check<4,4>(k);}
#elif HDGROUP_TOEPLITZ
        for(auto k:{8192U,24576U}){check<2,2>(k);check<4,2>(k);check<2,4>(k);check<4,4>(k);}
#else
        for(auto k:{8192U,24576U}){check<2,2>(k);check<4,2>(k);check<8,2>(k);check<16,2>(k);}
        check<2,4>(8192);check<4,4>(8192);check<8,4>(8192);check<16,4>(8192);
#endif
#if HDGROUP_WIDE_INNER
        for(auto k:{8192U,24576U}){
#if HDGROUP_LONG_RS
            check<1,16>(k);
#elif HDGROUP_XOR_MIXER==5
            check<16,16>(k);
#elif HDGROUP_XOR_MIXER>=3
            check<8,16>(k);
#elif HDGROUP_TOEPLITZ_COMPACT7
            check<4,16>(k);
#else
            check<2,16>(k);check<4,16>(k);
#if !HDGROUP_TOEPLITZ
            check<8,16>(k);check<16,16>(k);
#endif
#endif
        }
#endif
        return 0;
    }
    if(argc==2&&std::string(argv[1])=="--large-check"){
#if HDGROUP_LONG_RS
        check<1,4>(1U<<20);
#elif HDGROUP_XOR_MIXER==5
        check<16,4>(1U<<20);
#elif HDGROUP_XOR_MIXER>=3
        check<8,4>(1U<<20);
#elif HDGROUP_TOEPLITZ_COMPACT7
        check<4,4>(1U<<20);
#elif HDGROUP_TOEPLITZ
        if constexpr(HDGROUP_SUBFIELD){check<2,4>(1U<<20);check<4,4>(1U<<20);}
        else{check<2,2>(1U<<20);check<4,2>(1U<<20);}
#else
        if constexpr(HDGROUP_SUBFIELD){check<4,4>(1U<<20);check<8,4>(1U<<20);}
        else if constexpr(HDGROUP_FIELD32)check<2,2>(1U<<20);
        else{check<4,2>(1U<<20);check<8,2>(1U<<20);}
#endif
#if HDGROUP_WIDE_INNER
#if HDGROUP_LONG_RS
        check<1,16>(1U<<20);
#elif HDGROUP_XOR_MIXER==5
        check<16,16>(1U<<20);
#elif HDGROUP_XOR_MIXER>=3
        check<8,16>(1U<<20);
#elif HDGROUP_TOEPLITZ_COMPACT7
        check<4,16>(1U<<20);
#elif HDGROUP_TOEPLITZ
        check<2,16>(1U<<20);check<4,16>(1U<<20);
#else
        check<4,16>(1U<<20);check<8,16>(1U<<20);
#endif
#endif
        return 0;
    }
    if(argc<2||std::string(argv[1])!="--bench")throw std::invalid_argument("--check | --large-check | --bench [calls] [full|outer|inner] [logK] [nt|cached] "
#if HDGROUP_LONG_RS
        "[group:1] [rank:2|4"
#elif HDGROUP_XOR_MIXER==5
        "[group:16] [rank:2|4"
#elif HDGROUP_XOR_MIXER>=3
        "[group:8] [rank:2|4"
#elif HDGROUP_TOEPLITZ_COMPACT7
        "[group:4] [rank:2|4"
#else
        "[group:2|4|8|16] [rank:2|4"
#endif
#if HDGROUP_WIDE_INNER
        "|16"
#endif
        "]");
    const unsigned calls=argc>2?std::stoul(argv[2]):101,logK=argc>4?std::stoul(argv[4]):20;
    const std::string mode=argc>3?argv[3]:"full",stores=argc>5?argv[5]:"nt";
#if HDGROUP_LONG_RS
    const unsigned group=argc>6?std::stoul(argv[6]):1,rank=argc>7?std::stoul(argv[7]):4;
#elif HDGROUP_XOR_MIXER==5
    const unsigned group=argc>6?std::stoul(argv[6]):16,rank=argc>7?std::stoul(argv[7]):4;
#elif HDGROUP_XOR_MIXER>=3
    const unsigned group=argc>6?std::stoul(argv[6]):8,rank=argc>7?std::stoul(argv[7]):4;
#else
    const unsigned group=argc>6?std::stoul(argv[6]):4,rank=argc>7?std::stoul(argv[7]):2;
#endif
    if(!calls||logK<13||logK>22||(mode!="full"&&mode!="outer"&&mode!="inner")||(stores!="nt"&&stores!="cached"))throw std::invalid_argument("invalid benchmark parameters");
    if(rank==2)dispatch<2>(group,stores=="nt",calls,logK,mode);
    else if(rank==4)dispatch<4>(group,stores=="nt",calls,logK,mode);
#if HDGROUP_WIDE_INNER
    else if(rank==16)dispatch<16>(group,stores=="nt",calls,logK,mode);
#endif
    else throw std::invalid_argument("invalid inner rank");
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
