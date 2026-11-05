#pragma once
#include "SparseByteField128Simd.h"

// Select only the preserved paired-route fixed circuits. The new field keys
// and variable products below do not use the old byte-key implementations.
#ifndef HDGROUP_SUBFIELD
#define HDGROUP_SUBFIELD 1
#endif
#ifndef HDGROUP_FUSE_LAST
#define HDGROUP_FUSE_LAST 1
#endif
#ifndef HDGROUP_BYTE_KEYS
#define HDGROUP_BYTE_KEYS 1
#endif
#ifndef HD_SUMRANK_BYTE_KEYS
#define HD_SUMRANK_BYTE_KEYS 1
#endif
#ifndef HDGROUP_INNER
#define HDGROUP_INNER 5
#endif
#ifndef HDGROUP_FUSE_REFRESH
#define HDGROUP_FUSE_REFRESH 1
#endif
#ifndef HDGROUP_NATIVE_OUTPUT
#define HDGROUP_NATIVE_OUTPUT 1
#endif
#ifndef HDGROUP_TOEPLITZ
#define HDGROUP_TOEPLITZ 1
#endif
#ifndef HDGROUP_PACKET_ROUTE
#define HDGROUP_PACKET_ROUTE 1
#endif
#ifndef HDGROUP_ROUTE_WIDTH
#define HDGROUP_ROUTE_WIDTH 2
#endif
#ifndef HDGROUP_SHIFTED_PAIRS
#define HDGROUP_SHIFTED_PAIRS 2
#endif
#ifndef HDGROUP_PAIR_SWAPS
#define HDGROUP_PAIR_SWAPS 1
#endif
#include "../HdGroupedOuter.h"
#include "../HdGroupedInner.h"

namespace spin::research::char2::field128_lift::encoder {
namespace field=sparse_byte;
namespace old=spin::research::char2::hdgroup;
namespace sr=spin::research::char2::hdsumrank;
using Element=field::Element;
using Block=spin::research::char2::Block;
static_assert(HDGROUP_SUBFIELD==1&&HDGROUP_FUSE_LAST==1&&HDGROUP_BYTE_KEYS==1&&HD_SUMRANK_BYTE_KEYS==1);
static_assert(HDGROUP_INNER==5&&HDGROUP_FUSE_REFRESH==1&&HDGROUP_NATIVE_OUTPUT==1);
static_assert(HDGROUP_TOEPLITZ==1&&HDGROUP_PACKET_ROUTE==1&&HDGROUP_ROUTE_WIDTH==2);
static_assert(HDGROUP_SHIFTED_PAIRS==2&&HDGROUP_PAIR_SWAPS==1);
static_assert(HDGROUP_TOEPLITZ_COMPACT7==0&&HDGROUP_XOR_MIXER==0&&HDGROUP_LONG_RS==0&&HDGROUP_XOR_REFRESH==0);

#if defined(_MSC_VER)
#define LIFT_INLINE __forceinline
#define LIFT_NOINLINE __declspec(noinline)
#else
#define LIFT_INLINE inline __attribute__((always_inline))
#define LIFT_NOINLINE __attribute__((noinline))
#endif

inline bool nonzero(Element a) noexcept {return (a.lo|a.hi)!=0;}
// Untimed setup uses the checked arithmetic kernel rather than running a
// byte-polynomial oracle for millions of elimination products.
LIFT_INLINE Element setupMultiply(Element a,Element b) noexcept {
    const auto aa=_mm512_broadcast_i32x4(_mm_loadu_si128(reinterpret_cast<const __m128i*>(&a)));
    const auto bb=_mm512_broadcast_i32x4(_mm_loadu_si128(reinterpret_cast<const __m128i*>(&b)));
    Element result;
    _mm_storeu_si128(reinterpret_cast<__m128i*>(&result),_mm512_castsi512_si128(field::multiply4(aa,bb)));
    return result;
}
struct Refresh {std::array<Element,16> before{},after{};};
static_assert(sizeof(Refresh)==512);

// Deterministic research setup only: MT19937-64 streams are reproducible
// test/benchmark inputs, not ideal independent randomness and not a
// cryptographic seeded implementation of the uniform-setup ensemble.
struct Plan {
    std::size_t k,n,rows,stride;
    bool restrictedToGf256;
    std::size_t toeplitzRejected=0;
    std::vector<std::array<Element,7>> toeplitz;
    // Seven stored keys still evaluate the same nine-product Toeplitz4.
    std::vector<Element> outerKeys;
    std::vector<std::array<Element,16>> before,after;
    std::vector<Refresh> refresh;
    std::vector<std::uint16_t> route,packetRoute;

    static bool nonsingular(const std::array<Element,7>& t){
        Element a[4][4];
        for(unsigned i=0;i!=4;++i)for(unsigned j=0;j!=4;++j)a[i][j]=t[3+i-j];
        for(unsigned c=0;c!=4;++c){
            auto pivot=c;while(pivot!=4&&!nonzero(a[pivot][c]))++pivot;
            if(pivot==4)return false;
            for(unsigned j=c;j!=4;++j)std::swap(a[c][j],a[pivot][j]);
            for(unsigned i=c+1;i!=4;++i){
                const auto factor=a[i][c];if(!nonzero(factor))continue;
                for(unsigned j=c+1;j!=4;++j)
                    a[i][j]=setupMultiply(a[c][c],a[i][j])^setupMultiply(factor,a[c][j]);
                a[i][c]={};
            }
        }
        return true;
    }
    static std::array<Element,9> nineKeys(const std::array<Element,7>& t){
        std::array<Element,9> result;
        const auto pair=[&](unsigned i,Element a,Element b,Element c){result[i]=a;result[i+1]=a^b;result[i+2]=a^c;};
        pair(0,t[3],t[2],t[4]);
        pair(3,t[5]^t[3],t[4]^t[2],t[6]^t[4]);
        pair(6,t[1]^t[3],t[0]^t[2],t[2]^t[4]);
        return result;
    }
    explicit Plan(std::size_t size,std::uint64_t seed=113,bool gf256=false)
        :k(size),n(2*size),rows(size/128),stride(rows+4),restrictedToGf256(gf256){
        if(size<8192||size%8192||rows>32768)throw std::invalid_argument("literal GF128 dimensions: K multiple8192, K<=2^22");
        std::mt19937_64 rr(seed^0x91af823017ccb91aULL),dr(seed^0xe685639417f293b5ULL),tr(seed^0xb2d174839ae61713ULL);
        std::uniform_int_distribution<unsigned> byte(0,255),nzByte(1,255);
        const auto sample=[&](auto& random){return gf256?field::embedByte(static_cast<std::uint8_t>(byte(random))):Element{random(),random()};};
        const auto nonzeroSample=[&](auto& random){
            if(gf256)return field::embedByte(static_cast<std::uint8_t>(nzByte(random)));
            Element value;do{value={random(),random()};}while(!nonzero(value));return value;
        };
        toeplitz.resize(rows*64);outerKeys.resize(rows*64*7);
        for(auto& t:toeplitz){
            do{for(auto& c:t)c=sample(dr);if(nonsingular(t))break;++toeplitzRejected;}while(true);
        }
        constexpr unsigned selected[7]={0,1,2,3,5,6,7};
        for(std::size_t g=0;g!=rows;g+=4)for(unsigned c=0;c!=64;++c)for(unsigned lane=0;lane!=4;++lane){
            const auto keys=nineKeys(toeplitz[c*rows+g+lane]);
            for(unsigned i=0;i!=7;++i)outerKeys[g*448+4*(c*7+i)+lane]=keys[selected[i]];
        }
        route.resize(n);packetRoute.resize(n/2);
        std::mt19937_64 swapRandom(seed^0x28622e73287794afULL);
        std::uniform_int_distribution<unsigned> swapBit(0,1);
        const auto pairs=rows/2;
        for(unsigned c=0;c!=256;++c){
            auto* p=packetRoute.data()+c*pairs;
            for(std::size_t j=0;j!=pairs;++j)p[j]=std::uint16_t(2*j+(c&1));
            for(auto j=pairs-1;j;--j){std::uniform_int_distribution<std::size_t> d(0,j);std::swap(p[j],p[d(rr)]);}
            for(std::size_t j=0;j!=pairs;++j)p[j]|=std::uint16_t(swapBit(swapRandom)<<15);
            for(std::size_t first=0;first!=rows;first+=64)for(unsigned j=0;j!=16;++j)for(unsigned lane=0;lane!=4;++lane){
                const unsigned entry=p[first/2+2*j+lane/2];
                const auto record=std::size_t(entry&0x7fffU)+((lane%2)^(entry>>15));
                route[c*rows+first+16*lane+j]=std::uint16_t(record%rows);
            }
        }
        before.resize(n/64);after.resize(n/64);refresh.resize(n/64);
        for(std::size_t i=0;i!=n/64;++i){
            for(unsigned j=0;j!=16;++j){before[i][j]=nonzeroSample(tr);after[i][j]=nonzeroSample(tr);}
            for(unsigned j=0;j!=16;++j){
                refresh[i].before[j]=before[i][4*(j%4)+j/4];
                refresh[i].after[j]=after[i][4*(j%4)+j/4];
            }
        }
    }
    std::size_t scratchSize() const noexcept{return 256*stride;}
    std::size_t hotKeyBytes() const noexcept{return outerKeys.size()*sizeof(Element)+refresh.size()*sizeof(Refresh);}
    std::size_t setupBytes() const noexcept{
        return hotKeyBytes()+toeplitz.size()*sizeof(toeplitz[0])+(before.size()+after.size())*sizeof(before[0])+
               (route.size()+packetRoute.size())*sizeof(std::uint16_t);
    }
};

// Unchanged nine-product Toeplitz4 circuit. Keys are A0,A1,A2,V0,V2,W0,W1;
// V1=V0+(A1+A2), W2=W0+(A1+A2). XOR reconstruction is field-independent.
LIFT_INLINE void toeplitz4(__m512i* x,const Element* p){
    const auto a0=_mm512_loadu_si512(p),a1=_mm512_loadu_si512(p+4),a2=_mm512_loadu_si512(p+8);
    const auto keyXor=_mm512_xor_si512(a1,a2);
    const auto low=_mm512_xor_si512(x[0],x[2]),high=_mm512_xor_si512(x[1],x[3]);
    const auto sharedU=field::multiply4(_mm512_xor_si512(low,high),a0);
    const auto u0=_mm512_xor_si512(sharedU,field::multiply4(high,a1));
    const auto u1=_mm512_xor_si512(sharedU,field::multiply4(low,a2));
    const auto v0=_mm512_loadu_si512(p+12),v2=_mm512_loadu_si512(p+16);
    const auto sharedV=field::multiply4(_mm512_xor_si512(x[0],x[1]),v0);
    const auto outV0=_mm512_xor_si512(sharedV,field::multiply4(x[1],_mm512_xor_si512(v0,keyXor)));
    const auto outV1=_mm512_xor_si512(sharedV,field::multiply4(x[0],v2));
    const auto w0=_mm512_loadu_si512(p+20),w1=_mm512_loadu_si512(p+24);
    const auto sharedW=field::multiply4(_mm512_xor_si512(x[2],x[3]),w0);
    const auto outW0=_mm512_xor_si512(sharedW,field::multiply4(x[3],w1));
    const auto outW1=_mm512_xor_si512(sharedW,field::multiply4(x[2],_mm512_xor_si512(w0,keyXor)));
    x[0]=_mm512_xor_si512(u0,outW0);x[1]=_mm512_xor_si512(u1,outW1);
    x[2]=_mm512_xor_si512(u0,outV0);x[3]=_mm512_xor_si512(u1,outV1);
}
template<bool Stream>
LIFT_INLINE void mixStore(__m512i* column,unsigned j,Block* planes,std::size_t stride,const Element* keys){
    toeplitz4(column,keys+28*j);
    c8Unroll(std::make_index_sequence<4>{},[&](auto r){
        auto* out=planes+(4*j+r)*stride;
        if constexpr(Stream)_mm512_stream_si512(reinterpret_cast<__m512i*>(out),column[r]);
        else _mm512_storeu_si512(out,column[r]);
    });
}
template<unsigned Offset,bool Stream>
LIFT_NOINLINE inline void coset(const Block* source,Block* planes,std::size_t stride,const Element* keys){
    static_assert(Offset==0||Offset==32);
    alignas(64) __m512i values[128];
    for(unsigned r=0;r!=4;++r){
        c8Unroll(std::make_index_sequence<32>{},[&](auto j){values[32*r+j]=_mm512_loadu_si512(source+4*(32*r+j));});
        old::rsPrefix<32,Offset>(values+32*r);
    }
    for(unsigned j=0;j!=32;j+=2){
        alignas(64) __m512i even[4],odd[4];
        const auto coefficient=_mm512_set1_epi8(static_cast<char>(Offset+j));
        c8Unroll(std::make_index_sequence<4>{},[&](auto r){
            const auto a=values[32*r+j],b=values[32*r+j+1];
            even[r]=_mm512_xor_si512(a,_mm512_gf2p8mul_epi8(b,coefficient));
            odd[r]=_mm512_xor_si512(even[r],b);
        });
        mixStore<Stream>(even,Offset+j,planes,stride,keys);
        mixStore<Stream>(odd,Offset+j+1,planes,stride,keys);
    }
}
template<bool Stream>
LIFT_NOINLINE inline void outer(const Plan& p,const Block* input,Block* planes){
    if constexpr(Stream)if(reinterpret_cast<std::uintptr_t>(planes)&63U){outer<false>(p,input,planes);return;}
    constexpr std::size_t tile=64;
    for(std::size_t first=0;first<p.rows;first+=tile){
        const auto last=std::min(first+tile,p.rows);
        for(auto g=first;g<last;g+=4)coset<0,Stream>(input+g*128,planes+g,p.stride,p.outerKeys.data()+g*448);
        for(auto g=first;g<last;g+=4)coset<32,Stream>(input+g*128,planes+g,p.stride,p.outerKeys.data()+g*448);
    }
    if constexpr(Stream)_mm_sfence();
    for(unsigned c=0;c!=256;++c){
        auto* plane=planes+c*p.stride;
        plane[p.rows]=plane[0];
        std::memset(plane+p.rows+1,0,3*sizeof(Block));
    }
}

LIFT_INLINE void finishStriped(sr::State& state,sr::State feedback,const Refresh& t){
    c8Unroll(std::make_index_sequence<4>{},[&](auto i){
        feedback.values[i]=field::multiply4(_mm512_xor_si512(state.values[i],feedback.values[i]),t.before.data()+4*i);
    });
    old::stripedFft16<0,hdbw::detail::Kind::Inverse,true>(feedback);
    constexpr auto c=hdbw::detail::normalized(3,16);
    c8Unroll(std::make_index_sequence<4>{},[&](auto i){
        const auto x=feedback.values[i];
        const auto sum=_mm512_xor_si512(x,_mm512_shuffle_i32x4(x,x,0x4e));
        feedback.values[i]=_mm512_xor_si512(x,sr::byteProduct4<true,c>(sum));
    });
    old::stripedFft16<16,hdbw::detail::Kind::Forward,true>(feedback);
    c8Unroll(std::make_index_sequence<4>{},[&](auto i){state.values[i]=field::multiply4(feedback.values[i],t.after.data()+4*i);});
}
template<bool Stream>
LIFT_INLINE void innerPacket(const Block* input,const std::uint16_t* route,Block* output,sr::State& state,const Refresh& refresh){
    sr::State q=state,feedback;
    old::innerCoset<0,Stream>(input,route,output,q,feedback);
    old::innerCoset<4,Stream>(input,route,output,q,feedback);
    old::innerCoset<8,Stream>(input,route,output,q,feedback);
    old::innerCoset<12,Stream>(input,route,output,q,feedback);
    finishStriped(state,feedback,refresh);
}
template<bool Stream>
LIFT_NOINLINE inline void inner(const Plan& p,const Block* planes,Block* output){
    sr::State state{};
    for(unsigned region=0;region!=256;++region){
        const auto begin=region*p.rows;const auto* source=planes+region*p.stride;
        for(std::size_t j=0;j!=p.rows;j+=64){
            if(region+1<256){
                const auto* next=reinterpret_cast<const char*>(planes+(region+1)*p.stride+j);
                for(unsigned h=0;h!=16;++h)_mm_prefetch(next+64*h,_MM_HINT_T1);
            }
            innerPacket<Stream>(source,p.packetRoute.data()+(begin+j)/2,output+begin+j,state,p.refresh[(begin+j)/64]);
        }
    }
    if constexpr(Stream)_mm_sfence();
}
template<bool Stream=true>
inline void forward(const Plan& p,const Block* input,Block* output,Block* scratch){
    // Completing the outer first permits exact or partial input/output alias.
    // Scratch must remain disjoint from both payload ranges.
    outer<true>(p,input,scratch);
    if constexpr(Stream){
        if(reinterpret_cast<std::uintptr_t>(output)&63U)inner<false>(p,scratch,output);
        else inner<true>(p,scratch,output);
    }else inner<false>(p,scratch,output);
}

#undef LIFT_INLINE
#undef LIFT_NOINLINE
} // namespace spin::research::char2::field128_lift::encoder
