#pragma once
#include "RateSpeed4Code.h"
#include "RateSpeed4OuterNext.h"
#include "WideDiagonal16.h"
#include "WideDiagonal32.h"
#include <chrono>
#include <type_traits>

namespace spin::research::char2::hdwide {

// Full transposed-encoder COST PROTOTYPE, not a distance certificate.
// Embedded GF256 GRS[B,B*Keep/8], one regional permutation, full normalized
// four-input/four-state Cauchy packet, first Keep outputs retained. Both
// diagonals are uniform nonzero scalars in the actual selected wide field.
template<unsigned Bits>struct Field;
template<>struct Field<16>{
    using Scalar=std::uint16_t;using Coefficient=wide16::Coefficient;
    static Coefficient prepare(Scalar x){return wide16::prepareCoefficient(x);}
    static __m512i multiply(__m512i x,const Coefficient* c){return wide16::multiply512(x,c);}
};
template<>struct Field<32>{
    using Scalar=std::uint32_t;using Coefficient=wide32::Coefficient;
    static Coefficient prepare(Scalar x){return wide32::prepareCoefficient(x);}
    static __m512i multiply(__m512i x,const Coefficient* c){return wide32::multiply512(x,c);}
};

template<unsigned Bits,unsigned Keep,unsigned B=64>struct Plan {
    static_assert((Bits==16||Bits==32)&&(Keep==1||Keep==2||Keep==4));
    static_assert(B==64||(B==256&&Keep==2));
    using F=Field<Bits>;using Scalar=typename F::Scalar;using Coefficient=typename F::Coefficient;
    static constexpr unsigned bits=Bits,keep=Keep,dimension=B*Keep/8,outerLength=B;
    std::size_t requested,k,groups,inputN,outputN,planeStride;
    std::vector<std::uint32_t> route,foldedInverse;
    std::vector<std::uint16_t> localInverse;
    // Raw scalars are retained for independent validation, never read hot.
    std::vector<Scalar> dx,ds;
    std::vector<Coefficient> stateCoefficient,routedCoefficient;
    explicit Plan(std::size_t message,std::uint64_t seed):requested(message){
        constexpr std::size_t unit=dimension*64;
        constexpr std::size_t maxRows=(std::numeric_limits<std::uint32_t>::max()/B/64)*64;
        if(!message||message>maxRows*dimension)throw std::invalid_argument("wide cost route geometry");
        groups=((message+unit-1)/unit)*64;k=dimension*groups;
        inputN=B*groups;outputN=(inputN/4)*Keep;planeStride=groups+4;
        route.resize(inputN);foldedInverse.resize(inputN);dx.resize(inputN);ds.resize(inputN);
        stateCoefficient.resize(inputN);routedCoefficient.resize(inputN);
        std::mt19937_64 rr(seed^0xdf4a1b3285ec7096ULL),ir(seed^0x2e7d18c4f0b3a965ULL),sr(seed^0x918b745cf62d03eaULL);
        for(unsigned region=0;region!=B;++region){
            auto first=route.begin()+region*groups;
            std::iota(first,first+groups,static_cast<std::uint32_t>(region*groups));
            for(std::size_t j=groups-1;j;--j){std::uniform_int_distribution<std::size_t> draw(0,j);std::swap(first[j],first[draw(rr)]);}
        }
        // Truncate a uniform 64-bit draw and reject zero: exact uniformity
        // over all 2^Bits-1 nonzero field elements, not a restricted alphabet.
        for(std::size_t i=0;i!=inputN;++i){
            do{dx[i]=static_cast<Scalar>(ir());}while(!dx[i]);
            do{ds[i]=static_cast<Scalar>(sr());}while(!ds[i]);
            stateCoefficient[i]=F::prepare(ds[i]);routedCoefficient[route[i]]=F::prepare(dx[i]);
            foldedInverse[route[i]]=static_cast<std::uint32_t>(16*(i/16)+4*(i%4)+(i%16)/4);
        }
        if(groups<=65536){localInverse.resize(inputN);for(std::size_t i=0;i!=inputN;++i)localInverse[i]=static_cast<std::uint16_t>(foldedInverse[i]-(i/groups)*groups);}
    }
    std::size_t planeRecords()const noexcept{return B*planeStride;}
    static constexpr std::size_t sourceIndex(std::size_t row,unsigned j){return (row/4)*(4*dimension)+4*j+row%4;}
    static constexpr std::size_t outputIndex(std::size_t step,unsigned j){return 4*Keep*(step/4)+4*j+step%4;}
};

#if defined(_MSC_VER)
#define HDW_INLINE __forceinline
#define HDW_NOINLINE __declspec(noinline)
#else
#define HDW_INLINE inline __attribute__((always_inline))
#define HDW_NOINLINE __attribute__((noinline))
#endif

template<unsigned Keep>HDW_INLINE void prepareFour(const Block* input,__m512i& z0,__m512i& z1,__m512i& z2,__m512i& z3){
    z0=_mm512_loadu_si512(input);z1=_mm512_setzero_si512();z2=z1;z3=z1;
    if constexpr(Keep>=2)z1=_mm512_loadu_si512(input+4);
    if constexpr(Keep==4){z2=_mm512_loadu_si512(input+8);z3=_mm512_loadu_si512(input+12);}
    rs4Circulant<true,true>(z0,z1,z2,z3);
    c8TransposeRecords(z0,z1,z2,z3,z0,z1,z2,z3);
}
template<unsigned Bits>HDW_INLINE void historyStep(__m512i& moment,__m512i& history,const typename Field<Bits>::Coefficient* ds){
    const auto old=history;
    history=Field<Bits>::multiply(_mm512_xor_si512(old,moment),ds);
    moment=_mm512_xor_si512(old,rs8ConstantMultiply<true,rs4Scale>(moment));
}
template<unsigned Bits>HDW_INLINE void finishFour(__m512i z0,__m512i z1,__m512i z2,__m512i z3,Block* output,__m512i& history,const typename Field<Bits>::Coefficient* ds){
    historyStep<Bits>(z3,history,ds+12);historyStep<Bits>(z2,history,ds+8);
    historyStep<Bits>(z1,history,ds+4);historyStep<Bits>(z0,history,ds);
    c8TransposeRecords(z0,z1,z2,z3,z0,z1,z2,z3);
    rs4Circulant<false,true>(z0,z1,z2,z3);
    _mm512_storeu_si512(output,z0);_mm512_storeu_si512(output+4,z1);
    _mm512_storeu_si512(output+8,z2);_mm512_storeu_si512(output+12,z3);
}
// Same bounded two-live-tile schedule as the frozen 16-packet winner.
template<unsigned Bits,unsigned Keep>HDW_INLINE void inner16(const Block* input,Block* output,__m512i& history,const typename Field<Bits>::Coefficient* ds){
    __m512i z0,z1,z2,z3,w0,w1,w2,w3;
    prepareFour<Keep>(input+12*Keep,w0,w1,w2,w3);prepareFour<Keep>(input+8*Keep,z0,z1,z2,z3);
    finishFour<Bits>(w0,w1,w2,w3,output+48,history,ds+48);
    prepareFour<Keep>(input+4*Keep,w0,w1,w2,w3);finishFour<Bits>(z0,z1,z2,z3,output+32,history,ds+32);
    prepareFour<Keep>(input,z0,z1,z2,z3);finishFour<Bits>(w0,w1,w2,w3,output+16,history,ds+16);
    finishFour<Bits>(z0,z1,z2,z3,output,history,ds);
}
struct PhaseTimes{double inner=0,route=0,outer=0;};
using Clock=std::chrono::steady_clock;
template<bool Timed=false,unsigned Bits,unsigned Keep,unsigned B>
HDW_NOINLINE void regionInner(const Plan<Bits,Keep,B>& p,const Block* input,Block* planes,Block* temp,PhaseTimes* times=nullptr){
    auto history=_mm512_setzero_si512();
    for(unsigned region=B;region--;){
        Clock::time_point start;if constexpr(Timed)start=Clock::now();
        const auto offset=region*p.groups;
        for(auto begin=offset+p.groups;begin!=offset;){begin-=64;inner16<Bits,Keep>(input+(begin/4)*Keep,temp+begin-offset,history,p.stateCoefficient.data()+begin);}
        if constexpr(Timed){times->inner+=std::chrono::duration<double,std::micro>(Clock::now()-start).count();start=Clock::now();}
        auto* target=planes+region*p.planeStride;
        if(!p.localInverse.empty())for(std::size_t j=0;j!=p.groups;j+=4){
            auto value=c8GatherLocal4(temp,p.localInverse.data()+offset+j);
            _mm512_storeu_si512(target+j,Field<Bits>::multiply(value,p.routedCoefficient.data()+offset+j));
        }else for(std::size_t j=0;j!=p.groups;j+=4){
            auto value=sf16GatherRegion4(temp,p.foldedInverse.data()+offset+j,static_cast<std::uint32_t>(offset));
            _mm512_storeu_si512(target+j,Field<Bits>::multiply(value,p.routedCoefficient.data()+offset+j));
        }
        if constexpr(Timed)times->route+=std::chrono::duration<double,std::micro>(Clock::now()-start).count();
    }
}
// Frozen FFT prefix, generic retained coefficient stores. This deliberately
// does not claim the D39/40-only fused outer is tuned for these dimensions.
template<bool Stream,unsigned Bits,unsigned Keep>
HDW_NOINLINE void outer(const Plan<Bits,Keep,64>& p,const Block* planes,Block* output){
    constexpr unsigned D=Plan<Bits,Keep,64>::dimension;
    alignas(64)__m512i values[64];
    for(std::size_t row=0;row*D<p.requested;row+=4){
        rs4OuterPrefix(planes+row,p.planeStride*sizeof(Block),values);const auto base=row*D;
        if(base+4*D<=p.requested)c8Unroll(std::make_index_sequence<D>{},[&](auto j){sf16Store<Stream>(output+base+4*j,rs4OuterCoefficient<static_cast<unsigned>(decltype(j)::value)>(values));});
        else for(unsigned j=0;j!=D&&base+4*j<p.requested;++j){
            const auto value=_mm512_xor_si512(_mm512_load_si512(values+j),_mm512_load_si512(values+32+j));
            const auto left=p.requested-base-4*j;
            if(left>=4)sf16Store<Stream>(output+base+4*j,value);
            else{alignas(64)Block tail[4];_mm512_store_si512(tail,value);std::memcpy(output+base+4*j,tail,left*sizeof(Block));}
        }
    }
    if constexpr(Stream)_mm_sfence();
}
// Upper-half extension of the EXISTING coefficient FFT's sixteen-register
// strided-transpose leaf. The butterfly order and normalized additive basis
// are unchanged; only the compile-time evaluation-coset offset is new.
template<unsigned Offset>HDW_NOINLINE void upperLeaf(const unsigned char* source,std::size_t stride,__m512i* target){
    static_assert(Offset>=128&&Offset<256&&Offset%16==0);
#define HDW_LOAD(I) auto x##I=_mm512_loadu_si512(source+I*stride)
    HDW_LOAD(0);HDW_LOAD(1);HDW_LOAD(2);HDW_LOAD(3);HDW_LOAD(4);HDW_LOAD(5);HDW_LOAD(6);HDW_LOAD(7);
    HDW_LOAD(8);HDW_LOAD(9);HDW_LOAD(10);HDW_LOAD(11);HDW_LOAD(12);HDW_LOAD(13);HDW_LOAD(14);HDW_LOAD(15);
#undef HDW_LOAD
#define HDW_PAIR(A,B,L,D) rs4_outer_next_detail::pair<false,rs4_outer_next_detail::normalized(L,Offset^D)>(x##A,x##B)
    HDW_PAIR(14,15,0,14);HDW_PAIR(12,13,0,12);HDW_PAIR(13,15,1,12);HDW_PAIR(12,14,1,12);
    HDW_PAIR(10,11,0,10);HDW_PAIR(8,9,0,8);HDW_PAIR(9,11,1,8);HDW_PAIR(8,10,1,8);
    HDW_PAIR(11,15,2,8);HDW_PAIR(10,14,2,8);HDW_PAIR(9,13,2,8);HDW_PAIR(8,12,2,8);
    HDW_PAIR(6,7,0,6);HDW_PAIR(4,5,0,4);HDW_PAIR(5,7,1,4);HDW_PAIR(4,6,1,4);
    HDW_PAIR(2,3,0,2);HDW_PAIR(0,1,0,0);HDW_PAIR(1,3,1,0);HDW_PAIR(0,2,1,0);
    HDW_PAIR(3,7,2,0);HDW_PAIR(2,6,2,0);HDW_PAIR(1,5,2,0);HDW_PAIR(0,4,2,0);
    HDW_PAIR(7,15,3,0);HDW_PAIR(6,14,3,0);HDW_PAIR(5,13,3,0);HDW_PAIR(4,12,3,0);
    HDW_PAIR(3,11,3,0);HDW_PAIR(2,10,3,0);HDW_PAIR(1,9,3,0);HDW_PAIR(0,8,3,0);
#undef HDW_PAIR
#define HDW_STORE(I) _mm512_store_si512(target+I,x##I)
    HDW_STORE(0);HDW_STORE(1);HDW_STORE(2);HDW_STORE(3);HDW_STORE(4);HDW_STORE(5);HDW_STORE(6);HDW_STORE(7);
    HDW_STORE(8);HDW_STORE(9);HDW_STORE(10);HDW_STORE(11);HDW_STORE(12);HDW_STORE(13);HDW_STORE(14);HDW_STORE(15);
#undef HDW_STORE
}
template<unsigned Half,unsigned Coefficient>HDW_INLINE void fftStage(__m512i* values){
    c8Unroll(std::make_index_sequence<Half>{},[&](auto j){
        auto low=_mm512_load_si512(values+j),high=_mm512_load_si512(values+Half+j);
        rs4_outer_next_detail::pair<false,Coefficient>(low,high);
        _mm512_store_si512(values+j,low);_mm512_store_si512(values+Half+j,high);
    });
}
template<unsigned Offset>HDW_NOINLINE void coset64(const unsigned char* source,std::size_t stride,__m512i* values){
    using namespace c128_outer_coefficient_detail;
    if constexpr(Offset==0){
        leaf16_0_strided_transpose(source,stride,values);leaf16_16_strided_transpose(source+16*stride,stride,values+16);
        leaf16_32_strided_transpose(source+32*stride,stride,values+32);leaf16_48_strided_transpose(source+48*stride,stride,values+48);
    }else if constexpr(Offset==64){
        leaf16_64_strided_transpose(source,stride,values);leaf16_80_strided_transpose(source+16*stride,stride,values+16);
        leaf16_96_strided_transpose(source+32*stride,stride,values+32);leaf16_112_strided_transpose(source+48*stride,stride,values+48);
    }else{
        static_assert(Offset==128||Offset==192);
        upperLeaf<Offset>(source,stride,values);upperLeaf<Offset+16>(source+16*stride,stride,values+16);
        upperLeaf<Offset+32>(source+32*stride,stride,values+32);upperLeaf<Offset+48>(source+48*stride,stride,values+48);
    }
    fftStage<16,rs4_outer_next_detail::normalized(4,Offset)>(values);
    fftStage<16,rs4_outer_next_detail::normalized(4,Offset+32)>(values+32);
    fftStage<32,rs4_outer_next_detail::normalized(5,Offset)>(values);
}
// Pruned GRS[256,64] additive-basis transpose. Source degrees are <64,
// so the forward top128/top64 stages merely copy the same64 coefficients
// into four cosets. Their adjoints XOR four64-point transpose results:
// no full256 transform, systematic-basis conversion, or dense matrix action.
template<bool Stream,unsigned Bits>HDW_NOINLINE void outer(const Plan<Bits,2,256>& p,const Block* planes,Block* output){
    alignas(64)__m512i accumulator[64],values[64];const auto stride=p.planeStride*sizeof(Block);
    for(std::size_t row=0;row*64<p.requested;row+=4){
        const auto* source=reinterpret_cast<const unsigned char*>(planes+row);const auto base=row*64;
        coset64<0>(source,stride,accumulator);coset64<64>(source+64*stride,stride,values);
        c8Unroll(std::make_index_sequence<64>{},[&](auto j){_mm512_store_si512(accumulator+j,_mm512_xor_si512(_mm512_load_si512(accumulator+j),_mm512_load_si512(values+j)));});
        coset64<128>(source+128*stride,stride,values);
        c8Unroll(std::make_index_sequence<64>{},[&](auto j){_mm512_store_si512(accumulator+j,_mm512_xor_si512(_mm512_load_si512(accumulator+j),_mm512_load_si512(values+j)));});
        coset64<192>(source+192*stride,stride,values);
        if(base+256<=p.requested)c8Unroll(std::make_index_sequence<64>{},[&](auto j){sf16Store<Stream>(output+base+4*j,_mm512_xor_si512(_mm512_load_si512(accumulator+j),_mm512_load_si512(values+j)));});
        else for(unsigned j=0;j!=64&&base+4*j<p.requested;++j){
            const auto value=_mm512_xor_si512(_mm512_load_si512(accumulator+j),_mm512_load_si512(values+j));const auto left=p.requested-base-4*j;
            if(left>=4)sf16Store<Stream>(output+base+4*j,value);
            else{alignas(64)Block tail[4];_mm512_store_si512(tail,value);std::memcpy(output+base+4*j,tail,left*sizeof(Block));}
        }
    }
    if constexpr(Stream)_mm_sfence();
}
// Scratch is disjoint. Exact input/output alias is valid: consume all input
// before writing any final output. No allocation occurs in this entrypoint.
template<bool Timed=false,unsigned Bits,unsigned Keep,unsigned B>
inline void transpose(const Plan<Bits,Keep,B>& p,const Block* input,Block* output,Block* planes,Block* temp,PhaseTimes* times=nullptr){
    if constexpr(Timed)*times={};
    regionInner<Timed>(p,input,planes,temp,times);
    Clock::time_point start;if constexpr(Timed)start=Clock::now();
    if((reinterpret_cast<std::uintptr_t>(output)&63)==0)outer<true>(p,planes,output);else outer<false>(p,planes,output);
    if constexpr(Timed)times->outer=std::chrono::duration<double,std::micro>(Clock::now()-start).count();
}
#undef HDW_INLINE
#undef HDW_NOINLINE
} // namespace spin::research::char2::hdwide
