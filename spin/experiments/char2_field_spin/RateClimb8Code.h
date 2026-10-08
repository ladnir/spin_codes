#pragma once
#include "RateClimb8Inner.h"
#include "Cubic8Code.h"

namespace spin::research::char2 {

// Research-only coefficient-basis GRS[128,d] profiles. Setup fixtures are
// deterministic test sampling, outside the full precomputed transpose.
// planes and temp are caller-owned and disjoint from input/output. A custom
// planeStride must be >=groups and divisible by4 for the regional NT path.
template<class Profile>
struct RateClimb8Plan {
    using profile_type=Profile;
    static constexpr unsigned dimension=Profile::dimension,outer=128,width=8,state=8;
    static constexpr unsigned period=Profile::period,count=Profile::count,tileRecords=Profile::tileRecords;
    std::size_t requested,k,inputN,outputN,groups,planeStride;
    std::vector<std::uint32_t> route,inverseRoute;
    std::vector<std::uint16_t> inverseLocalRoute;
    std::vector<std::uint8_t> inputDiagonal,stateDiagonal;
    std::vector<std::uint32_t> foldedRoute,foldedInverseRoute;
    std::vector<std::uint16_t> foldedInverseLocalRoute;
    std::vector<std::uint8_t> packedInputDiagonal;
    explicit RateClimb8Plan(std::size_t message,std::uint64_t seed):requested(message){
        static_assert(dimension<128&&4%period==0&&64%(width*period)==0);
        static_assert(tileRecords*period==4*count);
        constexpr std::size_t unit=dimension*64;
        constexpr std::size_t maxRows=(std::numeric_limits<std::uint32_t>::max()/outer/64)*64;
        if(!message||message>maxRows*dimension)
            throw std::invalid_argument("RateClimb8 message size outside uint32 route geometry");
        groups=((message+unit-1)/unit)*64;k=groups*dimension;
        inputN=groups*outer;outputN=(inputN/(width*period))*count;planeStride=groups;
        // Same streams as the frozen plan permit a byte-for-byte Fixed7 control.
        std::mt19937_64 rr(seed^0xdf4a1b3285ec7096ULL),ir(seed^0x2e7d18c4f0b3a965ULL),
                        sr(seed^0x918b745cf62d03eaULL);
        route.resize(inputN);inverseRoute.resize(inputN);
        inputDiagonal.resize(inputN);stateDiagonal.resize(inputN);
        for(unsigned r=0;r!=outer;++r){
            auto first=route.begin()+r*groups;
            std::iota(first,first+groups,static_cast<std::uint32_t>(r*groups));
            for(std::size_t j=groups-1;j;--j){
                std::uniform_int_distribution<std::size_t> draw(0,j);
                std::swap(first[j],first[draw(rr)]);
            }
        }
        for(std::size_t i=0;i!=inputN;++i)inverseRoute[route[i]]=static_cast<std::uint32_t>(i);
        if(groups<=65536){
            inverseLocalRoute.resize(inputN);
            for(unsigned r=0;r!=outer;++r){
                const auto offset=r*groups;
                for(std::size_t j=0;j!=groups;++j)
                    inverseLocalRoute[offset+j]=static_cast<std::uint16_t>(inverseRoute[offset+j]-offset);
            }
        }
        for(auto& d:inputDiagonal)d=sfSampleNonzero(ir);
        for(auto& d:stateDiagonal)d=sfSampleNonzero(sr);
        // Absorb the final four-packet record transpose into setup. The
        // original tables remain intact for the baseline and scalar oracle.
        foldedRoute.resize(inputN);foldedInverseRoute.resize(inputN);
        packedInputDiagonal.resize(inputN);
        for(std::size_t i=0;i!=inputN;++i){
            const auto packed=foldedIndex(i);
            foldedRoute[packed]=route[i];
            foldedInverseRoute[route[i]]=static_cast<std::uint32_t>(packed);
            packedInputDiagonal[packed]=inputDiagonal[i];
        }
        if(groups<=65536){
            foldedInverseLocalRoute.resize(inputN);
            for(unsigned r=0;r!=outer;++r){
                const auto offset=r*groups;
                for(std::size_t j=0;j!=groups;++j)
                    foldedInverseLocalRoute[offset+j]=static_cast<std::uint16_t>(foldedInverseRoute[offset+j]-offset);
            }
        }
    }
    static constexpr std::size_t foldedIndex(std::size_t index)noexcept{
        return 32*(index/32)+4*(index%8)+(index%32)/8;
    }
    static constexpr unsigned retained(std::size_t step)noexcept{return Profile::retained(step);}
    static constexpr std::size_t sourceIndex(std::size_t row,unsigned coordinate)noexcept{
        return (row/4)*(4*dimension)+4*coordinate+row%4;
    }
    // coordinate must be strictly below retained(step).
    static constexpr std::size_t outputIndex(std::size_t step,unsigned coordinate)noexcept{
        const auto base=tileRecords*(step/4);
        const auto phase=step%4;
        if(coordinate<Profile::commonRetained)return base+4*coordinate+phase;
        if constexpr(Profile::tailLayout==R8TailLayout::Phase3)return base+24;
        else if constexpr(Profile::tailLayout==R8TailLayout::OddPhases)return base+20+phase/2;
        else return base+4*coordinate+phase;
    }
    std::size_t planeRecords()const noexcept{return 128*planeStride;}
};

template<bool Scatter,bool Folded,class Profile>
inline void r8InnerRange(const RateClimb8Plan<Profile>& p,const Block* input,Block* dest,
                          std::size_t first,std::size_t past,std::size_t base,
                          __m512i& h0,__m512i& h1){
    alignas(64) Block packet[32];
    const auto* diagonals=Folded?p.packedInputDiagonal.data():p.inputDiagonal.data();
    const auto* routing=Folded?p.foldedRoute.data():p.route.data();
    for(auto step=past;step!=first;){
        step-=4;const auto begin=8*step;
        r8InnerTransposeFour<Profile,Folded>(input+Profile::tileRecords*(step/4),packet,h0,h1,
                                            p.stateDiagonal.data()+begin);
        // Keep the frozen explicit/unrolled stores. Forced inlining allows
        // the compiler to eliminate the apparent 512-byte packet temporary.
        c8Unroll(std::make_index_sequence<8>{},[&](auto chunk){
            constexpr unsigned shift=4*unsigned(chunk);
            const auto value=_mm512_gf2p8mul_epi8(_mm512_load_si512(packet+shift),
                sfBroadcastFour(diagonals+begin+shift));
            if constexpr(Scatter){
                const auto* index=routing+begin+shift;
                dest[index[0]-base]=block(_mm512_castsi512_si128(value));
                dest[index[1]-base]=block(_mm512_extracti32x4_epi32(value,1));
                dest[index[2]-base]=block(_mm512_extracti32x4_epi32(value,2));
                dest[index[3]-base]=block(_mm512_extracti32x4_epi32(value,3));
            }else _mm512_storeu_si512(dest+begin-base+shift,value);
        });
    }
}
template<bool Stream,bool Compact,bool Folded,class Profile>
inline void r8RegionInner(const RateClimb8Plan<Profile>& p,const Block* input,Block* planes,Block* temp){
    auto h0=_mm512_setzero_si512(),h1=_mm512_setzero_si512();
    const auto* inverseLocal=Folded?p.foldedInverseLocalRoute.data():p.inverseLocalRoute.data();
    const auto* inverse=Folded?p.foldedInverseRoute.data():p.inverseRoute.data();
    for(unsigned r=128;r--;){
        const auto offset=r*p.groups;
        auto* region=planes+r*p.planeStride;
        r8InnerRange<false,Folded>(p,input,temp,offset/8,(offset+p.groups)/8,offset,h0,h1);
        if constexpr(Compact){
            for(std::size_t j=0;j!=p.groups;j+=4)
                sf16Store<Stream>(region+j,c8GatherLocal4(temp,inverseLocal+offset+j));
        }else{
            for(std::size_t j=0;j!=p.groups;j+=4)
                sf16Store<Stream>(region+j,sf16GatherRegion4(temp,inverse+offset+j,
                    static_cast<std::uint32_t>(offset)));
        }
    }
    if constexpr(Stream)_mm_sfence();
}
template<bool Folded,class Profile>
inline void r8InnerDispatchImpl(const RateClimb8Plan<Profile>& p,const Block* input,Block* planes,
                                Block* temp,Subfield16Kernel kernel){
    if(kernel==Subfield16Kernel::Scatter){
        auto h0=_mm512_setzero_si512(),h1=_mm512_setzero_si512();
        for(unsigned r=128;r--;){
            const auto offset=r*p.groups;
            r8InnerRange<true,Folded>(p,input,planes+r*p.planeStride,offset/8,(offset+p.groups)/8,offset,h0,h1);
        }
    }else{
        if(!temp)throw std::invalid_argument("RateClimb8 regional schedule requires temporary storage");
        const bool aligned=(reinterpret_cast<std::uintptr_t>(planes)&63)==0;
        const bool compact=Folded?!p.foldedInverseLocalRoute.empty():!p.inverseLocalRoute.empty();
        if(!compact){
            if(aligned)r8RegionInner<true,false,Folded>(p,input,planes,temp);
            else r8RegionInner<false,false,Folded>(p,input,planes,temp);
        }else{
            if(aligned)r8RegionInner<true,true,Folded>(p,input,planes,temp);
            else r8RegionInner<false,true,Folded>(p,input,planes,temp);
        }
    }
}
template<class Profile>
inline void r8InnerDispatch(const RateClimb8Plan<Profile>& p,const Block* input,Block* planes,
                            Block* temp,Subfield16Kernel kernel,bool folded=false){
    if(folded)r8InnerDispatchImpl<true>(p,input,planes,temp,kernel);
    else r8InnerDispatchImpl<false>(p,input,planes,temp,kernel);
}

#if defined(_MSC_VER)
#define R8_NOINLINE __declspec(noinline)
#else
#define R8_NOINLINE __attribute__((noinline))
#endif
template<bool Stream,bool Fused,class Profile>
R8_NOINLINE void r8OuterFinish(const RateClimb8Plan<Profile>& p,const Block* planes,Block* out){
    constexpr unsigned dimension=Profile::dimension;
    alignas(64) __m512i values[128];
    for(std::size_t row=0;row*dimension<p.requested;row+=4){
        if constexpr(Fused)c128OuterCoefficientTransposeStrided(planes+row,p.planeStride*sizeof(Block),values);
        else{
            c8Unroll(std::make_index_sequence<128>{},[&](auto j){
                values[j]=_mm512_loadu_si512(planes+std::size_t(j)*p.planeStride+row);
            });
            c128OuterCoefficientTransposeInplace(values);
        }
        const auto base=row*dimension;
        if(base+4*dimension<=p.requested){
            c8Unroll(std::make_index_sequence<dimension>{},[&](auto j){
                sf16Store<Stream>(out+base+4*j,values[j]);
            });
        }else{
            for(unsigned j=0;j!=dimension&&base+4*j<p.requested;++j){
                const auto left=p.requested-base-4*j;
                if(left>=4)sf16Store<Stream>(out+base+4*j,values[j]);
                else{alignas(64) Block tail[4];_mm512_store_si512(tail,values[j]);
                     std::memcpy(out+base+4*j,tail,left*sizeof(Block));}
            }
        }
    }
    if constexpr(Stream)_mm_sfence();
}
#undef R8_NOINLINE
template<class Profile>
inline void r8OuterDispatch(const RateClimb8Plan<Profile>& p,const Block* planes,Block* out,
                            Subfield16Kernel kernel,Cubic8OuterVariant variant){
    const bool stream=kernel==Subfield16Kernel::RegionAllNonTemporal&&(reinterpret_cast<std::uintptr_t>(out)&63)==0;
    if(variant==Cubic8OuterVariant::FusedLoad){
        if(stream)r8OuterFinish<true,true>(p,planes,out);else r8OuterFinish<false,true>(p,planes,out);
    }else{
        if(stream)r8OuterFinish<true,false>(p,planes,out);else r8OuterFinish<false,false>(p,planes,out);
    }
}
// Consume the entire compact input before writing the requested output
// prefix. Exact input/output alias preserves the unused input suffix.
template<class Profile>
inline void r8Transpose(const RateClimb8Plan<Profile>& p,const Block* input,Block* out,Block* planes,
                         Block* temp,Subfield16Kernel kernel,
                         Cubic8OuterVariant variant=Cubic8OuterVariant::Gather,bool folded=false){
    r8InnerDispatch(p,input,planes,temp,kernel,folded);r8OuterDispatch(p,planes,out,kernel,variant);
}

} // namespace spin::research::char2
