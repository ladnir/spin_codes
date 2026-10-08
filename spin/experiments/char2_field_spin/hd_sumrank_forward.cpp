#include "HdFastOuter.h"
#include "HdBlockwiseOuter256.h"
#include "HdBlockwiseRandomizedOuter.h"
#include "HdSumrankInner.h"
#include "HdSumrankOuter.h"
#include "BenchPayload.h"
#include "Punctured8Oracle.h"
#include <atomic>
#include <iomanip>
#include <iostream>

namespace c2=spin::research::char2;
namespace sr=c2::hdsumrank;
using Block=c2::Block;
using Clock=std::chrono::steady_clock;
#if defined(_MSC_VER)
#define SR_NOINLINE __declspec(noinline)
#else
#define SR_NOINLINE __attribute__((noinline))
#endif
#if defined(HD_SUMRANK_RESTRICT_PAYLOAD)
#if defined(_MSC_VER)
#define SR_PAYLOAD_RESTRICT __restrict
#else
#define SR_PAYLOAD_RESTRICT __restrict__
#endif
#else
#define SR_PAYLOAD_RESTRICT
#endif
static void require(bool b,const char* m){if(!b)throw std::runtime_error(m);}

// Setup-only log tables. Arithmetic validation below uses the independent
// bit-polynomial AES oracle instead. No table lookup occurs in hot encoding.
struct SetupField {
    std::array<std::uint16_t,131070> exp{};
    std::array<unsigned,65536> log{};
    SetupField(){
        std::uint16_t g=2;
        while(sr::power(g,65535/3)==1||sr::power(g,65535/5)==1||
              sr::power(g,65535/17)==1||sr::power(g,65535/257)==1)++g;
        std::uint16_t x=1;
        for(unsigned i=0;i!=65535;++i){exp[i]=exp[i+65535]=x;log[x]=i;x=sr::mul(x,g);}
        require(x==1,"setup field generator");
    }
    std::uint16_t mul(std::uint16_t a,std::uint16_t b)const{return a&&b?exp[log[a]+log[b]]:0;}
    bool invertible(sr::Matrix a)const{
        for(unsigned j=0;j!=16;++j){
            unsigned p=j;while(p!=16&&!a[p][j])++p;if(p==16)return false;
            std::swap(a[j],a[p]);const auto scale=exp[65535-log[a[j][j]]];
            for(unsigned k=j;k!=16;++k)a[j][k]=mul(a[j][k],scale);
            for(unsigned i=j+1;i!=16;++i){const auto c=a[i][j];for(unsigned k=j;k!=16;++k)a[i][k]^=mul(c,a[j][k]);}
        }return true;
    }
};
static const SetupField& setupField(){static const SetupField f;return f;}
struct Raw {std::array<std::uint16_t,16> u{},v{},u2{},v2{};sr::Matrix t,carry;};
#if defined(HD_SUMRANK_COMPACT_STEPS)
// Shared rank-one refresh never reads u. Keep its three used coefficient
// vectors on exactly three cache lines, without changing any setup coin.
struct alignas(64) CompactStep {
    std::array<c2::wide16::Coefficient,16> v;
    sr::DiagonalRefresh refresh;
};
static_assert(sizeof(CompactStep)==192);
#endif
struct Plan {
    std::size_t k,n,rows,stride;
    unsigned feedbackPeriod=2;
    sr::Fixed fixed;
    sr::Matrix alternateFeedback;
    std::vector<std::uint16_t> route,stripedRoute,diagonal;
#if defined(HD_SUMRANK_ROUTE_BYTES)||defined(HD_SUMRANK_ROUTE_UNITS8)
    std::vector<sr::PreparedOffset> offsetRoute;
#endif
    std::vector<c2::wide16::Coefficient> preparedDiagonal;
    std::vector<c2::wide16::Coefficient> paddedDiagonal;
    std::vector<c2::wide16::Coefficient> routedDiagonal;
    std::vector<sr::Mixer> mixers;
    std::vector<sr::Refresh> refreshes;
    std::vector<sr::DiagonalRefresh> diagonalRefreshes;
    std::vector<Raw> raw;
#if defined(HD_SUMRANK_COMPACT_STEPS)
    std::vector<CompactStep> compactSteps;
#endif
    template<unsigned Variant> const auto& stepMixer(std::size_t i)const{
#if defined(HD_SUMRANK_COMPACT_STEPS)
        if constexpr(Variant==42||Variant==45)return compactSteps[i];
        else
#endif
        return mixers[i];
    }
    template<unsigned Variant> const sr::DiagonalRefresh& stepRefresh(std::size_t i)const{
#if defined(HD_SUMRANK_COMPACT_STEPS)
        if constexpr(Variant==42||Variant==45)return compactSteps[i].refresh;
        else
#endif
        return diagonalRefreshes[i];
    }
    bool transformedCarry;
    Plan(std::size_t size,std::uint64_t seed,bool retainRaw=false,bool cauchy=false,bool towerBasis=false,bool fftWeights=false,bool systematic=false,bool globalMds=false,bool sortPackets=false,bool sharedDiagonal=false,unsigned localRank=16,bool alternating=false,bool globalNodes=false,bool blockEmission=false,unsigned carry=0,bool rankTwo=false):k(size),n(2*size),rows(size/128),stride(rows+4),fixed(towerBasis,fftWeights,systematic,globalMds),transformedCarry(carry){
        if(size<8192||size%8192||rows>65536)throw std::invalid_argument("K must be a multiple of 8192, at most 2^23");
        if(sharedDiagonal&&(!cauchy||!globalMds))throw std::invalid_argument("joint refresh requires MDS Cauchy maps");
        if(carry&&(!blockEmission||!sharedDiagonal))throw std::invalid_argument("carried-state research parameters");
        if(rankTwo&&(!sharedDiagonal||carry))throw std::invalid_argument("rank-two shared refresh parameters");
        if(localRank!=16){
            if(!sharedDiagonal||(localRank!=4&&localRank!=8))throw std::invalid_argument("block-feedback research parameters");
            fixed.maps.c=sr::blockFeedback(localRank,globalNodes);
        }
        if(alternating){if(localRank!=4)throw std::invalid_argument("alternating feedback starts with rank four");alternateFeedback=sr::blockFeedback(8,globalNodes);}
        if(blockEmission){
            if(localRank!=4&&localRank!=8)throw std::invalid_argument("block emission research parameters");
            const auto local=sr::blockFeedback(rankTwo?4:localRank);
            for(unsigned i=0;i!=64;++i)for(unsigned j=0;j!=16;++j)fixed.maps.a[i][j]=local[j][i];
        }
        if(carry==3){
            const auto shifted=sr::blockFeedback(4,false,16);
            for(unsigned r=0;r!=16;++r)fixed.maps.c[r]=shifted[4*(r%4)+r/4];
        }
        std::mt19937_64 rr(seed^0x91af823017ccb91aULL),mr(seed^0x14b2ec8399500213ULL),dr(seed^0xe685639417f293b5ULL),tr(seed^0xb2d174839ae61713ULL);
        std::uniform_int_distribution<unsigned> nz(1,65535);
        diagonal.resize(n);preparedDiagonal.resize(n);
        for(auto& d:diagonal)d=static_cast<std::uint16_t>(nz(dr));
        c2::wide16::prepare(diagonal.data(),preparedDiagonal.data(),n);
        route.resize(n);
        for(unsigned c=0;c!=256;++c){auto* p=route.data()+c*rows;std::iota(p,p+rows,std::uint16_t(0));
            for(auto j=rows-1;j;--j){std::uniform_int_distribution<std::size_t> d(0,j);std::swap(p[j],p[d(rr)]);}}
        if(sortPackets){
            if(!globalMds)throw std::invalid_argument("sorted packets require position-independent MDS moments");
            for(std::size_t i=0;i!=n;i+=64)std::sort(route.begin()+i,route.begin()+i+64);
        }
        routedDiagonal.resize(n);
        for(std::size_t i=0;i!=n;++i)routedDiagonal[i]=preparedDiagonal[(i/rows)*rows+route[i]];
        const auto& field=setupField();mixers.resize(n/64);
        if(cauchy)diagonalRefreshes.resize(n/64);else refreshes.resize(n/64);
        const auto h=cauchy&&retainRaw?sr::cauchyRefreshMatrix():sr::Matrix{};
        if(retainRaw)raw.reserve(n/64);
        for(std::size_t i=0;i!=n/64;++i){
            Raw r;bool nonzero;
            do{nonzero=false;for(auto& a:r.u){a=std::uint16_t(mr());nonzero|=a!=0;}}while(!nonzero);
            std::uint16_t dot;
            do{dot=0;for(unsigned j=0;j!=16;++j){r.v[j]=std::uint16_t(mr());dot^=field.mul(r.u[j],r.v[j]);}}while(dot==1);
            mixers[i]=c2::hdbw::prepareMixer(r.u.data(),r.v.data());
            if(cauchy){
                std::array<std::uint16_t,16>d1,d2;
                for(unsigned j=0;j!=16;++j){d1[j]=static_cast<std::uint16_t>(nz(tr));d2[j]=static_cast<std::uint16_t>(nz(tr));}
                diagonalRefreshes[i].prepare(d1.data(),d2.data());
                if(sharedDiagonal){
                    r.u=d2;
                    do{dot=0;for(unsigned j=0;j!=16;++j){r.v[j]=std::uint16_t(mr());dot^=field.mul(r.u[j],r.v[j]);}}while(dot==1);
                    mixers[i]=c2::hdbw::prepareMixer(r.u.data(),r.v.data());
                    if(rankTwo){
                        for(unsigned j=0;j!=16;++j)r.u2[j]=field.mul(d2[j],std::uint16_t(16+j));
                        std::uint16_t determinant;
                        do{
                            std::uint16_t aa=1,bb=0,cc=0,dd=1;
                            for(unsigned j=0;j!=16;++j){
                                r.v[j]=std::uint16_t(mr());r.v2[j]=std::uint16_t(mr());
                                aa^=field.mul(r.v[j],r.u[j]);bb^=field.mul(r.v[j],r.u2[j]);
                                cc^=field.mul(r.v2[j],r.u[j]);dd^=field.mul(r.v2[j],r.u2[j]);
                            }
                            determinant=field.mul(aa,dd)^field.mul(bb,cc);
                        }while(!determinant);
                        // The two stored rows are V, not rank-one (u,v).
                        // The joint 2x2 determinant was checked above;
                        // an extra rank-one rejection would change the law.
                        c2::wide16::prepare(r.v2.data(),mixers[i].u.data(),16);
                        c2::wide16::prepare(r.v.data(),mixers[i].v.data(),16);
                    }
                }
                if(retainRaw){
                    r.t=h;for(unsigned j=0;j!=16;++j)for(unsigned k=0;k!=16;++k)r.t[j][k]=field.mul(d2[j],field.mul(h[j][k],d1[k]));
                    if(carry>=2)r.carry=r.t;
                    else if(carry){r.carry=h;for(unsigned j=0;j!=16;++j)for(unsigned k=0;k!=16;++k)r.carry[j][k]=field.mul(d2[j],h[j][k]);}
                }
            }else{
                r.t=sr::Matrix(16,std::vector<std::uint16_t>(16));
                do{for(auto& row:r.t)for(auto& x:row)x=std::uint16_t(tr());}while(!field.invertible(r.t));
                refreshes[i].prepare(r.t);
            }
            if(retainRaw)raw.push_back(std::move(r));
        }
    }
    std::size_t planes()const{return 256*stride;}
    std::size_t scratch()const{return planes()+rows;}
};
// Implementation controls must keep their setup distribution identical to
// the construction they measure. This mapping is setup-only and constexpr.
template<unsigned Variant> inline constexpr unsigned ensemble=
    (Variant==42||Variant==44||Variant==45)?29:Variant>=38?38:Variant;
template<unsigned Variant> static Plan makePlan(std::size_t k,std::uint64_t seed,bool raw){
    constexpr unsigned R=ensemble<Variant>;
    Plan p(k,seed,raw,R>=2,R>=3,R>=7,R>=14,R>=16,R==18,R>=24,
        (R==26||R==32||R==38)?8:R>=27?4:16,
        (R>=28&&R<=31)||R==37,(R>=29&&R<=31)||R==37,
        R>=32,R==36?3:R==35?2:R==34?1:0,R>=37);
    if constexpr(Variant==45)p.feedbackPeriod=4;
#if defined(HD_SUMRANK_COMPACT_STEPS)
    if constexpr(Variant==42||Variant==45){
        p.compactSteps.resize(p.n/64);
        for(std::size_t i=0;i!=p.n/64;++i){
            p.compactSteps[i].v=p.mixers[i].v;
            p.compactSteps[i].refresh=p.diagonalRefreshes[i];
        }
    }
#endif
    if constexpr(Variant==43||Variant==44){
        p.stripedRoute.resize(p.n);
        for(std::size_t base=0;base!=p.n;base+=64)
            for(unsigned j=0;j!=16;++j)for(unsigned lane=0;lane!=4;++lane)
                p.stripedRoute[base+4*j+lane]=p.route[base+16*lane+j];
    }
#if defined(HD_SUMRANK_ROUTE_BYTES)||defined(HD_SUMRANK_ROUTE_UNITS8)
    if constexpr(Variant>=41&&Variant<=44){
        const auto& source=Variant>=43?p.stripedRoute:p.route;
        if constexpr(std::is_same_v<sr::PreparedOffset,sr::WordOffset>){
            // Larger regions retain the original record-index path.
            if(p.rows<=32768){
                p.offsetRoute.resize(p.n);
                for(std::size_t i=0;i!=p.n;++i)p.offsetRoute[i].value=std::uint16_t(2*source[i]);
            }
        }else{
            p.offsetRoute.resize(p.n);
            for(std::size_t i=0;i!=p.n;++i)p.offsetRoute[i].value=std::uint32_t(source[i])*sizeof(Block);
        }
    }
#endif
    return p;
}
SR_NOINLINE static void outer(const Plan& p,const Block* input,Block* planes){
    c2::hdbwouter::outerForward256<64,true,true>(p.k,p.stride,input,planes);
}
template<unsigned Refresh>
SR_NOINLINE static void selectedOuter(const Plan& p,const Block* input,Block* planes){
    if constexpr(Refresh==39)c2::hdsumrankouter::forward<64,true,4>(p.k,p.stride,input,planes,p.preparedDiagonal.data());
    else if constexpr(Refresh==40)c2::hdsumrankouter::forward<64,true,16>(p.k,p.stride,input,planes,p.preparedDiagonal.data());
    else if constexpr(Refresh==31)c2::hdsumrankouter::forward(p.k,p.stride,input,planes,p.paddedDiagonal.data(),p.rows+16);
    else if constexpr(Refresh==19)c2::hdsumrankouter::forward<32>(p.k,p.stride,input,planes,p.preparedDiagonal.data());
    else if constexpr(Refresh==20)c2::hdsumrankouter::forward<128>(p.k,p.stride,input,planes,p.preparedDiagonal.data());
    else if constexpr(Refresh==21)c2::hdsumrankouter::forward<64,false>(p.k,p.stride,input,planes,p.preparedDiagonal.data());
    else if constexpr(Refresh>=6)c2::hdsumrankouter::forward(p.k,p.stride,input,planes,p.preparedDiagonal.data());
    else outer(p,input,planes);
}

template<unsigned Rank>
SR_NOINLINE static void alternatingPacket(const Block* source,const std::uint16_t* route,Block* output,
                                          sr::State& state,const sr::Mixer& mixer,const sr::DiagonalRefresh& refresh){
    sr::mdsStriped<true,std::uint16_t,true,Rank>(source,route,output,state,mixer,refresh);
}

// Keep each fixed circuit in its own register-allocation unit. Combining
// both in one loop caused a measured code-generation regression; these
// direct calls include their state-transfer cost in the benchmark.
SR_NOINLINE static void alternatingRegion(const Plan& p,unsigned region,const Block* planes,Block* output,sr::State& state){
    const auto begin=region*p.rows;const auto* source=planes+region*p.stride;
    for(std::size_t j=0;j!=p.rows;j+=128){
        if(region+1<256){
            const auto* next=reinterpret_cast<const char*>(planes+(region+1)*p.stride+j);
            for(unsigned h=0;h!=32;++h)_mm_prefetch(next+64*h,_MM_HINT_T1);
        }
        const auto i=(begin+j)/64;
        alternatingPacket<4>(source,p.route.data()+begin+j,output+begin+j,state,p.mixers[i],p.diagonalRefreshes[i]);
        alternatingPacket<8>(source,p.route.data()+begin+j+64,output+begin+j+64,state,p.mixers[i+1],p.diagonalRefreshes[i+1]);
    }
}
// Same common circuit and region-prefetch schedule as mode42. Only the
// feedback phase changes. No new routing, temporary payload, or setup coin.
template<bool Routed,bool Stream>
SR_NOINLINE static void periodicInner(const Plan& p,const Block* SR_PAYLOAD_RESTRICT planes,Block* SR_PAYLOAD_RESTRICT output){
    sr::State state{};
    for(unsigned region=0;region!=256;++region){
        const auto begin=region*p.rows;
        const auto* source=planes+(Routed?region*p.stride:begin);
        for(std::size_t j=0;j!=p.rows;j+=64){
            if constexpr(Routed){
                if(region+1<256){
                    const auto* next=reinterpret_cast<const char*>(planes+(region+1)*p.stride+j);
                    for(unsigned h=0;h!=16;++h)_mm_prefetch(next+64*h,_MM_HINT_T1);
                }
            }
            const auto i=(begin+j)/64;
            sr::mdsStriped<Routed,std::uint16_t,true,0,16,0,false,Stream>(
                source+(Routed?0:j),Routed?p.route.data()+begin+j:nullptr,
                output+begin+j,state,p.stepMixer<45>(i),p.stepRefresh<45>(i),(i&3)==3);
        }
    }
    if constexpr(Stream)_mm_sfence();
}

template<unsigned Refresh=2>
// Full encoding consumes the complete input into disjoint scratch first.
// Only this internal boundary may promise disjoint source/output storage;
// packet helpers below still support exact input/output aliasing.
SR_NOINLINE static void inner(const Plan& p,const Block* SR_PAYLOAD_RESTRICT planes,Block* SR_PAYLOAD_RESTRICT output,Block* scratch){
    if constexpr(Refresh==45){
        if(reinterpret_cast<std::uintptr_t>(output)&63U)periodicInner<true,false>(p,planes,output);
        else periodicInner<true,true>(p,planes,output);
    }else{
#if defined(HD_SUMRANK_ROUTE_UNITS8)
    if constexpr(Refresh>=41)if(p.rows>32768){inner<ensemble<Refresh>>(p,planes,output,scratch);return;}
#endif
    if constexpr(Refresh==13||Refresh>=41){
        if(reinterpret_cast<std::uintptr_t>(output)&63U){inner<Refresh>=41?ensemble<Refresh>:10>(p,planes,output,scratch);return;}
    }
    sr::State state{};
    if constexpr(Refresh==28){
        if(p.rows%128==0){
            for(unsigned region=0;region!=256;++region)alternatingRegion(p,region,planes,output,state);
            return;
        }
    }
    for(unsigned region=0;region!=256;++region){
        const auto begin=region*p.rows;
        const Block* source;
        if constexpr(Refresh>=4)source=planes+region*p.stride;
        else{c2::hdbwrandom::stageRegion<true>(planes+region*p.stride,p.preparedDiagonal.data()+begin,scratch,p.rows);source=scratch;}
        for(std::size_t j=0;j!=p.rows;j+=64){
            if constexpr(Refresh>=5){
                // Stream the next region into L2 while working on this one.
                // It is only 128 KiB at K20, so both regions fit together.
                if(region+1<256){
                    const auto* next=reinterpret_cast<const char*>(planes+(region+1)*p.stride+j);
                    for(unsigned h=0;h!=16;++h)_mm_prefetch(next+64*h,_MM_HINT_T1);
                }
            }
            if constexpr(Refresh>=8){
                constexpr unsigned ahead=Refresh==9?128:256;
                if constexpr(Refresh<10)if(j+ahead<p.rows)for(unsigned h=0;h!=64;h+=8)
                    c2::c8Unroll(std::make_index_sequence<8>{},[&](auto u){
                        _mm_prefetch(reinterpret_cast<const char*>(source+p.route[begin+j+ahead+h+u]),_MM_HINT_T0);
                    });
            }else if(j+256<p.rows)for(unsigned h=0;h!=64;++h)
                _mm_prefetch(reinterpret_cast<const char*>(source+p.route[begin+j+256+h]),_MM_HINT_T0);
            const auto i=(begin+j)/64;
#if defined(HD_SUMRANK_ROUTE_BYTES)||defined(HD_SUMRANK_ROUTE_UNITS8)
            if constexpr(Refresh>=41)sr::mdsStriped<true,sr::PreparedOffset,true,(ensemble<Refresh>)==29?0:8,(ensemble<Refresh>)==29?16:4,0,(ensemble<Refresh>)==38,true,Refresh>=43>(source,p.offsetRoute.data()+begin+j,output+begin+j,state,p.stepMixer<Refresh>(i),p.stepRefresh<Refresh>(i),bool(i&1));
#else
            if constexpr(Refresh>=41)sr::mdsStriped<true,std::uint16_t,true,(ensemble<Refresh>)==29?0:8,(ensemble<Refresh>)==29?16:4,0,(ensemble<Refresh>)==38,true,Refresh>=43>(source,(Refresh>=43?p.stripedRoute:p.route).data()+begin+j,output+begin+j,state,p.stepMixer<Refresh>(i),p.stepRefresh<Refresh>(i),bool(i&1));
#endif
            else if constexpr(Refresh==23||Refresh>=25)sr::stripedByStep<Refresh,true,std::uint16_t>(source,p.route.data()+begin+j,output+begin+j,state,p.mixers[i],p.diagonalRefreshes[i],i);
            else if constexpr(Refresh>=16)sr::mdsPacket<true,std::uint16_t,Refresh>=17,Refresh==22,Refresh>=24>(source,p.route.data()+begin+j,output+begin+j,state,p.mixers[i],p.diagonalRefreshes[i]);
            else if constexpr(Refresh>=7)sr::fftPacket<true,std::uint16_t,Refresh==11?1:Refresh==12?2:0,Refresh==13,Refresh>=14,Refresh==15>(source,p.route.data()+begin+j,output+begin+j,state,p.mixers[i],p.diagonalRefreshes[i],p.fixed);
        else if constexpr(Refresh>=2)sr::packet<true,std::uint16_t,true,Refresh>=3,Refresh==4||Refresh==5>(source,p.route.data()+begin+j,output+begin+j,state,p.mixers[i],p.diagonalRefreshes[i],p.fixed,p.routedDiagonal.data()+begin+j);
            else sr::packet<true,std::uint16_t,Refresh!=0>(source,p.route.data()+begin+j,output+begin+j,state,p.mixers[i],p.refreshes[i],p.fixed);
        }
    }
    if constexpr(Refresh==13||Refresh>=41)_mm_sfence();
    }
}
template<unsigned Refresh=2>
SR_NOINLINE static void forward(const Plan& p,const Block* in,Block* out,Block* scratch){
    selectedOuter<Refresh>(p,in,scratch);inner<Refresh>(p,scratch,out,scratch+p.planes());
}
SR_NOINLINE static void routeOnly(const Plan& p,const Block* planes,Block* output,Block* scratch){
    for(unsigned c=0;c!=256;++c){
        const auto begin=c*p.rows;
        c2::hdbwrandom::stageRegion<true>(planes+c*p.stride,p.preparedDiagonal.data()+begin,scratch,p.rows);
        for(std::size_t j=0;j!=p.rows;j+=4)
            _mm512_storeu_si512(output+begin+j,c2::hdbw::detail::input4<true>(scratch,p.route.data()+begin+j,0));
    }
}
template<unsigned Refresh>
SR_NOINLINE static void plainInner(const Plan& p,const Block* input,Block* output){
    if constexpr(Refresh==45){
        if(reinterpret_cast<std::uintptr_t>(output)&63U)periodicInner<false,false>(p,input,output);
        else periodicInner<false,true>(p,input,output);
    }else{
    if constexpr(Refresh>=41){
        if(reinterpret_cast<std::uintptr_t>(output)&63U){plainInner<ensemble<Refresh>>(p,input,output);return;}
    }
    sr::State state{};
    for(std::size_t i=0;i!=p.n;i+=64){
        if constexpr(Refresh>=41)sr::mdsStriped<false,std::uint16_t,true,(ensemble<Refresh>)==29?0:8,(ensemble<Refresh>)==29?16:4,0,(ensemble<Refresh>)==38,true>(input+i,nullptr,output+i,state,p.stepMixer<Refresh>(i/64),p.stepRefresh<Refresh>(i/64),bool((i/64)&1));
        else if constexpr(Refresh==23||Refresh>=25)sr::stripedByStep<Refresh,false,std::uint16_t>(input+i,nullptr,output+i,state,p.mixers[i/64],p.diagonalRefreshes[i/64],i/64);
        else if constexpr(Refresh>=16)sr::mdsPacket<false,std::uint16_t,Refresh>=17,Refresh==22,Refresh>=24>(input+i,nullptr,output+i,state,p.mixers[i/64],p.diagonalRefreshes[i/64]);
        else if constexpr(Refresh>=7)sr::fftPacket<false,std::uint16_t,Refresh==11?1:Refresh==12?2:0,false,Refresh>=14,Refresh==15>(input+i,nullptr,output+i,state,p.mixers[i/64],p.diagonalRefreshes[i/64],p.fixed);
        else if constexpr(Refresh>=2)sr::packet<false,std::uint16_t,true,Refresh>=3>(input+i,nullptr,output+i,state,p.mixers[i/64],p.diagonalRefreshes[i/64],p.fixed);
        else sr::packet<false,std::uint16_t,Refresh!=0>(input+i,nullptr,output+i,state,p.mixers[i/64],p.refreshes[i/64],p.fixed);
    }
    if constexpr(Refresh>=41)_mm_sfence();
    }
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
static void checkScaledOuter(){
    std::mt19937_64 rng(813);const Block guard{0x93e13ULL,0xab812ULL};
    {
        // Boundary checks for both prepared address representations.
        std::vector<Block> input(65536);
        for(auto& x:input)x={rng(),rng()};
        const std::array<unsigned,4> small{0,1,32766,32767},large{32768,0,65534,65535};
        std::array<sr::WordOffset,4> words;
        std::array<sr::ByteOffset,4> bytes;
        for(unsigned j=0;j!=4;++j){words[j].value=std::uint16_t(2*small[j]);bytes[j].value=16*large[j];}
        std::array<Block,4> actual;
        _mm512_storeu_si512(actual.data(),sr::stripedInput4<true>(input.data(),words.data(),0));
        for(unsigned j=0;j!=4;++j)require(actual[j]==input[small[j]],"word route boundary");
        _mm512_storeu_si512(actual.data(),sr::stripedInput4<true>(input.data(),bytes.data(),0));
        for(unsigned j=0;j!=4;++j)require(actual[j]==input[large[j]],"byte route boundary");
    }
    for(std::size_t k:{512,1536,8192,24576}){
        const auto rows=k/128,stride=rows+4,count=256*stride;
        c2::BenchPayload input(k+4),actual(count+8);std::vector<Block> control(count);
        std::vector<std::uint16_t> d(256*rows);
        std::vector<c2::wide16::Coefficient> prepared(d.size());
        for(auto& x:d)x=std::uint16_t(1+rng()%65535);
        c2::wide16::prepare(d.data(),prepared.data(),d.size());
        std::vector<c2::wide16::Coefficient> padded(256*(rows+16));
        for(unsigned c=0;c!=256;++c)std::copy_n(prepared.data()+c*rows,rows,padded.data()+c*(rows+16));
        for(std::size_t i=0;i!=k+4;++i)input.p[i]={rng(),rng()};
        for(unsigned offset=0;offset!=4;++offset){
            std::fill_n(actual.p,count+8,guard);
            c2::hdbwouter::outerForward256(k,stride,input.p+offset,control.data());
            for(bool pad:{false,true}){
                std::fill_n(actual.p,count+8,guard);
                c2::hdsumrankouter::forward(k,stride,input.p+offset,actual.p+offset,pad?padded.data():prepared.data(),pad?rows+16:rows);
                for(unsigned c=0;c!=256;++c)for(std::size_t r=0;r!=stride;++r){
                    const auto expected=r<rows?scaled(control[c*stride+r],d[c*rows+r]):Block{};
                    require(actual.p[offset+c*stride+r]==expected,"fused outer scalar oracle and padding");
                }
                for(unsigned i=0;i!=offset;++i)require(actual.p[i]==guard,"fused outer prefix guard");
                for(auto i=count+offset;i!=count+8;++i)require(actual.p[i]==guard,"fused outer suffix guard");
            }
        }
    }
    std::cout<<"check,fused-outer-aligned-and-unaligned-scalar-oracle,pass\n";
}
static void scalarPacket(const Plan& p,const Block* x,Block* y,std::array<Block,16>& q,std::size_t step){
    const auto& m=p.fixed.maps;const auto& r=p.raw[step];
    const auto& feedbackMap=(step%p.feedbackPeriod==p.feedbackPeriod-1)&&!p.alternateFeedback.empty()?p.alternateFeedback:m.c;
    std::array<Block,16> feedback{},next=q;Block dot{},dot2{};
    for(unsigned i=0;i!=16;++i){add(dot,scaled(q[i],r.v[i]));add(dot2,scaled(q[i],r.v2[i]));}
    for(unsigned i=0;i!=16;++i){
        if(p.transformedCarry){next[i]={};for(unsigned j=0;j!=16;++j)add(next[i],scaled(q[j],r.carry[i][j]));}
        else{add(next[i],scaled(dot,r.u[i]));add(next[i],scaled(dot2,r.u2[i]));}
        for(unsigned j=0;j!=64;++j)add(feedback[i],scaled(x[j],feedbackMap[i][j]));
    }
    for(unsigned i=0;i!=64;++i){
        y[i]={};for(unsigned j=8*(i/8);j!=8*(i/8+1);++j)add(y[i],scaled(x[j],m.f[i][j]));
        for(unsigned j=0;j!=16;++j)add(y[i],scaled(q[j],m.a[i][j]));
    }
    for(unsigned i=0;i!=16;++i)for(unsigned j=0;j!=16;++j)add(next[i],scaled(feedback[j],r.t[i][j]));
    q=next;
}
static void preparePaddedDiagonal(Plan& p){
    // One cache line of padding rotates coefficient planes across L1 sets.
    // No setup distribution changes; the extra storage is never encoded.
    p.paddedDiagonal.resize(256*(p.rows+16));
    for(unsigned c=0;c!=256;++c)std::copy_n(p.preparedDiagonal.data()+c*p.rows,p.rows,p.paddedDiagonal.data()+c*(p.rows+16));
}
template<unsigned Refresh>static void check(){
    if constexpr(Refresh>=37){
        // Seed113 exercises V rows whose dot product equals one. That is
        // legal for the rank-two law and must not trigger rank-one checks.
        Plan setupRegression(1U<<20,113,false,true,true,true,true,true,false,true,4,true,true,true,0,true);
    }
    std::mt19937_64 rng(47);auto p=makePlan<Refresh>(8192,113,true);
    if constexpr(Refresh>=3){
        std::uint64_t h=0;
        for(const auto* m:{&p.fixed.maps.f,&p.fixed.maps.a,&p.fixed.maps.c})for(const auto& row:*m)
            for(auto x:row)h=(h^x)*0x100000001b3ULL;
        require(h==((ensemble<Refresh>)>=38?6742237189988151056ULL:(ensemble<Refresh>)==37?7680047191703725800ULL:(ensemble<Refresh>)==36?1225673016892679168ULL:(ensemble<Refresh>)>=33?5069833586872545664ULL:(ensemble<Refresh>)==32?9075152494480514616ULL:(ensemble<Refresh>)>=29?12097737318272255584ULL:(ensemble<Refresh>)==26?5242589667199767848ULL:(ensemble<Refresh>)>=27?16990099864024496792ULL:(ensemble<Refresh>)>=16?11358755321456897600ULL:(ensemble<Refresh>)>=14?1274355414671835662ULL:(ensemble<Refresh>)>=7?17847517350586406002ULL:14806900258954508744ULL),"independent Python exact-map digest");
    }
    for(unsigned i=0;i!=20000;++i){auto a=std::uint16_t(rng()),b=std::uint16_t(rng());require(setupField().mul(a,b)==oracleMul(a,b),"independent setup field");}
    for(unsigned trial=0;trial!=176;++trial){
        std::array<Block,64> x{},y{},want{};std::array<Block,16> q{},expectedState{},actualState{};
        if(trial<160){const unsigned v=trial%80,component=trial/80;const Block one{0x0001000100010001ULL<<(8*component),0x0001000100010001ULL<<(8*component)};
            if(v<64)x[v]=one;else q[v-64]=one;
        }else{for(auto& a:x)a={rng(),rng()};for(auto& a:q)a={rng(),rng()};}
        expectedState=q;scalarPacket(p,x.data(),want.data(),expectedState,trial);
        sr::State state;c2::hdbw::loadState(state,q.data());
        if constexpr(Refresh==23||Refresh>=25)sr::stripedByStep<Refresh,false,std::uint16_t>(x.data(),nullptr,y.data(),state,p.mixers[trial],p.diagonalRefreshes[trial],trial);
        else if constexpr(Refresh>=16)sr::mdsPacket<false,std::uint16_t,Refresh>=17,Refresh==22,Refresh>=24>(x.data(),nullptr,y.data(),state,p.mixers[trial],p.diagonalRefreshes[trial]);
        else if constexpr(Refresh>=7)sr::fftPacket<false,std::uint16_t,Refresh==11?1:Refresh==12?2:0,false,Refresh>=14,Refresh==15>(x.data(),nullptr,y.data(),state,p.mixers[trial],p.diagonalRefreshes[trial],p.fixed);
        else if constexpr(Refresh>=2)sr::packet<false,std::uint16_t,true,Refresh>=3>(x.data(),nullptr,y.data(),state,p.mixers[trial],p.diagonalRefreshes[trial],p.fixed);
        else sr::packet<false,std::uint16_t>(x.data(),nullptr,y.data(),state,p.mixers[trial],p.refreshes[trial],p.fixed);
        c2::hdbw::storeState(actualState.data(),state);
        require(y==want,"sumrank packet dense output oracle");require(actualState==expectedState,"sumrank packet dense state oracle");
        c2::hdbw::loadState(state,q.data());
        if constexpr(Refresh==23||Refresh>=25)sr::stripedByStep<Refresh,false,std::uint16_t>(x.data(),nullptr,x.data(),state,p.mixers[trial],p.diagonalRefreshes[trial],trial);
        else if constexpr(Refresh>=16)sr::mdsPacket<false,std::uint16_t,Refresh>=17,Refresh==22,Refresh>=24>(x.data(),nullptr,x.data(),state,p.mixers[trial],p.diagonalRefreshes[trial]);
        else if constexpr(Refresh>=7)sr::fftPacket<false,std::uint16_t,Refresh==11?1:Refresh==12?2:0,false,Refresh>=14,Refresh==15>(x.data(),nullptr,x.data(),state,p.mixers[trial],p.diagonalRefreshes[trial],p.fixed);
        else if constexpr(Refresh>=2)sr::packet<false,std::uint16_t,true,Refresh>=3>(x.data(),nullptr,x.data(),state,p.mixers[trial],p.diagonalRefreshes[trial],p.fixed);
        else sr::packet<false,std::uint16_t>(x.data(),nullptr,x.data(),state,p.mixers[trial],p.refreshes[trial],p.fixed);
        require(x==want,"sumrank packet alias");
        if(trial>=160){
            const Block guard{0x391234a89cULL,0x15e921bc7aULL};
            std::array<Block,72> source,result;source.fill(guard);result.fill(guard);
            std::array<std::uint32_t,64> route;
            std::iota(route.begin(),route.end(),0);std::shuffle(route.begin(),route.end(),rng);
            // Restore the original packet after the alias test.
            for(auto& b:x)b={rng(),rng()};expectedState=q;
            scalarPacket(p,x.data(),want.data(),expectedState,trial);
            for(unsigned offset=0;offset!=4;++offset){
                for(unsigned i=0;i!=64;++i)source[offset+route[i]]=x[i];
                c2::hdbw::loadState(state,q.data());
                if constexpr(Refresh==23||Refresh>=25)sr::stripedByStep<Refresh,true,std::uint32_t>(source.data()+offset,route.data(),result.data()+offset,state,p.mixers[trial],p.diagonalRefreshes[trial],trial);
                else if constexpr(Refresh>=16)sr::mdsPacket<true,std::uint32_t,Refresh>=17,Refresh==22,Refresh>=24>(source.data()+offset,route.data(),result.data()+offset,state,p.mixers[trial],p.diagonalRefreshes[trial]);
                else if constexpr(Refresh>=7)sr::fftPacket<true,std::uint32_t,Refresh==11?1:Refresh==12?2:0,false,Refresh>=14,Refresh==15>(source.data()+offset,route.data(),result.data()+offset,state,p.mixers[trial],p.diagonalRefreshes[trial],p.fixed);
        else if constexpr(Refresh>=2)sr::packet<true,std::uint32_t,true,Refresh>=3>(source.data()+offset,route.data(),result.data()+offset,state,p.mixers[trial],p.diagonalRefreshes[trial],p.fixed);
                else sr::packet<true,std::uint32_t>(source.data()+offset,route.data(),result.data()+offset,state,p.mixers[trial],p.refreshes[trial],p.fixed);
                c2::hdbw::storeState(actualState.data(),state);
                require(std::equal(want.begin(),want.end(),result.begin()+offset)&&actualState==expectedState,"routed offsets and nonzero-state equality");
                for(unsigned i=0;i<offset;++i)require(result[i]==guard,"packet prefix guard");
                for(unsigned i=offset+64;i!=72;++i)require(result[i]==guard,"packet suffix guard");
                result.fill(guard);
            }
        }
    }
    std::cout<<"check,sumrank-all80-bases-both-components-and-random-state,pass\n";
    for(std::size_t k:Refresh>=28?std::vector<std::size_t>{8192,16384,24576,32768}:std::vector<std::size_t>{8192,24576}){
        auto full=makePlan<Refresh>(k,919,true);std::vector<Block> x(k),planes(full.scratch()),out(full.n),want(full.n);for(auto& a:x)a={rng(),rng()};
        if constexpr(Refresh==31)preparePaddedDiagonal(full);
        outer(full,x.data(),planes.data());std::array<Block,16> q{};std::array<Block,64> packet{};
        for(std::size_t i=0;i!=full.n;i+=64){
            const auto c=i/full.rows;
            for(unsigned j=0;j!=64;++j){const auto row=full.route[i+j];packet[j]=scaled(planes[c*full.stride+row],full.diagonal[c*full.rows+row]);}
            scalarPacket(full,packet.data(),want.data()+i,q,i/64);
        }
        forward<Refresh>(full,x.data(),out.data(),planes.data());require(out==want,"full dense oracle");
        outer(full,x.data(),planes.data());
        routeOnly(full,planes.data(),out.data(),planes.data()+full.planes());
        std::vector<Block> plain(full.n);plainInner<Refresh>(full,out.data(),plain.data());require(plain==want,"plain arithmetic control");
        std::vector<Block> alias(full.n);std::copy(x.begin(),x.end(),alias.begin());forward<Refresh>(full,alias.data(),alias.data(),planes.data());require(alias==want,"full alias");
        if constexpr(Refresh==13||Refresh>=41){
            c2::BenchPayload aligned(full.n+8);const Block guard{18373,38183};
            for(unsigned offset=0;offset!=4;++offset){
                std::fill_n(aligned.p,full.n+8,guard);
                forward<Refresh>(full,x.data(),aligned.p+offset,planes.data());
                require(std::equal(want.begin(),want.end(),aligned.p+offset),"streaming full output alignment");
                for(unsigned j=0;j!=offset;++j)require(aligned.p[j]==guard,"stream prefix guard");
                for(auto j=full.n+offset;j!=full.n+8;++j)require(aligned.p[j]==guard,"stream suffix guard");
            }
            std::copy(x.begin(),x.end(),aligned.p);
            forward<Refresh>(full,aligned.p,aligned.p,planes.data());
            require(std::equal(want.begin(),want.end(),aligned.p),"streaming full alias");
        }
        std::cout<<"check,K,"<<k<<",dense-reference-route-and-full-alias,pass\n";
    }
}
static void produce(Block* p,std::size_t n){std::uint64_t x=913,y=113;
    for(std::size_t i=0;i!=n;++i){x^=x<<13;x^=x>>7;x^=x<<17;y+=0x9e3779b97f4a7c15ULL;_mm_stream_si128(reinterpret_cast<__m128i*>(p+i),_mm_set_epi64x(y,x));}_mm_sfence();}
static std::uint64_t hash(const Block* p,std::size_t n){std::uint64_t h=0;for(std::size_t i=0;i!=n;++i){h=(h^p[i].lo)*0x100000001b3ULL;h=(h^p[i].hi)*0x100000001b3ULL;}return h;}
static void stats(const char* label,std::vector<double> v){std::sort(v.begin(),v.end());std::cout<<label<<','<<v[v.size()/2]<<','<<v[v.size()/10]<<','<<v[9*v.size()/10]<<'\n';}
// Untimed target-size validation. Keep it separate from both the fast
// small-size suite and benchmark clocks. The reference performs dense
// packet arithmetic with an independent byte-polynomial field oracle.
template<unsigned Variant> static void checkLarge(){
    static_assert(Variant==41||Variant==42||Variant==45);
    auto p=makePlan<Variant>(std::size_t(1)<<20,113,true);
    c2::BenchPayload input(p.k),scratch(p.scratch()),actual(p.n+8),expected(p.n);
    produce(input.p,p.k);outer(p,input.p,scratch.p);
    std::array<Block,16> state{};std::array<Block,64> packet{};
    for(unsigned region=0;region!=256;++region){
        const auto begin=region*p.rows;
        for(std::size_t j=0;j!=p.rows;j+=64){
            const auto first=begin+j;
            for(unsigned h=0;h!=64;++h){
                const auto row=p.route[first+h];require(row<p.rows,"K20 route index");
                packet[h]=scaled(scratch.p[region*p.stride+row],p.diagonal[begin+row]);
            }
            scalarPacket(p,packet.data(),expected.p+first,state,first/64);
        }
        if((region&63)==63)std::cout<<"large-check,variant,"<<Variant
            <<",reference_regions,"<<region+1<<"/256\n"<<std::flush;
    }
    const Block guard{0x631943ULL,0x91837ULL};
    for(unsigned offset:{0U,1U}){
        std::fill_n(actual.p,p.n+8,guard);
        forward<Variant>(p,input.p,actual.p+offset,scratch.p);
        require(std::equal(expected.p,expected.p+p.n,actual.p+offset),"K20 dense oracle");
        for(unsigned j=0;j!=offset;++j)require(actual.p[j]==guard,"K20 prefix guard");
        for(auto j=p.n+offset;j!=p.n+8;++j)require(actual.p[j]==guard,"K20 suffix guard");
    }
    std::copy_n(input.p,p.k,actual.p);
    forward<Variant>(p,actual.p,actual.p,scratch.p);
    require(std::equal(expected.p,expected.p+p.n,actual.p),"K20 full alias");
    constexpr std::uint64_t golden=Variant==41?7812899659201065080ULL:
        Variant==42?16994121132285933612ULL:7582624837270167152ULL;
    require(hash(expected.p,p.n)==golden,"K20 frozen checksum");
    std::cout<<"large-check,K,"<<p.k<<",variant,"<<Variant
        <<",dense-reference-alignment-and-alias,pass,checksum,"<<golden<<'\n'<<std::flush;
}

template<unsigned Refresh=2>static void bench(unsigned calls,const std::string& mode,unsigned logK){
    auto p=makePlan<Refresh>(std::size_t(1)<<logK,113,false);c2::BenchPayload in(p.k),out(p.n),scratch(p.scratch());produce(in.p,p.k);outer(p,in.p,scratch.p);
    if constexpr(Refresh==31)preparePaddedDiagonal(p);
    if(mode=="plain-inner"){routeOnly(p,scratch.p,out.p,scratch.p+p.planes());std::memcpy(scratch.p,out.p,p.n*sizeof(Block));}
    else if(mode!="route-only")selectedOuter<Refresh>(p,in.p,scratch.p);
    auto run=[&]{if(mode=="inner")inner<Refresh>(p,scratch.p,out.p,scratch.p+p.planes());else if(mode=="outer")selectedOuter<Refresh>(p,in.p,scratch.p);else if(mode=="route-only")routeOnly(p,scratch.p,out.p,scratch.p+p.planes());else if(mode=="plain-inner")plainInner<Refresh>(p,scratch.p,out.p);else forward<Refresh>(p,in.p,out.p,scratch.p);};
    for(unsigned j=0;j!=3;++j)run();std::vector<double> times(calls);std::uint64_t consume=0;
    const Block* result=mode=="outer"?scratch.p:out.p;
    for(unsigned j=0;j!=calls;++j){produce(in.p,p.k);const auto a=Clock::now();run();times[j]=std::chrono::duration<double,std::micro>(Clock::now()-a).count();std::atomic_signal_fence(std::memory_order_seq_cst);consume^=result[(j*104729U)%p.n].lo;}
    std::cout<<"profile,sumrank-s16,K,"<<p.k<<",N,"<<p.n<<",refresh,"<<((ensemble<Refresh>)>=37?"rank-two-joint-diagonal-Cauchy":(ensemble<Refresh>)>=35?"common-diagonal-Cauchy":(ensemble<Refresh>)==34?"carried-diagonal-Cauchy":(ensemble<Refresh>)>=24?"joint-diagonal-Cauchy":(ensemble<Refresh>)>=2?"diagonal-Cauchy":(ensemble<Refresh>)==1?"uniform-GL16":"identity-UNCERTIFIED-CONTROL")<<",inner_family,"<<((ensemble<Refresh>)>=32?"block-RS":(ensemble<Refresh>)>=16?"systematic-RS":"layered-LRS")<<",sorted_packets,"<<((ensemble<Refresh>)==18)<<",tower_basis,"<<((ensemble<Refresh>)>=3&&(ensemble<Refresh>)<16)<<",direct_route,"<<((ensemble<Refresh>)>=4)<<",stream_ahead,"<<((ensemble<Refresh>)>=5)<<",fused_outer,"<<((ensemble<Refresh>)>=6)<<",fft_weights,"<<((ensemble<Refresh>)>=7&&(ensemble<Refresh>)<16)<<",variant,"<<Refresh<<",calls,"<<calls<<",checksum,"<<hash(result,p.n)<<",consume,"<<consume<<"\nname,median_us,p10_us,p90_us\n"<<std::fixed<<std::setprecision(3);stats(mode.c_str(),times);
}
int main(int argc,char** argv){try{
    if(!c2::gfni512Available())return 77;
    if(argc==2&&std::string(argv[1])=="--large-check"){checkLarge<41>();checkLarge<42>();checkLarge<45>();return 0;}
    if(argc==2&&std::string(argv[1])=="--check"){checkScaledOuter();check<1>();check<2>();check<3>();check<4>();check<5>();check<6>();check<7>();check<8>();check<9>();check<10>();check<11>();check<12>();check<13>();check<14>();check<15>();check<16>();check<17>();check<18>();check<19>();check<20>();check<21>();check<22>();check<23>();check<24>();check<25>();check<26>();check<27>();check<28>();check<29>();check<30>();check<31>();check<32>();check<33>();check<34>();check<35>();check<36>();check<37>();check<38>();check<39>();check<40>();check<41>();check<42>();check<43>();check<44>();check<45>();return 0;}
    if(argc<2||std::string(argv[1])!="--bench")throw std::invalid_argument("--check | --large-check | --bench [calls] [full|inner|outer|plain-inner|route-only] [logK] [bucket4448-fused-nt|bucket48-striped-route|bucket8-rank2-striped-route|bucket48-fused-nt|bucket8-rank2-nt|bucket8-rank2-prefetch|bucket8-rank2-prefetch-l2|bucket8-rank2|bucket48-rank2|bucket4-scatter|bucket4-common|bucket4-carry|bucket4-both|bucket8-both|bucket48-padded|bucket48-affine|bucket48-fused|bucket48|bucket8|bucket4|shared|shared-striped|striped|bounded|mds-sorted|outer32|outer128|outer-cached|mds-affine|mds|affine|systematic|nt|paired-load|gather|prefetch|near|l2only|fft|fused|stream|direct|tower|cauchy|refresh|identity]");
    unsigned calls=argc>2?std::stoul(argv[2]):31,logK=argc>4?std::stoul(argv[4]):20;std::string mode=argc>3?argv[3]:"full",refresh=argc>5?argv[5]:"cauchy";
    if(refresh=="bucket4448-fused-nt"){
        if(argc>6||!calls||calls>100000||logK<13||logK>23||(mode!="full"&&mode!="inner"&&mode!="outer"&&mode!="plain-inner"&&mode!="route-only"))throw std::invalid_argument("arguments");
        bench<45>(calls,mode,logK);return 0;
    }
    if(refresh=="bucket48-striped-route"||refresh=="bucket8-rank2-striped-route"||refresh=="bucket48-fused-nt"||refresh=="bucket8-rank2-nt"||refresh=="bucket8-rank2-prefetch"||refresh=="bucket8-rank2-prefetch-l2"||refresh=="bucket8-rank2"||refresh=="bucket48-rank2"||refresh=="bucket4-scatter"||refresh=="bucket4-common"||refresh=="bucket4-carry"||refresh=="bucket4-both"||refresh=="bucket8-both"||refresh=="bucket48-padded"||refresh=="bucket48-affine"||refresh=="bucket48-fused"||refresh=="bucket48"||refresh=="bucket8"||refresh=="bucket4"||refresh=="shared"||refresh=="shared-striped"||refresh=="striped"||refresh=="bounded"||refresh=="mds"||refresh=="mds-affine"||refresh=="mds-sorted"||refresh=="outer32"||refresh=="outer128"||refresh=="outer-cached"){
        if(argc>6||!calls||calls>100000||logK<13||logK>23||(mode!="full"&&mode!="inner"&&mode!="outer"&&mode!="plain-inner"&&mode!="route-only"))throw std::invalid_argument("arguments");
        if(refresh=="bucket48-striped-route")bench<44>(calls,mode,logK);else if(refresh=="bucket8-rank2-striped-route")bench<43>(calls,mode,logK);else if(refresh=="bucket48-fused-nt")bench<42>(calls,mode,logK);else if(refresh=="bucket8-rank2-nt")bench<41>(calls,mode,logK);else if(refresh=="bucket8-rank2-prefetch")bench<39>(calls,mode,logK);else if(refresh=="bucket8-rank2-prefetch-l2")bench<40>(calls,mode,logK);else if(refresh=="bucket8-rank2")bench<38>(calls,mode,logK);else if(refresh=="bucket48-rank2")bench<37>(calls,mode,logK);else if(refresh=="bucket4-scatter")bench<36>(calls,mode,logK);else if(refresh=="bucket4-common")bench<35>(calls,mode,logK);else if(refresh=="bucket4-carry")bench<34>(calls,mode,logK);else if(refresh=="bucket4-both")bench<33>(calls,mode,logK);else if(refresh=="bucket8-both")bench<32>(calls,mode,logK);else if(refresh=="bucket48-padded")bench<31>(calls,mode,logK);else if(refresh=="bucket48-affine")bench<30>(calls,mode,logK);else if(refresh=="bucket48-fused")bench<29>(calls,mode,logK);else if(refresh=="bucket48")bench<28>(calls,mode,logK);else if(refresh=="bucket8")bench<26>(calls,mode,logK);else if(refresh=="bucket4")bench<27>(calls,mode,logK);else if(refresh=="shared")bench<24>(calls,mode,logK);else if(refresh=="shared-striped")bench<25>(calls,mode,logK);else if(refresh=="striped")bench<23>(calls,mode,logK);else if(refresh=="bounded")bench<22>(calls,mode,logK);else if(refresh=="mds")bench<16>(calls,mode,logK);else if(refresh=="mds-sorted")bench<18>(calls,mode,logK);else if(refresh=="outer32")bench<19>(calls,mode,logK);else if(refresh=="outer128")bench<20>(calls,mode,logK);else if(refresh=="outer-cached")bench<21>(calls,mode,logK);else bench<17>(calls,mode,logK);return 0;
    }
    if(argc>6||!calls||calls>100000||logK<13||logK>23||(mode!="full"&&mode!="inner"&&mode!="outer"&&mode!="plain-inner"&&mode!="route-only")||(refresh!="refresh"&&refresh!="identity"&&refresh!="cauchy"&&refresh!="tower"&&refresh!="direct"&&refresh!="stream"&&refresh!="fused"&&refresh!="fft"&&refresh!="prefetch"&&refresh!="near"&&refresh!="l2only"&&refresh!="paired-load"&&refresh!="gather"&&refresh!="nt"&&refresh!="systematic"&&refresh!="affine"))throw std::invalid_argument("arguments");
    if(refresh=="identity")bench<0>(calls,mode,logK);else if(refresh=="refresh")bench<1>(calls,mode,logK);else if(refresh=="tower")bench<3>(calls,mode,logK);else if(refresh=="direct")bench<4>(calls,mode,logK);else if(refresh=="stream")bench<5>(calls,mode,logK);else if(refresh=="fused")bench<6>(calls,mode,logK);else if(refresh=="fft")bench<7>(calls,mode,logK);else if(refresh=="prefetch")bench<8>(calls,mode,logK);else if(refresh=="near")bench<9>(calls,mode,logK);else if(refresh=="l2only")bench<10>(calls,mode,logK);else if(refresh=="paired-load")bench<11>(calls,mode,logK);else if(refresh=="gather")bench<12>(calls,mode,logK);else if(refresh=="nt")bench<13>(calls,mode,logK);else if(refresh=="systematic")bench<14>(calls,mode,logK);else if(refresh=="affine")bench<15>(calls,mode,logK);else bench<2>(calls,mode,logK);return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
