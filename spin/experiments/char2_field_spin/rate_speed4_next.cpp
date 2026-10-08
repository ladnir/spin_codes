#include "RateSpeed4Code.h"
#include "RateSpeed4RouteNext.h"
#include "RateSpeed4OuterNext.h"
#include "RateSpeed4Oracle.h"
#include "BenchPayload.h"
#include <chrono>
#include <iostream>
#include <iomanip>
#include <string>
namespace c2=spin::research::char2;
using Clock=std::chrono::steady_clock;
static void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
static double us(Clock::time_point a,Clock::time_point b){return std::chrono::duration<double,std::micro>(b-a).count();}
static constexpr unsigned candidates=29;
static const char* names[candidates]={"cached-eight-page","cached-eight-tight","cached-eight-pad4","cached-eight-pad16","cached-eight-pad64","cached-eight-page-plus4","route-prefetch16","route-prefetch64","route-hardware-gather","route-diagonal","route-diagonal-prefetch16","outer-fused-multiply","outer-fused-affine","outer-fused-affine-prefetch","inner-gauge","inner-pairs","inner-precomputed-pairs","gauge-route-diagonal","combined-route-fused","combined-route-fused-pad4","combined-gauge-route-fused","route-batched-diagonal","combined-batched-pad4","combined-pairs-pad4","combined-gauge-pad4","combined-direct-temp","combined-batched-direct-temp","combined-pipeline16","combined-batched-pipeline16"};
static std::size_t strideFor(unsigned i,std::size_t rows){switch(i){case 1:return rows;case 2:case 19:case 22:case 23:case 24:case 25:case 26:case 27:case 28:return rows+4;case 3:return rows+16;case 4:return rows+64;case 5:return ((rows+255)/256)*256+4;default:return ((rows+255)/256)*256;}}
template<unsigned D>static void inner(unsigned i,c2::RateSpeed4SchedulePlan<D>& p,const c2::Block* in,c2::Block* planes,c2::Block* temp){
    p.planeStride=strideFor(i,p.groups);
    switch(i){
    case 6:c2::rs4RegionRouteNext<D,16>(p,in,planes,temp);break;
    case 7:c2::rs4RegionRouteNext<D,64>(p,in,planes,temp);break;
    case 8:c2::rs4RegionRouteNext<D,0,true>(p,in,planes,temp);break;
    case 9:c2::rs4RegionRouteNext<D,0,false,true>(p,in,planes,temp);break;
    case 10:c2::rs4RegionRouteNext<D,16,false,true>(p,in,planes,temp);break;
    case 14:c2::rs4RegionRouteNext<D,0,false,false,1>(p,in,planes,temp);break;
    case 15:c2::rs4RegionRouteNext<D,0,false,false,2>(p,in,planes,temp);break;
    case 16:c2::rs4RegionRouteNext<D,0,false,false,3>(p,in,planes,temp);break;
    case 17:case 20:case 24:c2::rs4RegionRouteNext<D,0,false,true,1>(p,in,planes,temp);break;
    case 18:case 19:c2::rs4RegionRouteNext<D,0,false,true>(p,in,planes,temp);break;
    case 21:case 22:c2::rs4RegionRouteNext<D,0,false,true,0,true>(p,in,planes,temp);break;
    case 23:c2::rs4RegionRouteNext<D,0,false,true,3>(p,in,planes,temp);break;
    case 25:c2::rs4RegionRouteNext<D,0,false,true,0,false,true>(p,in,planes,temp);break;
    case 26:c2::rs4RegionRouteNext<D,0,false,true,0,true,true>(p,in,planes,temp);break;
    case 27:c2::rs4RegionRouteNext<D,0,false,true,4,false,true>(p,in,planes,temp);break;
    case 28:c2::rs4RegionRouteNext<D,0,false,true,4,true,true>(p,in,planes,temp);break;
    default:c2::rs4RegionInner<true,D,false,true,false,true>(p,in,planes,temp);break;
    }
}
template<unsigned D,bool Stream>static void outerAligned(unsigned i,c2::RateSpeed4Plan<D>& p,const c2::Block* planes,c2::Block* out){
    using V=c2::Rs4OuterNextVariant;
    switch(i){
    case 11:case 18:case 19:case 20:case 22:case 23:case 24:case 25:case 26:case 27:case 28:c2::rs4OuterNext<D,Stream,V::FusedMultiply>(p,planes,out);break;
    case 12:c2::rs4OuterNext<D,Stream,V::FusedAffine>(p,planes,out);break;
    case 13:c2::rs4OuterNext<D,Stream,V::FusedAffinePrefetch>(p,planes,out);break;
    default:c2::rs4Outer<D,Stream>(p,planes,out);break;
    }
}
template<unsigned D>static void outer(unsigned i,c2::RateSpeed4Plan<D>& p,const c2::Block* planes,c2::Block* out){
    if((reinterpret_cast<std::uintptr_t>(out)&63)==0)outerAligned<D,true>(i,p,planes,out);else outerAligned<D,false>(i,p,planes,out);
}
template<unsigned D>static void execute(unsigned i,c2::RateSpeed4SchedulePlan<D>& p,const c2::Block* in,c2::Block* out,c2::Block* planes,c2::Block* temp){
    if(i==27){p.planeStride=strideFor(i,p.groups);c2::rs4TransposeRouted16<D>(p,in,out,planes,temp);}
    else{inner(i,p,in,planes,temp);outer(i,p,planes,out);}
}
static void checkLocalHistory(){
    std::mt19937_64 random(717);
    for(unsigned trial=0;trial!=256;++trial){
        alignas(64)c2::Block input[48],initial[4],expected[64],output[64],finalExpected[4],finalActual[4];
        for(auto& v:input)v={random(),random()};for(auto& v:initial)v={random(),random()};
        std::uint8_t ds[64],gauge[128],pairs[96];for(auto& d:ds)d=static_cast<std::uint8_t>(1+random()%255);
        const auto p0=c2::rs4PrepareGauge(ds,64,gauge);c2::rs4PreparePairs(ds,64,pairs);
        auto h=_mm512_load_si512(initial);
        c2::rs4InnerTransposeEight<true,true,true>(input+24,expected+32,h,ds+32);c2::rs4InnerTransposeEight<true,true,true>(input,expected,h,ds);_mm512_store_si512(finalExpected,h);
        for(unsigned variant=0;variant!=4;++variant){
            h=_mm512_load_si512(initial);
            if(variant==0){c2::rs4InnerTransposeEightGauge(input+24,output+32,h,gauge+64);c2::rs4InnerTransposeEightGauge(input,output,h,gauge);h=_mm512_gf2p8mul_epi8(h,c2::sfBroadcastFour(p0.data()));}
            else if(variant==1){c2::rs4InnerTransposeEightPairs<true,true,false>(input+24,output+32,h,ds+32);c2::rs4InnerTransposeEightPairs<true,true,false>(input,output,h,ds);}
            else if(variant==2){c2::rs4InnerTransposeEightPairs<true,true,true>(input+24,output+32,h,pairs+48);c2::rs4InnerTransposeEightPairs<true,true,true>(input,output,h,pairs);}
            else c2::rs4InnerTransposeSixteenPipeline(input,output,h,ds);
            _mm512_store_si512(finalActual,h);require(std::equal(expected,expected+64,output),"next: arbitrary-history output");require(std::equal(finalExpected,finalExpected+4,finalActual),"next: arbitrary-history state");
        }
    }
    std::cout<<"next gauge/pairs/pipeline arbitrary-history checks pass\n";
}
template<unsigned D>static void check(bool large){
    std::mt19937_64 random(113);const c2::Block guard{91,717};
    const std::vector<std::size_t> sizes=large?std::vector<std::size_t>{65536,262144,1048576}:std::vector<std::size_t>{1,D-1,D,D+1,129,130,131,4*D+1,4*D+2,4*D+3,64*D,64*D+1};
    for(auto k:sizes){
        c2::RateSpeed4SchedulePlan<D> p(k,991);std::size_t maxStride=0;for(unsigned i=0;i!=candidates;++i)maxStride=std::max(maxStride,strideFor(i,p.groups));
        c2::BenchPayload input(p.outputN+8),out(k+8),planes(64*maxStride+8),temp(p.groups+8);
        std::vector<c2::Block> original(p.outputN),expected(k);for(auto& v:original)v={random(),random()};
        if(large){c2::rs4TransposeCached(p,original.data(),out.p,planes.p,temp.p);std::copy(out.p,out.p+k,expected.begin());}else expected=c2::rs4OracleTranspose(p,original.data());
        for(unsigned i=0;i!=candidates;++i)for(unsigned offset:{0,1,2,3}){
            for(auto* b:{&input,&out,&planes,&temp})std::fill(b->p,b->p+b->count,guard);
            std::memcpy(input.p+offset,original.data(),original.size()*sizeof(c2::Block));
            execute(i,p,input.p+offset,out.p+offset,planes.p+offset,temp.p+offset);
            require(std::equal(expected.begin(),expected.end(),out.p+offset),"next: reference equality");
            require(std::equal(original.begin(),original.end(),input.p+offset),"next: input preservation");
            execute(i,p,input.p+offset,input.p+offset,planes.p+offset,temp.p+offset);
            require(std::equal(expected.begin(),expected.end(),input.p+offset),"next: alias prefix");
            require(std::equal(original.begin()+k,original.end(),input.p+offset+k),"next: alias suffix");
            for(auto* b:{&input,&out,&planes,&temp}){for(unsigned j=0;j!=offset;++j)require(b->p[j]==guard,"next: prefix guard");for(auto j=b->count-8+offset;j!=b->count;++j)require(b->p[j]==guard,"next: suffix guard");}
        }
        p.foldedInverseLocalRoute.clear();
        for(unsigned i=0;i!=candidates;++i){execute(i,p,original.data(),out.p,planes.p,temp.p);require(std::equal(expected.begin(),expected.end(),out.p),"next: wide inverse route");}
        std::cout<<"next D="<<D<<" K="<<k<<" exact/aliases/guards/wide-route pass\n";
    }
}
static std::uint64_t checksum(const c2::Block* data,std::size_t count){std::uint64_t h=0;for(std::size_t i=0;i!=count;++i){h=(h^data[i].lo)*0x100000001b3ULL;h=(h^data[i].hi)*0x100000001b3ULL;}return h;}
static void stats(const std::string& name,std::vector<double> values){std::sort(values.begin(),values.end());std::cout<<name<<','<<std::fixed<<std::setprecision(3)<<values[values.size()/2]<<','<<values[values.size()/10]<<','<<values[9*values.size()/10]<<'\n';}
static void bench(std::size_t k,unsigned calls,bool phases,bool fresh){
    c2::RateSpeed4SchedulePlan<40> p(k,113);std::size_t maxStride=0;for(unsigned i=0;i!=candidates;++i)maxStride=std::max(maxStride,strideFor(i,p.groups));
    c2::BenchPayload input(p.outputN),planes(64*maxStride),temp(p.groups);
    std::vector<c2::Block> original(p.outputN);std::mt19937_64 random(913);for(auto& v:original)v={random(),random()};
    const auto restore=[&]{std::memcpy(input.p,original.data(),original.size()*sizeof(c2::Block));};
    for(unsigned i=0;i!=candidates;++i){restore();for(unsigned j=0;j!=5;++j)execute(i,p,input.p,input.p,planes.p,temp.p);}
    std::array<std::vector<double>,candidates> total,inside,outside;std::array<std::uint64_t,candidates>sums{};
    for(unsigned i=0;i!=candidates;++i){total[i].resize(calls);inside[i].resize(calls);outside[i].resize(calls);}
    for(unsigned done=0,pass=0;done<calls;++pass){
        const auto batch=std::min(5U,calls-done);
        for(unsigned order=0;order!=candidates;++order){
            const auto i=(order+pass)%candidates;restore();
            for(unsigned j=0;j!=batch;++j){
                if(fresh&&j)restore();const auto a=Clock::now();
                if(phases){inner(i,p,input.p,planes.p,temp.p);const auto b=Clock::now();outer(i,p,planes.p,input.p);const auto c=Clock::now();inside[i][done+j]=us(a,b);outside[i][done+j]=us(b,c);total[i][done+j]=us(a,c);}
                else{execute(i,p,input.p,input.p,planes.p,temp.p);total[i][done+j]=us(a,Clock::now());}
            }
            sums[i]=(sums[i]^checksum(input.p,k)^batch)*0x100000001b3ULL;
        }
        done+=batch;
    }
    for(unsigned i=1;i!=candidates;++i)require(sums[0]==sums[i],"next: benchmark same-map checksum");
    std::cout<<"K,"<<k<<",N,"<<p.outputN<<",calls,"<<calls<<",restore_each_call,"<<fresh<<",setup_excluded,1\nname,median_us,p10_us,p90_us\n";
    for(unsigned i=0;i!=candidates;++i){stats(names[i],total[i]);if(phases){stats(std::string(names[i])+"-inner-route",inside[i]);stats(std::string(names[i])+"-outer",outside[i]);}}
}
// Final two-way benchmark uses only the selected plan, without unused
// gauge/pair setup. Both calls see identical payload/scratch addresses.
static void benchSelected(std::size_t k,unsigned calls,bool fresh){
    c2::RateSpeed4RoutedPlan<40> p(k,113);
    const auto page=((p.groups+255)/256)*256,selected=p.groups+4;
    c2::BenchPayload input(p.outputN),planes(64*std::max(page,selected)),temp(p.groups);
    std::vector<c2::Block> original(p.outputN);std::mt19937_64 random(913);for(auto& v:original)v={random(),random()};
    const auto restore=[&]{std::memcpy(input.p,original.data(),original.size()*sizeof(c2::Block));};
    const auto run=[&](unsigned i){if(i)c2::rs4TransposeRouted16(p,input.p,input.p,planes.p,temp.p);else c2::rs4TransposeCached(p,input.p,input.p,planes.p,temp.p);};
    for(unsigned i=0;i!=2;++i){p.planeStride=i?selected:page;restore();for(unsigned j=0;j!=5;++j)run(i);}
    std::array<std::vector<double>,2> times;std::array<std::uint64_t,2>sums{};for(auto& values:times)values.resize(calls);
    for(unsigned done=0,pass=0;done<calls;++pass){const auto batch=std::min(5U,calls-done);
        for(unsigned order=0;order!=2;++order){const auto i=(order+pass)%2;p.planeStride=i?selected:page;restore();
            for(unsigned j=0;j!=batch;++j){if(fresh&&j)restore();const auto start=Clock::now();run(i);times[i][done+j]=us(start,Clock::now());}
            sums[i]=(sums[i]^checksum(input.p,k)^batch)*0x100000001b3ULL;
        }done+=batch;
    }
    require(sums[0]==sums[1],"next: selected entry-point checksum");
    std::cout<<"K,"<<k<<",N,"<<p.outputN<<",calls,"<<calls<<",restore_each_call,"<<fresh<<",selected_only_plan,1\nname,median_us,p10_us,p90_us\n";
    stats("previous-cached-eight",times[0]);stats("selected-routed16",times[1]);
}
int main(int argc,char** argv){try{
    if(!c2::gfni512Available())return 77;
    if(argc==2&&std::string(argv[1])=="--check"){checkLocalHistory();check<40>(false);check<39>(false);return 0;}
    if(argc==2&&std::string(argv[1])=="--large-check"){check<40>(true);check<39>(true);return 0;}
    if(argc<2||argc>4)throw std::invalid_argument("usage: K [calls=51] [full|phases|fresh|selected|selected-fresh]");
    const auto k=std::stoull(argv[1]),calls=argc>2?std::stoull(argv[2]):51;const std::string mode=argc>3?argv[3]:"full";
    if(!calls||calls>1000000||(mode!="full"&&mode!="phases"&&mode!="fresh"&&mode!="selected"&&mode!="selected-fresh"))throw std::invalid_argument("calls/mode");
    if(mode=="selected"||mode=="selected-fresh")benchSelected(k,static_cast<unsigned>(calls),mode=="selected-fresh");else bench(k,static_cast<unsigned>(calls),mode=="phases",mode=="fresh");return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
