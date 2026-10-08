#include "RateSpeed4Profile.h"
#include "RateSpeed4Pipeline.h"
#include "RateSpeed4Oracle.h"
#include "BenchPayload.h"
#include <iostream>
#include <iomanip>
#include <string>
namespace c2=spin::research::char2;
using Clock=std::chrono::steady_clock;
static constexpr unsigned candidates=16;
static const char* names[candidates]={"frozen-routed16","profile-control","balanced-gather","dword-broadcast","balanced-dword","byte-offsets","byte-offsets-balanced-dword","compact-inner","rolling-inner","compact-rolling-inner","double-buffer-pipeline","qword-offsets","qword-offsets-balanced-dword","qword-compact-rolling","explicit-record-loads","explicit-balanced-dword"};
static void require(bool b,const char* what){if(!b)throw std::runtime_error(what);}
template<bool Timed=false,unsigned D>static void execute(unsigned variant,const c2::RateSpeed4ProfilePlan<D>& p,const c2::Block* input,c2::Block* output,c2::Block* planes,c2::Block* temp,c2::Rs4PhaseTimes* phases=nullptr){
    switch(variant){
    case 0:c2::rs4TransposeRouted16(p,input,output,planes,temp);break;
    case 1:c2::rs4TransposeProfile<0,Timed>(p,input,output,planes,temp,phases);break;
    case 2:c2::rs4TransposeProfile<1,Timed>(p,input,output,planes,temp,phases);break;
    case 3:c2::rs4TransposeProfile<2,Timed>(p,input,output,planes,temp,phases);break;
    case 4:c2::rs4TransposeProfile<3,Timed>(p,input,output,planes,temp,phases);break;
    case 5:c2::rs4TransposeProfile<4,Timed>(p,input,output,planes,temp,phases);break;
    case 6:c2::rs4TransposeProfile<7,Timed>(p,input,output,planes,temp,phases);break;
    case 7:c2::rs4TransposeProfile<8,Timed>(p,input,output,planes,temp,phases);break;
    case 8:c2::rs4TransposeProfile<16,Timed>(p,input,output,planes,temp,phases);break;
    case 9:c2::rs4TransposeProfile<24,Timed>(p,input,output,planes,temp,phases);break;
    case 10:
        if constexpr(Timed){const auto a=Clock::now();c2::rs4RegionPipeline<D>(p,input,planes,temp);const auto b=Clock::now();
            if((reinterpret_cast<std::uintptr_t>(output)&63)==0)c2::rs4OuterNext<D,true,c2::Rs4OuterNextVariant::FusedMultiply>(p,planes,output);else c2::rs4OuterNext<D,false,c2::Rs4OuterNextVariant::FusedMultiply>(p,planes,output);
            phases->compute=std::chrono::duration<double,std::micro>(b-a).count();phases->outer=std::chrono::duration<double,std::micro>(Clock::now()-b).count();
        }else c2::rs4TransposePipeline<D>(p,input,output,planes,temp);
        break;
    case 11:c2::rs4TransposeProfile<64,Timed>(p,input,output,planes,temp,phases);break;
    case 12:c2::rs4TransposeProfile<67,Timed>(p,input,output,planes,temp,phases);break;
    case 13:c2::rs4TransposeProfile<88,Timed>(p,input,output,planes,temp,phases);break;
    case 14:c2::rs4TransposeProfile<128,Timed>(p,input,output,planes,temp,phases);break;
    case 15:c2::rs4TransposeProfile<131,Timed>(p,input,output,planes,temp,phases);break;
    }
}
static void checkInner(){
    std::mt19937_64 random(819);
    for(std::size_t packets:{16,32,64,160})for(unsigned trial=0;trial!=32;++trial){
        std::vector<c2::Block> input(3*packets),expected(4*packets),out(4*packets);std::vector<std::uint8_t> ds(4*packets);
        alignas(64)c2::Block initial[4],wantState[4],gotState[4];for(auto& v:input)v={random(),random()};for(auto& v:initial)v={random(),random()};for(auto& d:ds)d=static_cast<std::uint8_t>(1+random()%255);
        auto h=_mm512_load_si512(initial);for(auto remaining=packets;remaining;){remaining-=16;c2::rs4InnerTransposeSixteenPipeline(input.data()+3*remaining,expected.data()+4*remaining,h,ds.data()+4*remaining);}_mm512_store_si512(wantState,h);
        for(unsigned v=0;v!=3;++v){h=_mm512_load_si512(initial);
            if(v==0)c2::rs4InnerTransposeCompact(input.data(),out.data(),h,ds.data(),packets);else if(v==1)c2::rs4InnerTransposeRolling<false>(input.data(),out.data(),h,ds.data(),packets);else c2::rs4InnerTransposeRolling<true>(input.data(),out.data(),h,ds.data(),packets);
            _mm512_store_si512(gotState,h);require(expected==out,"profile: inner arbitrary history");require(std::equal(wantState,wantState+4,gotState),"profile: inner final state");
        }
    }
    std::cout<<"profile compact/rolling arbitrary-history checks pass\n";
}
template<unsigned D>static void check(bool large){
    std::mt19937_64 rng(119);const c2::Block guard{771,997};
    const std::vector<std::size_t> sizes=large?std::vector<std::size_t>{65536,262144,1048576}:std::vector<std::size_t>{1,D-1,D,D+1,129,130,131,4*D+1,4*D+2,4*D+3,64*D+1};
    for(auto k:sizes){
        c2::RateSpeed4ProfilePlan<D> p(k,993);
        c2::BenchPayload input(p.outputN+8),output(k+8),planes(p.planeRecords()+8),temp(2*p.groups+8);
        std::vector<c2::Block> original(p.outputN),expected(k);for(auto& v:original)v={rng(),rng()};
        if(large){c2::rs4TransposeRouted16(p,original.data(),output.p,planes.p,temp.p);std::copy(output.p,output.p+k,expected.begin());}
        else expected=c2::rs4OracleTranspose(p,original.data());
        for(unsigned i=0;i!=candidates;++i)for(unsigned offset:{0,1,2,3}){
            for(auto* b:{&input,&output,&planes,&temp})std::fill(b->p,b->p+b->count,guard);
            std::memcpy(input.p+offset,original.data(),original.size()*sizeof(c2::Block));
            execute(i,p,input.p+offset,output.p+offset,planes.p+offset,temp.p+offset);
            require(std::equal(expected.begin(),expected.end(),output.p+offset),"profile: reference equality");
            require(std::equal(original.begin(),original.end(),input.p+offset),"profile: input preservation");
            execute(i,p,input.p+offset,input.p+offset,planes.p+offset,temp.p+offset);
            require(std::equal(expected.begin(),expected.end(),input.p+offset),"profile: alias prefix");
            require(std::equal(original.begin()+k,original.end(),input.p+offset+k),"profile: alias suffix");
            for(auto* b:{&input,&output,&planes,&temp}){for(unsigned j=0;j!=offset;++j)require(b->p[j]==guard,"profile: prefix guard");for(auto j=b->count-8+offset;j!=b->count;++j)require(b->p[j]==guard,"profile: suffix guard");}
            c2::Rs4PhaseTimes phases;execute<true>(i,p,original.data(),output.p,planes.p,temp.p,&phases);require(std::equal(expected.begin(),expected.end(),output.p),"profile: measured schedule");
        }
        p.foldedInverseLocalRoute.clear();p.byteOffset.clear();p.qwordOffset.clear();for(unsigned i=0;i!=candidates;++i){execute(i,p,original.data(),output.p,planes.p,temp.p);require(std::equal(expected.begin(),expected.end(),output.p),"profile: wide route");}
        std::cout<<"profile D="<<D<<" K="<<k<<" exact/alias/guards/wide/timed pass\n";
    }
}
static std::uint64_t checksum(const c2::Block* data,std::size_t n){std::uint64_t h=0;for(std::size_t i=0;i!=n;++i){h=(h^data[i].lo)*0x100000001b3ULL;h=(h^data[i].hi)*0x100000001b3ULL;}return h;}
static void stats(const std::string& name,std::vector<double> v){std::sort(v.begin(),v.end());std::cout<<name<<','<<std::fixed<<std::setprecision(3)<<v[v.size()/2]<<','<<v[v.size()/10]<<','<<v[9*v.size()/10]<<'\n';}
static void bench(std::size_t k,unsigned calls,bool timed,bool fresh){
    c2::RateSpeed4ProfilePlan<40> p(k,113);c2::BenchPayload input(p.outputN),planes(p.planeRecords()),temp(2*p.groups);
    std::vector<c2::Block> original(p.outputN);std::mt19937_64 random(913);for(auto& v:original)v={random(),random()};
    const auto restore=[&]{std::memcpy(input.p,original.data(),original.size()*sizeof(c2::Block));};
    std::array<std::vector<double>,candidates> total,compute,route,outer;std::array<std::uint64_t,candidates> sums{};
    for(unsigned i=0;i!=candidates;++i){total[i].resize(calls);compute[i].resize(calls);route[i].resize(calls);outer[i].resize(calls);restore();for(unsigned j=0;j!=5;++j)execute(i,p,input.p,input.p,planes.p,temp.p);}
    for(unsigned done=0,pass=0;done<calls;++pass){const auto batch=std::min(5U,calls-done);
        for(unsigned order=0;order!=candidates;++order){const auto i=(order+pass)%candidates;restore();
            for(unsigned j=0;j!=batch;++j){if(fresh&&j)restore();c2::Rs4PhaseTimes phases;const auto start=Clock::now();
                if(timed)execute<true>(i,p,input.p,input.p,planes.p,temp.p,&phases);else execute(i,p,input.p,input.p,planes.p,temp.p);
                total[i][done+j]=std::chrono::duration<double,std::micro>(Clock::now()-start).count();compute[i][done+j]=phases.compute;route[i][done+j]=phases.route;outer[i][done+j]=phases.outer;
            }sums[i]=(sums[i]^checksum(input.p,k)^batch)*0x100000001b3ULL;
        }done+=batch;
    }
    for(unsigned i=1;i!=candidates;++i)require(sums[0]==sums[i],"profile: benchmark checksum");
    std::cout<<"K,"<<k<<",N,"<<p.outputN<<",calls,"<<calls<<",region_clocks,"<<timed<<",restore_each_call,"<<fresh<<"\nname,median_us,p10_us,p90_us\n";
    for(unsigned i=0;i!=candidates;++i){stats(names[i],total[i]);if(timed&&i){stats(std::string(names[i])+(i==10?"-inner-route":"-inner"),compute[i]);if(i!=10)stats(std::string(names[i])+"-route",route[i]);stats(std::string(names[i])+"-outer",outer[i]);}}
}
int main(int argc,char** argv){try{
    if(!c2::gfni512Available())return 77;
    if(argc==2&&std::string(argv[1])=="--check"){checkInner();check<40>(false);check<39>(false);return 0;}
    if(argc==2&&std::string(argv[1])=="--large-check"){check<40>(true);check<39>(true);return 0;}
    if(argc<2||argc>4)throw std::invalid_argument("usage: K [calls=71] [full|phases|fresh]");
    const auto k=std::stoull(argv[1]),calls=argc>2?std::stoull(argv[2]):71;const std::string mode=argc>3?argv[3]:"full";
    if(!calls||calls>1000000||(mode!="full"&&mode!="phases"&&mode!="fresh"))throw std::invalid_argument("calls/mode");
    bench(k,static_cast<unsigned>(calls),mode=="phases",mode=="fresh");return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
