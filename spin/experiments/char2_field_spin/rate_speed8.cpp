#include "RateSpeed8Code.h"
#include "RateClimb8Oracle.h"
#include "BenchPayload.h"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>

namespace c2=spin::research::char2;
using Plan=c2::RateClimb8Plan<c2::R8Periodic84>;
using Clock=std::chrono::steady_clock;
constexpr unsigned Candidates=12;
static const char* names[Candidates]={"frozen84","column","tile16","tile64","tile256","tile1024","column-pf16","tile64-pf16","affine-maps","affine-all","fused-outer","affine-fused"};
static void require(bool good,const char* text){if(!good)throw std::runtime_error(text);}
static double elapsed(Clock::time_point a,Clock::time_point b){return std::chrono::duration<double,std::micro>(b-a).count();}

static void inner(unsigned candidate,const Plan& p,const c2::Block* input,c2::Block* planes,c2::Block* temp){
    switch(candidate){
    case 0:c2::r8InnerDispatch(p,input,planes,temp,c2::Subfield16Kernel::RegionAllNonTemporal,true);break;
    case 1:c2::s8RegionInner<0>(p,input,planes,temp);break;
    case 2:c2::s8RegionInner<16>(p,input,planes,temp);break;
    case 3:c2::s8RegionInner<64>(p,input,planes,temp);break;
    case 4:c2::s8RegionInner<256>(p,input,planes,temp);break;
    case 5:c2::s8RegionInner<1024>(p,input,planes,temp);break;
    case 6:c2::s8RegionInner<0,16>(p,input,planes,temp);break;
    case 7:c2::s8RegionInner<64,16>(p,input,planes,temp);break;
    case 8:c2::s8AffineRegionInner<true,false>(p,input,planes,temp);break;
    case 9:case 11:c2::s8AffineRegionInner<true,true>(p,input,planes,temp);break;
    default:c2::r8InnerDispatch(p,input,planes,temp,c2::Subfield16Kernel::RegionAllNonTemporal,true);break;
    }
}
static void outer(unsigned candidate,const Plan& p,const c2::Block* planes,c2::Block* output){
    switch(candidate){
    case 0:c2::r8OuterDispatchFinalStore(p,planes,output,c2::Subfield16Kernel::RegionAllNonTemporal,c2::Cubic8OuterVariant::FusedLoad);break;
    case 1:case 6:case 8:case 9:c2::s8Outer<0>(p,planes,output);break;
    case 2:c2::s8Outer<16>(p,planes,output);break;
    case 3:case 7:c2::s8Outer<64>(p,planes,output);break;
    case 4:c2::s8Outer<256>(p,planes,output);break;
    case 5:c2::s8Outer<1024>(p,planes,output);break;
    default:c2::r8SpeedOuterDispatch(p,planes,output,c2::Subfield16Kernel::RegionAllNonTemporal,c2::Cubic8OuterVariant::FusedLoad);break;
    }
}
static void encode(unsigned candidate,const Plan& p,const c2::Block* input,c2::Block* output,c2::Block* planes,c2::Block* temp){
    inner(candidate,p,input,planes,temp);outer(candidate,p,planes,output);
}

static void check(bool large){
    std::vector<std::size_t> sizes=large?std::vector<std::size_t>{65536,262144,1048576}:
        std::vector<std::size_t>{1,83,84,85,335,336,337,5376,5377};
    std::mt19937_64 random(137);const c2::Block guard{0x992271fe1ULL,0x3762471ULL};
    for(auto k:sizes){
        Plan p(k,113);p.planeStride=((p.groups+255)/256)*256;
        c2::BenchPayload input(p.outputN+8),out(k+8),planes(c2::s8PlaneRecords<1024>(p)+8),temp(p.groups+8);
        std::vector<c2::Block> original(p.outputN),expected(k);
        for(auto& v:original)v={random(),random()};
        if(large){std::memcpy(input.p,original.data(),p.outputN*sizeof(c2::Block));encode(0,p,input.p,out.p,planes.p,temp.p);std::copy(out.p,out.p+k,expected.begin());}
        else expected=c2::r8OracleTranspose(p,original.data());
        for(unsigned candidate=0;candidate!=Candidates;++candidate){
            for(auto* b:{&input,&out,&planes,&temp})std::fill(b->p,b->p+b->count,guard);
            std::memcpy(input.p,original.data(),p.outputN*sizeof(c2::Block));
            encode(candidate,p,input.p,out.p,planes.p,temp.p);
            require(std::equal(expected.begin(),expected.end(),out.p),"independent/schedule result");
            require(std::equal(original.begin(),original.end(),input.p),"input preservation");
            encode(candidate,p,input.p,input.p,planes.p,temp.p);
            require(std::equal(expected.begin(),expected.end(),input.p),"alias prefix");
            require(std::equal(original.begin()+k,original.end(),input.p+k),"alias suffix");
            for(auto* b:{&input,&out,&planes,&temp})for(auto j=b->count-8;j!=b->count;++j)require(b->p[j]==guard,"guard");
        }
        std::cout<<"K="<<k<<" all exact-map schedules/aliases/guards pass\n";
    }
}
static std::uint64_t checksum(const c2::Block* input,std::size_t count){
    std::uint64_t value=0;for(std::size_t i=0;i!=count;++i){value=(value^input[i].lo)*0x100000001b3ULL;value=(value^input[i].hi)*0x100000001b3ULL;}return value;
}
static void stats(const std::string& name,std::vector<double> samples){
    std::sort(samples.begin(),samples.end());
    std::cout<<name<<','<<std::fixed<<std::setprecision(3)<<samples[samples.size()/2]<<','<<samples[samples.size()/10]<<','<<samples[9*samples.size()/10]<<'\n';
}
static void bench(std::size_t k,unsigned calls,bool phases){
    Plan p(k,113);p.planeStride=((p.groups+255)/256)*256;
    c2::BenchPayload input(p.outputN),planes(c2::s8PlaneRecords<1024>(p)),temp(p.groups);
    std::vector<c2::Block> original(p.outputN);std::mt19937_64 random(913);for(auto& v:original)v={random(),random()};
    const auto restore=[&]{std::memcpy(input.p,original.data(),p.outputN*sizeof(c2::Block));};
    for(unsigned candidate=0;candidate!=Candidates;++candidate){restore();for(unsigned i=0;i!=5;++i)encode(candidate,p,input.p,input.p,planes.p,temp.p);}
    std::array<std::vector<double>,Candidates> full,inside,outside;std::array<std::uint64_t,Candidates> sums{};
    for(unsigned i=0;i!=Candidates;++i){full[i].resize(calls);inside[i].resize(calls);outside[i].resize(calls);}
    for(unsigned done=0,pass=0;done<calls;++pass){
        const auto count=std::min(5U,calls-done);
        for(unsigned order=0;order!=Candidates;++order){
            const auto candidate=(order+pass)%Candidates;restore();
            for(unsigned j=0;j!=count;++j){
                auto a=Clock::now();
                if(phases){inner(candidate,p,input.p,planes.p,temp.p);auto b=Clock::now();outer(candidate,p,planes.p,input.p);auto c=Clock::now();inside[candidate][done+j]=elapsed(a,b);outside[candidate][done+j]=elapsed(b,c);full[candidate][done+j]=elapsed(a,c);}
                else{encode(candidate,p,input.p,input.p,planes.p,temp.p);full[candidate][done+j]=elapsed(a,Clock::now());}
            }
            sums[candidate]=(sums[candidate]^checksum(input.p,k)^count)*0x100000001b3ULL;
        }
        done+=count;
    }
    for(auto sum:sums)require(sum==sums[0],"benchmark exact-map checksum");
    std::cout<<"K,"<<k<<",N,"<<p.outputN<<",same_input,"<<static_cast<void*>(input.p)<<",same_scratch,"<<static_cast<void*>(planes.p)<<",calls,"<<calls<<",setup_excluded,1\nname,median_us,p10_us,p90_us\n";
    for(unsigned i=0;i!=Candidates;++i){stats(names[i],full[i]);if(phases){stats(std::string(names[i])+"-inner-route",inside[i]);stats(std::string(names[i])+"-outer",outside[i]);}}
}
int main(int argc,char** argv){try{
    if(!c2::gfni512Available())return 77;
    if(argc==2&&std::string(argv[1])=="--check"){check(false);return 0;}
    if(argc==2&&std::string(argv[1])=="--large-check"){check(true);return 0;}
    if(argc<2||argc>4)throw std::invalid_argument("usage: K [calls=51] [full|phases]");
    const auto k=std::stoull(argv[1]),calls=argc>2?std::stoull(argv[2]):51;
    const std::string mode=argc>3?argv[3]:"full";
    if(!calls||calls>1000000||(mode!="full"&&mode!="phases"))throw std::invalid_argument("calls/mode");
    bench(k,static_cast<unsigned>(calls),mode=="phases");return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
