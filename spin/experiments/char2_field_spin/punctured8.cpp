#include "Punctured8Code.h"
#include "Punctured8Oracle.h"
#include "BenchPayload.h"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>

namespace c2=spin::research::char2;
using Clock=std::chrono::steady_clock;
static void require(bool good,const char* message){if(!good)throw std::runtime_error(message);}
static double micros(Clock::time_point a,Clock::time_point b){return std::chrono::duration<double,std::micro>(b-a).count();}
static constexpr auto gather=c2::Cubic8OuterVariant::Gather,fused=c2::Cubic8OuterVariant::FusedLoad;
static constexpr c2::Subfield16Kernel schedules[]={c2::Subfield16Kernel::Scatter,
    c2::Subfield16Kernel::RegionNonTemporal,c2::Subfield16Kernel::RegionAllNonTemporal};

static void checkPacket(){
    const auto& o=c2::p8Oracle();std::mt19937_64 rng(894);
    for(unsigned i=0;i!=8;++i)for(unsigned j=0;j!=8;++j)
        require(c2::p8A[i^j]==o.effectiveA[i][j]&&c2::p8C[i^j]==o.c[i][j]&&
                c2::p8F[i^j]==o.f[i][j],"native/oracle Cauchy matrix identity");
    for(unsigned trial=0;trial!=584;++trial){
        alignas(64) c2::Block input[28]{},state[8]{},expected[32]{},expectedState[8]{},actual[32],actualState[8];
        if(trial<576){
            const auto index=trial/16,base=16*(trial%16);auto* target=index<28?input+index:state+index-28;
            for(unsigned b=0;b!=16;++b)c2::p8OracleBytes(target)[b]=static_cast<std::uint8_t>(base+b);
        }else{
            for(auto& v:input)v={rng(),rng()};for(auto& v:state)v={rng(),rng()};
        }
        std::uint8_t dx[8],ds[32];std::fill(std::begin(dx),std::end(dx),1);
        for(auto& v:ds)v=static_cast<std::uint8_t>(1+rng()%255);
        std::copy(std::begin(state),std::end(state),expectedState);
        for(unsigned packet=4;packet--;){
            c2::Block in[8]{},next[8];
            for(unsigned j=0;j!=7;++j)in[j]=input[4*j+packet];
            c2::p8OraclePacket(in,expectedState,dx,ds+8*packet,expected+8*packet,next,true);
            std::copy(std::begin(next),std::end(next),expectedState);
        }
        auto h0=_mm512_load_si512(state),h1=_mm512_load_si512(state+4);
        c2::p8InnerTransposeFour(input,actual,h0,h1,ds);
        _mm512_store_si512(actualState,h0);_mm512_store_si512(actualState+4,h1);
        require(std::equal(std::begin(expected),std::end(expected),actual),"four-packet all-byte scalar output");
        require(std::equal(std::begin(expectedState),std::end(expectedState),actualState),"four-packet all-byte scalar state");
    }
}
static void check(){
    c2::p8OracleCheckLocal();checkPacket();std::mt19937_64 rng(173);
    const c2::Block guard{0x32456789ULL,0xfabc1234ULL};
    for(std::size_t k:{1,83,84,85,335,336,337,5376,5377,6000}){
        c2::Punctured8Plan p(k,113),full(p.k,113);
        require(p.outputN==full.outputN&&p.route==full.route&&p.inputDiagonal==full.inputDiagonal&&
            p.stateDiagonal==full.stateDiagonal,"padded setup identity");
        std::vector<c2::Block> original(p.outputN);for(auto& value:original)value={rng(),rng()};
        const auto expected=c2::p8OracleTranspose(p,original.data());
        c2::BenchPayload input(p.outputN+8,false),out(k+8,false),planes(128*(p.groups+4)+8,false),temp(p.groups+8,false);
        for(auto stride:{p.groups,p.groups+4})for(auto variant:{gather,fused})for(auto mode:schedules)for(unsigned offset:{0,1,2,3}){
            p.planeStride=stride;
            for(auto* buffer:{&input,&out,&planes,&temp})std::fill(buffer->p,buffer->p+buffer->count,guard);
            std::memcpy(input.p+offset,original.data(),p.outputN*sizeof(c2::Block));
            c2::p8Transpose(p,input.p+offset,out.p+offset,planes.p+offset,temp.p+offset,mode,variant);
            require(std::equal(expected.begin(),expected.end(),out.p+offset),"full independent scalar transpose");
            require(std::equal(original.begin(),original.end(),input.p+offset),"input preservation");
            c2::p8Transpose(p,input.p+offset,input.p+offset,planes.p+offset,temp.p+offset,mode,variant);
            require(std::equal(expected.begin(),expected.end(),input.p+offset),"alias prefix");
            require(std::equal(original.begin()+k,original.end(),input.p+offset+k),"alias suffix");
            for(auto* buffer:{&input,&out,&planes,&temp}){
                for(unsigned j=0;j!=offset;++j)require(buffer->p[j]==guard,"prefix guard");
                for(std::size_t j=buffer->count-8+offset;j!=buffer->count;++j)require(buffer->p[j]==guard,"suffix guard");
            }
        }
        p.inverseLocalRoute.clear();
        c2::p8Transpose(p,original.data(),out.p,planes.p,temp.p,c2::Subfield16Kernel::RegionNonTemporal);
        require(std::equal(expected.begin(),expected.end(),out.p),"uint32 route fallback");
        if(k==337||k==6000){
            std::vector<c2::Block> source(k);for(auto& value:source)value={rng(),rng()};
            const auto image=c2::p8OracleForward(p,source.data());
            require(c2::p8OracleDot(image.data(),original.data(),p.outputN)==
                    c2::p8OracleDot(source.data(),expected.data(),k),"full independent shortened adjoint");
        }
    }
    std::cout<<"punctured8 independent scalar/map/shortened-adjoint/schedule/alias/guard/tail checks pass\n";
}
static void checkLarge(){
    std::mt19937_64 rng(8439);
    for(unsigned exponent:{16,18,20}){
        c2::Punctured8Plan p(std::size_t(1)<<exponent,113);
        c2::BenchPayload input(p.outputN),out(p.requested),expected(p.requested),planes(128*(p.groups+4)),temp(p.groups),alias(p.outputN);
        for(std::size_t i=0;i!=p.outputN;++i)input.p[i]={rng(),rng()};
        c2::p8Transpose(p,input.p,expected.p,planes.p,temp.p,c2::Subfield16Kernel::Scatter,gather);
        for(auto stride:{p.groups,p.groups+4})for(auto variant:{gather,fused})for(auto mode:schedules){
            p.planeStride=stride;
            c2::p8Transpose(p,input.p,out.p,planes.p,temp.p,mode,variant);
            require(std::equal(expected.p,expected.p+p.requested,out.p),"large schedules");
            std::memcpy(alias.p,input.p,p.outputN*sizeof(c2::Block));
            c2::p8Transpose(p,alias.p,alias.p,planes.p,temp.p,mode,variant);
            require(std::equal(expected.p,expected.p+p.requested,alias.p),"large alias prefix");
            require(std::equal(input.p+p.requested,input.p+p.outputN,alias.p+p.requested),"large alias suffix");
        }
        std::cout<<"K"<<exponent<<" punctured8 schedules and aliases pass\n";
    }
}
static std::uint64_t checksum(const c2::Block* values,std::size_t n){
    std::uint64_t result=0;
    for(std::size_t i=0;i!=n;++i){result=(result^values[i].lo)*0x100000001b3ULL;result=(result^values[i].hi)*0x100000001b3ULL;}
    return result;
}
static void stats(const std::string& name,std::size_t k,std::size_t n,std::vector<double> values,std::uint64_t sum){
    std::sort(values.begin(),values.end());
    std::cout<<name<<','<<k<<','<<n<<','<<values.size()<<','<<std::fixed<<std::setprecision(4)
        <<values[values.size()/2]<<','<<values[values.size()/10]<<','<<values[9*values.size()/10]
        <<','<<std::hex<<sum<<std::dec<<'\n';
}
static void benchmark(std::size_t k,unsigned calls,bool phases){
    c2::Punctured8Plan p(k,113);c2::Cubic8Plan old(k,113);
    const auto mode=k>=1048576?c2::Subfield16Kernel::RegionAllNonTemporal:c2::Subfield16Kernel::Scatter;
    const auto oldVariant=k>=1048576?gather:fused;
    const auto records=std::max(p.outputN,old.n),planeRecords=std::max(128*(p.groups+4),old.n);
    c2::BenchPayload input(records+4,k>65536),planes(planeRecords,k>65536),temp(std::max(p.groups,old.groups),k>65536);
    std::vector<c2::Block> original(records);std::mt19937_64 rng(913);
    for(auto& value:original)value={rng(),rng()};
    const auto restore=[&]{std::memcpy(input.p,original.data(),records*sizeof(c2::Block));};
    const auto encode=[&](unsigned candidate){
        p.planeStride=candidate==1?p.groups:p.groups+4;
        if(candidate==0)c2::c8TransposePaired(old,input.p,input.p,planes.p,temp.p,mode,oldVariant);
        else c2::p8Transpose(p,input.p,input.p,planes.p,temp.p,mode,candidate==3?fused:gather);
    };
    for(unsigned candidate=0;candidate!=4;++candidate){restore();for(unsigned i=0;i!=5;++i)encode(candidate);}
    std::array<std::vector<double>,4> samples,inner,outer;std::array<std::uint64_t,4> sums{};
    for(unsigned candidate=0;candidate!=4;++candidate){samples[candidate].resize(calls);
        if(phases){inner[candidate].resize(calls);outer[candidate].resize(calls);}}
    for(unsigned done=0,pass=0;done<calls;++pass){
        const unsigned count=std::min(5U,calls-done);
        for(unsigned order=0;order!=4;++order){
            const unsigned candidate=(order+pass)%4;restore();
            p.planeStride=candidate==1?p.groups:p.groups+4;
            for(unsigned i=0;i!=count;++i){
                const auto a=Clock::now();
                if(phases){
                    if(candidate==0)c2::c8InnerDispatchPaired(old,input.p,planes.p,temp.p,mode);
                    else c2::p8InnerDispatch(p,input.p,planes.p,temp.p,mode);
                    const auto b=Clock::now();
                    if(candidate==0)c2::c8OuterDispatch(old,planes.p,input.p,mode,oldVariant);
                    else c2::p8OuterDispatch(p,planes.p,input.p,mode,candidate==3?fused:gather);
                    const auto c=Clock::now();inner[candidate][done+i]=micros(a,b);outer[candidate][done+i]=micros(b,c);
                    samples[candidate][done+i]=micros(a,c);
                }else{encode(candidate);samples[candidate][done+i]=micros(a,Clock::now());}
            }
            sums[candidate]=(sums[candidate]^checksum(input.p,k)^count)*0x100000001b3ULL;
        }
        done+=count;
    }
    require(sums[1]==sums[2]&&sums[1]==sums[3],"punctured layout/outer benchmark checksums");
    std::cout<<"same_input_address,"<<static_cast<void*>(input.p)<<",setup_excluded,1,restore_outside_timing,batch5,rotating_order"
        <<",padded_K,"<<p.k<<",input_N,"<<p.inputN<<",output_N,"<<p.outputN
        <<",different_codes_and_rates,1,schedule,"<<(k>=1048576?"region-nt-all":"scatter")
        <<"\nname,K,N,calls,median_us,p10_us,p90_us,checksum\n";
    const char* names[]={"cubic8-coefficients83-paired","punctured8-tight-gather","punctured8-padded-gather","punctured8-padded-fused"};
    for(unsigned candidate=0;candidate!=4;++candidate){
        const auto n=candidate?p.outputN:old.n;stats(names[candidate],k,n,samples[candidate],sums[candidate]);
        if(phases){stats(std::string(names[candidate])+"-inner-route",k,n,inner[candidate],0);
                   stats(std::string(names[candidate])+"-outer",k,n,outer[candidate],0);}
    }
}
int main(int argc,char** argv){try{
    if(!c2::gfni512Available())return 77;
    if(argc==2&&std::string(argv[1])=="--check"){check();return 0;}
    if(argc==2&&std::string(argv[1])=="--large-check"){checkLarge();return 0;}
    if(argc<2||argc>4)throw std::invalid_argument("usage: K [calls=51] [full|phases]");
    const auto k=std::stoull(argv[1]),calls=argc>2?std::stoull(argv[2]):51;
    const std::string mode=argc>3?argv[3]:"full";
    if(!calls||calls>1000000||(mode!="full"&&mode!="phases"))throw std::invalid_argument("calls or mode");
    benchmark(k,static_cast<unsigned>(calls),mode=="phases");return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
