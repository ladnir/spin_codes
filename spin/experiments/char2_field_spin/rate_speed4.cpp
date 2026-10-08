#include "RateSpeed4Code.h"
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
static void checkPacket(){
    std::mt19937_64 random(933);std::uint8_t dx[4]={1,1,1,1};
    for(unsigned trial=0;trial!=16*16+16;++trial){
        alignas(64)c2::Block input[12]{},state[4]{},expected[16]{},expectedState[4]{},actual[16]{},actualState[4]{};
        std::uint8_t ds[16];for(auto& d:ds)d=static_cast<std::uint8_t>(1+random()%255);
        if(trial<256){const unsigned index=trial/16,first=16*(trial%16);auto* ptr=index<12?input+index:state+index-12;for(unsigned byte=0;byte!=16;++byte)c2::p8OracleBytes(ptr)[byte]=static_cast<std::uint8_t>(first+byte);}
        else{for(auto& v:input)v={random(),random()};for(auto& v:state)v={random(),random()};}
        std::copy(state,state+4,expectedState);
        for(unsigned phase=4;phase--;){
            c2::Block packet[4]{},raw[4]{},next[4]{};for(unsigned j=0;j!=3;++j)packet[j]=input[4*j+phase];
            c2::rs4OraclePacket(packet,expectedState,dx,ds+4*phase,raw,next,true);
            for(unsigned j=0;j!=4;++j)expected[4*j+phase]=raw[j];std::copy(next,next+4,expectedState);
        }
        for(bool affine:{false,true})for(bool packed:{false,true}){
            auto h=_mm512_load_si512(state);
            if(packed){if(affine)c2::rs4InnerTransposeFour<true,true,true>(input,actual,h,ds);else c2::rs4InnerTransposeFour<false,false,true>(input,actual,h,ds);}
            else if(affine)c2::rs4InnerTransposeFour<true,true>(input,actual,h,ds);
            else c2::rs4InnerTransposeFour<false,false>(input,actual,h,ds);
            _mm512_store_si512(actualState,h);
            require(std::equal(expected,expected+16,actual),"dense packet output");require(std::equal(expectedState,expectedState+4,actualState),"dense packet state");
        }
    }
    std::cout<<"4state all-byte dense packet/state checks pass\n";
}
template<unsigned D>static void execute4(const c2::RateSpeed4Plan<D>& p,const c2::Block* input,c2::Block* output,c2::Block* planes,c2::Block* temp,bool affine,bool packed,unsigned extra=0){
    if(extra==3)c2::rs4TransposeCached(p,input,output,planes,temp);
    else if(extra==1){if(affine)c2::rs4Transpose<true,D,true,true>(p,input,output,planes,temp);else c2::rs4Transpose<false,D,true,true>(p,input,output,planes,temp);}
    else if(extra==2){if(affine)c2::rs4Transpose<true,D,true,false,true>(p,input,output,planes,temp);else c2::rs4Transpose<false,D,true,false,true>(p,input,output,planes,temp);}
    else if(packed){if(affine)c2::rs4Transpose<true,D,true>(p,input,output,planes,temp);else c2::rs4Transpose<false,D,true>(p,input,output,planes,temp);}
    else if(affine)c2::rs4Transpose<true>(p,input,output,planes,temp);else c2::rs4Transpose<false>(p,input,output,planes,temp);
}
template<unsigned D>static void checkProfile(bool large){
    const std::vector<std::size_t> sizes=large?std::vector<std::size_t>{65536,262144,1048576}:
        std::vector<std::size_t>{1,D-1,D,D+1,4*D-1,4*D,4*D+1,64*D,64*D+1};
    const c2::Block guard{0x235562166ULL,0x91afe8171ULL};std::mt19937_64 random(777);
    for(auto k:sizes){
        c2::RateSpeed4Plan<D> p(k,113);
        for(std::size_t i=0;i!=p.inputN;++i){
            require(p.route[i]/p.groups==i/p.groups&&p.inverseRoute[p.route[i]]==i,"route regional bijection");
            const auto folded=16*(i/16)+4*(i%4)+(i%16)/4;
            require(p.foldedInverseRoute[p.route[i]]==folded&&p.packedInputDiagonal[folded]==p.inputDiagonal[i],"folded setup map");
        }
        std::vector<c2::Block> original(p.outputN),expected(k);for(auto& v:original)v={random(),random()};
        c2::BenchPayload input(p.outputN+8),out(k+8),planes(p.planeRecords()+8),temp(p.groups+8);
        if(large){std::memcpy(input.p,original.data(),p.outputN*sizeof(c2::Block));c2::rs4Transpose(p,input.p,out.p,planes.p,temp.p);std::copy(out.p,out.p+k,expected.begin());}
        else expected=c2::rs4OracleTranspose(p,original.data());
        for(auto stride:{p.groups,((p.groups+255)/256)*256})for(bool affine:{false,true})for(bool packed:{false,true})for(unsigned extra:{0,1,2,3})for(unsigned offset:{0,1,2,3}){
            p.planeStride=stride;for(auto* b:{&input,&out,&planes,&temp})std::fill(b->p,b->p+b->count,guard);
            std::memcpy(input.p+offset,original.data(),p.outputN*sizeof(c2::Block));
            execute4(p,input.p+offset,out.p+offset,planes.p+offset,temp.p+offset,affine,packed,extra);
            require(std::equal(expected.begin(),expected.end(),out.p+offset),"independent full transpose");
            require(std::equal(original.begin(),original.end(),input.p+offset),"input preservation");
            execute4(p,input.p+offset,input.p+offset,planes.p+offset,temp.p+offset,affine,packed,extra);
            require(std::equal(expected.begin(),expected.end(),input.p+offset),"alias prefix");
            require(std::equal(original.begin()+k,original.end(),input.p+offset+k),"alias suffix");
            for(auto* b:{&input,&out,&planes,&temp}){for(unsigned j=0;j!=offset;++j)require(b->p[j]==guard,"prefix guard");for(auto j=b->count-8+offset;j!=b->count;++j)require(b->p[j]==guard,"suffix guard");}
        }
        p.foldedInverseLocalRoute.clear();c2::rs4Transpose(p,original.data(),out.p,planes.p,temp.p);
        require(std::equal(expected.begin(),expected.end(),out.p),"wide inverse index fallback");
        if(!large&&k==4*D+1){std::vector<c2::Block> source(k);for(auto& v:source)v={random(),random()};const auto encoded=c2::rs4OracleForward(p,source.data());require(c2::p8OracleDot(source.data(),expected.data(),k)==c2::p8OracleDot(encoded.data(),original.data(),p.outputN),"forward-transpose adjoint");}
        std::cout<<"D="<<D<<" K="<<k<<" oracle/schedules/aliases/guards pass\n";
    }
}
static std::uint64_t checksum(const c2::Block* data,std::size_t count){std::uint64_t h=0;for(std::size_t i=0;i!=count;++i){h=(h^data[i].lo)*0x100000001b3ULL;h=(h^data[i].hi)*0x100000001b3ULL;}return h;}
static void stats(const std::string& name,std::vector<double> values){std::sort(values.begin(),values.end());std::cout<<name<<','<<std::fixed<<std::setprecision(3)<<values[values.size()/2]<<','<<values[values.size()/10]<<','<<values[9*values.size()/10]<<'\n';}
static void bench(std::size_t k,unsigned calls,bool phases,bool fresh=false){
    c2::RateClimb8Plan<c2::R8Periodic84> old(k,113);old.planeStride=((old.groups+255)/256)*256;
    c2::RateSpeed4Plan<40> code(k,113);c2::RateSpeed4Plan<39> lower(k,113);
    const auto count=std::max({old.outputN,code.outputN,lower.outputN});
    c2::BenchPayload input(count),planes(std::max({old.planeRecords(),code.planeRecords(),lower.planeRecords()})),temp(std::max({old.groups,code.groups,lower.groups}));
    std::vector<c2::Block> original(count);std::mt19937_64 random(913);for(auto& v:original)v={random(),random()};
    constexpr unsigned candidates=12;
    const auto inner=[&](unsigned i){switch(i){case 0:c2::r8InnerDispatch(old,input.p,planes.p,temp.p,c2::Subfield16Kernel::RegionAllNonTemporal,true);break;case 1:c2::rs4RegionInner<false>(code,input.p,planes.p,temp.p);break;case 2:c2::rs4RegionInner<true>(code,input.p,planes.p,temp.p);break;case 3:c2::rs4RegionInner<false>(lower,input.p,planes.p,temp.p);break;case 4:c2::rs4RegionInner<true>(lower,input.p,planes.p,temp.p);break;case 5:c2::rs4RegionInner<false,40,true,true>(code,input.p,planes.p,temp.p);break;case 6:c2::rs4RegionInner<true,40,true,true>(code,input.p,planes.p,temp.p);break;case 7:c2::rs4RegionInner<false,40,true,true,false,true>(code,input.p,planes.p,temp.p);break;case 8:c2::rs4RegionInner<true,40,true,true,true>(code,input.p,planes.p,temp.p);break;case 9:c2::rs4RegionInner<true,40,true,true,false,true>(code,input.p,planes.p,temp.p);break;case 10:c2::rs4RegionInner<true,40,false,true>(code,input.p,planes.p,temp.p);break;default:c2::rs4RegionInner<true,40,false,true,false,true>(code,input.p,planes.p,temp.p);break;}};
    const auto outer=[&](unsigned i){if(i==0)c2::r8OuterDispatchFinalStore(old,planes.p,input.p,c2::Subfield16Kernel::RegionAllNonTemporal,c2::Cubic8OuterVariant::FusedLoad);else if(i==3||i==4)c2::rs4Outer(lower,planes.p,input.p);else c2::rs4Outer(code,planes.p,input.p);};
    const auto restore=[&]{std::memcpy(input.p,original.data(),count*sizeof(c2::Block));};
    for(unsigned i=0;i!=candidates;++i){restore();for(unsigned j=0;j!=5;++j){inner(i);outer(i);}}
    std::array<std::vector<double>,candidates> total,inside,outside;std::array<std::uint64_t,candidates>sums{};
    for(unsigned i=0;i!=candidates;++i){total[i].resize(calls);inside[i].resize(calls);outside[i].resize(calls);}
    for(unsigned done=0,pass=0;done<calls;++pass){
        const auto batch=std::min(5U,calls-done);
        for(unsigned order=0;order!=candidates;++order){
            const auto i=(order+pass)%candidates;restore();
            for(unsigned j=0;j!=batch;++j){
                if(fresh&&j)restore();
                const auto a=Clock::now();
                if(phases){
                    inner(i);auto b=Clock::now();outer(i);auto c=Clock::now();
                    inside[i][done+j]=us(a,b);outside[i][done+j]=us(b,c);total[i][done+j]=us(a,c);
                }else{
                    if(i==11)c2::rs4TransposeCached(code,input.p,input.p,planes.p,temp.p);
                    else{inner(i);outer(i);}
                    total[i][done+j]=us(a,Clock::now());
                }
            }
            sums[i]=(sums[i]^checksum(input.p,k)^batch)*0x100000001b3ULL;
        }
        done+=batch;
    }
    require(sums[1]==sums[2]&&sums[3]==sums[4],"same-map benchmark checksum");
    for(unsigned i=5;i!=candidates;++i)require(sums[1]==sums[i],"same-map batch/route checksum");
    std::cout<<"K,"<<k<<",N84,"<<old.outputN<<",N40,"<<code.outputN<<",N39,"<<lower.outputN<<",same_input,"<<static_cast<void*>(input.p)<<",same_scratch,"<<static_cast<void*>(planes.p)<<",calls,"<<calls<<",setup_excluded,1,restore_each_call,"<<fresh<<"\nname,median_us,p10_us,p90_us\n";
    const char* names[candidates]={"frozen84","outer64-d40","outer64-d40-affine","outer64-d39","outer64-d39-affine","outer64-d40-packed","outer64-d40-affine-packed","outer64-d40-eight","outer64-d40-scatter","outer64-d40-affine-eight","outer64-d40-cached","outer64-d40-cached-eight"};
    for(unsigned i=0;i!=candidates;++i){stats(names[i],total[i]);if(phases){stats(std::string(names[i])+"-inner-route",inside[i]);stats(std::string(names[i])+"-outer",outside[i]);}}
}
int main(int argc,char** argv){try{
    if(!c2::gfni512Available())return 77;
    if(argc==2&&std::string(argv[1])=="--check"){checkPacket();checkProfile<40>(false);checkProfile<39>(false);return 0;}
    if(argc==2&&std::string(argv[1])=="--large-check"){checkProfile<40>(true);checkProfile<39>(true);return 0;}
    if(argc<2||argc>4)throw std::invalid_argument("usage: K [calls=51] [full|phases|fresh]");
    const auto k=std::stoull(argv[1]),calls=argc>2?std::stoull(argv[2]):51;const std::string mode=argc>3?argv[3]:"full";
    if(!calls||calls>1000000||(mode!="full"&&mode!="phases"&&mode!="fresh"))throw std::invalid_argument("calls/mode");bench(k,static_cast<unsigned>(calls),mode=="phases",mode=="fresh");return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
