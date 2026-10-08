#include "RateClimb8Code.h"
#include "RateClimb8Oracle.h"
#include "RateClimb8Outer.h"
#include "Punctured8Code.h"
#include "BenchPayload.h"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>

namespace c2=spin::research::char2;
using Clock=std::chrono::steady_clock;
static constexpr auto gather=c2::Cubic8OuterVariant::Gather,fused=c2::Cubic8OuterVariant::FusedLoad;
static constexpr c2::Subfield16Kernel schedules[]={c2::Subfield16Kernel::Scatter,
    c2::Subfield16Kernel::RegionNonTemporal,c2::Subfield16Kernel::RegionAllNonTemporal};
static void require(bool good,const char* message){if(!good)throw std::runtime_error(message);}
static double micros(Clock::time_point a,Clock::time_point b){return std::chrono::duration<double,std::micro>(b-a).count();}

template<class Profile> static void finishCandidate(const c2::RateClimb8Plan<Profile>& p,
    const c2::Block* planes,c2::Block* output,c2::Subfield16Kernel mode,c2::Cubic8OuterVariant variant,bool finalStore){
    if(finalStore)c2::r8OuterDispatchFinalStore(p,planes,output,mode,variant);
    else c2::r8OuterDispatch(p,planes,output,mode,variant);
}
template<class Profile> static void encodeCandidate(const c2::RateClimb8Plan<Profile>& p,
    const c2::Block* input,c2::Block* output,c2::Block* planes,c2::Block* temp,
    c2::Subfield16Kernel mode,c2::Cubic8OuterVariant variant,bool folded,bool finalStore){
    c2::r8InnerDispatch(p,input,planes,temp,mode,folded);
    finishCandidate(p,planes,output,mode,variant,finalStore);
}

template<class Profile> static void checkPacket(){
    using Plan=c2::RateClimb8Plan<Profile>;
    constexpr unsigned records=Plan::tileRecords;
    std::mt19937_64 random(0x8988+Plan::dimension);
    for(unsigned trial=0;trial!=16*(records+8)+8;++trial){
        alignas(64) c2::Block input[records]{},state[8]{},expected[32]{},expectedState[8]{},actual[32],actualState[8];
        if(trial<16*(records+8)){
            const auto index=trial/16,base=16*(trial%16);
            auto* target=index<records?input+index:state+index-records;
            for(unsigned byte=0;byte!=16;++byte)c2::p8OracleBytes(target)[byte]=static_cast<std::uint8_t>(base+byte);
        }else{
            for(auto& v:input)v={random(),random()};
            for(auto& v:state)v={random(),random()};
        }
        std::uint8_t dx[8],ds[32];std::fill(std::begin(dx),std::end(dx),1);
        for(auto& v:ds)v=static_cast<std::uint8_t>(1+random()%255);
        std::copy(std::begin(state),std::end(state),expectedState);
        // Independent compact serialization: enumerate coordinates then phases.
        unsigned index[4][8]{};unsigned next=0;
        for(unsigned j=0;j!=8;++j)for(unsigned phase=0;phase!=4;++phase)
            if(j<Plan::retained(phase))index[phase][j]=next++;
        require(next==records,"compact packet layout size");
        for(unsigned phase=4;phase--;){
            c2::Block in[8]{},previous[8];
            for(unsigned j=0;j!=Plan::retained(phase);++j)in[j]=input[index[phase][j]];
            c2::p8OraclePacket(in,expectedState,dx,ds+8*phase,expected+8*phase,previous,true);
            std::copy(std::begin(previous),std::end(previous),expectedState);
        }
        auto h0=_mm512_load_si512(state),h1=_mm512_load_si512(state+4);
        c2::r8InnerTransposeFour<Profile>(input,actual,h0,h1,ds);
        _mm512_store_si512(actualState,h0);_mm512_store_si512(actualState+4,h1);
        require(std::equal(std::begin(expected),std::end(expected),actual),"all-byte packet output");
        require(std::equal(std::begin(expectedState),std::end(expectedState),actualState),"all-byte packet state");
        h0=_mm512_load_si512(state);h1=_mm512_load_si512(state+4);
        c2::r8InnerTransposeFour<Profile,true>(input,actual,h0,h1,ds);
        _mm512_store_si512(actualState,h0);_mm512_store_si512(actualState+4,h1);
        for(unsigned coordinate=0;coordinate!=8;++coordinate)for(unsigned phase=0;phase!=4;++phase)
            require(actual[4*coordinate+phase]==expected[8*phase+coordinate],"folded packet independent layout");
        require(std::equal(std::begin(expectedState),std::end(expectedState),actualState),"folded packet state");
    }
}

template<class Profile> static void checkProfile(){
    using Plan=c2::RateClimb8Plan<Profile>;
    checkPacket<Profile>();std::mt19937_64 random(193);
    const c2::Block guard{0x9278361,0xfad89712};
    const std::size_t d=Plan::dimension;
    for(std::size_t k:{std::size_t(1),d-1,d,d+1,4*d-1,4*d,4*d+1,64*d,64*d+1}){
        Plan p(k,113),padded(p.k,113);
        require(p.outputN==padded.outputN&&p.route==padded.route&&p.inputDiagonal==padded.inputDiagonal&&
            p.stateDiagonal==padded.stateDiagonal,"padded setup identity");
        c2::r8OracleValidatePlan(p);
        std::vector<c2::Block> original(p.outputN);for(auto& v:original)v={random(),random()};
        const auto expected=c2::r8OracleTranspose(p,original.data());
        c2::BenchPayload input(p.outputN+8,false),out(k+8,false),planes(128*(p.groups+4)+8,false),temp(p.groups+8,false);
        for(auto stride:{p.groups,p.groups+4})for(auto variant:{gather,fused})for(auto mode:schedules)
          for(bool folded:{false,true})for(bool finalStore:{false,true})for(unsigned offset:{0,1,2,3}){
            p.planeStride=stride;
            for(auto* b:{&input,&out,&planes,&temp})std::fill(b->p,b->p+b->count,guard);
            std::memcpy(input.p+offset,original.data(),p.outputN*sizeof(c2::Block));
            encodeCandidate(p,input.p+offset,out.p+offset,planes.p+offset,temp.p+offset,mode,variant,folded,finalStore);
            require(std::equal(expected.begin(),expected.end(),out.p+offset),"full independent scalar transpose");
            require(std::equal(original.begin(),original.end(),input.p+offset),"input preservation");
            encodeCandidate(p,input.p+offset,input.p+offset,planes.p+offset,temp.p+offset,mode,variant,folded,finalStore);
            require(std::equal(expected.begin(),expected.end(),input.p+offset),"alias prefix");
            require(std::equal(original.begin()+k,original.end(),input.p+offset+k),"alias suffix");
            for(auto* b:{&input,&out,&planes,&temp}){
                for(unsigned j=0;j!=offset;++j)require(b->p[j]==guard,"prefix guard");
                for(std::size_t j=b->count-8+offset;j!=b->count;++j)require(b->p[j]==guard,"suffix guard");
            }
        }
        p.inverseLocalRoute.clear();p.foldedInverseLocalRoute.clear();
        for(bool folded:{false,true}){
            c2::r8Transpose(p,original.data(),out.p,planes.p,temp.p,c2::Subfield16Kernel::RegionNonTemporal,gather,folded);
            require(std::equal(expected.begin(),expected.end(),out.p),"uint32 route fallback");
        }
        if(k==4*d+1||k==64*d+1){
            std::vector<c2::Block> source(k);for(auto& v:source)v={random(),random()};
            const auto image=c2::r8OracleForward(p,source.data());
            require(c2::p8OracleDot(image.data(),original.data(),p.outputN)==
                c2::p8OracleDot(source.data(),expected.data(),k),"shortened full field adjoint");
        }
    }
    std::cout<<"profile_d"<<d<<"_tile"<<Plan::tileRecords<<" independent scalar, packet, adjoint, alias, guard, unaligned checks pass\n";
}

template<class Profile> static void checkLargeProfile(){
    using Plan=c2::RateClimb8Plan<Profile>;std::mt19937_64 random(7767);
    for(unsigned exponent:{16,18,20}){
        Plan p(std::size_t(1)<<exponent,113);
        c2::BenchPayload input(p.outputN),out(p.requested),expected(p.requested),planes(128*(p.groups+256)),temp(p.groups),alias(p.outputN);
        for(std::size_t i=0;i!=p.outputN;++i)input.p[i]={random(),random()};
        c2::r8Transpose(p,input.p,expected.p,planes.p,temp.p,c2::Subfield16Kernel::Scatter,gather);
        for(auto stride:{p.groups,p.groups+4,((p.groups+255)/256)*256})for(auto variant:{gather,fused})for(auto mode:schedules)
          for(bool folded:{false,true})for(bool finalStore:{false,true}){
            p.planeStride=stride;
            encodeCandidate(p,input.p,out.p,planes.p,temp.p,mode,variant,folded,finalStore);
            require(std::equal(expected.p,expected.p+p.requested,out.p),"large schedule agreement");
            std::memcpy(alias.p,input.p,p.outputN*sizeof(c2::Block));
            encodeCandidate(p,alias.p,alias.p,planes.p,temp.p,mode,variant,folded,finalStore);
            require(std::equal(expected.p,expected.p+p.requested,alias.p),"large alias prefix");
            require(std::equal(input.p+p.requested,input.p+p.outputN,alias.p+p.requested),"large alias suffix");
        }
        if constexpr(Plan::dimension==84&&Plan::tileRecords==28){
            c2::Punctured8Plan frozen(p.requested,113);
            require(frozen.route==p.route&&frozen.inputDiagonal==p.inputDiagonal&&frozen.stateDiagonal==p.stateDiagonal,"frozen setup identity");
            c2::p8Transpose(frozen,input.p,out.p,planes.p,temp.p,c2::Subfield16Kernel::RegionAllNonTemporal,gather);
            require(std::equal(expected.p,expected.p+p.requested,out.p),"frozen75 full map identity");
        }
        std::cout<<"K"<<exponent<<" d"<<Plan::dimension<<" tile"<<Plan::tileRecords<<" schedule/alias checks pass\n";
    }
}

static std::uint64_t checksum(const c2::Block* values,std::size_t n){
    std::uint64_t result=0;
    for(std::size_t i=0;i!=n;++i){result=(result^values[i].lo)*0x100000001b3ULL;result=(result^values[i].hi)*0x100000001b3ULL;}
    return result;
}
static void stats(const std::string& name,std::size_t k,std::size_t padded,std::size_t raw,std::size_t n,std::vector<double> values,std::uint64_t sum){
    std::sort(values.begin(),values.end());
    std::cout<<name<<','<<k<<','<<padded<<','<<raw<<','<<n<<','<<values.size()<<','<<std::fixed<<std::setprecision(4)
        <<values[values.size()/2]<<','<<values[values.size()/10]<<','<<values[9*values.size()/10]
        <<','<<std::hex<<sum<<std::dec<<'\n';
}

static void benchmark(std::size_t k,unsigned calls,bool phases,c2::Subfield16Kernel mode,bool padded,bool pageAligned,bool folded,bool finalStore){
    c2::Punctured8Plan frozen(k,113);
    c2::RateClimb8Plan<c2::R8Fixed7> control(k,113);
    c2::RateClimb8Plan<c2::R8Fixed6> fixed(k,113);
    c2::RateClimb8Plan<c2::R8Periodic84> periodic84(k,113);
    c2::RateClimb8Plan<c2::R8Periodic74> periodic74(k,113);
    if(padded){control.planeStride+=4;fixed.planeStride+=4;periodic84.planeStride+=4;periodic74.planeStride+=4;}
    if(pageAligned){
        control.planeStride=((control.groups+255)/256)*256;fixed.planeStride=((fixed.groups+255)/256)*256;
        periodic84.planeStride=((periodic84.groups+255)/256)*256;periodic74.planeStride=((periodic74.groups+255)/256)*256;
    }
    const auto records=std::max({frozen.outputN,fixed.outputN,periodic84.outputN,periodic74.outputN});
    const auto rows=std::max({frozen.groups,fixed.groups,periodic84.groups,periodic74.groups});
    const auto stride=std::max({frozen.planeStride,control.planeStride,fixed.planeStride,periodic84.planeStride,periodic74.planeStride});
    c2::BenchPayload input(records,k>65536),planes(128*stride,k>65536),temp(rows,k>65536);
    std::vector<c2::Block> original(records);std::mt19937_64 random(913);
    for(auto& value:original)value={random(),random()};
    const auto restore=[&]{std::memcpy(input.p,original.data(),records*sizeof(c2::Block));};
    const auto runInner=[&](unsigned candidate){
        switch(candidate){
        case 0:c2::p8InnerDispatch(frozen,input.p,planes.p,temp.p,mode);break;
        case 1:c2::r8InnerDispatch(control,input.p,planes.p,temp.p,mode,folded);break;
        case 2:case 7:c2::r8InnerDispatch(fixed,input.p,planes.p,temp.p,mode,folded);break;
        case 3:case 5:c2::r8InnerDispatch(periodic84,input.p,planes.p,temp.p,mode,folded);break;
        default:c2::r8InnerDispatch(periodic74,input.p,planes.p,temp.p,mode,folded);break;
        }
    };
    const auto runOuter=[&](unsigned candidate){
        const auto variant=candidate>=5?fused:gather;
        switch(candidate){
        case 0:c2::p8OuterDispatch(frozen,planes.p,input.p,mode,gather);break;
        case 1:finishCandidate(control,planes.p,input.p,mode,gather,finalStore);break;
        case 2:case 7:finishCandidate(fixed,planes.p,input.p,mode,variant,finalStore);break;
        case 3:case 5:finishCandidate(periodic84,planes.p,input.p,mode,variant,finalStore);break;
        default:finishCandidate(periodic74,planes.p,input.p,mode,variant,finalStore);break;
        }
    };
    const auto encode=[&](unsigned candidate){
        const auto variant=candidate>=5?fused:gather;
        switch(candidate){
        case 0:c2::p8Transpose(frozen,input.p,input.p,planes.p,temp.p,mode,gather);break;
        case 1:encodeCandidate(control,input.p,input.p,planes.p,temp.p,mode,gather,folded,finalStore);break;
        case 2:case 7:encodeCandidate(fixed,input.p,input.p,planes.p,temp.p,mode,variant,folded,finalStore);break;
        case 3:case 5:encodeCandidate(periodic84,input.p,input.p,planes.p,temp.p,mode,variant,folded,finalStore);break;
        default:encodeCandidate(periodic74,input.p,input.p,planes.p,temp.p,mode,variant,folded,finalStore);break;
        }
    };
    for(unsigned candidate=0;candidate!=8;++candidate){restore();for(unsigned i=0;i!=5;++i)encode(candidate);}
    std::array<std::vector<double>,8> samples,inner,outer;std::array<std::uint64_t,8> sums{};
    for(unsigned candidate=0;candidate!=8;++candidate){samples[candidate].resize(calls);
        if(phases){inner[candidate].resize(calls);outer[candidate].resize(calls);}}
    for(unsigned done=0,pass=0;done<calls;++pass){
        const unsigned count=std::min(5U,calls-done);
        for(unsigned order=0;order!=8;++order){
            const unsigned candidate=(order+pass)%8;restore();
            for(unsigned i=0;i!=count;++i){
                const auto a=Clock::now();
                if(phases){
                    runInner(candidate);const auto b=Clock::now();
                    runOuter(candidate);const auto c=Clock::now();
                    inner[candidate][done+i]=micros(a,b);outer[candidate][done+i]=micros(b,c);samples[candidate][done+i]=micros(a,c);
                }else{encode(candidate);samples[candidate][done+i]=micros(a,Clock::now());}
            }
            sums[candidate]=(sums[candidate]^checksum(input.p,k)^count)*0x100000001b3ULL;
        }
        done+=count;
    }
    require(sums[0]==sums[1],"frozen75 benchmark map identity");
    require(sums[2]==sums[7]&&sums[3]==sums[5]&&sums[4]==sums[6],"exact outer schedule benchmark identity");
    const char* names[]={"frozen75","new75-control","fixed-six","periodic84","periodic74","periodic84-fused","periodic74-fused","fixed-six-fused"};
    const std::size_t dims[]={frozen.k,control.k,fixed.k,periodic84.k,periodic74.k,periodic84.k,periodic74.k,fixed.k};
    const std::size_t raw[]={frozen.inputN,control.inputN,fixed.inputN,periodic84.inputN,periodic74.inputN,periodic84.inputN,periodic74.inputN,fixed.inputN};
    const std::size_t n[]={frozen.outputN,control.outputN,fixed.outputN,periodic84.outputN,periodic74.outputN,periodic84.outputN,periodic74.outputN,fixed.outputN};
    std::cout<<"same_input_address,"<<static_cast<void*>(input.p)<<",same_scratch_address,"<<static_cast<void*>(planes.p)
        <<",setup_excluded,1,restore_outside_timing,batch5,rotating_order,different_codes_and_rates,1,new_column_padding,"<<padded
        <<",page_aligned_columns,"<<pageAligned<<",folded_route,"<<folded<<",final_store,"<<finalStore<<",schedule,"<<static_cast<unsigned>(mode)
        <<"\nname,K,Kpad,Nraw,Nout,calls,median_us,p10_us,p90_us,checksum\n";
    for(unsigned i=0;i!=8;++i){
        stats(names[i],k,dims[i],raw[i],n[i],samples[i],sums[i]);
        if(phases){stats(std::string(names[i])+"-inner-route",k,dims[i],raw[i],n[i],inner[i],0);
            stats(std::string(names[i])+"-outer",k,dims[i],raw[i],n[i],outer[i],0);}
    }
}

int main(int argc,char** argv){try{
    if(!c2::gfni512Available())return 77;
    if(argc==2&&std::string(argv[1])=="--check"){
        c2::r8OracleCheckLocal();checkProfile<c2::R8Fixed7>();checkProfile<c2::R8Fixed6>();
        checkProfile<c2::R8Periodic84>();checkProfile<c2::R8Periodic74>();return 0;
    }
    if(argc==2&&std::string(argv[1])=="--large-check"){
        checkLargeProfile<c2::R8Fixed7>();checkLargeProfile<c2::R8Fixed6>();
        checkLargeProfile<c2::R8Periodic84>();checkLargeProfile<c2::R8Periodic74>();return 0;
    }
    if(argc<2||argc>8)throw std::invalid_argument("usage: K [calls=51] [full|phases] [auto|scatter|region|region-all] [tight|padded|page] [packet|folded] [materialized|final-store]");
    const auto k=std::stoull(argv[1]),calls=argc>2?std::stoull(argv[2]):51;
    const std::string mode=argc>3?argv[3]:"full",route=argc>4?argv[4]:"auto";
    const std::string layout=argc>5?argv[5]:"tight";
    const std::string folding=argc>6?argv[6]:"packet",finish=argc>7?argv[7]:"materialized";
    if(!calls||calls>1000000||(mode!="full"&&mode!="phases"))throw std::invalid_argument("calls or mode");
    if(layout!="tight"&&layout!="padded"&&layout!="page")throw std::invalid_argument("column layout");
    if(folding!="packet"&&folding!="folded")throw std::invalid_argument("routing layout");
    if(finish!="materialized"&&finish!="final-store")throw std::invalid_argument("outer finish");
    auto kernel=k>=1048576?c2::Subfield16Kernel::RegionAllNonTemporal:c2::Subfield16Kernel::Scatter;
    if(route=="scatter")kernel=c2::Subfield16Kernel::Scatter;
    else if(route=="region")kernel=c2::Subfield16Kernel::RegionNonTemporal;
    else if(route=="region-all")kernel=c2::Subfield16Kernel::RegionAllNonTemporal;
    else if(route!="auto")throw std::invalid_argument("route schedule");
    benchmark(k,static_cast<unsigned>(calls),mode=="phases",kernel,layout=="padded",layout=="page",folding=="folded",finish=="final-store");return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
