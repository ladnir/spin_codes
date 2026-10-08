#include "Cubic8Code.h"
#include "Cauchy128CoefficientCode.h"
#include "BenchPayload.h"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>
#ifndef SPIN_CAUCHY_PUBLIC
#define SPIN_CAUCHY_PUBLIC 0
#endif
#if SPIN_CAUCHY_PUBLIC
#include <spin/Code.h>
#endif

namespace c2=spin::research::char2;
using Clock=std::chrono::steady_clock;
static void require(bool good,const char* message){if(!good)throw std::runtime_error(message);}
static double micros(Clock::time_point a,Clock::time_point b){return std::chrono::duration<double,std::micro>(b-a).count();}
static constexpr auto gather=c2::Cubic8OuterVariant::Gather;
static constexpr auto fused=c2::Cubic8OuterVariant::FusedLoad;
static constexpr auto compact=c2::Cauchy128Variant::CompactBatch;
static constexpr c2::Subfield16Kernel schedules[]={c2::Subfield16Kernel::Scatter,
    c2::Subfield16Kernel::RegionNonTemporal,c2::Subfield16Kernel::RegionAllNonTemporal};

static std::uint8_t scalarMultiply(unsigned a,unsigned b){
    unsigned r=0;
    while(b){if(b&1)r^=a;b>>=1;a<<=1;if(a&256)a^=0x11b;}
    return static_cast<std::uint8_t>(r);
}
struct CubicOracle {
    std::array<std::array<std::uint8_t,256>,256> mul{};
    std::array<std::array<std::uint8_t,83>,128> evaluation{};
    std::uint8_t a[8][8][4]{},c[8][4][8]{};
    CubicOracle(){
        for(unsigned x=0;x!=256;++x)for(unsigned y=0;y!=256;++y)mul[x][y]=scalarMultiply(x,y);
        const unsigned labels[8][8]={
            {0,1,188,189,92,93,224,225},{176,177,12,13,236,237,80,81},
            {2,3,190,191,94,95,226,227},{178,179,14,15,238,239,82,83},
            {6,7,186,187,90,91,230,231},{182,183,10,11,234,235,86,87},
            {4,5,184,185,88,89,228,229},{180,181,8,9,232,233,84,85}};
        const unsigned nodes[4]={16,17,172,173};
        for(unsigned p=0;p!=8;++p)for(unsigned i=0;i!=8;++i)for(unsigned j=0;j!=4;++j){
            const auto x=labels[p][i]^nodes[j];const auto cube=mul[mul[x][x]][x];
            a[p][i][j]=cube^1;c[p][j][i]=cube^3;
        }
        for(unsigned point=0;point!=128;++point)for(unsigned j=0;j!=83;++j){
            std::uint8_t value=1;
            for(unsigned bit=0;bit!=7;++bit)if(j&(1U<<bit)){
                std::uint8_t numerator=1,denominator=1;
                for(unsigned label=0;label!=(1U<<bit);++label){
                    numerator=mul[numerator][point^label];denominator=mul[denominator][(1U<<bit)^label];
                }
                unsigned inverse=1;while(mul[denominator][inverse]!=1)++inverse;
                value=mul[value][mul[numerator][inverse]];
            }
            evaluation[point][j]=value;
        }
    }
};
static const CubicOracle& oracle(){static const CubicOracle o;return o;}
static std::uint8_t* bytes(c2::Block* p){return reinterpret_cast<std::uint8_t*>(p);}
static const std::uint8_t* bytes(const c2::Block* p){return reinterpret_cast<const std::uint8_t*>(p);}

static void packetOracle(const c2::Block* input,const c2::Block* state,const std::uint8_t* diagonal,
                          c2::Block* output,c2::Block* following,unsigned phase,bool transpose){
    const auto& o=oracle();
    for(unsigned b=0;b!=16;++b){
        std::uint8_t r[4];
        for(unsigned j=0;j!=4;++j)r[j]=transpose?bytes(state+j)[b]:o.mul[diagonal[j]][bytes(state+j)[b]];
        for(unsigned i=0;i!=8;++i){
            auto value=bytes(input+i)[b];
            for(unsigned j=0;j!=4;++j)value^=o.mul[transpose?o.c[phase][j][i]:o.mul[2][o.a[phase][i][j]]][r[j]];
            bytes(output+i)[b]=value;
        }
        for(unsigned j=0;j!=4;++j){
            auto value=r[j];
            for(unsigned i=0;i!=8;++i)value^=o.mul[transpose?o.mul[2][o.a[phase][i][j]]:o.c[phase][j][i]][bytes(input+i)[b]];
            bytes(following+j)[b]=transpose?o.mul[diagonal[j]][value]:value;
        }
    }
}
template<bool Transpose>
static void packetNative(unsigned phase,const __m512i* in,__m512i* out,__m512i& state,__m512i diagonal){
#define C8_TEST_PHASE(P) case P: if constexpr(Transpose)c2::c8InnerTransposePacked<P>(in,out,state,diagonal); \
    else c2::c8InnerForwardPacked<P>(in,out,state,diagonal);break
    switch(phase){C8_TEST_PHASE(0);C8_TEST_PHASE(1);C8_TEST_PHASE(2);C8_TEST_PHASE(3);
                  C8_TEST_PHASE(4);C8_TEST_PHASE(5);C8_TEST_PHASE(6);C8_TEST_PHASE(7);}
#undef C8_TEST_PHASE
}
static void checkInner(){
    std::mt19937_64 rng(843);const auto& o=oracle();
    alignas(64) c2::Block input[8],state[4],expected[8],expectedState[4],actual[8],actualState[4];
    for(unsigned phase=0;phase!=8;++phase){
        for(unsigned i=0;i!=8;++i)for(unsigned j=0;j!=4;++j){
            require(c2::c8Coefficient(phase,i,j,false)==o.a[phase][i][j],"emission geometry");
            require(c2::c8Coefficient(phase,i,j,true)==o.c[phase][j][i],"injection geometry");
        }
        for(unsigned trial=0;trial!=200;++trial){
            std::fill(std::begin(input),std::end(input),c2::Block{});
            std::fill(std::begin(state),std::end(state),c2::Block{});
            if(trial<192){
                const auto column=trial/16,base=16*(trial%16);
                auto* target=column<8?input+column:state+column-8;
                for(unsigned b=0;b!=16;++b)bytes(target)[b]=static_cast<std::uint8_t>(base+b);
            }else{
                for(auto& value:input)value={rng(),rng()};for(auto& value:state)value={rng(),rng()};
            }
            std::uint8_t diagonal[4];for(auto& value:diagonal)value=static_cast<std::uint8_t>(1+rng()%255);
            for(bool transpose:{false,true})for(bool alias:{false,true}){
                packetOracle(input,state,diagonal,expected,expectedState,phase,transpose);
                alignas(64) __m512i vin[2]={_mm512_load_si512(input),_mm512_load_si512(input+4)},vout[2];
                auto h=_mm512_load_si512(state);auto* destination=alias?vin:vout;
                if(transpose)packetNative<true>(phase,vin,destination,h,c2::sfBroadcastFour(diagonal));
                else packetNative<false>(phase,vin,destination,h,c2::sfBroadcastFour(diagonal));
                _mm512_store_si512(actual,destination[0]);_mm512_store_si512(actual+4,destination[1]);
                _mm512_store_si512(actualState,h);
                require(std::equal(std::begin(expected),std::end(expected),actual),"packed packet scalar output");
                require(std::equal(std::begin(expectedState),std::end(expectedState),actualState),"packed packet scalar state");
                if(!alias)require(std::memcmp(vin,input,sizeof(input))==0,"packet input preservation");
            }
        }
    }
}
static void pairedPacketNative(unsigned firstPhase,const __m512i* in,__m512i* out,
                                 __m512i& state,const __m512i* diagonal){
    switch(firstPhase){
    case 7:c2::c8InnerTransposePairedPacked<7,6>(in,out,state,diagonal);break;
    case 5:c2::c8InnerTransposePairedPacked<5,4>(in,out,state,diagonal);break;
    case 3:c2::c8InnerTransposePairedPacked<3,2>(in,out,state,diagonal);break;
    case 1:c2::c8InnerTransposePairedPacked<1,0>(in,out,state,diagonal);break;
    default:throw std::invalid_argument("paired phase");
    }
}
static void checkPairedInner(){
    std::mt19937_64 rng(8476);
    alignas(64) c2::Block input[16],state[4],expected[16],middle[4],expectedState[4],actual[16],actualState[4];
    for(unsigned firstPhase:{7,5,3,1})for(unsigned trial=0;trial!=328;++trial){
        std::fill(std::begin(input),std::end(input),c2::Block{});
        std::fill(std::begin(state),std::end(state),c2::Block{});
        if(trial<320){
            const auto column=trial/16,base=16*(trial%16);
            auto* target=column<16?input+column:state+column-16;
            for(unsigned b=0;b!=16;++b)bytes(target)[b]=static_cast<std::uint8_t>(base+b);
        }else{
            for(auto& value:input)value={rng(),rng()};for(auto& value:state)value={rng(),rng()};
        }
        std::uint8_t diagonal[8];for(auto& value:diagonal)value=static_cast<std::uint8_t>(1+rng()%255);
        packetOracle(input,state,diagonal,expected,middle,firstPhase,true);
        packetOracle(input+8,middle,diagonal+4,expected+8,expectedState,firstPhase-1,true);
        const __m512i d[2]={c2::sfBroadcastFour(diagonal),c2::sfBroadcastFour(diagonal+4)};
        for(bool alias:{false,true}){
            alignas(64) __m512i vin[4]={_mm512_load_si512(input),_mm512_load_si512(input+4),
                                       _mm512_load_si512(input+8),_mm512_load_si512(input+12)},vout[4];
            auto h=_mm512_load_si512(state);auto* destination=alias?vin:vout;
            pairedPacketNative(firstPhase,vin,destination,h,d);
            for(unsigned j=0;j!=4;++j)_mm512_store_si512(actual+4*j,destination[j]);
            _mm512_store_si512(actualState,h);
            require(std::equal(std::begin(expected),std::end(expected),actual),"paired packed scalar output");
            require(std::equal(std::begin(expectedState),std::end(expectedState),actualState),"paired packed scalar state");
            if(!alias)require(std::memcmp(vin,input,sizeof(input))==0,"paired packet input preservation");
        }
    }
}
static void checkOuter(){
    const auto& o=oracle();std::mt19937_64 rng(12883);
    alignas(64) std::uint8_t input[128][64],expected[128][64];
    alignas(64) __m512i values[128];
    for(unsigned trial=0;trial!=132;++trial){
        std::memset(input,0,sizeof(input));
        if(trial<128)for(unsigned b=0;b!=64;++b)input[trial][b]=static_cast<std::uint8_t>(b+1);
        else for(auto& row:input)for(auto& value:row)value=static_cast<std::uint8_t>(rng());
        for(bool transpose:{false,true}){
            std::memset(expected,0,sizeof(expected));const unsigned outputs=transpose?83:128;
            for(unsigned j=0;j!=outputs;++j)for(unsigned b=0;b!=64;++b){
                auto& value=expected[j][b];
                if(transpose)for(unsigned i=0;i!=128;++i)value^=o.mul[o.evaluation[i][j]][input[i][b]];
                else for(unsigned i=0;i!=83;++i)value^=o.mul[o.evaluation[j][i]][input[i][b]];
            }
            std::memcpy(values,input,sizeof(input));
            if(transpose)c2::c128OuterCoefficientTransposeInplace(values);
            else{
                for(unsigned j=83;j!=128;++j)values[j]=_mm512_setzero_si512();
                c2::c128OuterCoefficientForwardInplace(values);
            }
            require(std::memcmp(values,expected,outputs*64)==0,"GRS12883 direct additive-basis outer");
            if(transpose)for(std::size_t stride:{64,83,192}){
                std::vector<std::uint8_t> strided(128*stride+67,0xb9);
                for(unsigned j=0;j!=128;++j)std::memcpy(strided.data()+3+j*stride,input[j],64);
                const auto original=strided;
                c2::c128OuterCoefficientTransposeStrided(strided.data()+3,stride,values);
                require(std::memcmp(values,expected,83*64)==0,"GRS12883 strided outer");
                require(strided==original,"strided outer source preservation");
            }
        }
    }
}

static std::vector<c2::Block> scalarTranspose(const c2::Cubic8Plan& p,const c2::Block* input){
    const auto& o=oracle();std::vector<c2::Block> planes(p.n),result(p.requested);
    c2::Block state[4]{},previous[4],packet[8];
    for(std::size_t step=p.n/8;step--;){
        packetOracle(input+8*step,state,p.stateDiagonal.data()+4*step,packet,previous,unsigned(step%8),true);
        for(unsigned j=0;j!=8;++j)for(unsigned b=0;b!=16;++b)
            bytes(planes.data()+p.route[8*step+j])[b]=o.mul[p.inputDiagonal[8*step+j]][bytes(packet+j)[b]];
        std::copy(std::begin(previous),std::end(previous),state);
    }
    for(std::size_t q=0;q!=p.requested;++q){
        const auto row=(q/332)*4+q%4;const auto coordinate=(q%332)/4;
        for(unsigned j=0;j!=128;++j)for(unsigned b=0;b!=16;++b)
            bytes(result.data()+q)[b]^=o.mul[o.evaluation[j][coordinate]][bytes(planes.data()+j*p.groups+row)[b]];
    }
    return result;
}
static std::vector<c2::Block> scalarForward(const c2::Cubic8Plan& p,const c2::Block* input){
    const auto& o=oracle();std::vector<c2::Block> planes(p.n),result(p.n);
    for(std::size_t row=0;row!=p.groups;++row)for(unsigned coordinate=0;coordinate!=83;++coordinate){
        const auto q=(row/4)*332+4*coordinate+row%4;if(q>=p.requested)continue;
        for(unsigned j=0;j!=128;++j)for(unsigned b=0;b!=16;++b)
            bytes(planes.data()+j*p.groups+row)[b]^=o.mul[o.evaluation[j][coordinate]][bytes(input+q)[b]];
    }
    c2::Block state[4]{},following[4],packet[8];
    for(std::size_t step=0;step!=p.n/8;++step){
        for(unsigned j=0;j!=8;++j)for(unsigned b=0;b!=16;++b)
            bytes(packet+j)[b]=o.mul[p.inputDiagonal[8*step+j]][bytes(planes.data()+p.route[8*step+j])[b]];
        packetOracle(packet,state,p.stateDiagonal.data()+4*step,result.data()+8*step,following,unsigned(step%8),false);
        std::copy(std::begin(following),std::end(following),state);
    }
    return result;
}
static std::array<std::uint8_t,16> scalarDot(const c2::Block* a,const c2::Block* b,std::size_t count){
    std::array<std::uint8_t,16> result{};const auto& o=oracle();
    for(std::size_t i=0;i!=count;++i)for(unsigned j=0;j!=16;++j)result[j]^=o.mul[bytes(a+i)[j]][bytes(b+i)[j]];
    return result;
}
static void check(){
    checkInner();checkPairedInner();checkOuter();std::mt19937_64 rng(173);
    const c2::Block guard{0x32456789ULL,0xfabc1234ULL};
    for(std::size_t k:{1,82,83,84,331,332,333,5312,5313,6000}){
        c2::Cubic8Plan p(k,113),full(p.k,113);
        require(p.n==full.n&&p.route==full.route&&p.inputDiagonal==full.inputDiagonal&&
                p.stateDiagonal==full.stateDiagonal,"full padded setup");
        std::vector<c2::Block> original(p.n);for(auto& value:original)value={rng(),rng()};
        const auto expected=scalarTranspose(p,original.data());
        c2::BenchPayload input(p.n+8,false),out(k+8,false),planes(p.n+8,false),temp(p.groups+8,false);
        for(bool paired:{false,true})for(auto selected:{gather,fused})for(auto mode:schedules)for(unsigned offset:{0,1,2,3}){
            for(auto* buffer:{&input,&out,&planes,&temp})std::fill(buffer->p,buffer->p+buffer->count,guard);
            std::memcpy(input.p+offset,original.data(),p.n*sizeof(c2::Block));
            if(paired)c2::c8TransposePaired(p,input.p+offset,out.p+offset,planes.p+offset,temp.p+offset,mode,selected);
            else c2::c8Transpose(p,input.p+offset,out.p+offset,planes.p+offset,temp.p+offset,mode,selected);
            require(std::equal(expected.begin(),expected.end(),out.p+offset),"full scalar transpose");
            require(std::equal(original.begin(),original.end(),input.p+offset),"full input preservation");
            if(paired)c2::c8TransposePaired(p,input.p+offset,input.p+offset,planes.p+offset,temp.p+offset,mode,selected);
            else c2::c8Transpose(p,input.p+offset,input.p+offset,planes.p+offset,temp.p+offset,mode,selected);
            require(std::equal(expected.begin(),expected.end(),input.p+offset),"full alias prefix");
            require(std::equal(original.begin()+k,original.end(),input.p+offset+k),"full alias tail");
            for(auto* buffer:{&input,&out,&planes,&temp}){
                for(unsigned j=0;j!=offset;++j)require(buffer->p[j]==guard,"full prefix guard");
                for(std::size_t j=buffer->count-8+offset;j!=buffer->count;++j)require(buffer->p[j]==guard,"full suffix guard");
            }
        }
        p.inverseLocalRoute.clear();
        c2::c8Transpose(p,original.data(),out.p,planes.p,temp.p,c2::Subfield16Kernel::RegionNonTemporal);
        require(std::equal(expected.begin(),expected.end(),out.p),"uint32 route fallback");
        c2::c8TransposePaired(p,original.data(),out.p,planes.p,temp.p,c2::Subfield16Kernel::RegionNonTemporal);
        require(std::equal(expected.begin(),expected.end(),out.p),"paired uint32 route fallback");
        if(k==333){
            std::vector<c2::Block> source(k);for(auto& value:source)value={rng(),rng()};
            const auto image=scalarForward(p,source.data());
            require(scalarDot(image.data(),original.data(),p.n)==scalarDot(source.data(),expected.data(),k),"full shortened adjoint");
        }
    }
    std::cout<<"cubic8 original/paired scalar/basis/all-byte/alias checks, direct GRS12883 outer, full scalar adjoint, schedules, guards and tails pass\n";
}
static void checkLarge(){
    std::mt19937_64 rng(9013);
    for(unsigned exponent:{16,18,20}){
        c2::Cubic8Plan p(std::size_t(1)<<exponent,113);
        c2::BenchPayload input(p.n),planes(p.n),temp(p.groups),expected(p.requested),out(p.requested),alias(p.n);
        for(std::size_t i=0;i!=p.n;++i)input.p[i]={rng(),rng()};
        c2::c8Transpose(p,input.p,expected.p,planes.p,temp.p,c2::Subfield16Kernel::Scatter,gather);
        for(bool paired:{false,true})for(auto selected:{gather,fused})for(auto mode:schedules){
            if(paired)c2::c8TransposePaired(p,input.p,out.p,planes.p,temp.p,mode,selected);
            else c2::c8Transpose(p,input.p,out.p,planes.p,temp.p,mode,selected);
            require(std::equal(expected.p,expected.p+p.requested,out.p),"large schedule");
            std::memcpy(alias.p,input.p,p.n*sizeof(c2::Block));
            if(paired)c2::c8TransposePaired(p,alias.p,alias.p,planes.p,temp.p,mode,selected);
            else c2::c8Transpose(p,alias.p,alias.p,planes.p,temp.p,mode,selected);
            require(std::equal(expected.p,expected.p+p.requested,alias.p),"large alias prefix");
            require(std::equal(input.p+p.requested,input.p+p.n,alias.p+p.requested),"large alias tail");
        }
        std::cout<<"K"<<exponent<<" cubic8 original/paired schedules and aliases pass\n";
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
static void benchmark(std::size_t k,std::uint64_t seed,unsigned calls,bool phases,const std::string& selection){
    c2::Cubic8Plan p(k,seed);c2::Cauchy128Plan old(k,seed);
    const auto mode=k>=1048576?c2::Subfield16Kernel::RegionAllNonTemporal:c2::Subfield16Kernel::Scatter;
    const auto oldVariant=k>=1048576?c2::Cauchy128CoefficientVariant::Gather:c2::Cauchy128CoefficientVariant::FusedLoad;
    const auto winner=k>=1048576?gather:fused;
    std::vector<unsigned> candidates{0,1};
    if(selection=="paired"){
        candidates.push_back(k>=1048576?2:3);candidates.push_back(5);
    }else{
        if(selection!="fused")candidates.push_back(2);if(selection!="gather")candidates.push_back(3);
    }
    const auto records=std::max({p.n,old.n,2*k}),planeRecords=std::max(p.n,old.n),tempRecords=std::max(p.groups,old.groups);
    c2::Block *input,*planes,*temp;
#if SPIN_CAUCHY_PUBLIC
    require(spin::valid_message_size(spin::Parameters::PacketRs8T64,k),"binary baseline K");
    spin::Code binary({k,spin::Parameters::PacketRs8T64,seed,seed});
    require(binary.backend()==spin::Backend::Avx512,"baseline AVX512 unavailable");
    const auto policy=k<=65536?spin::MemoryPolicy::Normal:spin::MemoryPolicy::PreferHugePages;
    spin::Buffer payload((records+12)*16,policy),planeBuffer(planeRecords*16,policy),tempBuffer(tempRecords*16,policy);
    auto work=binary.make_workspace(spin::Width::Bits128,policy);candidates.push_back(4);
    input=reinterpret_cast<c2::Block*>(payload.bytes().data()+64);
    planes=reinterpret_cast<c2::Block*>(planeBuffer.bytes().data());temp=reinterpret_cast<c2::Block*>(tempBuffer.bytes().data());
#else
    c2::BenchPayload payload(records+12,k>65536),planeBuffer(planeRecords,k>65536),tempBuffer(tempRecords,k>65536);
    input=payload.p+4;planes=planeBuffer.p;temp=tempBuffer.p;
#endif
    std::vector<c2::Block> initial(records);std::mt19937_64 rng(913);
    for(auto& value:initial)value={rng(),rng()};
    const auto restore=[&]{std::memcpy(input,initial.data(),records*sizeof(c2::Block));};
    const auto encode=[&](unsigned candidate){
        if(candidate==0)c2::c128Transpose(old,input,input,planes,temp,mode,compact);
        else if(candidate==1)c2::c128CoefficientTranspose(old,input,input,planes,temp,mode,oldVariant);
        else if(candidate<4)c2::c8Transpose(p,input,input,planes,temp,mode,candidate==2?gather:fused);
        else if(candidate==5)c2::c8TransposePaired(p,input,input,planes,temp,mode,winner);
#if SPIN_CAUCHY_PUBLIC
        else binary.transpose_inplace_bytes({reinterpret_cast<std::byte*>(input),2*k*16},work);
#endif
    };
    for(auto candidate:candidates){restore();for(unsigned i=0;i!=5;++i)encode(candidate);}
    std::array<std::vector<double>,6> samples,inner,outer;
    for(auto candidate:candidates){samples[candidate].resize(calls);
        if(phases&&candidate!=4){inner[candidate].resize(calls);outer[candidate].resize(calls);}}
    std::array<std::uint64_t,6> sums{};
    for(unsigned done=0,pass=0;done<calls;++pass){
        const unsigned count=std::min(5U,calls-done);
        for(unsigned order=0;order!=candidates.size();++order){
            const unsigned candidate=candidates[(order+pass)%candidates.size()];restore();
            for(unsigned i=0;i!=count;++i){
                const auto a=Clock::now();
                if(phases&&candidate!=4){
                    if(candidate<2)c2::c128InnerDispatch(old,input,planes,temp,mode,compact);
                    else if(candidate==5)c2::c8InnerDispatchPaired(p,input,planes,temp,mode);
                    else c2::c8InnerDispatch(p,input,planes,temp,mode);
                    const auto b=Clock::now();
                    if(candidate==0)c2::c128OuterDispatch(old,planes,input,mode);
                    else if(candidate==1)c2::c128CoefficientOuterDispatch(old,planes,input,mode,oldVariant);
                    else if(candidate==5)c2::c8OuterDispatch(p,planes,input,mode,winner);
                    else c2::c8OuterDispatch(p,planes,input,mode,candidate==2?gather:fused);
                    const auto c=Clock::now();inner[candidate][done+i]=micros(a,b);outer[candidate][done+i]=micros(b,c);
                    samples[candidate][done+i]=micros(a,c);
                }else{encode(candidate);samples[candidate][done+i]=micros(a,Clock::now());}
            }
            sums[candidate]=(sums[candidate]^checksum(input,k)^count)*0x100000001b3ULL;
        }
        done+=count;
    }
    if(selection=="all")require(sums[2]==sums[3],"cubic8 gather/fused timing checksums differ");
    if(selection=="paired")require(sums[5]==sums[k>=1048576?2:3],"cubic8 original/paired timing checksums differ");
    std::cout<<"same_buffer_address,"<<static_cast<void*>(input)<<",same_field_scratch,1,seed,"<<seed
        <<",cubic8_padded_K,"<<p.k<<",old_padded_K,"<<old.k
        <<",schedule,"<<(mode==c2::Subfield16Kernel::Scatter?"scatter":"region-nt-all")
        <<",pages,"<<(k<=65536?"normal":"huge-preferred")<<",public_allocator,"<<int(SPIN_CAUCHY_PUBLIC)
        <<",cubic8_outer_variants,"<<selection<<",different_codes_and_source_bases,1,restore_outside_timing,batch5,rotating_order\n"
        <<"name,K,N,calls,median_us,p10_us,p90_us,checksum\n";
    const char* names[]={"cauchy128-systematic87",k>=1048576?"cauchy128-coefficients87-gather":"cauchy128-coefficients87-fused",
                        "cubic8-coefficients83-gather","cubic8-coefficients83-fused","binary-rs8-public",
                        k>=1048576?"cubic8-coefficients83-paired-gather":"cubic8-coefficients83-paired-fused"};
    for(auto candidate:candidates){
        const auto n=candidate==4?2*k:(candidate<2?old.n:p.n);
        stats(names[candidate],k,n,samples[candidate],sums[candidate]);
        if(phases&&candidate!=4){stats(std::string(names[candidate])+"-inner-route",k,n,inner[candidate],0);
                               stats(std::string(names[candidate])+"-outer",k,n,outer[candidate],0);}
    }
}
int main(int argc,char** argv){try{
    if(!c2::gfni512Available())return 77;
    if(argc==2&&std::string(argv[1])=="--check"){check();return 0;}
    if(argc==2&&(std::string(argv[1])=="--large-check"||std::string(argv[1])=="--check-large")){checkLarge();return 0;}
    if(argc<2||argc>6)throw std::invalid_argument("usage: K [seed=113] [calls=101] [full|phases] [all|gather|fused|paired]");
    const auto k=std::stoull(argv[1]),seed=argc>2?std::stoull(argv[2]):113;
    const auto calls=argc>3?std::stoull(argv[3]):101;const std::string mode=argc>4?argv[4]:"full",selection=argc>5?argv[5]:"all";
    if(!calls||calls>1000000||(mode!="full"&&mode!="phases")||
       (selection!="all"&&selection!="gather"&&selection!="fused"&&selection!="paired"))throw std::invalid_argument("mode, variants, or calls");
    benchmark(k,seed,static_cast<unsigned>(calls),mode=="phases",selection);return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
