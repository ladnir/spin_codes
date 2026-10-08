#include "Cauchy128CoefficientCode.h"
#include "BenchPayload.h"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>
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
static double micros(Clock::time_point a,Clock::time_point b){
    return std::chrono::duration<double,std::micro>(b-a).count();
}
static constexpr auto compact=c2::Cauchy128Variant::CompactBatch;
static constexpr auto gather=c2::Cauchy128CoefficientVariant::Gather;
static constexpr auto fused=c2::Cauchy128CoefficientVariant::FusedLoad;
static constexpr c2::Subfield16Kernel schedules[]={c2::Subfield16Kernel::Scatter,
    c2::Subfield16Kernel::RegionNonTemporal,c2::Subfield16Kernel::RegionAllNonTemporal};

// Direct products of vanishing polynomials, independent of the FFT circuit.
struct CoefficientOracle {
    std::array<std::array<std::uint8_t,256>,256> mul{};
    std::array<std::array<std::uint8_t,87>,128> evaluation{};
    CoefficientOracle(){
        for(unsigned a=0;a!=256;++a)for(unsigned b=0;b!=256;++b)mul[a][b]=c2::sfMultiply(a,b);
        for(unsigned point=0;point!=128;++point)for(unsigned j=0;j!=87;++j){
            std::uint8_t value=1;
            for(unsigned bit=0;bit!=7;++bit)if(j&(1U<<bit)){
                std::uint8_t numerator=1,denominator=1;
                for(unsigned label=0;label!=(1U<<bit);++label){
                    numerator=mul[numerator][point^label];
                    denominator=mul[denominator][(1U<<bit)^label];
                }
                value=mul[value][mul[numerator][c2::sfInverse(denominator)]];
            }
            evaluation[point][j]=value;
        }
    }
};
static const CoefficientOracle& oracle(){static const CoefficientOracle value;return value;}

static void checkOuter(){
    const auto& o=oracle();std::mt19937_64 rng(12887321);
    alignas(64) std::uint8_t input[128][64],expected[128][64];
    alignas(64) std::uint8_t raw[128*64+128],source[128*64+128];
    // Every scalar coordinate is exercised on all64 SIMD byte lanes; random
    // trials then test sums, arbitrary byte values, and all output alignments.
    for(unsigned trial=0;trial!=132;++trial){
        std::memset(input,0,sizeof(input));
        if(trial<128)for(unsigned b=0;b!=64;++b)input[trial][b]=static_cast<std::uint8_t>(b+1);
        else for(auto& row:input)for(auto& value:row)value=static_cast<std::uint8_t>(rng());
        for(bool transpose:{false,true}){
            const unsigned outputs=transpose?87:128;
            std::memset(expected,0,sizeof(expected));
            for(unsigned j=0;j!=outputs;++j)for(unsigned b=0;b!=64;++b){
                auto& value=expected[j][b];
                if(trial<128){
                    if(transpose)value=o.mul[o.evaluation[trial][j]][input[trial][b]];
                    else if(trial<87)value=o.mul[o.evaluation[j][trial]][input[trial][b]];
                }else if(transpose){
                    for(unsigned i=0;i!=128;++i)value^=o.mul[o.evaluation[i][j]][input[i][b]];
                }else{
                    for(unsigned i=0;i!=87;++i)value^=o.mul[o.evaluation[j][i]][input[i][b]];
                }
            }
            for(unsigned offset=0;offset!=(trial<128?1U:4U);++offset)for(bool alias:{false,true}){
                std::memset(raw,0xa7,sizeof(raw));std::memset(source,0xb9,sizeof(source));
                auto* destination=raw+64+16*offset;
                auto* initial=alias?destination:source+64+16*offset;
                std::memcpy(initial,input,sizeof(input));
                const auto* from=reinterpret_cast<const __m512i*>(initial);
                auto* to=reinterpret_cast<__m512i*>(destination);
                if(transpose)c2::c128OuterCoefficientTranspose(from,to);
                else c2::c128OuterCoefficientForward(from,to);
                require(std::memcmp(destination,expected,outputs*64)==0,"coefficient outer dense oracle");
                if(!alias)require(std::memcmp(initial,input,sizeof(input))==0,"outer input preservation");
                for(unsigned i=0;i!=64+16*offset;++i)require(raw[i]==0xa7,"outer prefix guard");
                if(transpose&&alias)require(std::memcmp(destination+87*64,input[87],41*64)==0,"outer alias tail");
                else for(unsigned i=outputs*64;i!=128*64;++i)require(destination[i]==0xa7,"outer output tail");
                for(unsigned i=64+16*offset+128*64;i!=sizeof(raw);++i)require(raw[i]==0xa7,"outer suffix guard");
            }
            if(transpose)for(std::size_t stride:{64,83,192})
              for(unsigned offset=0;offset!=(trial<128?1U:4U);++offset){
                std::vector<std::uint8_t> packed(128*stride+64,0xb9);
                const auto inputOffset=3+16*offset; // Deliberately byte-unaligned.
                for(unsigned i=0;i!=128;++i)std::memcpy(packed.data()+inputOffset+i*stride,input[i],64);
                const auto original=packed;
                std::memset(raw,0xa7,sizeof(raw));auto* destination=raw+64+16*offset;
                c2::c128OuterCoefficientTransposeStrided(packed.data()+inputOffset,stride,
                                                        reinterpret_cast<__m512i*>(destination));
                require(std::memcmp(destination,expected,87*64)==0,"strided leaf direct polynomial oracle");
                require(packed==original,"strided leaf input preservation");
                for(unsigned i=0;i!=64+16*offset;++i)require(raw[i]==0xa7,"strided leaf prefix guard");
                for(unsigned i=64+16*offset+128*64;i!=sizeof(raw);++i)require(raw[i]==0xa7,"strided leaf suffix guard");
            }
        }
    }
}

static std::vector<c2::Block> basisTransposePrefix(const c2::Cauchy128Plan& p,const c2::Block* systematic){
    const auto& o=oracle();std::vector<c2::Block> result(p.requested);
    for(std::size_t q=0;q!=p.requested;++q){
        const auto row=(q/348)*4+q%4;const auto j=(q%348)/4;
        auto* destination=reinterpret_cast<std::uint8_t*>(&result[q]);
        for(unsigned i=0;i!=87;++i){
            const auto* source=reinterpret_cast<const std::uint8_t*>(systematic+c2::Cauchy128Plan::sourceIndex(row,i));
            const auto coefficient=o.evaluation[i][j];
            for(unsigned b=0;b!=16;++b)destination[b]^=o.mul[coefficient][source[b]];
        }
    }
    return result;
}

static void check(){
    checkOuter();std::mt19937_64 rng(173);const c2::Block guard{0x32456789ULL,0xfabc1234ULL};
    for(std::size_t k:{1,86,87,348,349,5568,6000}){
        c2::Cauchy128Plan p(k,113),full(p.k,113);
        require(p.n==full.n&&p.route==full.route&&p.inputDiagonal==full.inputDiagonal&&
                p.stateDiagonal==full.stateDiagonal,"full padded oracle setup");
        std::vector<c2::Block> original(p.n);for(auto& value:original)value={rng(),rng()};
        c2::BenchPayload input(p.n+8,false),out(k+8,false),planes(p.n+8,false),temp(p.groups+8,false),old(full.k,false);
        c2::c128Transpose(full,original.data(),old.p,planes.p,temp.p,c2::Subfield16Kernel::Scatter,compact);
        // Obtain EVERY padded systematic coordinate before the B^T map.
        const auto expected=basisTransposePrefix(p,old.p);
        for(auto selected:{gather,fused})for(auto mode:schedules)for(unsigned offset:{0,1,2,3}){
            for(auto* buffer:{&input,&out,&planes,&temp})std::fill(buffer->p,buffer->p+buffer->count,guard);
            std::memcpy(input.p+offset,original.data(),p.n*sizeof(c2::Block));
            c2::c128CoefficientTranspose(p,input.p+offset,out.p+offset,planes.p+offset,temp.p+offset,mode,selected);
            require(std::equal(expected.begin(),expected.end(),out.p+offset),"full B^T padded oracle");
            require(std::equal(original.begin(),original.end(),input.p+offset),"full input preservation");
            c2::c128CoefficientTranspose(p,input.p+offset,input.p+offset,planes.p+offset,temp.p+offset,mode,selected);
            require(std::equal(expected.begin(),expected.end(),input.p+offset),"full alias prefix");
            require(std::equal(original.begin()+k,original.end(),input.p+offset+k),"full alias tail");
            for(auto* buffer:{&input,&out,&planes,&temp}){
                for(unsigned j=0;j!=offset;++j)require(buffer->p[j]==guard,"full prefix guard");
                for(std::size_t j=buffer->count-8+offset;j!=buffer->count;++j)
                    require(buffer->p[j]==guard,"full suffix guard");
            }
        }
    }
    std::cout<<"coefficient outer direct polynomial/basis/random/alias/alignment checks and full B^T padded oracle, schedules, guards, tails pass\n";
}

static void checkLarge(){
    std::mt19937_64 rng(9013);
    for(unsigned exponent:{16,18,20}){
        c2::Cauchy128Plan p(std::size_t(1)<<exponent,113);
        c2::BenchPayload input(p.n),planes(p.n),temp(p.groups),expected(p.requested),out(p.requested),alias(p.n);
        for(std::size_t i=0;i!=p.n;++i)input.p[i]={rng(),rng()};
        c2::c128CoefficientTranspose(p,input.p,expected.p,planes.p,temp.p,c2::Subfield16Kernel::Scatter);
        for(auto selected:{gather,fused})for(auto mode:schedules){
            c2::c128CoefficientTranspose(p,input.p,out.p,planes.p,temp.p,mode,selected);
            require(std::equal(expected.p,expected.p+p.requested,out.p),"large coefficient schedule");
            std::memcpy(alias.p,input.p,p.n*sizeof(c2::Block));
            c2::c128CoefficientTranspose(p,alias.p,alias.p,planes.p,temp.p,mode,selected);
            require(std::equal(expected.p,expected.p+p.requested,alias.p),"large coefficient alias prefix");
            require(std::equal(input.p+p.requested,input.p+p.n,alias.p+p.requested),"large coefficient alias tail");
        }
        std::cout<<"K"<<exponent<<" coefficient schedules and aliases pass\n";
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
    c2::Cauchy128Plan p(k,seed);
    const auto mode=k>=1048576?c2::Subfield16Kernel::RegionAllNonTemporal:c2::Subfield16Kernel::Scatter;
    std::vector<unsigned> candidates{0};
    if(selection!="fused")candidates.push_back(1);
    if(selection!="gather")candidates.push_back(2);
    const std::size_t records=std::max(p.n,2*k);
    c2::Block *input,*planes,*temp;
#if SPIN_CAUCHY_PUBLIC
    require(spin::valid_message_size(spin::Parameters::PacketRs8T64,k),"binary baseline K");
    spin::Code binary({k,spin::Parameters::PacketRs8T64,seed,seed});
    require(binary.backend()==spin::Backend::Avx512,"baseline AVX512 unavailable");
    const auto policy=k<=65536?spin::MemoryPolicy::Normal:spin::MemoryPolicy::PreferHugePages;
    spin::Buffer payload((records+12)*16,policy),planeBuffer(p.n*16,policy),tempBuffer(p.groups*16,policy);
    auto work=binary.make_workspace(spin::Width::Bits128,policy);candidates.push_back(3);
    input=reinterpret_cast<c2::Block*>(payload.bytes().data()+64);
    planes=reinterpret_cast<c2::Block*>(planeBuffer.bytes().data());
    temp=reinterpret_cast<c2::Block*>(tempBuffer.bytes().data());
#else
    c2::BenchPayload payload(records+12,k>65536),planeBuffer(p.n,k>65536),tempBuffer(p.groups,k>65536);
    input=payload.p+4;planes=planeBuffer.p;temp=tempBuffer.p;
#endif
    std::vector<c2::Block> initial(records);std::mt19937_64 rng(913);
    for(auto& value:initial)value={rng(),rng()};
    const auto restore=[&]{std::memcpy(input,initial.data(),records*sizeof(c2::Block));};
    const auto encode=[&](unsigned candidate){
        if(candidate==0)c2::c128Transpose(p,input,input,planes,temp,mode,compact);
        else if(candidate==1)c2::c128CoefficientTranspose(p,input,input,planes,temp,mode);
        else if(candidate==2)c2::c128CoefficientTranspose(p,input,input,planes,temp,mode,fused);
#if SPIN_CAUCHY_PUBLIC
        else binary.transpose_inplace_bytes({reinterpret_cast<std::byte*>(input),2*k*16},work);
#endif
    };
    for(auto candidate:candidates){restore();for(unsigned i=0;i!=5;++i)encode(candidate);}
    std::array<std::vector<double>,4> samples,inner,outer;
    for(auto candidate:candidates){
        samples[candidate].resize(calls);if(phases&&candidate<3){inner[candidate].resize(calls);outer[candidate].resize(calls);}
    }
    std::array<std::uint64_t,4> sums{};
    for(unsigned done=0,pass=0;done<calls;++pass){
        const unsigned count=std::min(5U,calls-done);
        for(unsigned order=0;order!=candidates.size();++order){
            const unsigned candidate=candidates[(order+pass)%candidates.size()];restore();
            for(unsigned i=0;i!=count;++i){
                const auto a=Clock::now();
                if(phases&&candidate<3){
                    c2::c128InnerDispatch(p,input,planes,temp,mode,compact);const auto b=Clock::now();
                    if(candidate==0)c2::c128OuterDispatch(p,planes,input,mode);
                    else c2::c128CoefficientOuterDispatch(p,planes,input,mode,candidate==2?fused:gather);
                    const auto c=Clock::now();inner[candidate][done+i]=micros(a,b);outer[candidate][done+i]=micros(b,c);
                    samples[candidate][done+i]=micros(a,c);
                }else{encode(candidate);samples[candidate][done+i]=micros(a,Clock::now());}
            }
            sums[candidate]=(sums[candidate]^checksum(input,k)^count)*0x100000001b3ULL;
        }
        done+=count;
    }
    if(selection=="all")require(sums[1]==sums[2],"gather/fused timing checksums differ");
    std::cout<<"same_buffer_address,"<<static_cast<void*>(input)<<",same_field_scratch,1,seed,"<<seed
        <<",padded_K,"<<p.k<<",schedule,"<<(mode==c2::Subfield16Kernel::Scatter?"scatter":"region-nt-all")
        <<",pages,"<<(k<=65536?"normal":"huge-preferred")
        <<",public_allocator,"<<int(SPIN_CAUCHY_PUBLIC)
        <<",coefficient_variants,"<<selection
        <<",different_source_bases,1,restore_outside_timing,batch5,rotating_order\n"
        <<"name,K,N,calls,median_us,p10_us,p90_us,checksum\n";
    const char* names[]={"cauchy128-systematic","cauchy128-coefficients","cauchy128-coefficients-fused","binary-rs8-public"};
    for(auto candidate:candidates){
        stats(names[candidate],k,candidate==3?2*k:p.n,samples[candidate],sums[candidate]);
        if(phases&&candidate<3){
            stats(std::string(names[candidate])+"-inner-route",k,p.n,inner[candidate],0);
            stats(std::string(names[candidate])+"-outer",k,p.n,outer[candidate],0);
        }
    }
}

int main(int argc,char** argv){try{
    if(!c2::gfni512Available())return 77;
    if(argc==2&&std::string(argv[1])=="--check"){check();return 0;}
    if(argc==2&&(std::string(argv[1])=="--large-check"||std::string(argv[1])=="--check-large")){checkLarge();return 0;}
    if(argc<2||argc>6)throw std::invalid_argument("usage: K [seed=113] [calls=101] [full|phases] [all|gather|fused]");
    const auto k=std::stoull(argv[1]),seed=argc>2?std::stoull(argv[2]):113;
    const auto calls=argc>3?std::stoull(argv[3]):101;
    const std::string mode=argc>4?argv[4]:"full";
    const std::string selection=argc>5?argv[5]:"all";
    if(!calls||calls>1000000||(mode!="full"&&mode!="phases")||
       (selection!="all"&&selection!="gather"&&selection!="fused"))throw std::invalid_argument("mode, variants, or calls");
    benchmark(k,seed,static_cast<unsigned>(calls),mode=="phases",selection);return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
