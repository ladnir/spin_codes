#include "Gamma30RouteEncoder.h"
#include "../BenchPayload.h"
#include <atomic>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>

// Tagged copy of the frozen completion gamma30 timing driver: only engine
// namespace/name change. Setup, generator, samples, policy, and hashes match.
// Isolated timing driver, separate from the dense-oracle correctness driver.
// Run gamma30_route_encoder --check successfully on this build/CPU before timing.
// Only root/operator launches workloads; never run concurrent benchmarks.
namespace enc=spin::research::char2::field128_lift::gamma30_route;
namespace literal=enc::literal;
namespace c2=spin::research::char2;
using Block=enc::Block;
using Clock=std::chrono::steady_clock;

static void require(bool condition,const char* message){
    if(!condition)throw std::runtime_error(message);
}

// Deliberately identical to the frozen literal_encoder.cpp input generator.
static void produce(Block* output,std::size_t n,std::uint64_t seed=913){
    std::uint64_t x=seed,y=113;
    for(std::size_t i=0;i!=n;++i){x^=x<<13;x^=x>>7;x^=x<<17;y+=0x9e3779b97f4a7c15ULL;output[i]={x,y};}
}
static std::uint64_t hash(const Block* p,std::size_t n){
    std::uint64_t h=0;
    for(std::size_t i=0;i!=n;++i){h=(h^p[i].lo)*0x100000001b3ULL;h=(h^p[i].hi)*0x100000001b3ULL;}
    return h;
}

enum class Mode {Full,Outer,Inner};
struct SampleSummary {double median,mean,best,worst;};
static SampleSummary summarize(const std::vector<double>& samples){
    require(!samples.empty(),"empty timing samples");
    auto sorted=samples;std::sort(sorted.begin(),sorted.end());
    const auto count=sorted.size();
    const auto median=count%2?sorted[count/2]:0.5*(sorted[count/2-1]+sorted[count/2]);
    return {median,std::accumulate(sorted.begin(),sorted.end(),0.0)/count,sorted.front(),sorted.back()};
}

static void releaseOracleCopies(enc::Plan& plan){
    auto& p=plan.outerPlan;
    // The hot loops use compact outer keys and packed routes only. The dense
    // oracle retains these copies in its own separate correctness process.
    decltype(p.toeplitz){}.swap(p.toeplitz);
    decltype(p.route){}.swap(p.route);
}
static std::size_t setupPayloadBytes(const enc::Plan& plan){
    return plan.outerPlan.setupBytes()+plan.refresh.size()*sizeof(enc::Refresh);
}

template<bool Stream,Mode Which>
static inline void runOne(const enc::Plan& plan,const Block* input,Block* output,Block* scratch){
    // All dispatch is compile-time and outside the measured calls.
    if constexpr(Which==Mode::Full)enc::forward<Stream>(plan,input,output,scratch);
    else if constexpr(Which==Mode::Outer)enc::outer<true>(plan.outerPlan,input,scratch);
    else enc::inner<Stream>(plan,scratch,output);
}

template<bool Stream,Mode Which>
static void benchmark(std::size_t k,unsigned calls){
    const auto setupBegin=Clock::now();
    enc::Plan plan(k,113);const auto& p=plan.outerPlan;
    // This is the payload size after Plan construction, not peak allocator
    // usage during its temporary frozen-state16 setup or overall process RSS.
    const auto payloadBeforeRelease=setupPayloadBytes(plan);
    const auto coefficientBytes=plan.hotKeyBytes();
    releaseOracleCopies(plan);
    const auto retainedBytes=setupPayloadBytes(plan);
    c2::BenchPayload input(p.k),output(p.n),scratch(plan.scratchSize());
    produce(input.p,p.k);
    std::vector<double> samples(calls);
    // Match the frozen literal driver: populate scratch outside timing for
    // every mode, including the standalone inner workload.
    enc::outer<true>(p,input.p,scratch.p);
    const auto setupEnd=Clock::now();
    // Three (not five) warmups match literal_encoder.cpp exactly. Keep input,
    // allocation policy, store policy, and CPU affinity matched externally.
    constexpr unsigned warmups=3;
    for(unsigned i=0;i!=warmups;++i)runOne<Stream,Which>(plan,input.p,output.p,scratch.p);
    for(unsigned i=0;i!=calls;++i){
        std::atomic_signal_fence(std::memory_order_seq_cst);
        const auto begin=Clock::now();
        runOne<Stream,Which>(plan,input.p,output.p,scratch.p);
        std::atomic_signal_fence(std::memory_order_seq_cst);
        const auto end=Clock::now();
        samples[i]=std::chrono::duration<double,std::milli>(end-begin).count();
    }
    // No output scan, regeneration, validation, or sorting between samples.
    const auto checksum=Which==Mode::Outer?hash(scratch.p,plan.scratchSize()):hash(output.p,p.n);
    const auto result=summarize(samples);
    const char* name=Which==Mode::Full?"full":Which==Mode::Outer?"outer":"inner";
    std::cout<<"engine=gamma30_route_GF128 mode="<<name
             <<" logK="<<std::countr_zero(k)<<" samples="<<calls<<" warmups="<<warmups
             <<" scratch=NT output="<<(Which==Mode::Outer?"unused":Stream?"NT":"cached")
             <<" prefetch="<<(Which==Mode::Outer?"unused":"T1-next-region-16lines")
             <<" median_ms="<<std::fixed<<std::setprecision(6)<<result.median<<" mean_ms="<<result.mean
             <<" best_ms="<<result.best<<" worst_ms="<<result.worst
             <<" setup_ms="<<std::chrono::duration<double,std::milli>(setupEnd-setupBegin).count()
             <<" hot_coefficient_bytes="<<coefficientBytes
             <<" hot_route_bytes="<<p.packetRoute.size()*sizeof(std::uint16_t)
             <<" state_live_records=30 state_storage_records=32"
             <<" retained_setup_bytes="<<retainedBytes
             <<" setup_payload_before_release_bytes="<<payloadBeforeRelease
             <<" rejected_toeplitz="<<p.toeplitzRejected
             <<" final_checksum="<<std::hex<<checksum<<std::dec<<"\n";
}

template<bool Stream>
static void benchmarkMode(std::size_t k,unsigned calls,Mode mode){
    if(mode==Mode::Full)benchmark<Stream,Mode::Full>(k,calls);
    else if(mode==Mode::Outer)benchmark<Stream,Mode::Outer>(k,calls);
    else benchmark<Stream,Mode::Inner>(k,calls);
}
static unsigned integer(const char* text){
    const std::string value(text);std::size_t end=0;
    const auto result=std::stoul(value,&end);
    require(end==value.size()&&result<=1000000,"invalid or excessive integer argument");
    return static_cast<unsigned>(result);
}

int main(int argc,char** argv){
    try{
        // ISA-targeted research binary, not a runtime CPU dispatcher.
        // Defaults: 101 samples, K20, full path, NT output. Explicit flags use
        // the same ordering as the frozen literal benchmark driver.
        if(argc==2&&std::string(argv[1])=="--help"){
            std::cout<<"usage: gamma30_route_bench [--bench SAMPLES LOGK full|outer|inner [nt|cached]]\n"
                     <<"Default: --bench 101 20 full nt; warmups=3. CPU pinning is external.\n"
                     <<"Run correctness first. NEVER run concurrent benchmarks.\n";return 0;
        }
        unsigned calls=101,logK=20;Mode mode=Mode::Full;bool stream=true;
        if(argc!=1){
            require((argc==5||argc==6)&&std::string(argv[1])=="--bench","use --help for syntax");
            calls=integer(argv[2]);logK=integer(argv[3]);
            const std::string name=argv[4];
            require(name=="full"||name=="outer"||name=="inner","mode must be full, outer, or inner");
            mode=name=="full"?Mode::Full:name=="outer"?Mode::Outer:Mode::Inner;
            if(argc==6){const std::string store=argv[5];require(store=="nt"||store=="cached","store must be nt or cached");stream=store=="nt";}
        }
        require(calls>=1&&calls<=1000000,"samples must be in 1..1000000");
        require(logK>=13&&logK<=22,"logK must be in 13..22");
        std::cout<<"Research MT19937 setup; no ideal-uniform or cryptographic-seed guarantee. Run only one benchmark process at a time.\n";
        if(stream)benchmarkMode<true>(std::size_t(1)<<logK,calls,mode);
        else benchmarkMode<false>(std::size_t(1)<<logK,calls,mode);
        return 0;
    }catch(const std::exception& error){std::cerr<<"FAIL: "<<error.what()<<"\n";return 1;}
}
