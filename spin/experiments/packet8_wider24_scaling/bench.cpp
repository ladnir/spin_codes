// Retain the literal frozen scalar/oracle helpers without modifying them.
#define main packet8wide24_frozen_main
#include "../packet8_wider24/bench.cpp"
#undef main
#include "../packet8_wider24_opt/RandomizerVariants.h"
#include "StoreVariants.h"
#include <fstream>
#include <string>
namespace opt=spin::research::packet8wide24opt;
namespace scaling=spin::research::packet8wide24scaling;

// 0: original frozen kernel; 1: retained flat9; 2: flat9 + NT route;
// 3: cached route + NT output; 4: NT route/output; 5: copied flat9 cached.
static void route(const Block* x,Block* s,const candidate::Plan& p,unsigned mode) {
    if(mode==2||mode==4)scaling::reverseRoute(x,s,p,true);
    else candidate::reverseRouteFast(x,s,p);
}
static void outer(const Block* s,Block* y,const candidate::Plan& p,
                  const opt::RandomizerPlan& tables,unsigned mode) {
    if(mode==0)candidate::outerFast(s,y,p);
    else if(mode==1||mode==2)opt::outerRandomizer(s,y,p,tables,1);
    else scaling::outer(s,y,p,mode!=5);
}
static void checkModes(const candidate::Plan& p,const opt::RandomizerPlan& tables,unsigned mode) {
    checkSetup(p);
    if(!opt::randomizerTablesMatch(p,tables))throw std::runtime_error("table algebra mismatch");
    constexpr auto policy=spin::MemoryPolicy::Normal;
    spin::Buffer input((p.n+12)*16,policy),packed(p.scratchBlocks()*16,policy),raw(p.scratchBlocks()*16,policy),
        reference(p.scratchBlocks()*16,policy),expected(p.k*16,policy),output((p.k+12)*16,policy),
        message(p.k*16,policy),encoded(p.n*16,policy);
    fill(input,1701);fill(message,3949);
    const std::vector<std::byte> saved(input.bytes().begin(),input.bytes().end());
    candidate::transposeScalar(blocks(input)+4,blocks(expected),blocks(raw),p);
    candidate::packScratchScalar(blocks(raw),blocks(reference),p);
    candidate::forwardScalar(blocks(message),blocks(encoded),p);
    if(dotLanes(blocks(message),blocks(expected),p.k)!=dotLanes(blocks(encoded),blocks(input)+4,p.n))
        throw std::runtime_error("literal forward/transpose adjoint mismatch");
    for(unsigned v=0;v<6;++v)if(mode==99||mode==v)for(unsigned offset=0;offset<4;++offset) {
        auto* x=blocks(input)+4+offset;
        std::memcpy(x,saved.data()+64,p.n*16);
        std::fill(packed.bytes().begin(),packed.bytes().end(),std::byte{0xa5});
        route(x,blocks(packed),p,v);
        equalScratch(blocks(packed),blocks(reference),p,"route differs from scalar oracle");
        for(std::size_t g=0;g<p.groups;++g)
            for(std::size_t b=512*16;b<candidate::groupStride*16;++b)
                if(packed.bytes()[g*candidate::groupStride*16+b]!=std::byte{0xa5})
                    throw std::runtime_error("route padding overwritten");
        std::fill(output.bytes().begin(),output.bytes().end(),std::byte{0xa5});
        auto* y=blocks(output)+4+offset;
        outer(blocks(packed),y,p,tables,v);
        if(std::memcmp(y,blocks(expected),p.k*16))throw std::runtime_error("outer differs from scalar oracle");
        if(std::memcmp(x,saved.data()+64,p.n*16))throw std::runtime_error("input changed");
        for(std::size_t b=0;b<output.bytes().size();++b)
            if((b<(4+offset)*16||b>=(4+offset+p.k)*16)&&output.bytes()[b]!=std::byte{0xa5})
                throw std::runtime_error("output guard overwritten");
        outer(blocks(packed),x,p,tables,v);
        if(std::memcmp(x,blocks(expected),p.k*16)||std::memcmp(x+p.k,saved.data()+64+p.k*16,(p.n-p.k)*16))
            throw std::runtime_error("in-place prefix or preserved suffix mismatch");
    }
}
static std::size_t anonymousHugeKiB() {
#if defined(__linux__)
    std::ifstream in("/proc/self/smaps_rollup");std::string line;
    while(std::getline(in,line))if(line.rfind("AnonHugePages:",0)==0)
        return std::stoull(line.substr(line.find(':')+1));
#endif
    return 0;
}
int main(int argc,char** argv) {try {
    if(!spin::packet_fast_available()||!spin::capabilities().forward512)return 77;
    if(argc>7)throw std::invalid_argument("[K=65536] [seed=1] [calls=0] [mode=1] [normal|huge] [phases=0]");
    const std::size_t k=argc>1?std::stoull(argv[1]):65536;
    const std::uint64_t seed=argc>2?std::stoull(argv[2]):1;
    const auto requested=argc>3?std::stoull(argv[3]):0;
    const unsigned mode=argc>4?std::stoul(argv[4]):1;
    const std::string memory=argc>5?argv[5]:"normal";
    const bool phases=argc>6&&std::stoul(argv[6]);
    if(requested>1000000||(mode>5&&mode!=99)||(mode==99&&requested)
        ||(memory!="normal"&&memory!="huge"))throw std::invalid_argument("invalid calls/mode/page policy");
    candidate::Plan plan(k,seed);opt::RandomizerPlan tables(plan);checkModes(plan,tables,mode);
    if(!requested){std::cout<<"PASS packet8 scaling K="<<k<<" seed="<<seed<<" mode="<<mode<<'\n';return 0;}
    const auto policy=memory=="huge"?spin::MemoryPolicy::PreferHugePages:spin::MemoryPolicy::Normal;
    spin::Buffer input(plan.n*16,policy),scratch(plan.scratchBlocks()*16,policy);
    fill(input);auto* x=blocks(input);auto* s=blocks(scratch);
    auto doRoute=[&]{route(x,s,plan,mode);};
    auto doOuter=[&]{outer(s,x,plan,tables,mode);};
    for(unsigned i=0;i<5;++i){doRoute();doOuter();}
    const auto hugeKiB=anonymousHugeKiB();
    const auto calls=unsigned(requested);
    std::vector<double> total(calls),inner(calls),outerTimes(calls);
    for(unsigned i=0;i<calls;++i) {
        const auto a=Clock::now();doRoute();auto b=a;
        if(phases)b=Clock::now();
        doOuter();const auto c=Clock::now();
        total[i]=std::chrono::duration<double,std::micro>(c-a).count();
        inner[i]=phases?std::chrono::duration<double,std::micro>(b-a).count():0.;
        outerTimes[i]=phases?std::chrono::duration<double,std::micro>(c-b).count():0.;
    }
    for(auto* values:{&total,&inner,&outerTimes})std::sort(values->begin(),values->end());
    std::uint64_t checksum=0;
    for(std::size_t i=0;i<plan.k*16;i+=8){std::uint64_t v;std::memcpy(&v,input.bytes().data()+i,8);checksum=(checksum^v)*0x100000001b3ULL;}
    std::cout<<"mode,K,seed,calls,memory,phases,median_us,p10_us,p90_us,inner_us,outer_us,anon_huge_kib,checksum\n"
        <<mode<<','<<k<<','<<seed<<','<<calls<<','<<memory<<','<<phases<<','<<std::fixed<<std::setprecision(4)
        <<total[calls/2]<<','<<total[calls/10]<<','<<total[9*calls/10]<<','<<inner[calls/2]<<','<<outerTimes[calls/2]
        <<','<<hugeKiB<<','<<std::hex<<checksum<<'\n';
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
