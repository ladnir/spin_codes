// Reuse the frozen oracle/guard/adjoint checks without editing their source.
#define main packet8wide24_baseline_main
#include "../packet8_wider24/bench.cpp"
#undef main
#include "OuterTraffic.h"
#include "RandomizerVariants.h"
#include "FastTraffic.h"
#include "StructuredVariants.h"
namespace opt=spin::research::packet8wide24opt;

static bool validMode(unsigned mode) { return mode<=15||mode==20||mode==21; }
static void routeVariant(const Block* x,Block* s,const candidate::Plan& p,unsigned mode) {
    if(mode>=10&&mode<=15)opt::reverseRouteTraffic(x,s,p,mode-9);
    else candidate::reverseRouteFast(x,s,p);
}
static void outerVariant(const Block* s,Block* y,const candidate::Plan& p,
                         const opt::RandomizerPlan& tables,const opt::StructuredPlan& structured,unsigned mode) {
    if(mode>=1&&mode<=3)opt::outerTraffic(s,y,p,mode);
    else if(mode>=4&&mode<=6)opt::outerRandomizer(s,y,p,tables,mode-3);
    else if(mode==7)opt::outerRandomizer(s,y,p,tables,4);
    else if(mode==8||mode==9)opt::outerTraffic(s,y,p,mode-4);
    else if(mode==20||mode==21)opt::outerStructured(s,y,p,structured,mode-20);
    else candidate::outerFast(s,y,p);
}
static void checkVariant(const candidate::Plan& p,const opt::RandomizerPlan& tables,
                         const opt::StructuredPlan& structured,unsigned mode) {
    spin::Buffer input((p.n+12)*16),scratch(p.scratchBlocks()*16),reference(p.scratchBlocks()*16),
        expected(p.k*16),output((p.k+12)*16),raw(p.scratchBlocks()*16),message(p.k*16),encoded(p.n*16);
    fill(input,1701);
    const std::vector<std::byte> saved(input.bytes().begin(),input.bytes().end());
    candidate::transposeScalar(blocks(input)+4,blocks(expected),blocks(raw),p);
    candidate::packScratchScalar(blocks(raw),blocks(reference),p);
    fill(message,3949);
    candidate::forwardScalar(blocks(message),blocks(encoded),p);
    if(dotLanes(blocks(message),blocks(expected),p.k)!=dotLanes(blocks(encoded),blocks(input)+4,p.n))
        throw std::runtime_error("variant literal forward/transpose adjoint mismatch");
    for(unsigned offset=0;offset<4;++offset) {
        auto* x=blocks(input)+4+offset;
        std::memcpy(x,saved.data()+64,p.n*16);
        std::fill(scratch.bytes().begin(),scratch.bytes().end(),std::byte{0xa5});
        routeVariant(x,blocks(scratch),p,mode);
        equalScratch(blocks(scratch),blocks(reference),p,"optimized packed route mismatch");
        for(std::size_t g=0;g<p.groups;++g)
            for(std::size_t b=512*16;b<candidate::groupStride*16;++b)
                if(scratch.bytes()[g*candidate::groupStride*16+b]!=std::byte{0xa5})
                    throw std::runtime_error("optimized route padding overwritten");
        std::fill(output.bytes().begin(),output.bytes().end(),std::byte{0xa5});
        auto* y=blocks(output)+4+offset;
        outerVariant(blocks(scratch),y,p,tables,structured,mode);
        if(std::memcmp(y,blocks(expected),p.k*16))throw std::runtime_error("optimized outer mismatch");
        for(std::size_t b=0;b<output.bytes().size();++b)
            if((b<(4+offset)*16||b>=(4+offset+p.k)*16)&&output.bytes()[b]!=std::byte{0xa5})
                throw std::runtime_error("optimized output guard changed");
        if(std::memcmp(x,saved.data()+64,p.n*16))throw std::runtime_error("optimized input changed");
        outerVariant(blocks(scratch),x,p,tables,structured,mode);
        if(std::memcmp(x,blocks(expected),p.k*16)||std::memcmp(x+p.k,saved.data()+64+p.k*16,(p.n-p.k)*16))
            throw std::runtime_error("optimized in-place output or suffix mismatch");
    }
}

int main(int argc,char** argv) {try {
    if(!spin::packet_fast_available()||!spin::capabilities().forward512)return 77;
    if(argc>6)throw std::invalid_argument("[K=65536] [seed=1] [calls=0] [mode=0] [phases=0]");
    const std::size_t k=argc>1?std::stoull(argv[1]):65536;
    const std::uint64_t seed=argc>2?std::stoull(argv[2]):1;
    const auto requested=argc>3?std::stoull(argv[3]):0;
    const unsigned mode=argc>4?std::stoul(argv[4]):0;
    const bool phases=argc>5&&std::stoul(argv[5]);
    if(requested>1000000||(!validMode(mode)&&mode!=99)||(mode==99&&requested))
        throw std::invalid_argument("invalid calls/mode;99 means test every variant");
    const unsigned calls=unsigned(requested);
    candidate::Plan plan(k,seed);
    opt::RandomizerPlan tables(plan);
    opt::StructuredPlan structured(plan,seed);
    if(!opt::randomizerTablesMatch(plan,tables))throw std::runtime_error("randomizer setup algebra failed");
    checks(plan);
    if(!opt::structuredTablesMatch(structured))throw std::runtime_error("structured setup algebra failed");
    if(mode==99) {
        for(unsigned v:std::array<unsigned,13>{0,1,2,3,4,5,6,10,11,12,13,14,15})checkVariant(plan,tables,structured,v);
        auto omitted=plan;opt::applyOmittedRows(omitted);checkVariant(omitted,tables,structured,7);
        checkVariant(plan,tables,structured,8);checkVariant(omitted,tables,structured,9);
        auto changed=plan;opt::applyStructuredRows(changed,structured);
        checkVariant(changed,tables,structured,20);checkVariant(changed,tables,structured,21);
    } else {
        if(mode==7||mode==9)opt::applyOmittedRows(plan);
        if(mode==20||mode==21)opt::applyStructuredRows(plan,structured);
        checkVariant(plan,tables,structured,mode);
    }
    if(!calls){std::cout<<"PASS packet8wide24opt K="<<k<<" seed="<<seed<<" mode="<<mode<<'\n';return 0;}
    spin::Buffer input(plan.n*16,spin::MemoryPolicy::Normal),scratch(plan.scratchBlocks()*16,spin::MemoryPolicy::Normal);
    fill(input);auto* x=blocks(input);auto* s=blocks(scratch);
    auto route=[&]{routeVariant(x,s,plan,mode);};
    auto outer=[&]{outerVariant(s,x,plan,tables,structured,mode);};
    for(unsigned i=0;i<5;++i){route();outer();}
    std::vector<double> total(calls),inner(calls),outerTimes(calls);
    for(unsigned i=0;i<calls;++i) {
        const auto a=Clock::now();route();auto b=a;
        if(phases)b=Clock::now();
        outer();const auto c=Clock::now();
        total[i]=std::chrono::duration<double,std::micro>(c-a).count();
        inner[i]=phases?std::chrono::duration<double,std::micro>(b-a).count():0.;
        outerTimes[i]=phases?std::chrono::duration<double,std::micro>(c-b).count():0.;
    }
    for(auto* v:{&total,&inner,&outerTimes})std::sort(v->begin(),v->end());
    std::uint64_t checksum=0;
    for(std::size_t i=0;i<plan.k*16;i+=8){std::uint64_t v;std::memcpy(&v,input.bytes().data()+i,8);checksum=(checksum^v)*0x100000001b3ULL;}
    std::cout<<"mode,K,seed,calls,phases,median_us,p10_us,p90_us,inner_us,outer_us,checksum\n"
        <<mode<<','<<k<<','<<seed<<','<<calls<<','<<phases<<','<<std::fixed<<std::setprecision(4)
        <<total[calls/2]<<','<<total[calls/10]<<','<<total[9*calls/10]<<','
        <<inner[calls/2]<<','<<outerTimes[calls/2]<<','<<std::hex<<checksum<<'\n';
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
