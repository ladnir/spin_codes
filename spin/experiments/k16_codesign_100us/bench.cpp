#include "kernel/K16CodeSign.h"
#include "kernel/K16T128.h"
#include "kernel/K16NativeOuter.h"
#include "kernel/K16Monomial.h"
#include "kernel/K16Paired.h"
#include "kernel/K16PairedOptimized.h"
#include "kernel/K16PairedWide.h"
#include "kernel/K16PairedBasis.h"
#include "kernel/K16PairedBasisRaw.h"
#include "kernel/K16Paired15.h"
#include "kernel/K16PairedShuffle.h"
#include "kernel/K16Paired15Shuffle.h"
#include "kernel/K16PairedNoAlias.h"
#include "kernel/K16Paired15NoAlias.h"
#include "kernel/K16PairedFold.h"
#include "outer_variants/OuterVariants.h"
#include "outer_variants/OuterShuffleVariants.h"
#include "outer_variants/OuterNonTemporal.h"
#include <spin/Code.h>
#include <spin/PacketCode.h>
#include <algorithm>
#include <bit>
#include <chrono>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

namespace rs = spin::research::rs;
namespace candidate = spin::research::k16codesign;
namespace spin::research::k16codesign {
void outerScalar(const rs::Block*,rs::Block*,const rs::Plan&);
void outerForwardScalar(const rs::Block*,rs::Block*,const rs::Plan&);
}
using Block = rs::Block;
using Clock = std::chrono::steady_clock;
template<unsigned Mode> inline constexpr bool paired15Mode = Mode==30 || Mode==31 || Mode==36 || (Mode>=42 && Mode<=44) || Mode==46 || (Mode>=49 && Mode!=51);
static Block* blocks(spin::Buffer& b) { return reinterpret_cast<Block*>(b.bytes().data()); }
static void fill(spin::Buffer& b, std::uint64_t state=913) {
    for (std::size_t i=0;i<b.bytes().size();i+=8) {
        state+=0x9e3779b97f4a7c15ULL; auto v=state;
        v=(v^(v>>30))*0xbf58476d1ce4e5b9ULL;
        v=(v^(v>>27))*0x94d049bb133111ebULL; v^=v>>31;
        std::memcpy(b.bytes().data()+i,&v,8);
    }
}
static unsigned dot(const Block* a,const Block* b,std::size_t n) {
    unsigned p=0;
    for (std::size_t i=0;i<n;++i) {
        std::uint64_t x[2],y[2]; std::memcpy(x,a+i,16); std::memcpy(y,b+i,16);
        p^=std::popcount(x[0]&y[0])^std::popcount(x[1]&y[1]);
    }
    return p&1;
}
template<unsigned Mode>
static void route(const Block* x,Block* scratch,const rs::Plan& p,const candidate::Tables& tables,const candidate::T128Tables& large,const auto& simple,const candidate::PairedOptimizedTables& opt) {
    if constexpr (Mode==51) candidate::reverseRoutePairedFold(x,scratch,p,simple,opt);
    else if constexpr (Mode==52) candidate::reverseRoutePaired15Fold(x,scratch,p,simple,opt);
    else if constexpr (Mode>=49) candidate::reverseRoutePaired15NoAlias(x,scratch,p,simple,opt,Mode-49);
    else if constexpr (Mode>=47) candidate::reverseRoutePairedNoAlias(x,scratch,p,simple,opt,Mode-47);
    else if constexpr (Mode>=42 && Mode!=45) candidate::reverseRoutePaired15Shuffle(x,scratch,p,simple,opt,Mode==42?0:1);
    else if constexpr (Mode>=37) candidate::reverseRoutePairedShuffle(x,scratch,p,simple,opt,(Mode==41 || Mode==45)?2:Mode-37);
    else if constexpr (paired15Mode<Mode>) candidate::reverseRoutePaired15(x,scratch,p,simple,opt,Mode==30?0:1);
    else if constexpr (Mode>=28) candidate::reverseRoutePairedBasisRaw(x,scratch,p,simple,opt,Mode==29?1:0);
    else if constexpr (Mode>=22) candidate::reverseRoutePairedBasis(x,scratch,p,simple,opt,Mode==22?0:1);
    else if constexpr (Mode>=18) candidate::reverseRoutePairedWide(x,scratch,p,simple,opt,Mode-18);
    else if constexpr (Mode>=13) candidate::reverseRoutePairedOptimized(x,scratch,p,simple,opt,Mode-13);
    else if constexpr (Mode==12) candidate::reverseRoutePaired(x,scratch,p,simple);
    else if constexpr (Mode==11) candidate::reverseRouteMonomial(x,scratch,p,simple);
    else if constexpr (Mode==7 || Mode==8) candidate::reverseRouteT128(x,scratch,p,large);
    else if constexpr (Mode==0 || Mode==6 || Mode==9 || Mode==10) rs::reverseRouteFast(x,scratch,p);
    else if constexpr (Mode==4 || Mode==5) candidate::reverseRouteCachedField16(x,scratch,p,tables);
    else candidate::reverseRouteCachedGl16(x,scratch,p);
}
template<unsigned Mode>
static void outer(const Block* scratch,Block* x,const rs::Plan& p,const candidate::Tables& tables,const candidate::NativeOuterTables& native) {
    if constexpr (Mode==45 || Mode==46) candidate::outervariants::fieldLoopSharedNt(scratch,x,p,native);
    else if constexpr (Mode==35 || Mode==36 || Mode==41 || Mode==44 || Mode>=47) candidate::outervariants::fieldLoopSharedParity(scratch,x,p,native);
    else if constexpr (Mode==32) candidate::outervariants::fieldUnaryPack(scratch,x,p,native);
    else if constexpr (Mode==33) candidate::outervariants::fieldUnaryUnpack(scratch,x,p,native);
    else if constexpr (Mode==34) candidate::outervariants::fieldUnaryBoth(scratch,x,p,native);
    else if constexpr (Mode==24) candidate::outervariants::fieldInline(scratch,x,p,native);
    else if constexpr (Mode==25) candidate::outervariants::fieldUnrolled(scratch,x,p,native);
    else if constexpr (Mode==26) candidate::outervariants::fieldSharedParity(scratch,x,p,native);
    else if constexpr (Mode==27) candidate::outervariants::fieldBatch2(scratch,x,p,native);
    else if constexpr (Mode==9) candidate::outerNativeGl16(scratch,x,p,native);
    else if constexpr (Mode>=10) candidate::outerNativeField16(scratch,x,p,native);
    else if constexpr (Mode==0 || Mode==1 || Mode==7) rs::outerFast(scratch,x,p);
    else if constexpr (Mode==3 || Mode==5 || Mode==6 || Mode==8) candidate::outerFusedField16(scratch,x,p,tables);
    else candidate::outerFusedGl16(scratch,x,p);
}
template<unsigned Mode>
static void run(std::size_t k,std::uint64_t seed,unsigned calls,bool phases,bool huge) {
    rs::Plan plan(k,seed,seed,rs::Variant::Rs16Gf16,
        (Mode==0 || Mode==6 || Mode==9 || Mode==10)?rs::InnerKernel::Cached:rs::InnerKernel::RetainedStreaming);
    candidate::Tables tables;
    if constexpr (Mode==3 || Mode==5 || Mode==6 || Mode==8) candidate::customizeOuterField16(plan,seed,tables);
    if constexpr (Mode==4 || Mode==5) candidate::customizeInnerField16(plan,seed,tables);
    candidate::T128Tables large;
    if constexpr (Mode==7 || Mode==8) candidate::prepareT128(plan,large);
    std::conditional_t<(Mode>=12),candidate::PairedTables,candidate::MonomialTables> simple;
    candidate::PairedOptimizedTables opt;
    if constexpr (Mode==11) candidate::prepareMonomial(plan,simple);
    if constexpr (Mode==51) candidate::preparePairedShuffle(plan,simple,opt);
    else if constexpr (Mode>=49) candidate::customizePaired15Shuffle(plan,seed,simple,opt);
    else if constexpr (Mode==45 || Mode>=47) candidate::preparePairedShuffle(plan,simple,opt);
    else if constexpr (Mode>=42) candidate::customizePaired15Shuffle(plan,seed,simple,opt);
    else if constexpr (Mode==38 || Mode==40) candidate::customizePairedField16Shuffle(plan,seed,simple,opt);
    else if constexpr (Mode>=37) candidate::preparePairedShuffle(plan,simple,opt);
    else if constexpr (paired15Mode<Mode>) candidate::customizePaired15(plan,seed,simple,opt);
    else if constexpr (Mode>=23 && Mode<=29 && Mode!=28) candidate::customizePairedField16Basis(plan,seed,simple,opt);
    else if constexpr (Mode==22 || Mode==28 || Mode>=32) candidate::preparePairedBasis(plan,simple,opt);
    else if constexpr (Mode==16 || Mode==17 || Mode==19 || Mode==21)
        candidate::customizePairedField16(plan,seed,simple,opt);
    else if constexpr (Mode>=12) {
        candidate::preparePaired(plan,simple);
        if constexpr (Mode>=13) candidate::preparePairedOptimized(plan,opt);
    }
    candidate::NativeOuterTables native;
    if constexpr (Mode==9) candidate::prepareNativeGl16(plan,native);
    if constexpr (Mode>=10) candidate::customizeNativeField16(plan,seed,native);
    const auto policy=huge?spin::MemoryPolicy::PreferHugePages:spin::MemoryPolicy::Normal;
    spin::Buffer input((plan.n+12)*16,policy),output((plan.k+12)*16,policy),scratch(plan.scratchBlocks()*16,policy);
    spin::Buffer expected(plan.k*16),message(plan.k*16),encoded(plan.n*16);
    fill(input); fill(message,1949);
    const std::vector<std::byte> saved(input.bytes().begin(),input.bytes().end());
    const auto encode=[&](const Block* x,Block* y) { route<Mode>(x,blocks(scratch),plan,tables,large,simple,opt); outer<Mode>(blocks(scratch),y,plan,tables,native); };
    // Independent literal scalar transpose, separate buffers, alignment and
    // guards, in-place semantics, and forward/transpose dot-product identity.
    if constexpr (paired15Mode<Mode>) {
        candidate::reverseRoutePaired15Scalar(blocks(input)+4,blocks(scratch),plan);
        candidate::outerScalar(blocks(scratch),blocks(expected),plan);
    } else if constexpr (Mode>=12) {
        candidate::reverseRoutePairedScalar(blocks(input)+4,blocks(scratch),plan);
        candidate::outerScalar(blocks(scratch),blocks(expected),plan);
    } else if constexpr (Mode==11) {
        candidate::reverseRouteMonomialScalar(blocks(input)+4,blocks(scratch),plan);
        candidate::outerScalar(blocks(scratch),blocks(expected),plan);
    } else if constexpr (Mode==7 || Mode==8) {
        candidate::reverseRouteT128Scalar(blocks(input)+4,blocks(scratch),plan);
        candidate::outerScalar(blocks(scratch),blocks(expected),plan);
    } else rs::transposeScalar(blocks(input)+4,blocks(expected),blocks(scratch),plan);
    for (unsigned offset=0;offset<4;++offset) {
        std::fill(output.bytes().begin(),output.bytes().end(),std::byte{0xa5});
        std::memcpy(blocks(input)+4+offset,saved.data()+64,plan.n*16);
        auto* out=blocks(output)+4+offset;
        encode(blocks(input)+4+offset,out);
        if (std::memcmp(out,blocks(expected),plan.k*16)) throw std::runtime_error("scalar mismatch");
        if (!std::all_of(output.bytes().begin(),output.bytes().begin()+(4+offset)*16,[](std::byte x){return x==std::byte{0xa5};}) ||
            !std::all_of(output.bytes().begin()+(4+offset+plan.k)*16,output.bytes().end(),[](std::byte x){return x==std::byte{0xa5};}))
            throw std::runtime_error("output guard mismatch");
        if (std::memcmp(blocks(input)+4+offset,saved.data()+64,plan.n*16)) throw std::runtime_error("input modified");
    }
    std::memcpy(input.bytes().data(),saved.data(),saved.size());
    if constexpr (paired15Mode<Mode>) {
        candidate::outerForwardScalar(blocks(message),blocks(scratch),plan);
        candidate::forwardInnerPaired15Scalar(blocks(scratch),blocks(encoded),plan);
    } else if constexpr (Mode>=12) {
        candidate::outerForwardScalar(blocks(message),blocks(scratch),plan);
        candidate::forwardInnerPairedScalar(blocks(scratch),blocks(encoded),plan);
    } else if constexpr (Mode==11) {
        candidate::outerForwardScalar(blocks(message),blocks(scratch),plan);
        candidate::forwardInnerMonomialScalar(blocks(scratch),blocks(encoded),plan);
    } else if constexpr (Mode==7 || Mode==8) {
        candidate::outerForwardScalar(blocks(message),blocks(scratch),plan);
        candidate::forwardInnerT128Scalar(blocks(scratch),blocks(encoded),plan);
    } else rs::forwardScalar(blocks(message),blocks(encoded),plan);
    if (dot(blocks(encoded),blocks(input)+4,plan.n)!=dot(blocks(message),blocks(expected),plan.k)) throw std::runtime_error("adjoint mismatch");
    encode(blocks(input)+4,blocks(input)+4);
    if (std::memcmp(blocks(input)+4,blocks(expected),plan.k*16) ||
        std::memcmp(blocks(input)+4+plan.k,saved.data()+64+plan.k*16,(plan.n-plan.k)*16))
        throw std::runtime_error("in-place mismatch");
    if (!calls) { std::cout<<"PASS mode="<<Mode<<" K="<<k<<" seed="<<seed<<'\n'; return; }
    // Whole encoding is the headline; stage clocks are optional diagnostics.
    fill(input); auto* x=blocks(input)+4;
    for (unsigned i=0;i<5;++i) encode(x,x);
    std::vector<double> total(calls),inner(calls),out(calls);
    for (unsigned i=0;i<calls;++i) {
        const auto a=Clock::now(); route<Mode>(x,blocks(scratch),plan,tables,large,simple,opt);
        const auto b=phases?Clock::now():a; outer<Mode>(blocks(scratch),x,plan,tables,native);
        const auto c=Clock::now();
        total[i]=std::chrono::duration<double,std::micro>(c-a).count();
        if (phases) { inner[i]=std::chrono::duration<double,std::micro>(b-a).count();out[i]=std::chrono::duration<double,std::micro>(c-b).count(); }
    }
    std::sort(total.begin(),total.end());std::sort(inner.begin(),inner.end());std::sort(out.begin(),out.end());
    std::uint64_t hash=0;
    for(std::size_t i=0;i<plan.k*16;i+=8) {std::uint64_t v;std::memcpy(&v,reinterpret_cast<std::byte*>(x)+i,8);hash=(hash^v)*0x100000001b3ULL;}
    std::cout<<"mode,K,seed,calls,median_us,p10_us,p90_us,inner_us,outer_us,checksum,huge\n"
        <<Mode<<','<<k<<','<<seed<<','<<calls<<','<<std::fixed<<std::setprecision(4)
        <<total[calls/2]<<','<<total[calls/10]<<','<<total[9*calls/10]<<','<<inner[calls/2]<<','<<out[calls/2]<<','<<std::hex<<hash<<std::dec<<','<<huge<<'\n';
}
int main(int argc,char** argv) {try {
    if (!spin::packet_fast_available() || !spin::capabilities().forward512) return 77;
    if (argc<3 || argc>7) throw std::invalid_argument("mode K [seed=1] [calls=501;0=check] [phases=0] [huge=0]");
    auto m=std::stoul(argv[1]);auto k=std::stoull(argv[2]);auto seed=argc>3?std::stoull(argv[3]):1;
    auto calls=argc>4?std::stoul(argv[4]):501;bool phases=argc>5?std::stoul(argv[5])!=0:false;
    bool huge=argc>6?std::stoul(argv[6])!=0:false;
    if(calls>1000000)throw std::invalid_argument("too many calls");
#define CASE(M) case M:run<M>(k,seed,calls,phases,huge);break
    switch(m){CASE(0);CASE(1);CASE(2);CASE(3);CASE(4);CASE(5);CASE(6);CASE(7);CASE(8);CASE(9);CASE(10);CASE(11);CASE(12);CASE(13);CASE(14);CASE(15);CASE(16);CASE(17);CASE(18);CASE(19);CASE(20);CASE(21);CASE(22);CASE(23);CASE(24);CASE(25);CASE(26);CASE(27);CASE(28);CASE(29);CASE(30);CASE(31);CASE(32);CASE(33);CASE(34);CASE(35);CASE(36);CASE(37);CASE(38);CASE(39);CASE(40);CASE(41);CASE(42);CASE(43);CASE(44);CASE(45);CASE(46);CASE(47);CASE(48);CASE(49);CASE(50);CASE(51);CASE(52);default:throw std::invalid_argument("mode");}
#undef CASE
    return 0;
} catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
