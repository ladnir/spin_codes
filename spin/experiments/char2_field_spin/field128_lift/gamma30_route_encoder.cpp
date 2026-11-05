#include "Gamma30RouteEncoder.h"
#include "Gamma30TestSupport.h"
#include <type_traits>

// Correctness only. Frozen baseline and completion variants are compared on
// one identical Plan; independent dense maps do not use the fused route implementation.
namespace enc=spin::research::char2::field128_lift::gamma30_route;
namespace baseline=enc::baseline;
namespace previous=enc::previous;
namespace literal=enc::literal;
using namespace spin::research::char2::field128_lift::gamma30_test;
static_assert(std::is_same_v<enc::Plan,baseline::Plan>);
static_assert(std::is_same_v<enc::Plan,previous::Plan>);

static void checkOracleAndGeometry(const Oracle& oracle){
    require(oracle.mul[0x57][0x83]==0xc1,"AES known multiplication vector");
    for(unsigned i=1;i!=256;++i)require(oracle.mul[i][oracle.inverse[i]]==1,"AES inverse table");
    require(feedbackGamma==field::towerY,"gamma must be Y, not a byte-field constant");
    for(unsigned i=0;i!=64;++i){
        std::uint8_t sum=0;for(auto a:oracle.a[i])sum^=a;
        require(sum==1,"A preserves constant polynomials");
    }
    for(unsigned j=0;j!=30;++j)for(unsigned k=0;k!=30;++k){
        std::uint8_t dot=0;
        for(unsigned i=0;i!=64;++i)dot^=oracle.mul[oracle.a[i][j]][oracle.a[i][k]];
        require(dot==0,"dense A^T A must vanish");
    }
    auto h=oracle.h;
    for(unsigned c=0;c!=30;++c){
        unsigned pivot=c;while(pivot!=30&&!h[pivot][c])++pivot;
        require(pivot!=30,"dense H30 must be nonsingular");std::swap(h[c],h[pivot]);
        const auto inv=oracle.inverse[h[c][c]];
        for(unsigned j=c;j!=30;++j)h[c][j]=oracle.mul[h[c][j]][inv];
        for(unsigned i=c+1;i!=30;++i){
            const auto factor=h[i][c];
            for(unsigned j=c;j!=30;++j)h[i][j]^=oracle.mul[factor][h[c][j]];
        }
    }
    std::mt19937_64 random(0x2936aa150372dbf7ULL);
    for(unsigned sample=0;sample!=64;++sample){
        alignas(64) Element a[4],b[4],actual[4],actualY[4];
        for(unsigned lane=0;lane!=4;++lane){a[lane]=randomElement(random);b[lane]=randomElement(random);}
        const auto x=_mm512_load_si512(a);
        _mm512_store_si512(actual,field::multiply4(x,_mm512_load_si512(b)));
        _mm512_store_si512(actualY,field::multiplyY4(x));
        for(unsigned lane=0;lane!=4;++lane){
            require(actual[lane]==oracle.product(a[lane],b[lane]),"full product versus independent dense oracle");
            require(actualY[lane]==oracle.product(feedbackGamma,a[lane]),"Y shortcut versus generic polynomial division");
        }
    }
    std::cout<<"PASS independent AES/tower arithmetic, dense CA=0 and H30 rank\n";
}

static void checkFixedMaps(const Oracle& oracle){
    const auto checkState=[&](const State30& q){
        const auto state=packedState(q,true);
        auto y=packedPacket(Packet64{});enc::emitAdd(state,y);
        same(naturalPacket(y),oracle.emit(q),"A FFT map versus dense Lagrange");
        auto h=state;enc::refreshMap(h);
        same(naturalState(h),oracle.refresh(q),"H FFT map versus dense Lagrange");
        auto cy=y;const auto feedback=enc::feedback(cy);
        same(naturalState(feedback),State30{},"C*A zero including poisoned storage padding");
    };
    checkState(State30{});
    for(unsigned j=0;j!=30;++j)for(unsigned bit=0;bit!=8;++bit){
        State30 q{};q[j]={std::uint64_t(1)<<bit,0};checkState(q);
    }
    // Exercise every binary field coordinate as well as every state slot.
    for(unsigned bit=0;bit!=128;++bit){
        State30 q{};q[bit%30]=field::binaryBasis(bit);checkState(q);
    }
    for(unsigned i=0;i!=64;++i){
        Packet64 x{};x[i]=one;auto packed=packedPacket(x);
        same(naturalState(enc::feedback(packed)),oracle.feedback(x),"C transpose basis versus direct dense C");
    }
    std::mt19937_64 random(0x39018488cbf21356ULL);
    for(unsigned sample=0;sample!=8;++sample){
        State30 q{};Packet64 x{};
        for(auto& value:q)value=randomElement(random);
        for(auto& value:x)value=randomElement(random);
        checkState(q);
        auto packed=packedPacket(x);
        same(naturalState(enc::feedback(packed)),oracle.feedback(x),"C full-field random packet");
        auto y=packedPacket(x);enc::emitAdd(packedState(q,true),y);
        auto expected=oracle.emit(q);for(unsigned i=0;i!=64;++i)expected[i]^=x[i];
        same(naturalPacket(y),expected,"A add-to-payload semantics");
    }
    std::cout<<"PASS exact A/C/H dense basis/random maps and inert padding\n";
}

static void checkPacketRoutes(){
    c2::BenchPayload input(68,false),output(72,false);
    for(unsigned i=0;i!=64;++i)input.p[i]={0x91300000ULL+i,0xb6000000ULL+7*i};
    input.p[64]=input.p[0];
    for(unsigned parity=0;parity!=2;++parity)for(unsigned pattern=0;pattern!=4;++pattern){
        // For each J, both packed pairs independently take either orientation.
        std::array<std::uint16_t,32> route{};
        for(unsigned j=0;j!=32;++j)route[j]=std::uint16_t(2*j+parity+(((pattern>>(j%2))&1U)<<15));
        const auto packet=enc::loadPacket(input.p,route.data());
        const auto natural=naturalPacket(packet);
        for(unsigned j=0;j!=16;++j)for(unsigned lane=0;lane!=4;++lane){
            const unsigned entry=route[2*j+lane/2];
            const unsigned from=((entry&0x7fffU)+((lane%2)^(entry>>15)))%64;
            require(natural[16*lane+j]==element(input.p[from]),"paired route/parity/swap/wrap load mapping");
        }
        for(unsigned offset=0;offset!=4;++offset){
            std::fill_n(output.p,72,guard);enc::storePacket<false>(packet,output.p+offset);
            for(unsigned i=0;i!=64;++i)
                require(element(output.p[offset+4*(i%16)+i/16])==natural[i],"native output coordinate permutation");
            for(unsigned i=0;i!=offset;++i)require(output.p[i]==guard,"packet leading guard");
            for(unsigned i=offset+64;i!=72;++i)require(output.p[i]==guard,"packet trailing guard");
        }
        enc::storePacket<true>(packet,output.p);_mm_sfence();
        for(unsigned i=0;i!=64;++i)require(element(output.p[4*(i%16)+i/16])==natural[i],"NT packet store");
    }
    std::cout<<"PASS four swap combinations, both route parities, wrap and output layout\n";
}

static void checkFusedLoads(){
    c2::BenchPayload input(132,false);
    // Every eight-entry route group exercises all independent pair flips.
    // Compare both the unchanged completion loader and a scalar coordinate map.
    for(unsigned offset=0;offset!=4;++offset){
        auto* source=input.p+offset;produce(source,128,913+offset);source[128]=source[0];
        for(unsigned parity=0;parity!=2;++parity)for(unsigned repeated=0;repeated!=2;++repeated)
        for(unsigned flips=0;flips!=256;++flips){
            std::array<std::uint16_t,33> storage{};
            auto* route=storage.data()+1; // No vector alignment promise.
            for(unsigned i=0;i!=32;++i){
                const unsigned pair=repeated?((i/2)*7)%64:(13*i+17)%64;
                route[i]=std::uint16_t(2*pair+parity+(((flips>>(i%8))&1U)<<15));
            }
            route[31]=std::uint16_t(126+parity+(((flips>>7)&1U)<<15)); // Odd wrap.
            const auto actual=naturalPacket(enc::loadPacket(source,route));
            const auto control=naturalPacket(previous::loadPacket(source,route));
            same(actual,control,"fused route load equals frozen completion load");
            for(unsigned j=0;j!=16;++j)for(unsigned lane=0;lane!=4;++lane){
                const unsigned entry=route[2*j+lane/2];
                const unsigned from=(entry&0x7fffU)+((lane%2)^(entry>>15));
                require(actual[16*lane+j]==element(source[from]),"fused loader scalar route oracle");
            }
        }
    }
    std::cout<<"PASS fused loads: all256 swap patterns/group, both parities, wrap, permuted/repeated pairs and four alignments\n";
}

static void checkRecurrence(const Oracle& oracle){
    std::mt19937_64 random(0xe331d8cb934a6791ULL);
    State30 state{};
    auto actual=packedState(state,true);
    for(unsigned step=0;step!=17;++step){
        Packet64 x{};enc::Refresh diagonal{};
        for(auto& value:x)value=step==0?Element{}:randomElement(random);
        for(unsigned j=0;j!=30;++j){
            do{diagonal.before[j]=randomElement(random);}while(diagonal.before[j]==Element{});
            do{diagonal.after[j]=randomElement(random);}while(diagonal.after[j]==Element{});
        }
        diagonal.before[30]=diagonal.before[31]=diagonal.after[30]=diagonal.after[31]=one;
        if(step==1){for(auto& value:state)value=randomElement(random);actual=packedState(state,true);}
        const auto expected=oracle.step(x,state,diagonal);
        auto y=packedPacket(x);enc::emitAdd(actual,y);
        same(naturalPacket(y),expected,"stateful emission versus dense oracle");
        enc::finish(actual,enc::feedback(y),diagonal);
        same(naturalState(actual),state,"D2 H D1(Q+Y C X) exact recurrence");
    }
    std::cout<<"PASS persistent gamma30 recurrence, nonzero full-field diagonals and zero packet\n";
}

static void checkPlan(const enc::Plan& plan,const Oracle& oracle,std::uint64_t seed){
    const auto& p=plan.outerPlan;
    // Independent instance of the frozen plan: retain exact legacy samples,
    // transformed key order and route, not merely a statistically similar map.
    const literal::Plan control(p.k,seed,false);
    require(p.toeplitz==control.toeplitz&&p.outerKeys==control.outerKeys,"outer keys changed from frozen literal plan");
    require(p.route==control.route&&p.packetRoute==control.packetRoute,"route changed from frozen literal plan");
    require(p.before.empty()&&p.after.empty()&&p.refresh.empty(),"unused state16 diagonal storage retained");
    require(!p.restrictedToGf256,"outer setup unexpectedly restricted to GF256");
    require(plan.refresh.size()==p.n/64,"one independent diagonal pair per inner packet");
    bool fullField=false;
    for(const auto& d:plan.refresh){
        for(unsigned j=0;j!=30;++j){
            require(d.before[j]!=Element{}&&d.after[j]!=Element{},"zero gamma30 diagonal");
            fullField|=bool(d.before[j].hi||(d.before[j].lo>>8)||d.after[j].hi||(d.after[j].lo>>8));
        }
        require(d.before[30]==one&&d.before[31]==one&&d.after[30]==one&&d.after[31]==one,"padding keys must be one");
    }
    require(fullField,"gamma30 diagonal samples unexpectedly byte-only");
    for(unsigned region=0;region!=256;++region){
        std::vector<bool> seen(p.rows,false);
        for(std::size_t i=0;i!=p.rows;++i){
            const auto index=p.route[region*p.rows+i];
            require(index<p.rows&&!seen[index],"regional route is not a permutation");seen[index]=true;
        }
        std::vector<bool> paired(p.rows/2,false);
        for(std::size_t j=0;j!=p.rows/2;++j){
            const unsigned entry=p.packetRoute[region*(p.rows/2)+j],base=entry&0x7fffU;
            require(base<p.rows&&(base&1U)==(region&1U),"regional pairing parity");
            require(!paired[base/2],"regional pair permutation repeats a pair");paired[base/2]=true;
        }
    }
    // A sample of accepted matrices is independently checked by a full
    // permutation determinant, not the setup's fraction-free elimination.
    for(unsigned i=0;i!=16;++i)require(oracle.determinant(p.toeplitz[(i*257)%p.toeplitz.size()])!=Element{},"accepted Toeplitz singular");
}

static void checkWhole(const Oracle& oracle,std::size_t k,std::uint64_t seed){
    const enc::Plan plan(k,seed);const auto& p=plan.outerPlan;checkPlan(plan,oracle,seed);
    c2::BenchPayload input(p.k,false),output(p.n+8,false),expected(p.n,false);
    c2::BenchPayload scratch(p.scratchSize()+8,false),referencePlanes(p.scratchSize(),false);
    produce(input.p,p.k,seed^0x985211c971324fffULL);
    oracle.outerReference(plan,input.p,referencePlanes.p);
    oracle.innerReference(plan,referencePlanes.p,expected.p);
    baseline::forward(plan,input.p,output.p,scratch.p);
    same(output.p,expected.p,p.n,"baseline whole encoder on identical Plan");
    same(scratch.p,referencePlanes.p,p.scratchSize(),"baseline outer on identical Plan");
    previous::forward(plan,input.p,output.p,scratch.p);
    same(output.p,expected.p,p.n,"completion predecessor on identical Plan");
    enc::outer<false>(p,input.p,scratch.p);
    same(scratch.p,referencePlanes.p,p.scratchSize(),"unchanged outer dense oracle and wrap padding");
    enc::inner<false>(plan,scratch.p,output.p);
    same(output.p,expected.p,p.n,"inner-only dense oracle, continuous state across regions");
    for(unsigned offset=0;offset!=4;++offset){
        std::fill_n(output.p,p.n+8,guard);std::fill_n(scratch.p,p.scratchSize()+8,guard);
        enc::forward(plan,input.p,output.p+offset,scratch.p+offset);
        same(output.p+offset,expected.p,p.n,"route-fused whole encoder dense oracle / alignment fallback");
        same(scratch.p+offset,referencePlanes.p,p.scratchSize(),"deferred NT outer / alignment / padding");
        for(unsigned i=0;i!=offset;++i){require(output.p[i]==guard,"whole output leading guard");require(scratch.p[i]==guard,"scratch leading guard");}
        for(std::size_t i=p.n+offset;i!=p.n+8;++i)require(output.p[i]==guard,"whole output trailing guard");
        for(std::size_t i=p.scratchSize()+offset;i!=p.scratchSize()+8;++i)require(scratch.p[i]==guard,"scratch trailing guard");
    }
    enc::forward<false>(plan,input.p,output.p,scratch.p);
    same(output.p,expected.p,p.n,"cached output whole encoder");
    std::copy_n(input.p,p.k,output.p);enc::forward(plan,output.p,output.p,scratch.p);
    same(output.p,expected.p,p.n,"exact input/output alias");
    std::copy_n(input.p,p.k,output.p+3);enc::forward(plan,output.p+3,output.p,scratch.p);
    same(output.p,expected.p,p.n,"backward partial alias and unaligned input");
    std::copy_n(input.p,p.k,output.p);enc::forward(plan,output.p,output.p+1,scratch.p);
    same(output.p+1,expected.p,p.n,"forward partial input/output alias");
    std::fill_n(input.p,p.k,Block{});enc::forward(plan,input.p,output.p,scratch.p);
    for(std::size_t i=0;i!=p.n;++i)require(output.p[i]==Block{},"whole zero input");
    std::cout<<"PASS dense whole encoder K="<<k<<" seed="<<seed
             <<" / aligned+unaligned / cached+NT / guards / exact+partial alias / zero\n";
}

static void checkLinearity(const Oracle& oracle){
    const enc::Plan plan(8192,819113);const auto& p=plan.outerPlan;
    c2::BenchPayload u(p.k,false),v(p.k,false),mix(p.k,false),eu(p.n,false),ev(p.n,false),emix(p.n,false),scratch(p.scratchSize(),false);
    produce(u.p,p.k,591);produce(v.p,p.k,8173);
    const Element lambda{0x0305070b0d111317ULL,0x1d1f25292b2f353bULL};
    for(std::size_t i=0;i!=p.k;++i)mix.p[i]=block(oracle.product(lambda,element(u.p[i]))^element(v.p[i]));
    enc::forward(plan,u.p,eu.p,scratch.p);enc::forward(plan,v.p,ev.p,scratch.p);
    enc::forward(plan,mix.p,emix.p,scratch.p);
    for(std::size_t i=0;i!=p.n;++i)
        require(element(emix.p[i])==(oracle.product(lambda,element(eu.p[i]))^element(ev.p[i])),"full-field linearity");
    std::cout<<"PASS full-field linearity for lambda outside GF256\n";
}

int main(int argc,char** argv){
    try{
        if(argc>2||(argc==2&&std::string(argv[1])!="--check")){
            std::cerr<<"usage: gamma30_route_encoder [--check]\nCorrectness only; no benchmark mode.\n";return 2;
        }
        std::cout<<"Gamma30 route-load A/B correctness checks; deterministic MT research setup, not a cryptographic sampler.\n";
        const Oracle oracle;
        checkOracleAndGeometry(oracle);checkFixedMaps(oracle);checkPacketRoutes();checkFusedLoads();checkRecurrence(oracle);
        // One and two packets per region respectively; K16384 also crosses
        // an outer tile boundary. Neither check is a K20 performance run.
        checkWhole(oracle,8192,113);checkWhole(oracle,16384,7919);checkLinearity(oracle);
        std::cout<<"PASS all gamma30 route encoder checks and same-Plan baseline/completion equivalence\n";return 0;
    }catch(const std::exception& error){std::cerr<<"FAIL: "<<error.what()<<"\n";return 1;}
}
