#include "Gamma30Encoder.h"
#include "../BenchPayload.h"
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

// Correctness driver only: no timing path and no hidden benchmark workload.
// The dense reference below shares neither the encoder's FFT maps nor its
// SIMD field arithmetic. BenchPayload is used only as an aligned allocator.
namespace enc=spin::research::char2::field128_lift::gamma30;
namespace literal=enc::literal;
namespace field=enc::field;
namespace c2=spin::research::char2;
using Element=enc::Element;
using Block=enc::Block;
using State30=std::array<Element,30>;
using Packet64=std::array<Element,64>;

static void require(bool condition,const char* label){
    if(!condition)throw std::runtime_error(label);
}
static Element element(Block x){return {x.lo,x.hi};}
static Block block(Element x){return {x.lo,x.hi};}
static Element randomElement(std::mt19937_64& random){return {random(),random()};}
static constexpr Element one{1,0},feedbackGamma{0x100,0};
static constexpr Block guard{0x2e71d6935a90bf24ULL,0xc031674f27e18ac5ULL};

static void same(const Block* actual,const Block* expected,std::size_t n,const char* label){
    for(std::size_t i=0;i!=n;++i)if(actual[i]!=expected[i])
        throw std::runtime_error(std::string(label)+": record "+std::to_string(i)+" differs");
}
template<std::size_t N>
static void same(const std::array<Element,N>& actual,const std::array<Element,N>& expected,const char* label){
    for(std::size_t i=0;i!=N;++i)if(actual[i]!=expected[i])
        throw std::runtime_error(std::string(label)+": coordinate "+std::to_string(i)+" differs");
}
static void produce(Block* output,std::size_t n,std::uint64_t seed){
    std::mt19937_64 random(seed);
    for(std::size_t i=0;i!=n;++i)output[i]=block(randomElement(random));
}

struct Oracle {
    std::array<std::array<std::uint8_t,256>,256> mul{};
    std::array<std::uint8_t,256> inverse{};
    std::array<std::array<std::uint8_t,32>,64> outer{};
    std::array<std::array<std::uint8_t,30>,64> a{};
    std::array<std::array<std::uint8_t,30>,30> h{};

    // Independent bit-serial AES multiplication, not field::byteMultiply.
    static std::uint8_t aes(unsigned left,unsigned right){
        unsigned result=0;
        for(unsigned bit=0;bit!=8;++bit){
            if(right&1U)result^=left;
            right>>=1;left=(left<<1)^((left&128U)?0x11bU:0U);
        }
        return static_cast<std::uint8_t>(result);
    }
    static std::uint8_t byte(Element x,unsigned i){
        return static_cast<std::uint8_t>((i<8?x.lo:x.hi)>>(8*(i%8)));
    }
    static Element pack(const std::array<std::uint8_t,16>& bytes){
        Element result{};
        for(unsigned i=0;i!=8;++i){
            result.lo|=std::uint64_t(bytes[i])<<(8*i);
            result.hi|=std::uint64_t(bytes[8+i])<<(8*i);
        }
        return result;
    }
    Oracle(){
        for(unsigned i=0;i!=256;++i)for(unsigned j=0;j!=256;++j){
            mul[i][j]=aes(i,j);
            if(mul[i][j]==1)inverse[i]=static_cast<std::uint8_t>(j);
        }
        for(unsigned i=0;i!=64;++i){
            for(unsigned j=0;j!=32;++j)outer[i][j]=novelBasis(j,i);
            for(unsigned j=0;j!=30;++j)a[i][j]=lagrange(128,30,j,i);
        }
        for(unsigned i=0;i!=30;++i)for(unsigned j=0;j!=30;++j)
            h[i][j]=lagrange(0,30,j,30+i);
    }
    std::uint8_t lagrange(unsigned first,unsigned count,unsigned column,unsigned point) const{
        std::uint8_t numerator=1,denominator=1;
        for(unsigned j=0;j!=count;++j)if(j!=column){
            numerator=mul[numerator][point^(first+j)];
            denominator=mul[denominator][(first+column)^(first+j)];
        }
        require(denominator!=0,"oracle Lagrange nodes must be distinct");
        return mul[numerator][inverse[denominator]];
    }
    std::uint8_t novelBasis(unsigned degree,unsigned point) const{
        std::uint8_t result=1;
        for(unsigned bit=0;bit!=5;++bit)if((degree>>bit)&1U){
            std::uint8_t value=1,normalizer=1;
            for(unsigned j=0;j!=(1U<<bit);++j){
                value=mul[value][point^j];normalizer=mul[normalizer][(1U<<bit)^j];
            }
            result=mul[result][mul[value][inverse[normalizer]]];
        }
        return result;
    }
    Element scale(Element x,std::uint8_t coefficient) const{
        std::array<std::uint8_t,16> bytes{};
        for(unsigned j=0;j!=16;++j)bytes[j]=mul[coefficient][byte(x,j)];
        return pack(bytes);
    }
    Element product(Element left,Element right) const{
        std::array<std::uint8_t,31> coefficients{};
        for(unsigned i=0;i!=16;++i)for(unsigned j=0;j!=16;++j)
            coefficients[i+j]^=mul[byte(left,i)][byte(right,j)];
        // Independent generic long division by the explicit tower modulus.
        // Do not call either the SIMD reduction or the multiply-by-Y helper.
        constexpr std::array<std::uint8_t,17> modulus{
            0x32,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,1};
        for(unsigned degree=31;degree-- >16;){
            const auto leading=coefficients[degree];
            for(unsigned j=0;j!=17;++j)coefficients[degree-16+j]^=mul[leading][modulus[j]];
        }
        std::array<std::uint8_t,16> bytes{};
        std::copy_n(coefficients.begin(),16,bytes.begin());return pack(bytes);
    }
    Packet64 emit(const State30& state) const{
        Packet64 result{};
        for(unsigned i=0;i!=64;++i)for(unsigned j=0;j!=30;++j)result[i]^=scale(state[j],a[i][j]);
        return result;
    }
    State30 feedback(const Packet64& input) const{
        State30 result{};
        for(unsigned j=0;j!=30;++j)for(unsigned i=0;i!=64;++i)result[j]^=scale(input[i],a[i][j]);
        return result;
    }
    State30 refresh(const State30& state) const{
        State30 result{};
        for(unsigned i=0;i!=30;++i)for(unsigned j=0;j!=30;++j)result[i]^=scale(state[j],h[i][j]);
        return result;
    }
    State30 finish(const State30& state,const State30& syndrome,const enc::Refresh& diagonal) const{
        State30 before{};
        for(unsigned j=0;j!=30;++j)
            before[j]=product(state[j]^product(feedbackGamma,syndrome[j]),diagonal.before[j]);
        auto result=refresh(before);
        for(unsigned j=0;j!=30;++j)result[j]=product(result[j],diagonal.after[j]);
        return result;
    }
    Packet64 step(const Packet64& input,State30& state,const enc::Refresh& diagonal) const{
        auto output=emit(state);
        for(unsigned i=0;i!=64;++i)output[i]^=input[i];
        // Consume C*x, deliberately independent of the hot path's C*y trick.
        state=finish(state,feedback(input),diagonal);
        return output;
    }
    Element determinant(const std::array<Element,7>& t) const{
        // Characteristic two: determinant equals the unsigned permutation sum.
        std::array<unsigned,4> permutation{0,1,2,3};Element result{};
        do{
            auto term=one;
            for(unsigned i=0;i!=4;++i)term=product(term,t[3+i-permutation[i]]);
            result^=term;
        }while(std::next_permutation(permutation.begin(),permutation.end()));
        return result;
    }
    void outerReference(const enc::Plan& plan,const Block* input,Block* planes) const{
        const auto& p=plan.outerPlan;
        for(std::size_t group=0;group!=p.rows;++group)for(unsigned c=0;c!=64;++c){
            std::array<Element,4> column{};
            for(unsigned r=0;r!=4;++r)for(unsigned d=0;d!=32;++d)
                column[r]^=scale(element(input[(group/4)*512+4*(r*32+d)+group%4]),outer[c][d]);
            const auto& t=p.toeplitz[c*p.rows+group];
            for(unsigned r=0;r!=4;++r){
                Element value{};
                for(unsigned j=0;j!=4;++j)value^=product(column[j],t[3+r-j]);
                planes[(4*c+r)*p.stride+group]=block(value);
            }
        }
        for(unsigned region=0;region!=256;++region){
            auto* plane=planes+region*p.stride;
            plane[p.rows]=plane[0];std::fill_n(plane+p.rows+1,3,Block{});
        }
    }
    void innerReference(const enc::Plan& plan,const Block* planes,Block* output) const{
        const auto& p=plan.outerPlan;
        State30 state{}; // One state for all packets and all 256 regions.
        for(std::size_t first=0;first!=p.n;first+=64){
            Packet64 input{};
            for(unsigned i=0;i!=64;++i)
                input[i]=element(planes[(first/p.rows)*p.stride+p.route[first+i]]);
            const auto y=step(input,state,plan.refresh[first/64]);
            for(unsigned i=0;i!=64;++i)output[first+4*(i%16)+i/16]=block(y[i]);
        }
    }
};

static enc::State packedState(const State30& state,bool poisonPadding=false){
    alignas(64) std::array<Element,32> storage{};
    std::copy(state.begin(),state.end(),storage.begin());
    if(poisonPadding){storage[30]={~0ULL,0x701239ca45b8ddedULL};storage[31]={0x528921a132456789ULL,~0ULL};}
    enc::State result;
    for(unsigned i=0;i!=8;++i)result.values[i]=_mm512_load_si512(storage.data()+4*i);
    return result;
}
static State30 naturalState(const enc::State& state,bool requirePadding=true){
    alignas(64) std::array<Element,32> storage{};
    for(unsigned i=0;i!=8;++i)_mm512_store_si512(storage.data()+4*i,state.values[i]);
    if(requirePadding)require(storage[30]==Element{}&&storage[31]==Element{},"state padding must be zero");
    State30 result{};std::copy_n(storage.begin(),30,result.begin());return result;
}
static enc::Packet packedPacket(const Packet64& input){
    enc::Packet result;
    for(unsigned i=0;i!=16;++i)result.values[i]=_mm512_loadu_si512(input.data()+4*i);
    return result;
}
static Packet64 naturalPacket(const enc::Packet& packet){
    Packet64 result{};
    for(unsigned i=0;i!=16;++i)_mm512_storeu_si512(result.data()+4*i,packet.values[i]);
    return result;
}

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
    literal::outer<false>(p,input.p,scratch.p);
    same(scratch.p,referencePlanes.p,p.scratchSize(),"unchanged outer dense oracle and wrap padding");
    enc::inner<false>(plan,scratch.p,output.p);
    same(output.p,expected.p,p.n,"inner-only dense oracle, continuous state across regions");
    for(unsigned offset=0;offset!=4;++offset){
        std::fill_n(output.p,p.n+8,guard);std::fill_n(scratch.p,p.scratchSize()+8,guard);
        enc::forward(plan,input.p,output.p+offset,scratch.p+offset);
        same(output.p+offset,expected.p,p.n,"whole encoder dense oracle / alignment fallback");
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
            std::cerr<<"usage: gamma30_encoder [--check]\nCorrectness only; no benchmark mode.\n";return 2;
        }
        std::cout<<"Gamma30 isolated correctness checks; deterministic MT research setup, not a cryptographic sampler.\n";
        const Oracle oracle;
        checkOracleAndGeometry(oracle);checkFixedMaps(oracle);checkPacketRoutes();checkRecurrence(oracle);
        // One and two packets per region respectively; K16384 also crosses
        // an outer tile boundary. Neither check is a K20 performance run.
        checkWhole(oracle,8192,113);checkWhole(oracle,16384,7919);checkLinearity(oracle);
        std::cout<<"PASS all gamma30 encoder checks\n";return 0;
    }catch(const std::exception& error){std::cerr<<"FAIL: "<<error.what()<<"\n";return 1;}
}
