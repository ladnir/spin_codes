#include "LiteralEncoder.h"
#include "../BenchPayload.h"
#include <atomic>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>

namespace enc=spin::research::char2::field128_lift::encoder;
namespace field=enc::field;
namespace old=enc::old;
namespace sr=enc::sr;
namespace c2=spin::research::char2;
using Element=field::Element;
using Block=enc::Block;
using Clock=std::chrono::steady_clock;

static void require(bool condition,const char* message){if(!condition)throw std::runtime_error(message);}
static Element element(Block a){return {a.lo,a.hi};}
static Block block(Element a){return {a.lo,a.hi};}
static void produce(Block* output,std::size_t n,std::uint64_t seed=913){
    std::uint64_t x=seed,y=113;
    for(std::size_t i=0;i!=n;++i){x^=x<<13;x^=x>>7;x^=x<<17;y+=0x9e3779b97f4a7c15ULL;output[i]={x,y};}
}
static std::uint64_t hash(const Block* p,std::size_t n){
    std::uint64_t h=0;
    for(std::size_t i=0;i!=n;++i){h=(h^p[i].lo)*0x100000001b3ULL;h=(h^p[i].hi)*0x100000001b3ULL;}
    return h;
}
static void same(const Block* actual,const Block* expected,std::size_t n,const char* label){
    for(std::size_t i=0;i!=n;++i)if(actual[i]!=expected[i])
        throw std::runtime_error(std::string(label)+": record "+std::to_string(i)+" differs");
}

// All oracle paths below use scalar AES-polynomial arithmetic and dense
// matrices. They share neither FFT/transpose schedules nor CLMUL reduction
// with the encoder. The small multiplication table only speeds up checking.
struct Oracle {
    std::array<std::array<std::uint8_t,256>,256> mul{};
    std::array<std::uint8_t,256> inverse{};
    std::array<std::array<std::uint8_t,32>,64> outer{};
    std::array<std::array<std::uint8_t,4>,16> emission{};
    std::array<std::array<std::uint8_t,16>,16> refresh{};
    Oracle(){
        for(unsigned a=0;a!=256;++a)for(unsigned b=0;b!=256;++b){
            mul[a][b]=field::byteMultiply(std::uint8_t(a),std::uint8_t(b));
            if(mul[a][b]==1)inverse[a]=std::uint8_t(b);
        }
        for(unsigned c=0;c!=64;++c)for(unsigned d=0;d!=32;++d)outer[c][d]=basis(d,c);
        for(unsigned i=0;i!=16;++i)for(unsigned j=0;j!=4;++j)emission[i][j]=basis(j,i);
        for(unsigned i=0;i!=16;++i)for(unsigned j=0;j!=16;++j){
            std::uint8_t numerator=1,denominator=1;
            for(unsigned k=0;k!=16;++k)if(k!=j){
                numerator=mul[numerator][(16+i)^k];denominator=mul[denominator][j^k];
            }
            refresh[i][j]=mul[numerator][inverse[denominator]];
        }
    }
    std::uint8_t basis(unsigned degree,unsigned point) const{
        std::uint8_t result=1;
        for(unsigned bit=0;bit!=7;++bit)if(degree>>bit&1){
            std::uint8_t value=1,normalizer=1;
            for(unsigned j=0;j!=(1U<<bit);++j){value=mul[value][point^j];normalizer=mul[normalizer][(1U<<bit)^j];}
            result=mul[result][mul[value][inverse[normalizer]]];
        }
        return result;
    }
    Element byteScale(Element x,std::uint8_t a) const{
        std::array<std::uint8_t,16> result{};
        for(unsigned j=0;j!=16;++j)result[j]=mul[a][field::coefficient(x,j)];
        return field::fromCoefficients(result);
    }
    Element product(Element a,Element b) const{
        std::array<std::uint8_t,31> coefficients{};
        for(unsigned i=0;i!=16;++i)for(unsigned j=0;j!=16;++j)
            coefficients[i+j]^=mul[field::coefficient(a,i)][field::coefficient(b,j)];
        // Generic monic polynomial division, not the two-fold SIMD reduction.
        for(unsigned degree=31;degree-- >16;){
            const auto leading=coefficients[degree];
            for(unsigned j=0;j!=17;++j)coefficients[degree-16+j]^=mul[leading][field::modulus[j]];
        }
        std::array<std::uint8_t,16> result{};
        for(unsigned j=0;j!=16;++j)result[j]=coefficients[j];
        return field::fromCoefficients(result);
    }
    Element determinant(const std::array<Element,7>& t) const{
        // Characteristic two removes all permutation signs. This is
        // deliberately independent of setup's fraction-free elimination.
        std::array<unsigned,4> permutation{0,1,2,3};Element result{};
        do{
            Element term=field::embedByte(1);
            for(unsigned row=0;row!=4;++row)term=product(term,t[3+row-permutation[row]]);
            result^=term;
        }while(std::next_permutation(permutation.begin(),permutation.end()));
        return result;
    }
    void reference(const enc::Plan& p,const Block* input,Block* output) const{
        std::vector<Element> planes(p.n);
        for(std::size_t group=0;group!=p.rows;++group)for(unsigned c=0;c!=64;++c){
            std::array<Element,4> column{};
            // Four interleaved independent groups; each group owns four
            // 32-coefficient RS inputs, not 128 consecutive input records.
            for(unsigned r=0;r!=4;++r)for(unsigned d=0;d!=32;++d)
                column[r]^=byteScale(element(input[(group/4)*512+4*(r*32+d)+group%4]),outer[c][d]);
            const auto& t=p.toeplitz[c*p.rows+group];
            for(unsigned r=0;r!=4;++r){
                Element value{};
                for(unsigned a=0;a!=4;++a)value^=product(column[a],t[3+r-a]);
                planes[(4*c+r)*p.rows+group]=value;
            }
        }
        // One persistent natural-order state across every region. The dense
        // oracle consumes C*x, whereas the optimized C*y path relies on CA=0.
        std::array<Element,16> state{};
        for(std::size_t first=0;first!=p.n;first+=64){
            const auto region=first/p.rows,step=first/64;
            std::array<Element,16> feedback{},next{};
            for(unsigned i=0;i!=64;++i){
                const auto x=planes[region*p.rows+p.route[first+i]];
                auto y=x;
                for(unsigned j=0;j!=4;++j){
                    y^=byteScale(state[4*(i/16)+j],emission[i%16][j]);
                    feedback[4*(i/16)+j]^=byteScale(x,emission[i%16][j]);
                }
                output[first+4*(i%16)+i/16]=block(y);
            }
            for(unsigned j=0;j!=16;++j)feedback[j]=product(feedback[j]^state[j],p.before[step][j]);
            for(unsigned i=0;i!=16;++i){
                for(unsigned j=0;j!=16;++j)next[i]^=byteScale(feedback[j],refresh[i][j]);
                next[i]=product(next[i],p.after[step][i]);
            }
            state=next;
        }
    }
};

struct NoMixer {};
struct ByteControl {
    std::vector<std::uint8_t> keys;
    std::vector<sr::ByteDiagonalRefresh> refresh;
    explicit ByteControl(const enc::Plan& p):keys(p.rows*64*9),refresh(p.n/64){
        require(p.restrictedToGf256,"byte control requires GF256-restricted setup");
        for(std::size_t g=0;g!=p.rows;g+=4)for(unsigned c=0;c!=64;++c)for(unsigned lane=0;lane!=4;++lane){
            const auto nine=enc::Plan::nineKeys(p.toeplitz[c*p.rows+g+lane]);
            for(unsigned j=0;j!=9;++j){
                require(!nine[j].hi&&nine[j].lo<256,"byte control outer coefficient range");
                keys[g*576+4*(c*9+j)+lane]=std::uint8_t(nine[j].lo);
            }
        }
        for(std::size_t step=0;step!=refresh.size();++step){
            std::array<std::uint16_t,16> a,b;
            for(unsigned j=0;j!=16;++j){
                const auto x=p.refresh[step].before[j],y=p.refresh[step].after[j];
                require(!x.hi&&!y.hi&&x.lo<256&&y.lo<256,"byte control refresh coefficient range");
                a[j]=std::uint16_t(x.lo);b[j]=std::uint16_t(y.lo);
            }
            refresh[step].prepare(a.data(),b.data());
        }
    }
    std::size_t hotKeyBytes() const noexcept{return keys.size()+refresh.size()*sizeof(refresh[0]);}
    void outer(const enc::Plan& p,const Block* input,Block* scratch) const{
        old::forward<4,true,64>(p.k,p.stride,input,scratch,keys.data(),nullptr);
        for(unsigned c=0;c!=256;++c)scratch[c*p.stride+p.rows]=scratch[c*p.stride];
    }
    template<bool Stream>
    HD35_NOINLINE void inner(const enc::Plan& p,const Block* planes,Block* output) const{
        // Match encoder::inner exactly: one persistent striped state, the
        // same next-region T1 prefetches, and the same selected store policy.
        sr::State state{};
        for(unsigned region=0;region!=256;++region){
            const auto begin=region*p.rows;const auto* source=planes+region*p.stride;
            for(std::size_t j=0;j!=p.rows;j+=64){
                if(region+1<256){
                    const auto* next=reinterpret_cast<const char*>(planes+(region+1)*p.stride+j);
                    for(unsigned h=0;h!=16;++h)_mm_prefetch(next+64*h,_MM_HINT_T1);
                }
                old::innerPacket4<Stream>(source,p.packetRoute.data()+(begin+j)/2,
                    output+begin+j,state,NoMixer{},refresh[(begin+j)/64]);
            }
        }
        if constexpr(Stream)_mm_sfence();
    }
    template<bool Stream=true>
    void forward(const enc::Plan& p,const Block* input,Block* output,Block* scratch) const{
        outer(p,input,scratch);
        if constexpr(Stream){
            if(reinterpret_cast<std::uintptr_t>(output)&63U)inner<false>(p,scratch,output);
            else inner<true>(p,scratch,output);
        }else inner<false>(p,scratch,output);
    }
};

static void checkScalarSetup(const Oracle& oracle){
    std::mt19937_64 random(0xa526c671a2144cf5ULL);
    for(unsigned i=0;i!=64;++i){
        const Element a{random(),random()},b{random(),random()};
        require(oracle.product(a,b)==field::multiplyOracle(a,b),"table scalar oracle disagrees with polynomial arithmetic");
        require(enc::setupMultiply(a,b)==oracle.product(a,b),"fast setup arithmetic disagrees with oracle");
    }
    for(unsigned i=0;i!=32;++i){
        std::array<Element,7> t{};
        if(i==1)t[3]=field::embedByte(1);
        else if(i>1)for(auto& a:t)a=i<10?field::embedByte(std::uint8_t(random()%3)):Element{random(),random()};
        require(enc::Plan::nonsingular(t)==enc::nonzero(oracle.determinant(t)),"Toeplitz rejection test disagrees with Leibniz determinant");
    }
    // Explicitly check the self-orthogonal RS emission used by C=A^T.
    for(unsigned i=0;i!=4;++i)for(unsigned j=0;j!=4;++j){
        std::uint8_t dot=0;
        for(unsigned row=0;row!=16;++row)dot^=oracle.mul[oracle.emission[row][i]][oracle.emission[row][j]];
        require(dot==0,"fixed C*A is not zero");
    }
    std::cout<<"PASS scalar arithmetic, setup rejection, dense C*A=0\n";
}

static void checkPlan(const enc::Plan& p){
    require(p.outerKeys.size()==p.rows*64*7,"seven field coefficients per Toeplitz matrix");
    unsigned swapPatterns=0;
    std::vector<unsigned> seen(p.rows);
    for(unsigned region=0;region!=256;++region){
        std::fill(seen.begin(),seen.end(),0);
        for(std::size_t j=0;j!=p.rows;++j){
            const auto r=p.route[region*p.rows+j];require(r<p.rows,"logical route range");++seen[r];
        }
        for(auto n:seen)require(n==1,"logical route must be a permutation");
        for(std::size_t j=0;j!=p.rows/2;j+=2){
            const auto a=p.packetRoute[region*p.rows/2+j],b=p.packetRoute[region*p.rows/2+j+1];
            require((a&1)==(region&1)&&(b&1)==(region&1),"alternating pair parity");
            swapPatterns|=1U<<((a>>15)|((b>>15)<<1));
        }
    }
    require(swapPatterns==15,"all four independent pair swap masks must occur");
    for(std::size_t step=0;step!=p.refresh.size();++step)for(unsigned j=0;j!=16;++j){
        require(enc::nonzero(p.before[step][j])&&enc::nonzero(p.after[step][j]),"refresh diagonal is nonzero");
        require(p.refresh[step].before[j]==p.before[step][4*(j%4)+j/4],"striped before diagonal order");
        require(p.refresh[step].after[j]==p.after[step][4*(j%4)+j/4],"striped after diagonal order");
    }
}

static void checkToeplitz(const Oracle& oracle){
    std::mt19937_64 random(0xc217f02a7d610613ULL);
    constexpr unsigned selected[7]={0,1,2,3,5,6,7};
    for(unsigned sample=0;sample!=16;++sample){
        std::array<std::array<Element,7>,4> t;
        std::array<Element,28> keys;
        alignas(64) Element input[4][4],output[4][4];
        alignas(64) __m512i values[4];
        for(unsigned lane=0;lane!=4;++lane){
            for(auto& a:t[lane])a=sample==0?Element{}:Element{random(),random()};
            const auto nine=enc::Plan::nineKeys(t[lane]);
            for(unsigned j=0;j!=7;++j)keys[4*j+lane]=nine[selected[j]];
            for(unsigned row=0;row!=4;++row)input[row][lane]={random(),random()};
        }
        for(unsigned row=0;row!=4;++row)values[row]=_mm512_load_si512(input[row]);
        enc::toeplitz4(values,keys.data());
        for(unsigned row=0;row!=4;++row)_mm512_store_si512(output[row],values[row]);
        for(unsigned row=0;row!=4;++row)for(unsigned lane=0;lane!=4;++lane){
            Element expected{};
            for(unsigned j=0;j!=4;++j)expected^=oracle.product(input[j][lane],t[lane][3+row-j]);
            require(output[row][lane]==expected,"compact Toeplitz circuit disagrees with dense full-field matrix");
        }
    }
    std::cout<<"PASS compact seven-key/nine-product full-field Toeplitz circuit\n";
}

template<bool Byte>
static void releaseOracleCopies(enc::Plan& p){
    decltype(p.toeplitz){}.swap(p.toeplitz);decltype(p.before){}.swap(p.before);
    decltype(p.after){}.swap(p.after);decltype(p.route){}.swap(p.route);
    if constexpr(Byte){decltype(p.outerKeys){}.swap(p.outerKeys);decltype(p.refresh){}.swap(p.refresh);}
}

static void checkOne(const Oracle& oracle,std::size_t k,std::uint64_t seed,bool restricted){
    enc::Plan p(k,seed,restricted);checkPlan(p);
    c2::BenchPayload input(p.k),actual(p.n+8),expected(p.n),scratch(p.scratchSize()+8);
    produce(input.p,p.k);
    oracle.reference(p,input.p,expected.p);
    constexpr Block guard{0x2e71d6935a90bf24ULL,0xc031674f27e18ac5ULL};
    for(unsigned offset=0;offset!=4;++offset){
        std::fill_n(actual.p,p.n+8,guard);std::fill_n(scratch.p,p.scratchSize()+8,guard);
        enc::forward(p,input.p,actual.p+offset,scratch.p+offset);
        same(actual.p+offset,expected.p,p.n,"dense whole-encoder oracle/alignment");
        for(unsigned i=0;i!=offset;++i){require(actual.p[i]==guard,"output leading guard");require(scratch.p[i]==guard,"scratch leading guard");}
        for(std::size_t i=p.n+offset;i!=p.n+8;++i)require(actual.p[i]==guard,"output trailing guard");
        for(std::size_t i=p.scratchSize()+offset;i!=p.scratchSize()+8;++i)require(scratch.p[i]==guard,"scratch trailing guard");
    }
    enc::forward<false>(p,input.p,actual.p,scratch.p);
    same(actual.p,expected.p,p.n,"cached output");
    std::copy_n(input.p,p.k,actual.p);enc::forward(p,actual.p,actual.p,scratch.p);
    same(actual.p,expected.p,p.n,"exact input/output alias");
    std::copy_n(input.p,p.k,actual.p+3);enc::forward(p,actual.p+3,actual.p,scratch.p);
    same(actual.p,expected.p,p.n,"partial input/output alias and unaligned input");
    std::copy_n(input.p,p.k,actual.p);enc::forward(p,actual.p,actual.p+1,scratch.p);
    same(actual.p+1,expected.p,p.n,"forward partial input/output alias");
    if(restricted){
        ByteControl control(p);
        for(unsigned offset=0;offset!=4;++offset){
            control.forward<true>(p,input.p,actual.p+offset,scratch.p+offset);
            same(actual.p+offset,expected.p,p.n,"retained GF256 kernel NT/unaligned restriction");
        }
        control.forward<false>(p,input.p,actual.p,scratch.p);
        same(actual.p,expected.p,p.n,"retained GF256 kernel cached restriction");
        // Exercise the exact split entry points selected by component timing.
        control.outer(p,input.p,scratch.p);control.inner<true>(p,scratch.p,actual.p);
        same(actual.p,expected.p,p.n,"retained GF256 split outer/inner restriction");
        releaseOracleCopies<true>(p);control.forward<true>(p,input.p,actual.p,scratch.p);
        same(actual.p,expected.p,p.n,"retained GF256 hot-only setup");
    }else{
        require(std::any_of(p.toeplitz.begin(),p.toeplitz.end(),[](const auto& t){return std::any_of(t.begin(),t.end(),[](Element a){return a.hi||(a.lo>>8);});}),"full-field outer samples unexpectedly byte-only");
        releaseOracleCopies<false>(p);enc::forward(p,input.p,actual.p,scratch.p);
        same(actual.p,expected.p,p.n,"literal GF128 hot-only setup");
    }
    std::cout<<"PASS whole encoder K="<<p.k<<" seed="<<seed<<" coefficients="<<(restricted?"GF256-restricted":"GF128")
             <<" oracle/alignment/guards/alias checksum="<<std::hex<<hash(expected.p,p.n)<<std::dec<<"\n";
}

static void checkLinearity(const Oracle& oracle){
    const enc::Plan p(8192,20261009,false);
    c2::BenchPayload u(p.k),v(p.k),mixed(p.k),eu(p.n),ev(p.n),actual(p.n),scratch(p.scratchSize());
    produce(u.p,p.k,591);produce(v.p,p.k,8173);
    const Element lambda{0x0305070b0d111317ULL,0x1d1f25292b2f353bULL};
    for(std::size_t j=0;j!=p.k;++j)mixed.p[j]=block(oracle.product(lambda,element(u.p[j]))^element(v.p[j]));
    enc::forward(p,u.p,eu.p,scratch.p);enc::forward(p,v.p,ev.p,scratch.p);
    enc::forward(p,mixed.p,actual.p,scratch.p);
    for(std::size_t j=0;j!=p.n;++j)
        require(element(actual.p[j])==(oracle.product(lambda,element(eu.p[j]))^element(ev.p[j])),"whole encoder is not full-field linear");
    std::fill_n(u.p,p.k,Block{});enc::forward(p,u.p,actual.p,scratch.p);
    for(std::size_t j=0;j!=p.n;++j)require(actual.p[j]==Block{},"zero input must give zero output");
    std::cout<<"PASS E(lambda*u+v)=lambda*E(u)+E(v) for non-subfield lambda, zero input\n";
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
static void checkSummary(){
    const auto odd=summarize({5,1,3}),even=summarize({9,1,3,5});
    require(odd.median==3&&odd.mean==3&&odd.best==1&&odd.worst==5,"odd sample summary");
    require(even.median==4&&even.mean==4.5&&even.best==1&&even.worst==9,"even sample summary");
}

struct NoByteControl {explicit NoByteControl(const enc::Plan&) {}};
template<bool Byte,bool Stream,Mode Which,class Control>
static inline void runOne(const enc::Plan& p,const Control& control,const Block* input,Block* output,Block* scratch){
    // Engine, component, and store policy dispatch are all outside timing.
    // No function pointer, virtual dispatch, or heap callback in this path.
    if constexpr(Byte){
        if constexpr(Which==Mode::Full)control.template forward<Stream>(p,input,output,scratch);
        else if constexpr(Which==Mode::Outer)control.outer(p,input,scratch);
        else control.template inner<Stream>(p,scratch,output);
    }else{
        if constexpr(Which==Mode::Full)enc::forward<Stream>(p,input,output,scratch);
        else if constexpr(Which==Mode::Outer)enc::outer<true>(p,input,scratch);
        else enc::inner<Stream>(p,scratch,output);
    }
}
template<bool Byte,bool Stream,Mode Which>
static void benchmark(std::size_t k,unsigned calls){
    const auto setupBegin=Clock::now();enc::Plan p(k,113,Byte);
    const std::conditional_t<Byte,ByteControl,NoByteControl> control(p);
    const auto setupPeakBytes=p.setupBytes()+[&]{if constexpr(Byte)return control.hotKeyBytes();else return std::size_t(0);}();
    const auto hotCoefficientBytes=[&]{if constexpr(Byte)return control.hotKeyBytes();else return p.hotKeyBytes();}();
    // Dense-oracle copies are useful to --check but never read by a hot
    // kernel. Release them before allocation/warmup, in both timing modes.
    releaseOracleCopies<Byte>(p);
    const auto retainedSetupBytes=p.setupBytes()+[&]{if constexpr(Byte)return control.hotKeyBytes();else return std::size_t(0);}();
    c2::BenchPayload input(p.k),output(p.n),scratch(p.scratchSize());produce(input.p,p.k);
    std::vector<double> samples(calls);
    // Populate scratch outside timing for the standalone inner benchmark;
    // identical preparation is used for all engine/component selections.
    if constexpr(Byte)control.outer(p,input.p,scratch.p);
    else enc::outer<true>(p,input.p,scratch.p);
    const auto setupEnd=Clock::now();
    // One explicitly requested, serial workload. No background jobs, threads,
    // throughput groups, or recursive benchmark invocations are launched.
    constexpr unsigned warmups=3;
    for(unsigned i=0;i!=warmups;++i)runOne<Byte,Stream,Which>(p,control,input.p,output.p,scratch.p);
    for(unsigned i=0;i!=calls;++i){
        std::atomic_signal_fence(std::memory_order_seq_cst);
        const auto begin=Clock::now();
        runOne<Byte,Stream,Which>(p,control,input.p,output.p,scratch.p);
        std::atomic_signal_fence(std::memory_order_seq_cst);
        const auto end=Clock::now();
        samples[i]=std::chrono::duration<double,std::milli>(end-begin).count();
    }
    // No payload hashing, regeneration, validation, or statistics touches
    // the cache between timed calls. Hash the final result only after batch.
    const auto checksum=Which==Mode::Outer?hash(scratch.p,p.scratchSize()):hash(output.p,p.n);
    const auto result=summarize(samples);
    const char* name=Which==Mode::Full?"full":Which==Mode::Outer?"outer":"inner";
    std::cout<<"engine="<<(Byte?"retained_GF256":"literal_GF128")<<" mode="<<name
             <<" logK="<<std::countr_zero(k)<<" samples="<<calls<<" warmups="<<warmups
             <<" scratch=NT output="<<(Which==Mode::Outer?"unused":Stream?"NT":"cached")
             <<" prefetch="<<(Which==Mode::Outer?"unused":"T1-next-region-16lines")
             <<" median_ms="<<std::fixed<<std::setprecision(6)<<result.median<<" mean_ms="<<result.mean
             <<" best_ms="<<result.best<<" worst_ms="<<result.worst
             <<" setup_ms="<<std::chrono::duration<double,std::milli>(setupEnd-setupBegin).count()
             <<" hot_coefficient_bytes="<<hotCoefficientBytes<<" hot_route_bytes="<<p.packetRoute.size()*sizeof(std::uint16_t)
             <<" retained_setup_bytes="<<retainedSetupBytes<<" peak_setup_payload_bytes="<<setupPeakBytes
             <<" rejected_toeplitz="<<p.toeplitzRejected<<" final_checksum="<<std::hex<<checksum<<std::dec<<"\n";
}
template<bool Byte,bool Stream>
static void benchmarkMode(std::size_t k,unsigned calls,Mode mode){
    if(mode==Mode::Full)benchmark<Byte,Stream,Mode::Full>(k,calls);
    else if(mode==Mode::Outer)benchmark<Byte,Stream,Mode::Outer>(k,calls);
    else benchmark<Byte,Stream,Mode::Inner>(k,calls);
}
template<bool Byte>
static void benchmarkStore(std::size_t k,unsigned calls,Mode mode,bool stream){
    if(stream)benchmarkMode<Byte,true>(k,calls,mode);
    else benchmarkMode<Byte,false>(k,calls,mode);
}

int main(int argc,char** argv){
    try{
        // Compile this isolated research driver for AVX512F/BW/DQ/VL, GFNI,
        // and VPCLMULQDQ (for example GCC -O3 -std=c++20 -march=znver4).
        // A globally ISA-targeted binary is not a portable CPU dispatcher.
        if(argc==1||(argc==2&&std::string(argv[1])=="--check")){
            const Oracle oracle;checkSummary();checkScalarSetup(oracle);checkToeplitz(oracle);
            // K16384 crosses outer tiles and contains two inner steps in
            // every region; K8192 alone would not exercise those boundaries.
            checkOne(oracle,8192,113,false);checkOne(oracle,16384,7919,false);
            checkOne(oracle,8192,113,true);checkLinearity(oracle);
            std::cout<<"PASS all literal GF128 encoder checks\n";return 0;
        }
        if(argc>=5&&argc<=7&&std::string(argv[1])=="--bench"){
            const auto calls=std::stoul(argv[2]),logK=std::stoul(argv[3]);
            require(calls>=1&&calls<=1000000,"calls must be in 1..1000000");
            require(logK>=13&&logK<=22,"logK must be in 13..22");
            const std::string name=argv[4];
            require(name=="full"||name=="outer"||name=="inner","mode must be full, outer, or inner");
            const auto mode=name=="full"?Mode::Full:name=="outer"?Mode::Outer:Mode::Inner;
            const std::string store=argc>=6?argv[5]:"nt";require(store=="nt"||store=="cached","store must be nt or cached");
            const std::string engine=argc==7?argv[6]:"gf128";require(engine=="gf128"||engine=="gf256","engine must be gf128 or gf256");
            std::cout<<"Research MT19937 setup; no ideal-uniform or cryptographic-seed guarantee. Run only one benchmark process at a time.\n";
            if(engine=="gf256")benchmarkStore<true>(std::size_t(1)<<logK,unsigned(calls),mode,store=="nt");
            else benchmarkStore<false>(std::size_t(1)<<logK,unsigned(calls),mode,store=="nt");
            return 0;
        }
        std::cerr<<"usage: literal_encoder [--check]\n       literal_encoder --bench SAMPLES LOGK full|outer|inner [nt|cached] [gf128|gf256]\n";
        return 2;
    }catch(const std::exception& error){std::cerr<<"FAIL: "<<error.what()<<"\n";return 1;}
}
