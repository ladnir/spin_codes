#include "HdWideMdsInner.h"
#include "HdSumrankOuter.h"
#include "BenchPayload.h"
#include "Punctured8Oracle.h"
#include <atomic>
#include <iomanip>
#include <iostream>

namespace c2=spin::research::char2;
namespace sr=c2::hdsumrank;
namespace wide=c2::hdwide;
using Block=c2::Block;
using Clock=std::chrono::steady_clock;
#if defined(_MSC_VER)
#define HDW_NOINLINE __declspec(noinline)
#else
#define HDW_NOINLINE __attribute__((noinline))
#endif
static void require(bool b,const char* why){if(!b)throw std::runtime_error(why);}
static std::uint16_t scalar(c2::wide16::Coefficient x){return std::uint16_t(x.a|(unsigned(x.b)<<8));}

template<unsigned S> struct Plan {
    std::size_t k,n,rows,stride;
    std::vector<std::uint16_t> route;
    std::vector<c2::wide16::Coefficient> diagonal;
    std::vector<wide::Setup<S>> steps;
    explicit Plan(std::size_t size):k(size),n(2*size),rows(size/128),stride(rows+4){
        if(size<16384||size%16384||rows>65536)throw std::invalid_argument("K must be a multiple of 16384, at most 2^23");
        std::mt19937_64 rr(113^0x91af823017ccb91aULL),mr(113^0x14b2ec8399500213ULL),dr(113^0xe685639417f293b5ULL),tr(113^0xb2d174839ae61713ULL);
        std::uniform_int_distribution<unsigned> nz(1,65535);
        diagonal.resize(n);for(auto& d:diagonal)d=c2::wide16::prepareCoefficient(std::uint16_t(nz(dr)));
        route.resize(n);
        for(unsigned c=0;c!=256;++c){auto* p=route.data()+c*rows;std::iota(p,p+rows,std::uint16_t(0));
            for(auto j=rows-1;j;--j){std::uniform_int_distribution<std::size_t> d(0,j);std::swap(p[j],p[d(rr)]);}}
        steps.resize(n/128);
        for(auto& step:steps){
            for(unsigned i=0;i!=S;++i){step.before[i]=c2::wide16::prepareCoefficient(std::uint16_t(nz(tr)));step.after[i]=c2::wide16::prepareCoefficient(std::uint16_t(nz(tr)));}
            std::uint16_t dot;
            do{dot=0;for(unsigned i=0;i!=S;++i){const auto v=std::uint16_t(mr());step.v[i]=c2::wide16::prepareCoefficient(v);dot^=sr::mul(v,scalar(step.after[i]));}}while(dot==1);
        }
    }
    std::size_t scratch()const{return 256*stride;}
};

template<unsigned S> HDW_NOINLINE static void outer(const Plan<S>& p,const Block* input,Block* planes){
    c2::hdsumrankouter::forward(p.k,p.stride,input,planes,p.diagonal.data());
}
template<unsigned S> HDW_NOINLINE static void inner(const Plan<S>& p,const Block* planes,Block* output){
    wide::State<S> state;
    for(unsigned region=0;region!=256;++region){
        const auto begin=region*p.rows;const auto* source=planes+region*p.stride;
        for(std::size_t j=0;j!=p.rows;j+=128){
            if(region+1<256){
                const auto* next=reinterpret_cast<const char*>(planes+(region+1)*p.stride+j);
                for(unsigned h=0;h!=32;++h)_mm_prefetch(next+64*h,_MM_HINT_T1);
            }
            wide::packet<S,true>(source,p.route.data()+begin+j,output+begin+j,state,p.steps[(begin+j)/128]);
        }
    }
}
template<unsigned S> HDW_NOINLINE static void forward(const Plan<S>& p,const Block* input,Block* output,Block* scratch){
    outer(p,input,scratch);inner(p,scratch,output);
}

static std::uint16_t oracleMul(std::uint16_t a,std::uint16_t b){
    const auto& f=c2::p8Oracle().mul;const auto h=f[a>>8][b>>8];
    return std::uint16_t((f[a&255][b&255]^f[0x20][h])|((f[a&255][b>>8]^f[a>>8][b&255]^h)<<8));
}
static Block scaled(Block x,std::uint16_t a){
    std::array<std::uint16_t,8> v;std::memcpy(v.data(),&x,16);
    for(auto& c:v)c=oracleMul(c,a);std::memcpy(&x,v.data(),16);return x;
}
static void add(Block& a,Block b){a.lo^=b.lo;a.hi^=b.hi;}
template<unsigned S> struct Oracle {
    sr::Matrix a,h;
    Oracle():a(128,std::vector<std::uint16_t>(S)){
        for(unsigned i=0;i!=128;++i)for(unsigned j=0;j!=S;++j){
            a[i][j]=1;for(unsigned b=0;b!=5;++b)if(j>>b&1)a[i][j]=sr::mul(a[i][j],c2::hdbw::detail::normalized(b,i));
        }
        sr::Matrix x(32,std::vector<std::uint16_t>(32)),y=x;
        for(unsigned i=0;i!=32;++i)for(unsigned j=0;j!=32;++j){x[i][j]=sr::power(i,j);y[i][j]=sr::power(i+32,j);}
        h=sr::product(y,sr::inverse(x));h.resize(S);for(auto& row:h)row.resize(S);
    }
    void packet(const Block* x,Block* y,std::array<Block,S>& q,const wide::Setup<S>& setup)const{
        Block dot{};std::array<Block,S> feedback{},next=q;
        for(unsigned i=0;i!=S;++i){
            add(dot,scaled(q[i],scalar(setup.v[i])));
            for(unsigned j=0;j!=128;++j)add(feedback[i],scaled(x[j],a[j][i]));
            feedback[i]=scaled(feedback[i],scalar(setup.before[i]));
        }
        for(unsigned i=0;i!=128;++i){y[i]=x[i];for(unsigned j=0;j!=S;++j)add(y[i],scaled(q[j],a[i][j]));}
        for(unsigned i=0;i!=S;++i){Block value=dot;for(unsigned j=0;j!=S;++j)add(value,scaled(feedback[j],h[i][j]));add(next[i],scaled(value,scalar(setup.after[i])));}
        q=next;
    }
};
template<unsigned S> static void check(){
    Plan<S> p(16384);Oracle<S> oracle;std::mt19937_64 rng(41731);
    const Block guard{0x837263821abcULL,0x15ecda785abULL};
    for(unsigned trial=0;trial!=2*(128+S)+16;++trial){
        std::array<Block,128> x{},want{};std::array<Block,S> q{},expected{},actual{};
        if(trial<2*(128+S)){
            const unsigned i=trial%(128+S),component=trial/(128+S);
            const Block one{0x0001000100010001ULL<<(8*component),0x0001000100010001ULL<<(8*component)};
            if(i<128)x[i]=one;else q[i-128]=one;
        }else{for(auto& a:x)a={rng(),rng()};for(auto& a:q)a={rng(),rng()};}
        const auto& setup=p.steps[trial%p.steps.size()];expected=q;oracle.packet(x.data(),want.data(),expected,setup);
        wide::State<S> state;std::memcpy(&state,q.data(),sizeof(q));
        std::array<Block,128> result{};
        wide::packet<S,false,std::uint16_t>(x.data(),nullptr,result.data(),state,setup);std::memcpy(actual.data(),&state,sizeof(actual));
        require(result==want&&actual==expected,"wide packet basis/random dense oracle");
        std::memcpy(&state,q.data(),sizeof(q));result=x;
        wide::packet<S,false,std::uint16_t>(result.data(),nullptr,result.data(),state,setup);
        require(result==want,"wide packet exact alias");
        if(trial>=2*(128+S)){
            std::array<std::uint32_t,128> route;std::iota(route.begin(),route.end(),0);std::shuffle(route.begin(),route.end(),rng);
            for(unsigned offset=0;offset!=4;++offset){
                std::array<Block,136> source,out;source.fill(guard);out.fill(guard);
                for(unsigned i=0;i!=128;++i)source[offset+route[i]]=x[i];
                std::memcpy(&state,q.data(),sizeof(q));
                wide::packet<S,true>(source.data()+offset,route.data(),out.data()+offset,state,setup);std::memcpy(actual.data(),&state,sizeof(actual));
                require(std::equal(want.begin(),want.end(),out.begin()+offset)&&actual==expected,"wide routed unaligned oracle");
                for(unsigned i=0;i!=offset;++i)require(out[i]==guard,"prefix guard");
                for(unsigned i=offset+128;i!=136;++i)require(out[i]==guard,"suffix guard");
            }
        }
    }
    for(std::size_t k:{16384,49152}){
        Plan<S> full(k);std::vector<Block> x(k),out(full.n),want(full.n),scratch(full.scratch());for(auto& a:x)a={rng(),rng()};
        c2::hdbwouter::outerForward256(k,full.stride,x.data(),scratch.data());
        std::array<Block,128> packet{};std::array<Block,S> q{};
        for(std::size_t i=0;i!=full.n;i+=128){
            const auto region=i/full.rows;
            for(unsigned j=0;j!=128;++j){const auto row=full.route[i+j];packet[j]=scaled(scratch[region*full.stride+row],scalar(full.diagonal[region*full.rows+row]));}
            oracle.packet(packet.data(),want.data()+i,q,full.steps[i/128]);
        }
        forward(full,x.data(),out.data(),scratch.data());require(out==want,"wide full dense oracle");
        std::copy(x.begin(),x.end(),out.begin());forward(full,out.data(),out.data(),scratch.data());require(out==want,"wide full alias");
    }
    std::cout<<"check,wide-MDS,state,"<<S<<",all-bases-random-routed-full-alias,pass\n";
}
static void produce(Block* p,std::size_t n){std::uint64_t x=913,y=113;
    for(std::size_t i=0;i!=n;++i){x^=x<<13;x^=x>>7;x^=x<<17;y+=0x9e3779b97f4a7c15ULL;_mm_stream_si128(reinterpret_cast<__m128i*>(p+i),_mm_set_epi64x(y,x));}_mm_sfence();}
static std::uint64_t hash(const Block* p,std::size_t n){std::uint64_t h=0;for(std::size_t i=0;i!=n;++i){h=(h^p[i].lo)*0x100000001b3ULL;h=(h^p[i].hi)*0x100000001b3ULL;}return h;}
template<unsigned S> static void bench(unsigned calls,const std::string& mode,unsigned logK){
    Plan<S> p(std::size_t(1)<<logK);c2::BenchPayload in(p.k),out(p.n),scratch(p.scratch());produce(in.p,p.k);outer(p,in.p,scratch.p);
    auto run=[&]{if(mode=="inner")inner(p,scratch.p,out.p);else if(mode=="outer")outer(p,in.p,scratch.p);else forward(p,in.p,out.p,scratch.p);};
    for(unsigned j=0;j!=3;++j)run();std::vector<double> times(calls);std::uint64_t consume=0;const Block* result=mode=="outer"?scratch.p:out.p;
    for(unsigned j=0;j!=calls;++j){produce(in.p,p.k);const auto a=Clock::now();run();times[j]=std::chrono::duration<double,std::micro>(Clock::now()-a).count();std::atomic_signal_fence(std::memory_order_seq_cst);consume^=result[(j*104729U)%p.n].lo;}
    std::sort(times.begin(),times.end());
    std::cout<<"profile,wide-MDS,K,"<<p.k<<",state,"<<S<<",packet_width,128,shared_diagonal,1,calls,"<<calls<<",checksum,"<<hash(result,p.n)<<",consume,"<<consume<<"\nname,median_us,p10_us,p90_us\n"<<std::fixed<<std::setprecision(3)<<mode<<','<<times[calls/2]<<','<<times[calls/10]<<','<<times[9*calls/10]<<'\n';
}
int main(int argc,char** argv){try{
    if(!c2::gfni512Available())return 77;
    if(argc==2&&std::string(argv[1])=="--check"){check<20>();check<24>();check<32>();return 0;}
    if(argc<2||std::string(argv[1])!="--bench")throw std::invalid_argument("--check | --bench [calls] [full|inner|outer] [state20|24|32] [logK]");
    unsigned calls=argc>2?std::stoul(argv[2]):31,state=argc>4?std::stoul(argv[4]):20,logK=argc>5?std::stoul(argv[5]):20;std::string mode=argc>3?argv[3]:"full";
    if(argc>6||!calls||calls>100000||logK<14||logK>23||(state!=20&&state!=24&&state!=32)||(mode!="full"&&mode!="inner"&&mode!="outer"))throw std::invalid_argument("arguments");
    if(state==20)bench<20>(calls,mode,logK);else if(state==24)bench<24>(calls,mode,logK);else bench<32>(calls,mode,logK);return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
