#include "Hd35Outer.h"
#include "HdSumrankInner.h"
#include "BenchPayload.h"
#include "Punctured8Oracle.h"
#include <atomic>
#include <iomanip>
#include <iostream>

namespace c2=spin::research::char2;
namespace sr=c2::hdsumrank;
using Block=c2::Block;
using Clock=std::chrono::steady_clock;
static void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}

// Isolated 35% research candidate; no production/default profile changes.
// MT is reproducible benchmark setup, not the ideal independent randomness
// in the probability theorem. Setup and all allocations precede encoding.
template<unsigned Outer=128>
struct Plan {
    static_assert(Outer==128||Outer==256);
    std::size_t k,n,rows,stride;
    std::vector<std::uint16_t> route,diagonal;
    std::vector<c2::wide16::Coefficient> prepared;
    std::vector<sr::DiagonalRefresh> refresh;
    std::vector<std::array<std::uint16_t,16>> before,after;
    explicit Plan(std::size_t size,std::uint64_t seed):k(size),n(2*size),rows(size/(Outer/2)),stride(rows+4){
        if(size<4096||size%4096||rows%64||rows>65536)throw std::invalid_argument("common-refresh research dimensions");
        std::mt19937_64 rr(seed^0x91af823017ccb91aULL),dr(seed^0xe685639417f293b5ULL),tr(seed^0xb2d174839ae61713ULL);
        std::uniform_int_distribution<unsigned> nz(1,65535);
        diagonal.resize(n);prepared.resize(n);route.resize(n);
        for(auto& d:diagonal)d=std::uint16_t(nz(dr));
        c2::wide16::prepare(diagonal.data(),prepared.data(),n);
        for(unsigned c=0;c!=Outer;++c){
            auto* p=route.data()+c*rows;std::iota(p,p+rows,std::uint16_t(0));
            for(auto j=rows-1;j;--j){std::uniform_int_distribution<std::size_t> d(0,j);std::swap(p[j],p[d(rr)]);}
        }
        refresh.resize(n/64);before.resize(n/64);after.resize(n/64);
        for(std::size_t i=0;i!=n/64;++i){
            for(unsigned j=0;j!=16;++j){before[i][j]=std::uint16_t(nz(tr));after[i][j]=std::uint16_t(nz(tr));}
            refresh[i].prepare(before[i].data(),after[i].data());
        }
    }
    std::size_t scratch()const{return Outer*stride;}
};
// The common-refresh instantiation never reads a mixer. An empty type makes
// that compile-time property explicit without allocating unused coefficients.
struct NoMixer {};

template<bool Stream,unsigned Rank=4,unsigned Outer>
HD35_NOINLINE static void inner(const Plan<Outer>& p,const Block* planes,Block* output){
    sr::State state{};
    for(unsigned region=0;region!=Outer;++region){
        const auto begin=region*p.rows;
        const auto* source=planes+region*p.stride;
        for(std::size_t j=0;j!=p.rows;j+=64){
            if(region+1<Outer){
                const auto* next=reinterpret_cast<const char*>(planes+(region+1)*p.stride+j);
                for(unsigned h=0;h!=16;++h)_mm_prefetch(next+64*h,_MM_HINT_T1);
            }
            sr::mdsStriped<true,std::uint16_t,true,Rank,Rank,2,false,Stream>(source,
                p.route.data()+begin+j,output+begin+j,state,NoMixer{},p.refresh[(begin+j)/64]);
        }
    }
    if constexpr(Stream)_mm_sfence();
}

template<unsigned Outer>
HD35_NOINLINE static void outer(const Plan<Outer>& p,const Block* input,Block* scratch){
    if constexpr(Outer==128)c2::hd35::outer128(p.k,p.stride,input,scratch,p.prepared.data());
    else c2::hdsumrankouter::forward(p.k,p.stride,input,scratch,p.prepared.data());
}
template<bool Stream,unsigned Rank=4,unsigned Outer>
HD35_NOINLINE static void forward(const Plan<Outer>& p,const Block* input,Block* output,Block* scratch){
    outer(p,input,scratch);
    if constexpr(Stream){
        if(reinterpret_cast<std::uintptr_t>(output)&63U)inner<false,Rank>(p,scratch,output);
        else inner<true,Rank>(p,scratch,output);
    }else inner<false,Rank>(p,scratch,output);
}

static std::uint16_t oracleMul(std::uint16_t a,std::uint16_t b){
    const auto& f=c2::p8Oracle().mul;
    const auto h=f[a>>8][b>>8];
    return std::uint16_t((f[a&255][b&255]^f[0x20][h])|((f[a&255][b>>8]^f[a>>8][b&255]^h)<<8));
}
static Block scaled(Block x,std::uint16_t a){
    std::array<std::uint16_t,8> v;std::memcpy(v.data(),&x,16);
    for(auto& c:v)c=oracleMul(c,a);std::memcpy(&x,v.data(),16);return x;
}
static void add(Block& a,Block b){a.lo^=b.lo;a.hi^=b.hi;}
static std::uint8_t byteInverse(std::uint8_t a){
    const auto& f=c2::p8Oracle().mul;
    for(unsigned b=1;b!=256;++b)if(f[a][b]==1)return std::uint8_t(b);
    throw std::runtime_error("zero inverse");
}
static std::uint8_t basis(unsigned degree,unsigned point){
    const auto& f=c2::p8Oracle().mul;std::uint8_t result=1;
    for(unsigned bit=0;bit!=7;++bit)if(degree>>bit&1){
        std::uint8_t value=1,normalizer=1;
        for(unsigned j=0;j!=(1U<<bit);++j){value=f[value][point^j];normalizer=f[normalizer][(1U<<bit)^j];}
        result=f[result][f[value][byteInverse(normalizer)]];
    }
    return result;
}

template<unsigned Rank=4,unsigned Outer>
static void reference(const Plan<Outer>& p,const Block* input,Block* output){
    constexpr unsigned dimension=Outer/2;
    std::array<std::array<std::uint8_t,dimension>,Outer> coefficients;
    for(unsigned c=0;c!=Outer;++c)for(unsigned d=0;d!=dimension;++d)coefficients[c][d]=basis(d,c);
    std::vector<Block> planes(p.n);
    for(unsigned c=0;c!=Outer;++c)for(std::size_t r=0;r!=p.rows;++r){
        Block v{};
        for(unsigned d=0;d!=dimension;++d)add(v,scaled(input[(r/4)*(4*dimension)+4*d+r%4],coefficients[c][d]));
        planes[c*p.rows+r]=scaled(v,p.diagonal[c*p.rows+r]);
    }
    // Independent dense interpolation: Lagrange weights from nodes0..15
    // to nodes16..31, with polynomial-byte arithmetic, not SIMD FFT calls.
    std::array<std::array<std::uint8_t,16>,16> h;
    const auto& f=c2::p8Oracle().mul;
    for(unsigned i=0;i!=16;++i)for(unsigned j=0;j!=16;++j){
        std::uint8_t a=1,b=1;
        for(unsigned k=0;k!=16;++k)if(k!=j){a=f[a][(16+i)^k];b=f[b][j^k];}
        h[i][j]=f[a][byteInverse(b)];
    }
    std::array<Block,16> state{};
    for(std::size_t first=0;first!=p.n;first+=64){
        std::array<Block,64> x;
        std::array<Block,16> feedback{},next{};
        const auto region=first/p.rows,step=first/64;
        for(unsigned i=0;i!=64;++i){
            x[i]=planes[region*p.rows+p.route[first+i]];
            Block y=x[i];
            for(unsigned j=0;j!=Rank;++j){
                const auto a=basis(j,i%(4*Rank));
                add(y,scaled(state[Rank*(i/(4*Rank))+j],a));
                add(feedback[Rank*(i/(4*Rank))+j],scaled(x[i],a));
            }
            output[first+i]=y;
        }
        for(unsigned j=0;j!=16;++j){add(feedback[j],state[j]);feedback[j]=scaled(feedback[j],p.before[step][j]);}
        for(unsigned i=0;i!=16;++i){
            for(unsigned j=0;j!=16;++j)add(next[i],scaled(feedback[j],h[i][j]));
            next[i]=scaled(next[i],p.after[step][i]);
        }
        state=next;
    }
}
static void produce(Block* p,std::size_t n){
    std::uint64_t x=913,y=113;
    for(std::size_t i=0;i!=n;++i){x^=x<<13;x^=x>>7;x^=x<<17;y+=0x9e3779b97f4a7c15ULL;_mm_stream_si128(reinterpret_cast<__m128i*>(p+i),_mm_set_epi64x(y,x));}_mm_sfence();
}
static std::uint64_t hash(const Block* p,std::size_t n){std::uint64_t h=0;for(std::size_t i=0;i!=n;++i){h=(h^p[i].lo)*0x100000001b3ULL;h=(h^p[i].hi)*0x100000001b3ULL;}return h;}

template<unsigned Outer=128,unsigned Rank=4>
static void check(std::size_t k){
    Plan<Outer> p(k,113);c2::BenchPayload input(k+4),actual(p.n+8),expected(p.n),scratch(p.scratch());
    produce(input.p,k);reference<Rank>(p,input.p,expected.p);
    const Block guard{0xa1369,0x73495};
    for(unsigned offset=0;offset!=4;++offset){
        std::fill_n(actual.p,p.n+8,guard);
        forward<true,Rank>(p,input.p,actual.p+offset,scratch.p);
        require(std::equal(expected.p,expected.p+p.n,actual.p+offset),"full independent dense oracle");
        for(unsigned j=0;j!=offset;++j)require(actual.p[j]==guard,"prefix guard");
        for(auto j=p.n+offset;j!=p.n+8;++j)require(actual.p[j]==guard,"suffix guard");
    }
    forward<false,Rank>(p,input.p,actual.p,scratch.p);
    require(std::equal(expected.p,expected.p+p.n,actual.p),"cached output oracle");
    std::copy_n(input.p,k,actual.p);forward<true,Rank>(p,actual.p,actual.p,scratch.p);
    require(std::equal(expected.p,expected.p+p.n,actual.p),"full alias");
    // Four interleaved outer rows must not impose input alignment.
    std::move_backward(input.p,input.p+k,input.p+k+1);
    forward<true,Rank>(p,input.p+1,actual.p,scratch.p);
    require(std::equal(expected.p,expected.p+p.n,actual.p),"unaligned input oracle");
    std::cout<<"check,common,outer,"<<Outer<<",rank,"<<Rank<<",K,"<<k<<",pass,checksum,"<<hash(expected.p,p.n)<<'\n'<<std::flush;
}
template<bool Stream,unsigned Outer=128,unsigned Rank=4>
static void bench(unsigned calls,unsigned logK,const std::string& mode){
    Plan<Outer> p(std::size_t(1)<<logK,113);c2::BenchPayload input(p.k),out(p.n),scratch(p.scratch());
    produce(input.p,p.k);outer(p,input.p,scratch.p);
    auto run=[&]{
        if(mode=="outer")outer(p,input.p,scratch.p);
        else if(mode=="inner")inner<Stream,Rank>(p,scratch.p,out.p);
        else forward<Stream,Rank>(p,input.p,out.p,scratch.p);
    };
    for(unsigned j=0;j!=3;++j)run();
    std::vector<double> times(calls);std::uint64_t consume=0;
    const auto* result=mode=="outer"?scratch.p:out.p;
    for(unsigned j=0;j!=calls;++j){
        produce(input.p,p.k);const auto start=Clock::now();run();
        times[j]=std::chrono::duration<double,std::micro>(Clock::now()-start).count();
        std::atomic_signal_fence(std::memory_order_seq_cst);consume^=result[(j*104729U)%p.n].lo;
    }
    std::sort(times.begin(),times.end());
    std::cout<<"profile,common,outer,"<<Outer<<",rank,"<<Rank<<",K,"<<p.k<<",N,"<<p.n<<",stream,"<<Stream<<",calls,"<<calls
        <<",checksum,"<<hash(result,p.n)<<",consume,"<<consume<<"\nname,median_us,p10_us,p90_us\n"
        <<std::fixed<<std::setprecision(3)<<mode<<','<<times[calls/2]<<','<<times[calls/10]<<','<<times[9*calls/10]<<'\n';
}
template<unsigned Rank>
static void dispatch(unsigned outer,bool nt,unsigned calls,unsigned logK,const std::string& mode){
    if(outer==128){if(nt)bench<true,128,Rank>(calls,logK,mode);else bench<false,128,Rank>(calls,logK,mode);}
    else{if(nt)bench<true,256,Rank>(calls,logK,mode);else bench<false,256,Rank>(calls,logK,mode);}
}
int main(int argc,char** argv){try{
    if(!c2::gfni512Available())return 77;
    if(argc==2&&std::string(argv[1])=="--check"){
        for(auto k:{4096U,8192U,24576U}){check(k);check<128,2>(k);check<128,1>(k);}
        for(auto k:{8192U,24576U}){check<256>(k);check<256,2>(k);check<256,1>(k);}return 0;
    }
    if(argc==2&&std::string(argv[1])=="--large-check"){check(1U<<20);check<256>(1U<<20);check<128,2>(1U<<20);check<256,2>(1U<<20);return 0;}
    if(argc<2||std::string(argv[1])!="--bench")throw std::invalid_argument("--check | --large-check | --bench [calls] [full|outer|inner] [logK] [nt|cached] [128|256] [1|2|4]");
    const unsigned calls=argc>2?std::stoul(argv[2]):101,logK=argc>4?std::stoul(argv[4]):20;
    const std::string mode=argc>3?argv[3]:"full",stores=argc>5?argv[5]:"nt";
    const unsigned outer=argc>6?std::stoul(argv[6]):128;
    const unsigned rank=argc>7?std::stoul(argv[7]):4;
    if(!calls||logK<12||logK>22||(mode!="full"&&mode!="outer"&&mode!="inner")||(stores!="nt"&&stores!="cached")||(outer!=128&&outer!=256))throw std::invalid_argument("invalid benchmark parameters");
    if(rank==4)dispatch<4>(outer,stores=="nt",calls,logK,mode);
    else if(rank==2)dispatch<2>(outer,stores=="nt",calls,logK,mode);
    else if(rank==1)dispatch<1>(outer,stores=="nt",calls,logK,mode);
    else throw std::invalid_argument("invalid local rank");
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
