#include "HdFastForward.h"
#include "RateSpeed4Oracle.h"
#include "BenchPayload.h"
#include <atomic>
#include <iomanip>
#include <iostream>
#include <string>

namespace c2=spin::research::char2;
namespace hd=c2::hdwide;
namespace fast=c2::hdfast;
static void require(bool v,const char* message){if(!v)throw std::runtime_error(message);}
static std::uint16_t mul16(std::uint16_t a,std::uint16_t b){
    const auto& t=c2::p8Oracle().mul;
    const auto high=t[a>>8][b>>8];
    return std::uint16_t(t[a&255][b&255]^t[0x20][high])
        |(std::uint16_t(t[a&255][b>>8]^t[a>>8][b&255]^high)<<8);
}
static c2::Block scale(c2::Block x,std::uint16_t a){
    std::uint16_t lanes[8];std::memcpy(lanes,&x,16);
    for(auto& v:lanes)v=mul16(v,a);
    std::memcpy(&x,lanes,16);return x;
}
static void add(c2::Block& a,c2::Block b){a.lo^=b.lo;a.hi^=b.hi;}
static const std::array<std::array<std::uint8_t,64>,256>& evaluation(){
    static const auto table=[] {
        std::array<std::array<std::uint8_t,64>,256> r{};const auto& o=c2::p8Oracle();
        for(unsigned p=0;p!=256;++p){
            std::uint8_t normalized[6]{};
            for(unsigned bit=0;bit!=6;++bit){std::uint8_t num=1,den=1;
                for(unsigned j=0;j!=(1U<<bit);++j){num=o.mul[num][p^j];den=o.mul[den][(1U<<bit)^j];}
                normalized[bit]=o.mul[num][o.inverse[den]];
            }
            for(unsigned degree=0;degree!=64;++degree){std::uint8_t v=1;
                for(unsigned bit=0;bit!=6;++bit)if(degree&(1U<<bit))v=o.mul[v][normalized[bit]];
                r[p][degree]=v;
            }
        }return r;
    }();return table;
}
static std::vector<c2::Block> scalarForward(const fast::ForwardPlan& p,const c2::Block* input){
    std::vector<c2::Block> planes(p.inputN),result(p.outputN);
    const auto& e=evaluation();const auto& matrix=c2::rs4Oracle().matrix;
    for(std::size_t i=0;i!=p.requested;++i){const auto row=(i/256)*4+i%4;const auto d=(i%256)/4;
        for(unsigned j=0;j!=256;++j)add(planes[j*p.groups+row],scale(input[i],e[j][d]));
    }
    std::array<c2::Block,4> state{};
    for(std::size_t step=0;step!=p.inputN/4;++step){std::array<c2::Block,8> x{},y{};
        for(unsigned j=0;j!=4;++j){x[j]=scale(planes[p.route[4*step+j]],p.dx[4*step+j]);x[j+4]=scale(state[j],p.ds[4*step+j]);}
        for(unsigned i=0;i!=8;++i)for(unsigned j=0;j!=8;++j)add(y[i],scale(x[j],matrix[i][j]));
        result[p.outputIndex(step,0)]=y[0];result[p.outputIndex(step,1)]=y[1];
        for(unsigned j=0;j!=4;++j)state[j]=y[j+4];
    }return result;
}
static std::array<std::uint16_t,8> dot(const c2::Block* a,const c2::Block* b,std::size_t n){
    std::array<std::uint16_t,8> sum{};
    for(std::size_t i=0;i!=n;++i){std::uint16_t x[8],y[8];std::memcpy(x,a+i,16);std::memcpy(y,b+i,16);
        for(unsigned j=0;j!=8;++j)sum[j]^=mul16(x[j],y[j]);
    }return sum;
}
static void check(){
    const c2::Block guard{771,997};std::mt19937_64 rng(917);
    for(std::size_t k:{1U,63U,64U,65U,259U}){
        fast::ForwardPlan p(k,113);std::vector<c2::Block> source(k);
        for(auto& x:source)x={rng(),rng()};const auto want=scalarForward(p,source.data());
        c2::BenchPayload input(p.outputN+8),output(p.outputN+8),planes(p.planeRecords()+8),temp(2*p.groups+8);
        for(unsigned offset=0;offset!=4;++offset){
            for(auto* a:{&input,&output,&planes,&temp})std::fill(a->p,a->p+a->count,guard);
            std::copy(source.begin(),source.end(),input.p+offset);
            fast::forward(p,input.p+offset,output.p+offset,planes.p+offset,nullptr,temp.p+offset);
            require(std::equal(want.begin(),want.end(),output.p+offset),"forward scalar equality");
            require(std::equal(source.begin(),source.end(),input.p+offset),"forward preserved input");
            fast::forward(p,input.p+offset,input.p+offset,planes.p+offset,nullptr,temp.p+offset);
            require(std::equal(want.begin(),want.end(),input.p+offset),"forward exact alias");
            for(auto* a:{&input,&output,&planes,&temp}){
                for(unsigned j=0;j!=offset;++j)require(a->p[j]==guard,"prefix guard");
                for(auto j=a->count-8+offset;j!=a->count;++j)require(a->p[j]==guard,"suffix guard");
            }
        }
        p.localRoute.clear();hd::PhaseTimes times;
        fast::forward<true>(p,source.data(),output.p,planes.p,&times,temp.p);
        require(std::equal(want.begin(),want.end(),output.p),"full route index/timed equality");
        std::cout<<"check,K,"<<k<<",scalar/alias/offset/guards/full-index,pass\n";
    }
    fast::ForwardPlan p(1U<<16,7919);
    c2::BenchPayload input(p.requested),dual(p.outputN),output(p.outputN),adjoint(p.requested),planes(p.planeRecords()),temp(2*p.groups);
    for(auto* a:{&input,&dual})for(std::size_t i=0;i!=a->count;++i)a->p[i]={rng(),rng()};
    fast::forward(p,input.p,output.p,planes.p,nullptr,temp.p);
    hd::transpose(static_cast<const hd::Plan<16,2,256>&>(p),dual.p,adjoint.p,planes.p,temp.p);
    require(dot(input.p,adjoint.p,p.requested)==dot(output.p,dual.p,p.outputN),"K16 extension-field adjoint");
    std::cout<<"check,K,65536,extension-field-adjoint,pass\n";
}
static std::uint64_t checksum(const c2::Block* p,std::size_t n){
    std::uint64_t h=0;for(std::size_t i=0;i!=n;++i){h=(h^p[i].lo)*0x100000001b3ULL;h=(h^p[i].hi)*0x100000001b3ULL;}return h;
}
static void stats(const char* name,std::vector<double> v){
    std::sort(v.begin(),v.end());std::cout<<name<<','<<std::fixed<<std::setprecision(3)<<v[v.size()/2]<<','<<v[v.size()/10]<<','<<v[9*v.size()/10]<<'\n';
}
static void produce(c2::Block* p,std::size_t n){
    std::uint64_t x=913,y=113;
    for(std::size_t i=0;i!=n;++i){x^=x<<13;x^=x>>7;x^=x<<17;y+=0x9e3779b97f4a7c15ULL;
        _mm_stream_si128(reinterpret_cast<__m128i*>(p+i),_mm_set_epi64x(y,x));
    }_mm_sfence();
}
static void bench(unsigned calls,bool phases,bool produced,unsigned klog){
    fast::ForwardPlan p(std::size_t(1)<<klog,113);
    c2::BenchPayload input(p.requested),output(p.outputN),planes(p.planeRecords()),temp(2*p.groups);produce(input.p,input.count);
    const auto inputHash=checksum(input.p,input.count);
    for(unsigned i=0;i!=3;++i)fast::forward(p,input.p,output.p,planes.p,nullptr,temp.p);
    std::vector<double> total(calls),outer(calls),inner(calls),route(calls);std::uint64_t consumed=0;
    for(unsigned i=0;i!=calls;++i){if(produced)produce(input.p,input.count);hd::PhaseTimes t;
        const auto start=hd::Clock::now();
        if(phases)fast::forward<true>(p,input.p,output.p,planes.p,&t,temp.p);else fast::forward(p,input.p,output.p,planes.p,nullptr,temp.p);
        total[i]=std::chrono::duration<double,std::micro>(hd::Clock::now()-start).count();outer[i]=t.outer;inner[i]=t.inner;route[i]=t.route;
        std::atomic_signal_fence(std::memory_order_seq_cst);consumed^=output.p[(i*104729U)%output.count].lo;
    }
    require(inputHash==checksum(input.p,input.count),"benchmark input modified");
    std::cout<<"profile,forward-grs256-64-gf16-keep2,K,"<<p.requested<<",N,"<<p.outputN<<",raw,"<<p.inputN<<",calls,"<<calls
        <<",produced_input,"<<produced<<",phase_clocks,"<<phases<<",output_checksum,"<<checksum(output.p,output.count)<<",consume,"<<consumed
        <<",source_preserved,1,outer_tile,"<<HDF_OUTER_TILE<<",split_route,"<<HDF_SPLIT_ROUTE<<",prefetch_region,"<<HDF_PREFETCH_REGION
        <<",hardware_gather,"<<HDF_HARDWARE_GATHER<<",route_copy,"<<HDF_ROUTE_COPY<<"\nname,median_us,p10_us,p90_us\n";
    stats("full",total);if(phases){stats("outer",outer);stats(HDF_SPLIT_ROUTE?"inner":"route_inner",inner);if(HDF_SPLIT_ROUTE)stats("route",route);}
}
int main(int argc,char** argv){try{
    if(!c2::gfni512Available())return 77;
    if(argc==2&&std::string(argv[1])=="--check"){check();return 0;}
    if(argc<2||std::string(argv[1])!="--bench")throw std::invalid_argument("--check | --bench [calls=31] [full|phases] [produced|reused] [klog=20]");
    const unsigned calls=argc>2?std::stoul(argv[2]):31,klog=argc>5?std::stoul(argv[5]):20;
    const std::string mode=argc>3?argv[3]:"full",policy=argc>4?argv[4]:"produced";
    if(!calls||calls>1000000||klog<12||klog>24||(mode!="full"&&mode!="phases")||(policy!="produced"&&policy!="reused"))throw std::invalid_argument("benchmark arguments");
    bench(calls,mode=="phases",policy=="produced",klog);return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
