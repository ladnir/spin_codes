#include "HighDistanceWide.h"
#include "RateSpeed4Oracle.h"
#include "BenchPayload.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <atomic>

namespace c2=spin::research::char2;
namespace hd=c2::hdwide;
static void require(bool value,const char* text){if(!value)throw std::runtime_error(text);}

// Independent extension-field arithmetic; the fixed-map oracle below uses
// literal normalized Cauchy matrices, not the optimized circulant circuit.
template<unsigned Bits>static std::uint32_t multiply(std::uint32_t a,std::uint32_t b){
    if constexpr(Bits==32)return c2::medium32::towerMultiplyOracle(a,b);
    else{
        const auto a0=std::uint8_t(a),a1=std::uint8_t(a>>8),b0=std::uint8_t(b),b1=std::uint8_t(b>>8);
        const auto high=c2::p8OracleMultiply(a1,b1);
        return (c2::p8OracleMultiply(a0,b0)^c2::p8OracleMultiply(0x20,high))|
            (std::uint32_t(c2::p8OracleMultiply(a0,b1)^c2::p8OracleMultiply(a1,b0)^high)<<8);
    }
}
template<unsigned Bits>static void scale(c2::Block& record,std::uint32_t coefficient){
    auto* bytes=reinterpret_cast<unsigned char*>(&record);
    for(unsigned lane=0;lane!=128/Bits;++lane){std::uint32_t value=0;std::memcpy(&value,bytes+lane*(Bits/8),Bits/8);value=multiply<Bits>(value,coefficient);std::memcpy(bytes+lane*(Bits/8),&value,Bits/8);}
}
static void fixedAdd(c2::Block& target,const c2::Block& source,std::uint8_t coefficient){
    auto* out=c2::p8OracleBytes(&target);const auto* in=c2::p8OracleBytes(&source);const auto& table=c2::p8Oracle().mul;
    for(unsigned b=0;b!=16;++b)out[b]^=table[coefficient][in[b]];
}
template<unsigned Bits>static void packet(const c2::Block* input,const c2::Block* state,const typename hd::Field<Bits>::Scalar* dx,const typename hd::Field<Bits>::Scalar* ds,c2::Block* output,c2::Block* next,bool transpose){
    std::array<c2::Block,8> values{},image{};const auto& matrix=c2::rs4Oracle().matrix;
    for(unsigned j=0;j!=4;++j){values[j]=input[j];values[4+j]=state[j];if(!transpose){scale<Bits>(values[j],dx[j]);scale<Bits>(values[4+j],ds[j]);}}
    for(unsigned i=0;i!=8;++i)for(unsigned j=0;j!=8;++j)fixedAdd(image[i],values[j],transpose?matrix[j][i]:matrix[i][j]);
    for(unsigned j=0;j!=4;++j){output[j]=image[j];next[j]=image[4+j];if(transpose){scale<Bits>(output[j],dx[j]);scale<Bits>(next[j],ds[j]);}}
}
// Direct normalized additive-polynomial evaluation at ALL256 AES points.
// This oracle does not call an FFT, butterfly, or optimized basis helper.
static const std::array<std::array<std::uint8_t,64>,256>& evaluation256(){
    static const auto table=[] {
        std::array<std::array<std::uint8_t,64>,256> result{};const auto& o=c2::p8Oracle();
        for(unsigned point=0;point!=256;++point){
            std::uint8_t normalized[6]{};
            for(unsigned bit=0;bit!=6;++bit){std::uint8_t numerator=1,denominator=1;
                for(unsigned label=0;label!=(1U<<bit);++label){numerator=o.mul[numerator][point^label];denominator=o.mul[denominator][(1U<<bit)^label];}
                normalized[bit]=o.mul[numerator][o.inverse[denominator]];
            }
            for(unsigned degree=0;degree!=64;++degree){std::uint8_t value=1;for(unsigned bit=0;bit!=6;++bit)if(degree&(1U<<bit))value=o.mul[value][normalized[bit]];result[point][degree]=value;}
        }
        return result;
    }();return table;
}
static void checkOuter256(){
    hd::Plan<16,2,256> p(256,113);c2::BenchPayload planes(p.planeRecords()),output(256);
    const std::array<c2::Block,4> basis{{{0x0123456789abcdefULL,0xfedcba9876543210ULL},{1,2},{3,5},{7,11}}};
    const auto& evaluation=evaluation256();
    for(unsigned point=0;point!=256;++point){
        if(point)std::fill(planes.p+(point-1)*p.planeStride,planes.p+(point-1)*p.planeStride+4,c2::Block{});
        std::copy(basis.begin(),basis.end(),planes.p+point*p.planeStride);hd::outer<false>(p,planes.p,output.p);
        for(unsigned degree=0;degree!=64;++degree)for(unsigned row=0;row!=4;++row){c2::Block want{};fixedAdd(want,basis[row],evaluation[point][degree]);require(want==output.p[4*degree+row],"wide256 FFT direct polynomial basis");}
    }
    std::cout<<"check,B,256,D,64,all256-point/64-coefficient/four-row-basis,pass\n";
}
template<unsigned Bits,unsigned Keep,unsigned B>static std::vector<c2::Block> oracle(const hd::Plan<Bits,Keep,B>& p,const c2::Block* input,bool transpose){
    using P=hd::Plan<Bits,Keep,B>;std::vector<c2::Block> planes(p.inputN),result(transpose?p.requested:p.outputN);
    std::array<c2::Block,4> state{},next{},in{},out{};const auto& evaluation=evaluation256();
    if(transpose){
        for(auto step=p.inputN/4;step--;){
            in.fill({});for(unsigned j=0;j!=Keep;++j)in[j]=input[P::outputIndex(step,j)];
            packet<Bits>(in.data(),state.data(),p.dx.data()+4*step,p.ds.data()+4*step,out.data(),next.data(),true);
            for(unsigned j=0;j!=4;++j)planes[p.route[4*step+j]]=out[j];state=next;
        }
        for(std::size_t source=0;source!=p.requested;++source){const auto row=(source/(4*P::dimension))*4+source%4;const auto j=(source%(4*P::dimension))/4;
            for(unsigned point=0;point!=B;++point)fixedAdd(result[source],planes[point*p.groups+row],evaluation[point][j]);}
    }else{
        for(std::size_t source=0;source!=p.requested;++source){const auto row=(source/(4*P::dimension))*4+source%4;const auto j=(source%(4*P::dimension))/4;
            for(unsigned point=0;point!=B;++point)fixedAdd(planes[point*p.groups+row],input[source],evaluation[point][j]);}
        for(std::size_t step=0;step!=p.inputN/4;++step){
            for(unsigned j=0;j!=4;++j)in[j]=planes[p.route[4*step+j]];
            packet<Bits>(in.data(),state.data(),p.dx.data()+4*step,p.ds.data()+4*step,out.data(),next.data(),false);
            for(unsigned j=0;j!=Keep;++j)result[P::outputIndex(step,j)]=out[j];state=next;
        }
    }
    return result;
}
template<unsigned Bits>static std::array<std::uint32_t,128/Bits> dot(const c2::Block* a,const c2::Block* b,std::size_t count){
    std::array<std::uint32_t,128/Bits> sum{};
    for(std::size_t i=0;i!=count;++i)for(unsigned lane=0;lane!=128/Bits;++lane){std::uint32_t x=0,y=0;
        std::memcpy(&x,c2::p8OracleBytes(a+i)+lane*(Bits/8),Bits/8);std::memcpy(&y,c2::p8OracleBytes(b+i)+lane*(Bits/8),Bits/8);sum[lane]^=multiply<Bits>(x,y);}
    return sum;
}
template<unsigned Bits,unsigned Keep,unsigned B=64>static void check(){
    using P=hd::Plan<Bits,Keep,B>;constexpr unsigned D=P::dimension;std::mt19937_64 random(913);const c2::Block guard{771,997};
    for(std::size_t k:{std::size_t(1),std::size_t(D-1),std::size_t(D+1),std::size_t(4*D+3),std::size_t(64*D+1)}){
        P p(k,113);std::vector<c2::Block> input(p.outputN),source(k);for(auto& v:input)v={random(),random()};for(auto& v:source)v={random(),random()};
        const auto want=oracle(p,input.data(),true),forward=oracle(p,source.data(),false);
        require(dot<Bits>(source.data(),want.data(),k)==dot<Bits>(forward.data(),input.data(),p.outputN),"wide scalar adjoint");
        c2::BenchPayload in(p.outputN+8),out(k+8),planes(p.planeRecords()+8),temp(p.groups+8);
        for(unsigned offset=0;offset!=4;++offset){
            for(auto* a:{&in,&out,&planes,&temp})std::fill(a->p,a->p+a->count,guard);
            std::copy(input.begin(),input.end(),in.p+offset);
            hd::transpose(p,in.p+offset,out.p+offset,planes.p+offset,temp.p+offset);
            require(std::equal(want.begin(),want.end(),out.p+offset),"wide scalar transpose equality");
            require(std::equal(input.begin(),input.end(),in.p+offset),"wide input preservation");
            hd::transpose(p,in.p+offset,in.p+offset,planes.p+offset,temp.p+offset);
            require(std::equal(want.begin(),want.end(),in.p+offset),"wide exact alias");
            require(std::equal(input.begin()+k,input.end(),in.p+offset+k),"wide alias suffix");
            for(auto* a:{&in,&out,&planes,&temp}){for(unsigned j=0;j!=offset;++j)require(a->p[j]==guard,"wide prefix guard");for(auto j=a->count-8+offset;j!=a->count;++j)require(a->p[j]==guard,"wide suffix guard");}
        }
        p.localInverse.clear();hd::PhaseTimes phases;hd::transpose<true>(p,input.data(),out.p,planes.p,temp.p,&phases);
        require(std::equal(want.begin(),want.end(),out.p),"wide full-index/timed equality");
        std::cout<<"check,bits,"<<Bits<<",B,"<<B<<",keep,"<<Keep<<",K,"<<k<<",scalar/adjoint/alias/offset/guards/wide-index,pass\n";
    }
}
static std::uint64_t checksum(const c2::Block* input,std::size_t count){std::uint64_t h=0;for(std::size_t j=0;j!=count;++j){h=(h^input[j].lo)*0x100000001b3ULL;h=(h^input[j].hi)*0x100000001b3ULL;}return h;}
static void stats(const char* name,std::vector<double> value){std::sort(value.begin(),value.end());std::cout<<name<<','<<std::fixed<<std::setprecision(3)<<value[value.size()/2]<<','<<value[value.size()/10]<<','<<value[9*value.size()/10]<<'\n';}
template<unsigned Bits,unsigned Keep,unsigned B=64>static void bench(unsigned calls,bool timed,bool alias=false){
    constexpr std::size_t k=1U<<20;using P=hd::Plan<Bits,Keep,B>;P p(k,113);
    c2::BenchPayload input(p.outputN),output(alias?4:k),planes(p.planeRecords()),temp(p.groups);auto* destination=alias?input.p:output.p;
    std::mt19937_64 random(913);for(std::size_t i=0;i!=input.count;++i)input.p[i]={random(),random()};
    for(unsigned i=0;i!=3;++i)hd::transpose(p,input.p,destination,planes.p,temp.p);
    std::vector<double> total(calls),inner(calls),route(calls),outer(calls);std::uint64_t consumed=0;
    for(unsigned i=0;i!=calls;++i){hd::PhaseTimes phases;const auto start=hd::Clock::now();
        if(timed)hd::transpose<true>(p,input.p,destination,planes.p,temp.p,&phases);else hd::transpose(p,input.p,destination,planes.p,temp.p);
        total[i]=std::chrono::duration<double,std::micro>(hd::Clock::now()-start).count();inner[i]=phases.inner;route[i]=phases.route;outer[i]=phases.outer;
        std::atomic_signal_fence(std::memory_order_seq_cst);consumed^=destination[(i*104729U)%k].lo;
    }
    consumed^=checksum(destination,k);
    std::cout<<"prototype,bits,"<<Bits<<",B,"<<B<<",D,"<<P::dimension<<",keep,"<<Keep<<",K,"<<k<<",N,"<<p.outputN<<",raw,"<<p.inputN<<",rate,"<<double(k)/p.outputN
        <<",coefficient_bytes_per_raw_record,"<<2*sizeof(typename P::Coefficient)<<",coefficient_read_bytes,"<<2*p.inputN*sizeof(typename P::Coefficient)
        <<",route_index_bytes_per_raw_record,"<<(p.localInverse.empty()?4:2)<<",planes_bytes,"<<p.planeRecords()*sizeof(c2::Block)<<",temp_bytes,"<<p.groups*sizeof(c2::Block)
        <<",calls,"<<calls<<",phase_clocks,"<<timed<<",exact_alias,"<<alias<<",restore_each_call,0,checksum,"<<consumed<<"\nname,median_us,p10_us,p90_us\n";
    stats("full",total);if(timed){stats("inner",inner);stats("route",route);stats("outer",outer);}
}
template<unsigned Bits>static void choose(unsigned keep,unsigned calls,bool timed,bool alias){if(keep==1)bench<Bits,1>(calls,timed,alias);else if(keep==2)bench<Bits,2>(calls,timed,alias);else if(keep==4)bench<Bits,4>(calls,timed,alias);else throw std::invalid_argument("keep must be 1,2,4");}
int main(int argc,char** argv){try{
    if(!c2::gfni512Available())return 77;
    if(argc==2&&std::string(argv[1])=="--check"){check<16,1>();check<16,2>();check<16,4>();check<32,1>();check<32,2>();check<32,4>();checkOuter256();check<16,2,256>();return 0;}
    if(argc==2&&std::string(argv[1])=="--check256"){checkOuter256();check<16,2,256>();return 0;}
    if(argc<4||argc>7||std::string(argv[1])!="--bench")throw std::invalid_argument("usage: --check | --check256 | --bench bits16or32 keep1or2or4 [calls=101] [full|phases|alias|alias-phases] [B=64|256]; B256 requires bits16/keep2; benchmark K is fixed at 1048576");
    const auto bits=std::stoul(argv[2]),keep=std::stoul(argv[3]),calls=argc>4?std::stoul(argv[4]):101;const std::string mode=argc>5?argv[5]:"full";
    const auto B=argc>6?std::stoul(argv[6]):64;
    const bool timed=mode=="phases"||mode=="alias-phases",alias=mode=="alias"||mode=="alias-phases";
    if(!calls||calls>1000000||(mode!="full"&&mode!="phases"&&mode!="alias"&&mode!="alias-phases"))throw std::invalid_argument("calls/mode");
    if(B==256){if(bits!=16||keep!=2)throw std::invalid_argument("B256 requires bits16/keep2");bench<16,2,256>(unsigned(calls),timed,alias);return 0;}
    if(B!=64)throw std::invalid_argument("B must be64 or256");
    if(bits==16)choose<16>(unsigned(keep),unsigned(calls),timed,alias);else if(bits==32)choose<32>(unsigned(keep),unsigned(calls),timed,alias);else throw std::invalid_argument("field bits must be 16 or 32");return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
