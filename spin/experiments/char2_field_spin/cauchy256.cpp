#include "Cauchy256.h"
#include "BenchPayload.h"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <span>
#include <string>
#if SPIN_CAUCHY_PUBLIC
#include <spin/Code.h>
#endif

namespace c2=spin::research::char2;
using Clock=std::chrono::steady_clock;
static void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
static double us(Clock::time_point a,Clock::time_point b){return std::chrono::duration<double,std::micro>(b-a).count();}
static c2::Subfield16Kernel kernel(const std::string& name){
    if(name=="scatter")return c2::Subfield16Kernel::Scatter;
    if(name=="region-nt")return c2::Subfield16Kernel::RegionNonTemporal;
    if(name=="region-nt-all")return c2::Subfield16Kernel::RegionAllNonTemporal;
    throw std::invalid_argument("expected scatter|region-nt|region-nt-all");
}
static c2::Cauchy256Variant variant(const std::string& name){
    if(name=="basis")return c2::Cauchy256Variant::Basis;
    if(name=="dense")return c2::Cauchy256Variant::Dense;
    if(name=="fft")return c2::Cauchy256Variant::Fft;
    if(name=="dense-fft")return c2::Cauchy256Variant::DenseFft;
    if(name=="dense-fft-fused")return c2::Cauchy256Variant::DenseFftFused;
    if(name=="fft-inner-fused")return c2::Cauchy256Variant::FftInnerFused;
    throw std::invalid_argument("expected basis|dense|fft|dense-fft|dense-fft-fused|fft-inner-fused");
}
static constexpr c2::Cauchy256Variant variants[]={c2::Cauchy256Variant::Basis,c2::Cauchy256Variant::Dense,
    c2::Cauchy256Variant::Fft,c2::Cauchy256Variant::DenseFft,c2::Cauchy256Variant::DenseFftFused,
    c2::Cauchy256Variant::FftInnerFused};

// Independent dense scalar interpolation and Cauchy formulas. These do not use
// the generated native circuit or its factorization.
struct Oracle {
    std::array<std::array<std::uint8_t,256>,256> mul{};
    std::array<std::array<std::uint8_t,176>,80> parity{};
    std::array<std::array<std::array<std::uint8_t,8>,16>,4> a{};
    Oracle(){
        for(unsigned x=0;x!=256;++x)for(unsigned y=0;y!=256;++y)mul[x][y]=c2::sfMultiply(x,y);
        std::uint8_t denominator[176];
        for(unsigned j=0;j!=176;++j){
            std::uint8_t d=1;
            for(unsigned h=0;h!=176;++h)if(h!=j)d=mul[d][h^j];
            denominator[j]=c2::sfInverse(d);
        }
        for(unsigned r=0;r!=80;++r){
            const unsigned x=176+r;
            std::uint8_t product=1;
            for(unsigned j=0;j!=176;++j)product=mul[product][x^j];
            for(unsigned j=0;j!=176;++j)parity[r][j]=mul[mul[product][c2::sfInverse(x^j)]][denominator[j]];
        }
        constexpr unsigned w[8]={0,1,92,93,224,225,188,189},phase[4]={0,2,6,4};
        for(unsigned p=0;p!=4;++p)for(unsigned i=0;i!=16;++i)for(unsigned j=0;j!=8;++j){
            const auto x=w[i%8]^(i>=8?80:0)^phase[p]^16^w[j];
            auto power=[&](unsigned n){std::uint8_t z=1;for(unsigned t=0;t!=n;++t)z=mul[z][x];return z;};
            a[p][i][j]=power(7)^power(3)^x^1;
        }
    }
};
static const Oracle& oracle(){static const Oracle value;return value;}

static std::vector<std::uint8_t> forward(const c2::Cauchy256Plan& p,const std::vector<std::uint8_t>& message,unsigned bytes){
    require(c2::sfValidWidth(bytes)&&message.size()==p.requested*bytes,"forward oracle dimensions");
    const auto& o=oracle();
    std::vector<std::uint8_t> planes(p.n*bytes),out(p.n*bytes);
    for(std::size_t row=0;row<(p.requested+175)/176;++row)
      for(unsigned j=0;j!=176 && 176*row+j<p.requested;++j)
        for(unsigned b=0;b!=bytes;++b){
            auto v=message[(176*row+j)*bytes+b];
            planes[(j*p.groups+row)*bytes+b]=v;
            for(unsigned t=0;t!=80;++t)planes[((176+t)*p.groups+row)*bytes+b]^=o.mul[o.parity[t][j]][v];
        }
    std::uint8_t state[8][16]{};
    for(std::size_t step=0;step!=p.n/16;++step){
        std::uint8_t next[8][16]{};
        for(unsigned j=0;j!=8;++j)for(unsigned b=0;b!=bytes;++b)
            next[j][b]=state[j][b]=o.mul[p.stateDiagonal[8*step+j]][state[j][b]];
        for(unsigned i=0;i!=16;++i)for(unsigned b=0;b!=bytes;++b){
            const auto index=step*16+i;
            const auto x=o.mul[p.inputDiagonal[index]][planes[std::size_t(p.route[index])*bytes+b]];
            auto y=x;
            for(unsigned j=0;j!=8;++j){
                y^=o.mul[2][o.mul[o.a[step%4][i][j]][state[j][b]]];
                next[j][b]^=o.mul[o.a[step%4][i][j]][x];
            }
            out[index*bytes+b]=y;
        }
        std::memcpy(state,next,sizeof(state));
    }
    return out;
}
static std::vector<std::uint8_t> transpose(const c2::Cauchy256Plan& p,const std::vector<std::uint8_t>& input,unsigned bytes){
    require(c2::sfValidWidth(bytes)&&input.size()==p.n*bytes,"transpose oracle dimensions");
    const auto& o=oracle();
    std::vector<std::uint8_t> planes(p.n*bytes),out(p.requested*bytes);
    std::uint8_t state[8][16]{};
    for(std::size_t step=p.n/16;step--;){
        std::uint8_t next[8][16];std::memcpy(next,state,sizeof(state));
        for(unsigned i=0;i!=16;++i)for(unsigned b=0;b!=bytes;++b){
            const auto index=16*step+i;const auto v=input[index*bytes+b];auto x=v;
            for(unsigned j=0;j!=8;++j){
                x^=o.mul[o.a[step%4][i][j]][state[j][b]];
                next[j][b]^=o.mul[2][o.mul[o.a[step%4][i][j]][v]];
            }
            planes[std::size_t(p.route[index])*bytes+b]=o.mul[p.inputDiagonal[index]][x];
        }
        for(unsigned j=0;j!=8;++j)for(unsigned b=0;b!=bytes;++b)
            state[j][b]=o.mul[p.stateDiagonal[8*step+j]][next[j][b]];
    }
    for(std::size_t row=0;row<(p.requested+175)/176;++row)
      for(unsigned j=0;j!=176 && 176*row+j<p.requested;++j)
        for(unsigned b=0;b!=bytes;++b){
            auto v=planes[(j*p.groups+row)*bytes+b];
            for(unsigned t=0;t!=80;++t)v^=o.mul[o.parity[t][j]][planes[((176+t)*p.groups+row)*bytes+b]];
            out[(176*row+j)*bytes+b]=v;
        }
    return out;
}
static std::uint8_t dot(const std::vector<std::uint8_t>& a,const std::vector<std::uint8_t>& b){
    require(a.size()==b.size(),"dot size");std::uint8_t v=0;
    for(std::size_t j=0;j!=a.size();++j)v^=oracle().mul[a[j]][b[j]];
    return v;
}
static void checkOuter(){
    std::mt19937_64 rng(671);const auto& o=oracle();
    alignas(64) std::uint8_t in[256][64],out[256][64];
    for(unsigned test=0;test!=3;++test){
        for(auto& row:in)for(auto& v:row)v=static_cast<std::uint8_t>(rng());
        for(bool fft:{false,true}){
        if(fft)c2::c256OuterForwardFft(reinterpret_cast<__m512i*>(in),reinterpret_cast<__m512i*>(out));
        else c2::c256OuterForward(reinterpret_cast<__m512i*>(in),reinterpret_cast<__m512i*>(out));
        for(unsigned j=0;j!=256;++j)for(unsigned b=0;b!=64;++b){
            std::uint8_t v=j<176?in[j][b]:0;
            if(j>=176)for(unsigned h=0;h!=176;++h)v^=o.mul[o.parity[j-176][h]][in[h][b]];
            require(v==out[j][b],"native outer forward/interpolation mismatch");
        }
        if(fft)c2::c256OuterTransposeFft(reinterpret_cast<__m512i*>(in),reinterpret_cast<__m512i*>(out));
        else c2::c256OuterTranspose(reinterpret_cast<__m512i*>(in),reinterpret_cast<__m512i*>(out));
        for(unsigned j=0;j!=176;++j)for(unsigned b=0;b!=64;++b){
            auto v=in[j][b];for(unsigned h=0;h!=80;++h)v^=o.mul[o.parity[h][j]][in[h+176][b]];
            require(v==out[j][b],"native outer transpose/interpolation mismatch");
        }
        }
        std::memset(out+176,0xa7,80*64);
        c2::c256OuterParityTransposeFft(reinterpret_cast<const __m512i*>(in+176),reinterpret_cast<__m512i*>(out));
        for(unsigned j=0;j!=176;++j)for(unsigned b=0;b!=64;++b){
            std::uint8_t v=0;
            for(unsigned h=0;h!=80;++h)v^=o.mul[o.parity[h][j]][in[h+176][b]];
            require(v==out[j][b],"native parity-only FFT/interpolation mismatch");
        }
        for(unsigned j=176;j!=256;++j)for(auto v:out[j])require(v==0xa7,"parity-only FFT suffix overwritten");
    }
}
template<unsigned Phase,bool Dense,bool Fft=false>
static void checkInnerPacket(const std::uint8_t* input,std::uint8_t* output,std::uint8_t* state,
                             const std::uint8_t* diagonal,bool transpose){
    const __m512i packedDiagonal[2]={c2::sfBroadcastFour(diagonal),c2::sfBroadcastFour(diagonal+4)};
    const auto* in=reinterpret_cast<const __m512i*>(input);
    auto* out=reinterpret_cast<__m512i*>(output);
    auto* persistent=reinterpret_cast<__m512i*>(state);
    static_assert(!Dense || !Fft);
    if constexpr(Fft){
        if(transpose)c2::c256InnerTransposeFftPacked<Phase>(in,out,persistent,packedDiagonal);
        else c2::c256InnerForwardFftPacked<Phase>(in,out,persistent,packedDiagonal);
    }else if constexpr(Dense){
        if(transpose)c2::c256InnerTransposeDensePacked<Phase>(in,out,persistent,packedDiagonal);
        else c2::c256InnerForwardDensePacked<Phase>(in,out,persistent,packedDiagonal);
    }else{
        if(transpose)c2::c256InnerTransposePacked<Phase>(in,out,persistent,packedDiagonal);
        else c2::c256InnerForwardPacked<Phase>(in,out,persistent,packedDiagonal);
    }
}
template<bool Dense,bool Fft=false>
static void checkInnerPacketPhase(const std::uint8_t* input,std::uint8_t* output,std::uint8_t* state,
                                  const std::uint8_t* diagonal,unsigned phase,bool transpose){
    switch(phase){
    case 0:checkInnerPacket<0,Dense,Fft>(input,output,state,diagonal,transpose);break;
    case 1:checkInnerPacket<1,Dense,Fft>(input,output,state,diagonal,transpose);break;
    case 2:checkInnerPacket<2,Dense,Fft>(input,output,state,diagonal,transpose);break;
    case 3:checkInnerPacket<3,Dense,Fft>(input,output,state,diagonal,transpose);break;
    default:throw std::invalid_argument("inner test phase");
    }
}
static void checkInner(){
    std::mt19937_64 rng(195);const auto& o=oracle();
    alignas(64) std::uint8_t in[16][16],working[16][16],out[16][16],state[8][16],saved[8][16],diagonal[8];
    for(unsigned phase=0;phase!=4;++phase)for(unsigned trial=0;trial!=8;++trial)
      for(bool trans:{false,true}){
        for(auto& row:in)for(auto& x:row)x=static_cast<std::uint8_t>(rng());
        for(auto& row:state)for(auto& x:row)x=static_cast<std::uint8_t>(rng());
        for(auto& d:diagonal)d=c2::sfSampleNonzero(rng);
        std::memcpy(saved,state,sizeof(saved));
        for(unsigned inner:{0,1,2})for(bool alias:{false,true}){
        std::memcpy(state,saved,sizeof(saved));std::memcpy(working,in,sizeof(in));
        auto* result=alias?working:out;
        if(inner==2)checkInnerPacketPhase<false,true>(&working[0][0],&result[0][0],&state[0][0],diagonal,phase,trans);
        else if(inner==1)checkInnerPacketPhase<true>(&working[0][0],&result[0][0],&state[0][0],diagonal,phase,trans);
        else checkInnerPacketPhase<false>(&working[0][0],&result[0][0],&state[0][0],diagonal,phase,trans);
        if(!alias)require(std::memcmp(working,in,sizeof(in))==0,"out-of-place packet input modified");
        for(unsigned b=0;b!=16;++b){
            for(unsigned i=0;i!=16;++i){auto v=in[i][b];
                for(unsigned j=0;j!=8;++j){auto h=saved[j][b];
                    if(!trans)h=o.mul[2][o.mul[diagonal[j]][h]];
                    v^=o.mul[o.a[phase][i][j]][h];
                }require(result[i][b]==v,"native inner output/scalar formula mismatch");
            }
            for(unsigned j=0;j!=8;++j){auto h=saved[j][b];if(!trans)h=o.mul[diagonal[j]][h];
                for(unsigned i=0;i!=16;++i){auto v=o.mul[o.a[phase][i][j]][in[i][b]];
                    h^=trans?o.mul[2][v]:v;
                }if(trans)h=o.mul[diagonal[j]][h];
                require(state[j][b]==h,"native inner state/scalar formula mismatch");
            }
        }
        }
    }
}
static void check(){
    checkOuter();checkInner();std::mt19937_64 rng(711);
    const c2::Block guard{0x914fab34f21f3aceULL,0x7126334455667788ULL};
    for(std::size_t k:{1,175,176,177,703,704,705,11263,11264,11265}){
        c2::Cauchy256Plan p(k,113);
        for(std::size_t i=0;i!=p.n;++i){
            require(p.route[i]/p.groups==i/p.groups&&p.inverseRoute[p.route[i]]==i,"regional route mismatch");
            require(p.inputDiagonal[i]!=0,"zero input diagonal");
        }
        for(auto d:p.stateDiagonal)require(d!=0,"zero state diagonal");
        for(unsigned bytes:{1,4,8,16}){
            std::vector<std::uint8_t> m(k*bytes),v(p.n*bytes);
            for(auto& x:m)x=static_cast<std::uint8_t>(rng());
            for(auto& x:v)x=static_cast<std::uint8_t>(rng());
            auto y=forward(p,m,bytes),expected=transpose(p,v,bytes);
            require(dot(y,v)==dot(m,expected),"scalar field adjoint mismatch");
            if(bytes!=16)continue;
            c2::BenchPayload input(p.n+4,false),output(k+4,false),planes(p.n+4,false),temp(p.groups+4,false);
            for(auto schedule:{c2::Subfield16Kernel::Scatter,c2::Subfield16Kernel::RegionNonTemporal,c2::Subfield16Kernel::RegionAllNonTemporal})
              for(auto selected:variants)for(unsigned offset:{0,1,2,3}){
                for(auto* buffer:{&input,&output,&planes,&temp})std::fill(buffer->p,buffer->p+buffer->count,guard);
                std::memcpy(input.p+offset,v.data(),v.size());
                c2::c256Transpose(p,input.p+offset,output.p+offset,planes.p+offset,temp.p+offset,schedule,selected);
                require(std::memcmp(output.p+offset,expected.data(),expected.size())==0,"SIMD/scalar transpose mismatch");
                require(std::memcmp(input.p+offset,v.data(),v.size())==0,"out-of-place input modified");
                c2::c256Transpose(p,input.p+offset,input.p+offset,planes.p+offset,temp.p+offset,schedule,selected);
                require(std::memcmp(input.p+offset,expected.data(),expected.size())==0,"in-place prefix mismatch");
                require(std::memcmp(input.p+offset+k,v.data()+16*k,16*(p.n-k))==0,"in-place requested suffix modified");
                if(k==705&&offset==0&&schedule==c2::Subfield16Kernel::Scatter){
                    auto updated=v;std::memcpy(updated.data(),expected.data(),expected.size());
                    const auto expectedTwice=transpose(p,updated,16);
                    c2::c256Transpose(p,input.p,input.p,planes.p,temp.p,schedule,selected);
                    require(std::memcmp(input.p,expectedTwice.data(),expectedTwice.size())==0,"state leaked between calls");
                    require(std::memcmp(input.p+k,v.data()+16*k,16*(p.n-k))==0,"repeat modified suffix");
                }
                for(auto* buffer:{&input,&output,&planes,&temp}){
                    for(unsigned j=0;j!=offset;++j)require(buffer->p[j]==guard,"prefix guard modified");
                    require(buffer->p[buffer->count-4+offset]==guard,"suffix guard modified");
                }
            }
        }
    }
    for(auto bad:{std::size_t(0),std::numeric_limits<std::size_t>::max()}){
        bool rejected=false;try{c2::Cauchy256Plan p(bad,1);}catch(const std::invalid_argument&){rejected=true;}
        require(rejected,"invalid size accepted");
    }
    std::cout<<"Cauchy256 all variants: independent interpolation, scalar widths1/4/8/16, adjoints, packet aliases, route, schedules, padding, alignments, guards and aliases pass\n";
}
static void checkLarge(){
    std::mt19937_64 rng(197);for(unsigned lg:{16,18,20}){
        c2::Cauchy256Plan p(std::size_t(1)<<lg,113);
        c2::BenchPayload input(p.n),planes(p.n),temp(p.groups),expected(p.requested),out(p.requested),alias(p.n);
        for(std::size_t i=0;i!=p.n;++i)input.p[i]={rng(),rng()};
        c2::c256Transpose(p,input.p,expected.p,planes.p,temp.p,c2::Subfield16Kernel::Scatter);
        for(auto selected:variants){
        for(auto schedule:{c2::Subfield16Kernel::Scatter,c2::Subfield16Kernel::RegionNonTemporal,c2::Subfield16Kernel::RegionAllNonTemporal}){
            c2::c256Transpose(p,input.p,out.p,planes.p,temp.p,schedule,selected);
            require(std::equal(expected.p,expected.p+p.requested,out.p),"large schedule mismatch");
        }
        std::memcpy(alias.p,input.p,p.n*sizeof(c2::Block));
        c2::c256Transpose(p,alias.p,alias.p,planes.p,temp.p,c2::Subfield16Kernel::RegionAllNonTemporal,selected);
        require(std::equal(expected.p,expected.p+p.requested,alias.p),"large alias prefix mismatch");
        require(std::equal(input.p+p.requested,input.p+p.n,alias.p+p.requested),"large alias suffix mismatch");
        }
        std::cout<<"K"<<lg<<" all variants, schedules and in-place tail pass\n";
    }
}
static std::uint64_t checksum(const c2::Block* p,std::size_t n){
    std::uint64_t v=0;for(std::size_t i=0;i!=n;++i){v=(v^p[i].lo)*0x100000001b3ULL;v=(v^p[i].hi)*0x100000001b3ULL;}return v;
}
static void stats(const char* name,std::size_t k,std::size_t n,std::vector<double> values,std::uint64_t sum){
    std::sort(values.begin(),values.end());
    std::cout<<name<<','<<k<<','<<n<<','<<values.size()<<','<<std::fixed<<std::setprecision(4)
        <<values[values.size()/2]<<','<<values[values.size()/10]<<','<<values[9*values.size()/10]<<','<<std::hex<<sum<<std::dec<<'\n';
}
static void benchmark(std::size_t k,std::uint64_t seed,unsigned calls,const std::string& schedule,bool phases,
                       const std::string& variantName){
    const auto start=Clock::now();c2::Cauchy256Plan p(k,seed);const auto setup=us(start,Clock::now());
    c2::BenchPayload input(p.n),planes(p.n),temp(p.groups);std::mt19937_64 rng(913);
    for(std::size_t i=0;i!=p.n;++i)input.p[i]={rng(),rng()};
    const auto mode=kernel(schedule);const auto selected=variant(variantName);
    auto encode=[&]{c2::c256Transpose(p,input.p,input.p,planes.p,temp.p,mode,selected);};
    for(unsigned i=0;i!=5;++i)encode();
    std::vector<double> total(calls),inner(calls),outer(calls);
    for(unsigned i=0;i!=calls;++i){
        const auto a=Clock::now();
        if(phases){c2::c256InnerDispatch(p,input.p,planes.p,temp.p,mode,selected);const auto b=Clock::now();
            c2::c256OuterDispatch(p,planes.p,input.p,mode,selected);const auto c=Clock::now();
            inner[i]=us(a,b);outer[i]=us(b,c);total[i]=us(a,c);
        }else{encode();total[i]=us(a,Clock::now());}
    }
    std::cout<<"setup_us,"<<setup<<",setup_bytes,"<<p.setupBytes()<<",scratch_bytes,"<<(p.n+p.groups)*16<<",padded_K,"<<p.k
        <<",seed,"<<seed<<",schedule,"<<schedule<<",mode,"<<(phases?"phases":"full")<<",variant,"<<variantName<<"\nname,K,N,calls,median_us,p10_us,p90_us,checksum\n";
    stats("cauchy256",k,p.n,total,checksum(input.p,k));
    if(phases){stats("inner-route",k,p.n,inner,0);stats("outer",k,p.n,outer,0);}
}
#if SPIN_CAUCHY_PUBLIC
static void compare(std::size_t k,std::uint64_t seed,unsigned calls,const std::string& schedule,
                     const std::string& variantName){
    require(spin::valid_message_size(spin::Parameters::PacketRs8T64,k),"comparison needs valid binary baseline K");
    c2::Cauchy256Plan p(k,seed);spin::Code baseline({k,spin::Parameters::PacketRs8T64,seed,seed});
    require(baseline.backend()==spin::Backend::Avx512,"baseline AVX512 unavailable");
    const auto policy=k<=65536?spin::MemoryPolicy::Normal:spin::MemoryPolicy::PreferHugePages;
    const auto records=std::max(p.n,2*k);
    spin::Buffer payload((records+12)*16,policy),planeStorage(p.n*16,policy),tempStorage(p.groups*16,policy);
    auto work=baseline.make_workspace(spin::Width::Bits128,policy);
    auto* input=reinterpret_cast<c2::Block*>(payload.bytes().data()+64);
    auto* planes=reinterpret_cast<c2::Block*>(planeStorage.bytes().data());
    auto* temp=reinterpret_cast<c2::Block*>(tempStorage.bytes().data());
    std::vector<c2::Block> initial(records);std::mt19937_64 rng(913);for(auto& x:initial)x={rng(),rng()};
    const auto mode=kernel(schedule);const auto selected=variant(variantName);
    auto encode=[&](unsigned v){if(v)c2::c256Transpose(p,input,input,planes,temp,mode,selected);
        else baseline.transpose_inplace_bytes({reinterpret_cast<std::byte*>(input),2*k*16},work);};
    auto restore=[&]{std::memcpy(input,initial.data(),records*16);};
    for(unsigned v=0;v!=2;++v){restore();for(unsigned j=0;j!=5;++j)encode(v);}
    std::array<std::vector<double>,2> samples;for(auto& s:samples)s.resize(calls);
    std::array<std::uint64_t,2> sums{};
    for(unsigned completed=0,pass=0;completed<calls;++pass){
        const auto count=std::min(5U,calls-completed);
        for(unsigned order=0;order!=2;++order){const auto v=order^(pass&1U);restore();
            for(unsigned j=0;j!=count;++j){const auto a=Clock::now();encode(v);samples[v][completed+j]=us(a,Clock::now());}
            sums[v]=(sums[v]^checksum(input,k)^count)*0x100000001b3ULL;
        }
        completed+=count;
    }
    std::cout<<"same_buffer_address,"<<static_cast<void*>(input)<<",seed,"<<seed<<",schedule,"<<schedule<<",variant,"<<variantName
        <<",pages,"<<(k<=65536?"normal":"huge-preferred")<<",restore_outside_timing,batch5,alternating_order\n"
        <<"name,K,N,calls,median_us,p10_us,p90_us,checksum\n";
    stats("binary-rs8-public",k,2*k,samples[0],sums[0]);stats("cauchy256",k,p.n,samples[1],sums[1]);
}
#endif
int main(int argc,char** argv){try{
    if(!c2::gfni512Available())return 77;
    if(argc==2&&std::string(argv[1])=="--check"){check();return 0;}
    if(argc==2&&std::string(argv[1])=="--check-large"){checkLarge();return 0;}
    const bool ab=argc>1&&std::string(argv[1])=="--compare";const int shift=ab?1:0;
    if(argc<2+shift||argc>7+shift)throw std::invalid_argument("usage: [--compare] K [seed=113] [calls=51] [scatter|region-nt|region-nt-all] [full|phases] [basis|dense|fft|dense-fft|dense-fft-fused|fft-inner-fused]");
    const auto k=std::stoull(argv[1+shift]),seed=argc>2+shift?std::stoull(argv[2+shift]):113;
    const auto count=argc>3+shift?std::stoull(argv[3+shift]):51;
    const std::string schedule=argc>4+shift?argv[4+shift]:"scatter",mode=argc>5+shift?argv[5+shift]:"full";
    const std::string variantName=argc>6+shift?argv[6+shift]:"basis";
    (void)variant(variantName);
    if(!count||count>1000000||(mode!="full"&&mode!="phases"))throw std::invalid_argument("invalid benchmark mode/call count");
    if(ab){
#if SPIN_CAUCHY_PUBLIC
        if(mode!="full")throw std::invalid_argument("comparison requires full mode");
        compare(k,seed,static_cast<unsigned>(count),schedule,variantName);
#else
        throw std::invalid_argument("comparison requires SPIN_CHAR2_PUBLIC_ROOT at configure time");
#endif
    }else benchmark(k,seed,static_cast<unsigned>(count),schedule,mode=="phases",variantName);
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
