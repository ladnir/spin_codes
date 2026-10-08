#include "Cauchy128.h"
#include "Cauchy256.h"
#include "BenchPayload.h"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>
#if SPIN_CAUCHY_PUBLIC
#include <spin/Code.h>
#endif

namespace c2=spin::research::char2;
using Clock=std::chrono::steady_clock;
static void require(bool x,const char* s){if(!x)throw std::runtime_error(s);}
static double micros(Clock::time_point a,Clock::time_point b){return std::chrono::duration<double,std::micro>(b-a).count();}
static c2::Subfield16Kernel schedule(const std::string& s){
    if(s=="scatter")return c2::Subfield16Kernel::Scatter;
    if(s=="region-nt")return c2::Subfield16Kernel::RegionNonTemporal;
    if(s=="region-nt-all")return c2::Subfield16Kernel::RegionAllNonTemporal;
    throw std::invalid_argument("expected scatter|region-nt|region-nt-all");
}
static c2::Cauchy128Variant variant(const std::string& s){
    if(s=="baseline")return c2::Cauchy128Variant::Baseline;
    if(s=="compact")return c2::Cauchy128Variant::Compact;
    if(s=="batch")return c2::Cauchy128Variant::BatchDiagonal;
    if(s=="compact-batch")return c2::Cauchy128Variant::CompactBatch;
    if(s=="lane-fft")return c2::Cauchy128Variant::LaneFft;
    if(s=="compact-lane-fft")return c2::Cauchy128Variant::CompactLaneFft;
    throw std::invalid_argument("expected baseline|compact|batch|compact-batch|lane-fft|compact-lane-fft");
}
static constexpr c2::Cauchy128Variant variants[]={c2::Cauchy128Variant::Baseline,c2::Cauchy128Variant::Compact,
    c2::Cauchy128Variant::BatchDiagonal,c2::Cauchy128Variant::CompactBatch,
    c2::Cauchy128Variant::LaneFft,c2::Cauchy128Variant::CompactLaneFft};
static constexpr const char* variantNames[]={"baseline","compact","batch","compact-batch","lane-fft","compact-lane-fft"};
static constexpr unsigned variantCount=sizeof(variants)/sizeof(variants[0]);
static_assert(variantCount==sizeof(variantNames)/sizeof(variantNames[0]));
// Independent dense Lagrange interpolation, not the native FFT factorization.
struct Oracle128 {
    std::array<std::array<std::uint8_t,256>,256> mul{};
    std::array<std::array<std::uint8_t,87>,41> parity{};
    std::array<std::array<std::array<std::uint8_t,8>,16>,4> a{};
    Oracle128(){
        for(unsigned x=0;x!=256;++x)for(unsigned y=0;y!=256;++y)mul[x][y]=c2::sfMultiply(x,y);
        std::uint8_t invDenom[87];
        for(unsigned j=0;j!=87;++j){std::uint8_t v=1;
            for(unsigned h=0;h!=87;++h)if(h!=j)v=mul[v][j^h];
            invDenom[j]=c2::sfInverse(v);
        }
        for(unsigned r=0;r!=41;++r){const auto x=87+r;std::uint8_t v=1;
            for(unsigned h=0;h!=87;++h)v=mul[v][x^h];
            for(unsigned j=0;j!=87;++j)parity[r][j]=mul[mul[v][c2::sfInverse(x^j)]][invDenom[j]];
        }
        constexpr unsigned w[8]={0,1,92,93,224,225,188,189},phase[4]={0,2,6,4};
        for(unsigned p=0;p!=4;++p)for(unsigned i=0;i!=16;++i)for(unsigned j=0;j!=8;++j){
            const auto x=w[i%8]^(i>=8?80:0)^phase[p]^16^w[j];
            auto power=[&](unsigned n){std::uint8_t v=1;for(unsigned h=0;h!=n;++h)v=mul[v][x];return v;};
            a[p][i][j]=power(7)^power(3)^x^1;
        }
    }
};
static const Oracle128& oracle(){static const Oracle128 o;return o;}
template<unsigned Phase>
static void checkInnerLaneFft(){
    const auto& o=oracle();std::mt19937_64 rng(922+Phase);
    alignas(64) std::uint8_t input[16][16],incoming[8][16],output[16][16],observed[8][16];
    alignas(64) std::uint8_t expected[16][16],next[8][16];
    std::uint8_t d[8];
    for(unsigned trial=0;trial!=416;++trial){
        std::memset(input,0,sizeof(input));std::memset(incoming,0,sizeof(incoming));
        if(trial<384){
            if(trial<256)input[trial/16][trial%16]=1;
            else incoming[(trial-256)/16][trial%16]=1;
        }else{
            for(auto& row:input)for(auto& x:row)x=static_cast<std::uint8_t>(rng());
            for(auto& row:incoming)for(auto& x:row)x=static_cast<std::uint8_t>(rng());
        }
        for(auto& x:d)x=1+rng()%255;
        const __m512i diagonal[2]={c2::c128LaneDiagonalFour(d),c2::c128LaneDiagonalFour(d+4)};
        for(bool forward:{false,true})for(bool alias:{false,true}){
            auto* initial=reinterpret_cast<const __m512i*>(incoming);
            __m512i state[2]={c2::c128LanePermutation(initial[0]),c2::c128LanePermutation(initial[1])};
            std::memcpy(output,input,sizeof(input));
            auto* from=reinterpret_cast<const __m512i*>(alias?output:input);
            auto* to=reinterpret_cast<__m512i*>(output);
            if(forward)c2::c128InnerForwardLaneFftPacked<Phase>(from,to,state,diagonal);
            else c2::c128InnerTransposeLaneFftPacked<Phase>(from,to,state,diagonal);
            _mm512_store_si512(observed,c2::c128LanePermutation(state[0]));
            _mm512_store_si512(observed+4,c2::c128LanePermutation(state[1]));
            for(unsigned j=0;j!=8;++j)for(unsigned b=0;b!=16;++b)
                next[j][b]=forward?o.mul[d[j]][incoming[j][b]]:incoming[j][b];
            for(unsigned i=0;i!=16;++i)for(unsigned b=0;b!=16;++b){
                auto v=input[i][b];
                for(unsigned j=0;j!=8;++j){
                    const auto s=forward?o.mul[d[j]][incoming[j][b]]:incoming[j][b];
                    const auto a=o.mul[o.a[Phase][i][j]][s];
                    v^=forward?o.mul[2][a]:a;
                    const auto m=o.mul[o.a[Phase][i][j]][input[i][b]];
                    next[j][b]^=forward?m:o.mul[2][m];
                }
                expected[i][b]=v;
            }
            if(!forward)for(unsigned j=0;j!=8;++j)for(unsigned b=0;b!=16;++b)
                next[j][b]=o.mul[d[j]][next[j][b]];
            require(std::memcmp(output,expected,sizeof(output))==0,"lane FFT packet output");
            require(std::memcmp(observed,next,sizeof(next))==0,"lane FFT persistent state");
        }
    }
}
static std::vector<std::uint8_t> forward(const c2::Cauchy128Plan& p,const std::vector<std::uint8_t>& m,unsigned bytes){
    require(c2::sfValidWidth(bytes)&&m.size()==p.requested*bytes,"forward dimensions");
    const auto& o=oracle();std::vector<std::uint8_t> planes(p.n*bytes),out(p.n*bytes);
    for(std::size_t q=0;q!=p.requested;++q){
        const auto row=(q/348)*4+q%4;const auto j=(q%348)/4;
        for(unsigned b=0;b!=bytes;++b){auto v=m[q*bytes+b];planes[(j*p.groups+row)*bytes+b]=v;
            for(unsigned h=0;h!=41;++h)planes[((87+h)*p.groups+row)*bytes+b]^=o.mul[o.parity[h][j]][v];}
    }
    std::uint8_t state[8][16]{};
    for(std::size_t step=0;step!=p.n/16;++step){
        std::uint8_t next[8][16]{};
        for(unsigned j=0;j!=8;++j)for(unsigned b=0;b!=bytes;++b)
            next[j][b]=state[j][b]=o.mul[p.stateDiagonal[8*step+j]][state[j][b]];
        for(unsigned i=0;i!=16;++i)for(unsigned b=0;b!=bytes;++b){const auto index=16*step+i;
            const auto x=o.mul[p.inputDiagonal[index]][planes[std::size_t(p.route[index])*bytes+b]];auto y=x;
            for(unsigned j=0;j!=8;++j){y^=o.mul[2][o.mul[o.a[step%4][i][j]][state[j][b]]];
                                      next[j][b]^=o.mul[o.a[step%4][i][j]][x];}
            out[index*bytes+b]=y;
        }
        std::memcpy(state,next,sizeof(state));
    }
    return out;
}
static std::vector<std::uint8_t> transpose(const c2::Cauchy128Plan& p,const std::vector<std::uint8_t>& in,unsigned bytes){
    require(c2::sfValidWidth(bytes)&&in.size()==p.n*bytes,"transpose dimensions");
    const auto& o=oracle();std::vector<std::uint8_t> planes(p.n*bytes),out(p.requested*bytes);
    std::uint8_t state[8][16]{};
    for(std::size_t step=p.n/16;step--;){std::uint8_t next[8][16];std::memcpy(next,state,sizeof(state));
        for(unsigned i=0;i!=16;++i)for(unsigned b=0;b!=bytes;++b){const auto index=16*step+i;auto x=in[index*bytes+b];
            for(unsigned j=0;j!=8;++j){x^=o.mul[o.a[step%4][i][j]][state[j][b]];
                                      next[j][b]^=o.mul[2][o.mul[o.a[step%4][i][j]][in[index*bytes+b]]];}
            planes[std::size_t(p.route[index])*bytes+b]=o.mul[p.inputDiagonal[index]][x];
        }
        for(unsigned j=0;j!=8;++j)for(unsigned b=0;b!=bytes;++b)
            state[j][b]=o.mul[p.stateDiagonal[8*step+j]][next[j][b]];
    }
    for(std::size_t q=0;q!=p.requested;++q){const auto row=(q/348)*4+q%4,j=(q%348)/4;
        for(unsigned b=0;b!=bytes;++b){auto v=planes[(j*p.groups+row)*bytes+b];
            for(unsigned h=0;h!=41;++h)v^=o.mul[o.parity[h][j]][planes[((87+h)*p.groups+row)*bytes+b]];
            out[q*bytes+b]=v;}
    }
    return out;
}
static std::uint8_t dot(const std::vector<std::uint8_t>& a,const std::vector<std::uint8_t>& b){
    require(a.size()==b.size(),"dot dimensions");std::uint8_t v=0;
    for(std::size_t i=0;i!=a.size();++i)v^=oracle().mul[a[i]][b[i]];return v;
}
static void checkOuter(){
    std::mt19937_64 rng(610);const auto& o=oracle();
    alignas(64) std::uint8_t input[128][64],output[128][64];
    for(unsigned trial=0;trial!=4;++trial){
        for(auto& row:input)for(auto& x:row)x=static_cast<std::uint8_t>(rng());
        c2::c128OuterForwardFft(reinterpret_cast<const __m512i*>(input),reinterpret_cast<__m512i*>(output));
        for(unsigned i=0;i!=128;++i)for(unsigned b=0;b!=64;++b){std::uint8_t v=i<87?input[i][b]:0;
            if(i>=87)for(unsigned j=0;j!=87;++j)v^=o.mul[o.parity[i-87][j]][input[j][b]];
            require(output[i][b]==v,"native outer forward/interpolation");}
        c2::c128OuterTransposeFft(reinterpret_cast<const __m512i*>(input),reinterpret_cast<__m512i*>(output));
        for(unsigned j=0;j!=87;++j)for(unsigned b=0;b!=64;++b){auto v=input[j][b];
            for(unsigned h=0;h!=41;++h)v^=o.mul[o.parity[h][j]][input[87+h][b]];
            require(output[j][b]==v,"native outer transpose/interpolation");}
        std::memset(output,0xa7,sizeof(output));
        c2::c128OuterParityTransposeFft(reinterpret_cast<const __m512i*>(input+87),reinterpret_cast<__m512i*>(output));
        for(unsigned j=0;j!=87;++j)for(unsigned b=0;b!=64;++b){std::uint8_t v=0;
            for(unsigned h=0;h!=41;++h)v^=o.mul[o.parity[h][j]][input[87+h][b]];
            require(output[j][b]==v,"native parity/interpolation");}
        for(unsigned j=87;j!=128;++j)for(auto v:output[j])require(v==0xa7,"outer suffix guard");
    }
}
static void check(){
    checkOuter();checkInnerLaneFft<0>();checkInnerLaneFft<1>();checkInnerLaneFft<2>();checkInnerLaneFft<3>();
    std::mt19937_64 rng(618);const c2::Block guard{0x32456789ULL,0xfabc1234ULL};
    for(std::size_t k:{1,86,87,88,347,348,349,5567,5568,5569}){
        c2::Cauchy128Plan p(k,113);
        for(std::size_t i=0;i!=p.n;++i){
            require(p.route[i]/p.groups==i/p.groups&&p.inverseRoute[p.route[i]]==i&&p.inputDiagonal[i],"route or diagonal");
            require(p.inverseLocalRoute[p.route[i]]==i%p.groups,"compact inverse route");
        }
        for(auto d:p.stateDiagonal)require(d!=0,"state diagonal");
        for(unsigned bytes:{1,4,8,16}){
            std::vector<std::uint8_t> m(k*bytes),in(p.n*bytes);
            for(auto& x:m)x=static_cast<std::uint8_t>(rng());for(auto& x:in)x=static_cast<std::uint8_t>(rng());
            const auto expected=transpose(p,in,bytes);
            require(dot(forward(p,m,bytes),in)==dot(m,expected),"scalar adjoint");
            if(bytes!=16)continue;
            c2::BenchPayload input(p.n+4,false),out(k+4,false),planes(p.n+4,false),temp(p.groups+4,false);
            for(auto mode:{c2::Subfield16Kernel::Scatter,c2::Subfield16Kernel::RegionNonTemporal,c2::Subfield16Kernel::RegionAllNonTemporal})
              for(auto selected:variants)for(unsigned offset:{0,1,2,3}){
                for(auto* b:{&input,&out,&planes,&temp})std::fill(b->p,b->p+b->count,guard);
                std::memcpy(input.p+offset,in.data(),in.size());
                c2::c128Transpose(p,input.p+offset,out.p+offset,planes.p+offset,temp.p+offset,mode,selected);
                require(std::memcmp(out.p+offset,expected.data(),expected.size())==0,"native/scalar");
                require(std::memcmp(input.p+offset,in.data(),in.size())==0,"input preservation");
                c2::c128Transpose(p,input.p+offset,input.p+offset,planes.p+offset,temp.p+offset,mode,selected);
                require(std::memcmp(input.p+offset,expected.data(),expected.size())==0,"in-place prefix");
                require(std::memcmp(input.p+offset+k,in.data()+16*k,16*(p.n-k))==0,"in-place suffix");
                if(k==349&&offset==0&&mode==c2::Subfield16Kernel::Scatter){
                    auto again=in;std::memcpy(again.data(),expected.data(),expected.size());const auto twice=transpose(p,again,16);
                    c2::c128Transpose(p,input.p,input.p,planes.p,temp.p,mode,selected);
                    require(std::memcmp(input.p,twice.data(),twice.size())==0,"persistent call state");
                }
                for(auto* b:{&input,&out,&planes,&temp}){for(unsigned j=0;j!=offset;++j)require(b->p[j]==guard,"prefix guard");
                    require(b->p[b->count-4+offset]==guard,"suffix guard");}
            }
            if(k==349){
                // Exercise the same 32-bit fallback used beyond uint16 regions,
                // without requiring a many-million-record allocation in this test.
                auto local=std::move(p.inverseLocalRoute);
                std::memcpy(input.p,in.data(),in.size());
                for(auto selected:{c2::Cauchy128Variant::CompactBatch,c2::Cauchy128Variant::CompactLaneFft}){
                    c2::c128Transpose(p,input.p,out.p,planes.p,temp.p,c2::Subfield16Kernel::RegionAllNonTemporal,selected);
                    require(std::memcmp(out.p,expected.data(),expected.size())==0,"compact fallback");
                }
                p.inverseLocalRoute=std::move(local);
            }
        }
    }
    for(auto n:{std::size_t(0),std::numeric_limits<std::size_t>::max()}){bool bad=false;
        try{c2::Cauchy128Plan p(n,1);}catch(const std::invalid_argument&){bad=true;}require(bad,"size rejection");}
    std::cout<<"Cauchy128: dense interpolation, scalar widths1/4/8/16, adjoints, four-row source ordering, routes, padding, alignment, guards and aliases pass\n";
}
static void checkLarge(){
    std::mt19937_64 rng(197);for(unsigned lg:{16,18,20}){
        c2::Cauchy128Plan p(std::size_t(1)<<lg,113);
        c2::BenchPayload input(p.n),planes(p.n),temp(p.groups),expected(p.requested),out(p.requested),alias(p.n);
        for(std::size_t i=0;i!=p.n;++i)input.p[i]={rng(),rng()};
        c2::c128Transpose(p,input.p,expected.p,planes.p,temp.p,c2::Subfield16Kernel::Scatter);
        for(auto selected:variants)for(auto mode:{c2::Subfield16Kernel::Scatter,c2::Subfield16Kernel::RegionNonTemporal,c2::Subfield16Kernel::RegionAllNonTemporal}){
            c2::c128Transpose(p,input.p,out.p,planes.p,temp.p,mode,selected);
            require(std::equal(expected.p,expected.p+p.requested,out.p),"large schedule");
            std::memcpy(alias.p,input.p,p.n*16);c2::c128Transpose(p,alias.p,alias.p,planes.p,temp.p,mode,selected);
            require(std::equal(expected.p,expected.p+p.requested,alias.p),"large alias prefix");
            require(std::equal(input.p+p.requested,input.p+p.n,alias.p+p.requested),"large alias suffix");
        }
        std::cout<<"K"<<lg<<" schedules and aliases pass\n";
    }
}
static std::uint64_t checksum(const c2::Block* p,std::size_t n){std::uint64_t v=0;
    for(std::size_t j=0;j!=n;++j){v=(v^p[j].lo)*0x100000001b3ULL;v=(v^p[j].hi)*0x100000001b3ULL;}return v;}
static void stats(const char* name,std::size_t k,std::size_t n,std::vector<double> v,std::uint64_t sum){
    std::sort(v.begin(),v.end());std::cout<<name<<','<<k<<','<<n<<','<<v.size()<<','<<std::fixed<<std::setprecision(4)
        <<v[v.size()/2]<<','<<v[v.size()/10]<<','<<v[9*v.size()/10]<<','<<std::hex<<sum<<std::dec<<'\n';
}
static void benchmark(std::size_t k,std::uint64_t seed,unsigned calls,c2::Subfield16Kernel mode,bool phases,
                       c2::Cauchy128Variant selected,const std::string& name){
    const auto start=Clock::now();c2::Cauchy128Plan p(k,seed);const auto setup=micros(start,Clock::now());
    c2::BenchPayload input(p.n),planes(p.n),temp(p.groups);std::mt19937_64 rng(913);
    for(std::size_t i=0;i!=p.n;++i)input.p[i]={rng(),rng()};
    auto encode=[&]{c2::c128Transpose(p,input.p,input.p,planes.p,temp.p,mode,selected);};
    for(unsigned j=0;j!=5;++j)encode();std::vector<double> total(calls),inner(calls),outer(calls);
    for(unsigned j=0;j!=calls;++j){const auto a=Clock::now();if(phases){
        c2::c128InnerDispatch(p,input.p,planes.p,temp.p,mode,selected);const auto b=Clock::now();
        c2::c128OuterDispatch(p,planes.p,input.p,mode);const auto c=Clock::now();
        inner[j]=micros(a,b);outer[j]=micros(b,c);total[j]=micros(a,c);
    }else{encode();total[j]=micros(a,Clock::now());}}
    std::cout<<"setup_us,"<<setup<<",setup_bytes,"<<p.setupBytes()<<",scratch_bytes,"<<16*(p.n+p.groups)
        <<",padded_K,"<<p.k<<",seed,"<<seed<<",variant,"<<name<<",source_order,four-row-interleaved\nname,K,N,calls,median_us,p10_us,p90_us,checksum\n";
    stats("cauchy128",k,p.n,total,checksum(input.p,k));
    if(phases){stats("inner-route",k,p.n,inner,0);stats("outer",k,p.n,outer,0);}
}
static void compareVariants(std::size_t k,std::uint64_t seed,unsigned calls,c2::Subfield16Kernel mode){
    c2::Cauchy128Plan p(k,seed);
    c2::BenchPayload input(p.n),planes(p.n),temp(p.groups);
    std::vector<c2::Block> initial(p.n);std::mt19937_64 rng(913);for(auto& x:initial)x={rng(),rng()};
    auto restore=[&]{std::memcpy(input.p,initial.data(),p.n*16);};
    auto encode=[&](unsigned v){c2::c128Transpose(p,input.p,input.p,planes.p,temp.p,mode,variants[v]);};
    for(unsigned v=0;v!=variantCount;++v){restore();for(unsigned i=0;i!=5;++i)encode(v);}
    std::array<std::vector<double>,variantCount> samples;for(auto& s:samples)s.resize(calls);
    std::array<std::uint64_t,variantCount> sums{};
    for(unsigned done=0,pass=0;done<calls;++pass){const auto count=std::min(5U,calls-done);
        for(unsigned order=0;order!=variantCount;++order){const auto v=(order+pass)%variantCount;restore();
            for(unsigned j=0;j!=count;++j){const auto a=Clock::now();encode(v);samples[v][done+j]=micros(a,Clock::now());}
            sums[v]=(sums[v]^checksum(input.p,k)^count)*0x100000001b3ULL;}
        done+=count;
    }
    for(unsigned v=1;v!=variantCount;++v)require(sums[v]==sums[0],"variant timing checksums differ");
    std::cout<<"same_buffer_variants,"<<static_cast<void*>(input.p)<<",seed,"<<seed
        <<",restore_outside_timing,batch5,rotating_order,diagnostic_allocator\nname,K,N,calls,median_us,p10_us,p90_us,checksum\n";
    for(unsigned v=0;v!=variantCount;++v)stats(variantNames[v],k,p.n,samples[v],sums[v]);
}
#if SPIN_CAUCHY_PUBLIC
static void compare(std::size_t k,std::uint64_t seed,unsigned calls,c2::Subfield16Kernel mode,
                    c2::Cauchy128Variant selected,const std::string& name){
    require(spin::valid_message_size(spin::Parameters::PacketRs8T64,k),"binary baseline K");
    c2::Cauchy128Plan p(k,seed);c2::Cauchy256Plan older(k,seed);
    spin::Code binary({k,spin::Parameters::PacketRs8T64,seed,seed});
    require(binary.backend()==spin::Backend::Avx512,"baseline AVX512 unavailable");
    const auto policy=k<=65536?spin::MemoryPolicy::Normal:spin::MemoryPolicy::PreferHugePages;
    const auto records=std::max({p.n,older.n,2*k});
    spin::Buffer payload((records+12)*16,policy),planeBuffer(p.n*16,policy),tempBuffer(p.groups*16,policy),
        oldPlaneBuffer(older.n*16,policy),oldTempBuffer(older.groups*16,policy);
    auto work=binary.make_workspace(spin::Width::Bits128,policy);
    auto* in=reinterpret_cast<c2::Block*>(payload.bytes().data()+64);
    auto* planes=reinterpret_cast<c2::Block*>(planeBuffer.bytes().data());
    auto* temp=reinterpret_cast<c2::Block*>(tempBuffer.bytes().data());
    auto* oldPlanes=reinterpret_cast<c2::Block*>(oldPlaneBuffer.bytes().data());
    auto* oldTemp=reinterpret_cast<c2::Block*>(oldTempBuffer.bytes().data());
    std::vector<c2::Block> initial(records);std::mt19937_64 rng(913);for(auto& x:initial)x={rng(),rng()};
    auto restore=[&]{std::memcpy(in,initial.data(),records*16);};
    auto encode=[&](unsigned v){if(v==2)c2::c128Transpose(p,in,in,planes,temp,mode,selected);
        else if(v==1)c2::c256Transpose(older,in,in,oldPlanes,oldTemp,mode,c2::Cauchy256Variant::DenseFftFused);
        else binary.transpose_inplace_bytes({reinterpret_cast<std::byte*>(in),2*k*16},work);};
    // Three candidates use the same payload address and initial bytes, but
    // distinct workspaces. Source ordering differs; checksums need not match.
    for(unsigned v=0;v!=3;++v){restore();for(unsigned j=0;j!=5;++j)encode(v);}
    std::array<std::vector<double>,3> samples;for(auto& s:samples)s.resize(calls);
    std::array<std::uint64_t,3> sums{};
    for(unsigned done=0,pass=0;done<calls;++pass){const auto count=std::min(5U,calls-done);
        for(unsigned order=0;order!=3;++order){const auto v=(order+pass)%3;restore();
            for(unsigned j=0;j!=count;++j){const auto a=Clock::now();encode(v);samples[v][done+j]=micros(a,Clock::now());}
            sums[v]=(sums[v]^checksum(in,k)^count)*0x100000001b3ULL;}
        done+=count;
    }
    std::cout<<"same_buffer_address,"<<static_cast<void*>(in)<<",seed,"<<seed<<",variant,"<<name<<",pages,"
        <<(k<=65536?"normal":"huge-preferred")<<",restore_outside_timing,batch5,rotating_order\nname,K,N,calls,median_us,p10_us,p90_us,checksum\n";
    stats("binary-rs8-public",k,2*k,samples[0],sums[0]);stats("cauchy256-dense-fft-fused",k,older.n,samples[1],sums[1]);
    stats("cauchy128",k,p.n,samples[2],sums[2]);
}
#endif
int main(int argc,char** argv){try{
    if(!c2::gfni512Available())return 77;
    if(argc==2&&std::string(argv[1])=="--check"){check();return 0;}
    if(argc==2&&std::string(argv[1])=="--check-large"){checkLarge();return 0;}
    const bool ab=argc>1&&std::string(argv[1])=="--compare";
    const bool allVariants=argc>1&&std::string(argv[1])=="--variants";const int shift=(ab||allVariants)?1:0;
    if(argc<2+shift||argc>7+shift)throw std::invalid_argument("usage: [--compare|--variants] K [seed=113] [calls=51] [scatter|region-nt|region-nt-all] [full|phases] [baseline|compact|batch|compact-batch|lane-fft|compact-lane-fft]");
    const auto k=std::stoull(argv[1+shift]),seed=argc>2+shift?std::stoull(argv[2+shift]):113;
    const auto calls=argc>3+shift?std::stoull(argv[3+shift]):51;
    const auto mode=schedule(argc>4+shift?argv[4+shift]:"scatter");
    const std::string phase=argc>5+shift?argv[5+shift]:"full";
    const std::string name=argc>6+shift?argv[6+shift]:"baseline";const auto selected=variant(name);
    if(!calls||calls>1000000||(phase!="full"&&phase!="phases"))throw std::invalid_argument("mode or calls");
    if(allVariants){require(phase=="full","variant comparison requires full mode");compareVariants(k,seed,static_cast<unsigned>(calls),mode);}
    else if(ab){
#if SPIN_CAUCHY_PUBLIC
        require(phase=="full","comparison requires full mode");compare(k,seed,static_cast<unsigned>(calls),mode,selected,name);
#else
        throw std::invalid_argument("comparison requires SPIN_CHAR2_PUBLIC_ROOT");
#endif
    }else benchmark(k,seed,static_cast<unsigned>(calls),mode,phase=="phases",selected,name);
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
