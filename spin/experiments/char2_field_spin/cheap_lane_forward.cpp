// Isolated research experiment. Preserve the certified encoder and its setup.
#define HDGROUP_SUBFIELD 1
#define HDGROUP_FUSE_LAST 1
#define HDGROUP_BYTE_KEYS 1
#define HD_SUMRANK_BYTE_KEYS 1
#define HDGROUP_INNER 5
#define HDGROUP_FUSE_REFRESH 1
#define HDGROUP_NATIVE_OUTPUT 1
#define HDGROUP_TOEPLITZ 1
#define HDGROUP_PACKET_ROUTE 1
#define HDGROUP_ROUTE_WIDTH 2
#define HDGROUP_SHIFTED_PAIRS 2
#define HDGROUP_PAIR_SWAPS 1
#define main preserved_grouped_main
#include "hd_grouped_forward.cpp"
#undef main
#include "field128_lift/SparseByteField128MulYSimd.h"

namespace lane_trial {
namespace tower=c2::field128_lift::sparse_byte;
struct RefreshKeys {
    std::vector<std::array<std::uint8_t,16>> a;
    explicit RefreshKeys(std::size_t n,std::uint64_t seed):a(n/64){
        std::mt19937_64 rng(seed);
        for(auto& step:a)for(auto& x:step)x=static_cast<std::uint8_t>(rng());
    }
    std::size_t bytes()const{return a.size()*16;}
};

// Extend the existing final nonzero byte diagonal b by (Y+a). No general
// extension-field product: one fixed companion map plus one byte product.
template<bool Scalar=false>
inline __attribute__((always_inline)) void refresh(sr::State& state,const std::array<std::uint8_t,16>& a){
    c2::c8Unroll(std::make_index_sequence<4>{},[&](auto i){
        if constexpr(Scalar){
            alignas(64) tower::Element values[4];_mm512_store_si512(values,state.values[i]);
            for(unsigned j=0;j<4;++j)values[j]=tower::multiplyOracle(values[j],{std::uint64_t(a[4*i+j])|0x100,0});
            state.values[i]=_mm512_load_si512(values);
        }else{
            const auto x=state.values[i];
            state.values[i]=_mm512_xor_si512(tower::multiplyY4(x),sr::byteRefreshProduct<i>(x,a));
        }
    });
}

template<bool Scalar=false>
HD35_NOINLINE void innerLifted(const Plan<4>& p,const RefreshKeys& keys,const Block* planes,Block* output){
    sr::State state{};
    for(unsigned region=0;region<256;++region){
        const auto begin=region*p.rows;const auto* source=planes+region*p.stride;
        for(std::size_t j=0;j<p.rows;j+=64){
            if(HDGROUP_PREFETCH&&region+1<256){
                const auto* next=reinterpret_cast<const char*>(planes+(region+1)*p.stride+j);
                for(unsigned h=0;h<16;++h)_mm_prefetch(next+64*h,_MM_HINT_T1);
            }
            const auto step=(begin+j)/64;
            c2::hdgroup::innerPacket4<true>(source,p.packetRoute.data()+(begin+j)/2,
                output+begin+j,state,NoMixer{},p.stripedRefresh[step]);
            refresh<Scalar>(state,keys.a[step]);
        }
    }
    _mm_sfence();
}

template<bool Scalar=false>
HD35_NOINLINE void encodeLifted(const Plan<4>& p,const RefreshKeys& keys,const Block* input,Block* output,Block* scratch){
    ::outer(p,input,scratch);innerLifted<Scalar>(p,keys,scratch,output);
}
// Each record independently samples two fair bits a,b; L=T^b P^a.
// P rotates byte i to i+1 mod16; T adds byte0 to byte1.
// Separate packed arrays avoid 16-byte struct padding. Setup expands masks.
struct Keys {
    std::vector<std::uint8_t> rotate;
    std::vector<std::uint64_t> shear;
    explicit Keys(std::size_t n,std::uint64_t seed):rotate(n/4),shear(n/4){
        std::mt19937_64 rng(seed); // reproducible benchmark, not ideal setup
        for(std::size_t i=0;i<n/4;++i){
            const auto bits=rng();
            for(unsigned j=0;j<4;++j){
                if((bits>>(2*j))&1)rotate[i]|=std::uint8_t(3U<<(2*j));
                if((bits>>(2*j+1))&1)shear[i]|=std::uint64_t(1)<<(16*j+1);
            }
        }
    }
    std::size_t bytes()const{return rotate.size()+8*shear.size();}
};

inline __attribute__((always_inline)) __m512i apply(__m512i x,std::uint8_t rotate,std::uint64_t shear){
    const auto rotation=_mm512_setr_epi64(
        0x060504030201000fLL,0x0e0d0c0b0a090807LL,
        0x060504030201000fLL,0x0e0d0c0b0a090807LL,
        0x060504030201000fLL,0x0e0d0c0b0a090807LL,
        0x060504030201000fLL,0x0e0d0c0b0a090807LL);
    const auto rotated=_mm512_shuffle_epi8(x,rotation);
    x=_mm512_mask_blend_epi64(rotate,x,rotated);
    const auto extra=_mm512_maskz_shuffle_epi8(shear,x,_mm512_setzero_si512());
    return _mm512_xor_si512(x,extra);
}

static Block scalar(Block x,bool a,bool b){
    std::array<std::uint8_t,16> in,out;std::memcpy(in.data(),&x,16);
    for(unsigned i=0;i<16;++i)out[i]=in[a?(i+15)%16:i];
    if(b)out[1]^=out[0];
    std::memcpy(&x,out.data(),16);return x;
}

template<bool Mix>
inline __attribute__((always_inline)) void store(__m512i* x,unsigned column,Block* planes,std::size_t stride,
    const c2::hdgroup::Coefficient* coeff,const Keys& keys,std::size_t batch){
    c2::hdgroup::toeplitzColumn<4>(x,coeff+column*36);
    c2::c8Unroll(std::make_index_sequence<4>{},[&](auto r){
        auto y=x[r];
        if constexpr(Mix){const auto i=batch*256+4*column+r;y=apply(y,keys.rotate[i],keys.shear[i]);}
        _mm512_stream_si512(reinterpret_cast<__m512i*>(planes+(4*column+r)*stride),y);
    });
}

template<bool Mix,unsigned Offset>
HD35_NOINLINE void coset(const Block* source,Block* planes,std::size_t stride,
    const c2::hdgroup::Coefficient* coeff,const Keys& keys,std::size_t batch){
    alignas(64) __m512i values[128];
    for(unsigned r=0;r<4;++r){
        c2::c8Unroll(std::make_index_sequence<32>{},[&](auto j){values[r*32+j]=_mm512_loadu_si512(source+4*(r*32+j));});
        c2::hdgroup::rsPrefix<32,Offset>(values+r*32);
    }
    for(unsigned j=0;j<32;j+=2){
        alignas(64) __m512i even[4],odd[4];
        const auto c=_mm512_set1_epi8(static_cast<char>(Offset+j));
        c2::c8Unroll(std::make_index_sequence<4>{},[&](auto r){
            const auto a=values[r*32+j],b=values[r*32+j+1];
            even[r]=_mm512_xor_si512(a,_mm512_gf2p8mul_epi8(b,c));
            odd[r]=_mm512_xor_si512(even[r],b);
        });
        store<Mix>(even,Offset+j,planes,stride,coeff,keys,batch);
        store<Mix>(odd,Offset+j+1,planes,stride,coeff,keys,batch);
    }
}

template<bool Mix>
HD35_NOINLINE void outer(const Plan<4>& p,const Keys& keys,const Block* input,Block* scratch){
    for(std::size_t first=0;first<p.rows;first+=64){
        const auto last=std::min(first+64,p.rows);
        for(auto g=first;g<last;g+=4)coset<Mix,0>(input+g*128,scratch+g,p.stride,p.preparedBefore.data()+g*576,keys,g/4);
        for(auto g=first;g<last;g+=4)coset<Mix,32>(input+g*128,scratch+g,p.stride,p.preparedBefore.data()+g*576,keys,g/4);
    }
    _mm_sfence();
    for(unsigned c=0;c<256;++c){auto* row=scratch+c*p.stride;row[p.rows]=row[0];std::memset(row+p.rows+1,0,3*sizeof(Block));}
}

template<bool Mix>
HD35_NOINLINE void encode(const Plan<4>& p,const Keys& keys,const Block* input,Block* output,Block* scratch){
    outer<Mix>(p,keys,input,scratch);::inner<true,4>(p,scratch,output);
}

static void checkLocal(){
    alignas(64) Block input[4],actual[4];
    std::mt19937_64 rng(58);
    for(unsigned trial=0;trial<256;++trial){
        std::uint8_t r=0;std::uint64_t s=0;
        for(unsigned j=0;j<4;++j){input[j]={rng(),rng()};if((trial>>(2*j))&1)r|=3U<<(2*j);if((trial>>(2*j+1))&1)s|=1ULL<<(16*j+1);}
        _mm512_store_si512(actual,apply(_mm512_load_si512(input),r,s));
        for(unsigned j=0;j<4;++j)require(actual[j]==scalar(input[j],(trial>>(2*j))&1,(trial>>(2*j+1))&1),"independent scalar lane oracle");
    }
    // Every binary vector: exact max target atom for the four-map family.
    unsigned maximum=0;
    for(unsigned v=1;v<65536;++v){
        std::uint16_t images[4];
        for(unsigned a=0;a<2;++a)for(unsigned b=0;b<2;++b){
            auto w=std::uint16_t(a?((v<<1)|(v>>15)):v);
            if(b&&(w&1))w^=2;images[2*a+b]=w;
        }
        for(auto w:images){unsigned count=0;for(auto u:images)count+=u==w;maximum=std::max(maximum,count);}
    }
    require(maximum==2,"four-map atom");
    std::cout<<"local_checks,pass,max_atom,1/2\n";
    for(unsigned bit=0;bit<128;++bit)for(unsigned a=0;a<256;++a){
        sr::State fast{},slow{};std::array<std::uint8_t,16> coeff;
        alignas(64) tower::Element basis[4];for(auto& x:basis)x=tower::binaryBasis(bit);
        for(auto& x:fast.values)x=_mm512_load_si512(basis);slow=fast;
        for(unsigned j=0;j<16;++j)coeff[j]=static_cast<std::uint8_t>(a+j);
        refresh(fast,coeff);refresh<true>(slow,coeff);
        for(unsigned j=0;j<4;++j)require(_mm512_cmpeq_epi8_mask(fast.values[j],slow.values[j])==~std::uint64_t(0),"companion plus byte product oracle");
    }
    std::cout<<"field_refresh_checks,pass,basis_vectors,128,byte_coefficients,256\n";
}

static void checkFull(){
    ::check<4,4>(8192);checkLocal();
    for(auto k:{8192U,16384U}){
        Plan<4> p(k,113);Keys keys(p.n,409);
        c2::BenchPayload input(p.k),a(p.scratch()),b(p.scratch()),out(p.n),expected(p.n);
        produce(input.p,p.k);::outer(p,input.p,a.p);outer<false>(p,keys,input.p,b.p);
        require(std::equal(a.p,a.p+p.scratch(),b.p),"unchanged outer clone");
        for(std::size_t batch=0;batch<p.rows/4;++batch)for(unsigned c=0;c<256;++c)for(unsigned j=0;j<4;++j){
            const auto i=batch*256+c;auto& v=a.p[c*p.stride+4*batch+j];
            v=scalar(v,(keys.rotate[i]>>(2*j))&1,(keys.shear[i]>>(16*j+1))&1);
        }
        for(unsigned c=0;c<256;++c)a.p[c*p.stride+p.rows]=a.p[c*p.stride];
        outer<true>(p,keys,input.p,b.p);
        require(std::equal(a.p,a.p+p.scratch(),b.p),"fused versus scalar outer");
        ::inner<true,4>(p,a.p,expected.p);encode<true>(p,keys,input.p,out.p,b.p);
        require(std::equal(out.p,out.p+p.n,expected.p),"whole encode composition");
        std::copy_n(input.p,p.k,out.p);encode<true>(p,keys,out.p,out.p,b.p);
        require(std::equal(out.p,out.p+p.n,expected.p),"in-place encode");
        std::cout<<"full_checks,pass,K,"<<k<<",checksum,"<<hash(out.p,p.n)<<'\n';
        RefreshKeys lift(p.n,501);
        encodeLifted<true>(p,lift,input.p,expected.p,a.p);encodeLifted(p,lift,input.p,out.p,b.p);
        require(std::equal(out.p,out.p+p.n,expected.p),"full lifted scalar oracle");
        std::copy_n(input.p,p.k,out.p);encodeLifted(p,lift,out.p,out.p,b.p);
        require(std::equal(out.p,out.p+p.n,expected.p),"lifted in-place encode");
        std::cout<<"full_lifted_checks,pass,K,"<<k<<",checksum,"<<hash(out.p,p.n)<<'\n';
        if(k==8192){
            for(std::size_t j=0;j<p.k;++j){tower::Element v;std::memcpy(&v,input.p+j,16);v=tower::multiplyOracle(v,tower::towerY);std::memcpy(input.p+j,&v,16);}
            encodeLifted(p,lift,input.p,out.p,b.p);
            for(std::size_t j=0;j<p.n;++j){tower::Element v;std::memcpy(&v,expected.p+j,16);v=tower::multiplyOracle(v,tower::towerY);require(std::memcmp(out.p+j,&v,16)==0,"whole encoder commutes with Y");}
            std::cout<<"full_field_linearity_check,pass,K,"<<k<<'\n';
        }
    }
    // Exercise the degree-growth lemma on rank-one input streams. This is
    // a regression check for the proved statement, not a distance search.
    Plan<4> p(8192,113);RefreshKeys lift(p.n,501);std::mt19937_64 rng(3911);
    alignas(64) Block input[65]{},output[64];std::array<std::uint16_t,32> route;
    for(unsigned j=0;j<32;++j)route[j]=2*j;
    for(unsigned trial=0;trial<4;++trial){
        sr::State state{};
        for(unsigned step=0;step<16;++step){
            for(auto& x:input)x={step?std::uint64_t(rng()&255):0,0};
            if(!step)input[0]={1,0};
            c2::hdgroup::innerPacket4<true>(input,route.data(),output,state,NoMixer{},p.stripedRefresh[16*trial+step]);
            _mm_sfence();
            if(step){
                unsigned weight=0;for(auto& x:output){tower::Element v;std::memcpy(&v,&x,16);weight+=tower::coefficient(v,step)!=0;}
                require(weight>=13,"rank-one emission floor during unreduced degree growth");
            }
            refresh(state,lift.a[16*trial+step]);
            if(step<15){
                alignas(64) tower::Element values[16];std::memcpy(values,&state,256);unsigned top=0;
                for(auto& x:values){top+=tower::coefficient(x,step+1)!=0;for(unsigned j=step+2;j<16;++j)require(!tower::coefficient(x,j),"degree growth upper bound");}
                require(top!=0,"degree growth nonzero leading coefficient");
            }
        }
    }
    std::cout<<"rank_one_degree_growth_checks,pass,trials,4,protected_output_steps,15\n";
}

template<int Mode>
static void benchOne(const Plan<4>& p,const Keys& keys,const RefreshKeys& lift,const Block* input,Block* output,Block* scratch,unsigned calls){
    auto run=[&]{if constexpr(Mode==0)::forward<true,4>(p,input,output,scratch);else if constexpr(Mode==3)encodeLifted(p,lift,input,output,scratch);else encode<Mode==2>(p,keys,input,output,scratch);};
    for(unsigned i=0;i<5;++i)run();
    std::vector<double> times(calls);
    for(auto& t:times){const auto start=Clock::now();run();t=std::chrono::duration<double,std::milli>(Clock::now()-start).count();std::atomic_signal_fence(std::memory_order_seq_cst);}
    std::sort(times.begin(),times.end());
    std::cout<<"timing,mode,"<<Mode<<",K,"<<p.k<<",calls,"<<calls<<",median_ms,"<<times[calls/2]<<",p10_ms,"<<times[calls/10]<<",p90_ms,"<<times[calls*9/10]<<",checksum,"<<hash(output,p.n)<<'\n'<<std::flush;
}
}

int main(int argc,char** argv){try{
    if(argc==2&&std::string(argv[1])=="--check"){lane_trial::checkFull();return 0;}
    const bool refreshOnly=argc==4&&std::string(argv[1])=="--bench-refresh";
    if(argc!=4||(!refreshOnly&&std::string(argv[1])!="--bench"))throw std::invalid_argument("--check | --bench[-refresh] logK calls");
    const unsigned logK=std::stoul(argv[2]),calls=std::stoul(argv[3]);
    if(logK<13||logK>20||calls<11||calls>101)throw std::invalid_argument("bounded benchmark");
    Plan<4> p(std::size_t(1)<<logK,113);lane_trial::Keys keys(p.n,409);
    lane_trial::RefreshKeys lift(p.n,501);
    c2::BenchPayload input(p.k),out(p.n),scratch(p.scratch());produce(input.p,p.k);
    std::cout<<std::setprecision(9)<<"metadata,lane_keys_bytes,"<<keys.bytes()<<",refresh_keys_bytes,"<<lift.bytes()<<",setup_excluded,1,fresh_input_between_calls,0\n";
    if(refreshOnly){
        lane_trial::benchOne<0>(p,keys,lift,input.p,out.p,scratch.p,calls);
        lane_trial::benchOne<3>(p,keys,lift,input.p,out.p,scratch.p,calls);
        lane_trial::benchOne<3>(p,keys,lift,input.p,out.p,scratch.p,calls);
        lane_trial::benchOne<0>(p,keys,lift,input.p,out.p,scratch.p,calls);return 0;
    }
    lane_trial::benchOne<0>(p,keys,lift,input.p,out.p,scratch.p,calls);
    lane_trial::benchOne<1>(p,keys,lift,input.p,out.p,scratch.p,calls);
    lane_trial::benchOne<2>(p,keys,lift,input.p,out.p,scratch.p,calls);
    lane_trial::benchOne<3>(p,keys,lift,input.p,out.p,scratch.p,calls);
    lane_trial::benchOne<3>(p,keys,lift,input.p,out.p,scratch.p,calls);
    lane_trial::benchOne<2>(p,keys,lift,input.p,out.p,scratch.p,calls);
    lane_trial::benchOne<1>(p,keys,lift,input.p,out.p,scratch.p,calls);
    lane_trial::benchOne<0>(p,keys,lift,input.p,out.p,scratch.p,calls);
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
