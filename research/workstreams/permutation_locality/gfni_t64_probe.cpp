// Isolated t64/s16 uniform-GL16 encoder. No production or t128 source changes.
#define main t64_unused_retained_driver_main
#include "packed_driver.cpp"
#undef main
#include "GfniT64.h"
#include "generated/BchCircuit.h"
namespace spin::detail::kernel {
void bchPackedCoeffCompact(const block*,block*,const std::uint64_t*);
unsigned packedCoeffTileMode();
}
namespace t64probe {
namespace gf=spin::research::gfni_r4;
using Fixed=spin::research::t64::Fixed;
using Matrix=std::array<unsigned,16>;
static_assert(Fixed::T==64 && Fixed::S==16);
static_assert(Fixed::columns==k::Map64S16::columns);
static_assert(Fixed::groupOrder==k::Map64S16::groupOrder);
static_assert(Fixed::groupedColumns==k::Map64S16::groupedColumns);
static_assert(Fixed::feedbackColumns==Fixed::columns);

static unsigned rank(Matrix rows) {
    unsigned count=0;
    for(unsigned c=0;c<16;++c) {
        unsigned pivot=count;while(pivot<16 && !(rows[pivot]&(1U<<c)))++pivot;
        if(pivot==16)continue;
        std::swap(rows[pivot],rows[count]);
        for(unsigned j=count+1;j<16;++j)if(rows[j]&(1U<<c))rows[j]^=rows[count];
        ++count;
    }
    return count;
}
static Matrix transpose(const Matrix& rows) {
    Matrix result{};
    for(unsigned i=0;i<16;++i)for(unsigned j=0;j<16;++j)result[j]|=((rows[i]>>j)&1U)<<i;
    return result;
}
struct Updates {
    std::vector<Matrix> forward,reverse;
    std::vector<gf::Dense16Row> packed;
    std::size_t draws=0;
    Updates(std::size_t epochs,std::uint64_t seed):forward(epochs),reverse(epochs),packed(epochs) {
        k::setup::Words rng(seed^0x3f625a92ULL);
        for(std::size_t epoch=0;epoch<epochs;++epoch) {
            // Sixteen independent uniform row masks, rejected unless full
            // rank: exactly uniform GL16 under ideal independent Words.
            // Sample one matrix at EVERY 64-symbol physical step.
            do {++draws;for(auto& row:forward[epoch])row=unsigned(rng())&65535U;}
            while(rank(forward[epoch])!=16);
            reverse[epoch]=transpose(forward[epoch]);
            packed[epoch]=gf::makeDense16(reverse[epoch].data());
        }
    }
};
struct Coefficients {
    std::vector<std::uint64_t> storage;
    std::uint64_t* compact;
    explicit Coefficients(const pd::Gl32& gl):storage(gl.coeff.size()/2+8),compact(reinterpret_cast<std::uint64_t*>(
        (reinterpret_cast<std::uintptr_t>(storage.data())+63)&~std::uintptr_t(63))) {
        for(std::size_t i=0;i<gl.coeff.size()/2;++i) {
            if(gl.coeff[2*i]!=gl.coeff[2*i+1])throw std::runtime_error("GL32 coefficient pair differs");
            compact[i]=gl.coeff[2*i];
        }
    }
};
static void equal(const void* a,const void* b,std::size_t size,const char* why) {
    if(std::memcmp(a,b,size))throw std::runtime_error(why);
}
template<bool Transpose> static void updateOracle(__m128i* state,const Updates& setup,std::size_t epoch) {
    alignas(32) __m128i next[16];
    // Only the original accepted matrix is read, never its transpose table
    // or the packed GFNI coefficients used by the measured path.
    for(unsigned j=0;j<16;++j) {
        auto value=_mm_setzero_si128();
        for(unsigned i=0;i<16;++i) {
            const bool bit=Transpose?((setup.forward[epoch][i]>>j)&1):((setup.forward[epoch][j]>>i)&1);
            if(bit)value=_mm_xor_si128(value,state[i]);
        }
        next[j]=value;
    }
    std::memcpy(state,next,sizeof(next));
}
template<bool Reverse> static void innerOracle(
    const k::block* input,k::block* out,std::size_t n,const Updates& setup) {
    alignas(32) __m128i state[16]{},syndrome[16];
    for(std::size_t at=0;at<n/64;++at) {
        const auto epoch=Reverse?n/64-1-at:at,base=64*epoch;
        std::fill(std::begin(syndrome),std::end(syndrome),_mm_setzero_si128());
        for(unsigned p=0;p<64;++p) {
            auto value=input[base+p].mData;
            for(unsigned j=0;j<16;++j) {
                // A is symmetric to C by declaration, not by an assumed
                // shortcut on running state. Feedback always uses raw x.
                const bool a=(Fixed::expansionRows[j]>>p)&1;
                const bool c=(Fixed::feedbackColumns[p]>>j)&1;
                if(Reverse?c:a)value=_mm_xor_si128(value,state[j]);
                if(Reverse?a:c)syndrome[j]=_mm_xor_si128(syndrome[j],input[base+p].mData);
            }
            out[base+p]=k::block(value);
        }
        updateOracle<Reverse>(state,setup,epoch);
        for(unsigned j=0;j<16;++j)state[j]=_mm_xor_si128(state[j],syndrome[j]);
    }
}
template<class Emit> SPIN_FORCEINLINE void reverseInner(
    const k::block* input,std::size_t n,const gf::Dense16Row* rows,Emit&& emit) {
    // Non-coroutine SIMD scratch, all fixed-size. Packed state, unpacked
    // words and packed feedback are separate owned objects (no aliasing).
    alignas(32) __m128i words[16],values[64],syndrome[16],table[4][16];
    gf::Packed<16> state{},feedback;
    for(std::size_t epoch=n/64;epoch-->0;) {
        const auto base=64*epoch;
        if(epoch+1==n/64) {
            for(unsigned p=64;p-->0;){values[p]=input[base+p].mData;emit(base+p,input[base+p]);}
        } else {
            gf::unpack(state,words);
            Fixed::emissionTable(words,table);
            Fixed::emitReverse(input+base,values,table,base,emit,std::make_index_sequence<64>{});
        }
        if(!epoch)break; // No flush and no unused state update after output.
        k::zeta<64>(values);Fixed::finish(values,syndrome);
        if(epoch+1==n/64)gf::pack<16>(syndrome,state);
        else {gf::pack<16>(syndrome,feedback);gf::denseStep16<true>(state,rows[epoch],&feedback);}
    }
}
static SPIN_NOINLINE void encode(k::block* input,k::block* scratch,std::size_t n,
    const Packets<4>& packets,const Coefficients& coeff,const gf::Dense16Row* rows,double* phases) {
    const auto start=phases?Clock::now():Clock::time_point{};
    pd::Emitter<false> emit{scratch,packets.bases.data(),packets.controls.data()};
    reverseInner(input,n,rows,emit);_mm_sfence();
    const auto middle=phases?Clock::now():Clock::time_point{};
    for(std::size_t tile=0;tile<n/1024;++tile)
        k::bchPackedCoeffCompact(scratch+tile*tileStride,input+tile*512,coeff.compact+512*tile);
    if(phases) {
        phases[0]+=std::chrono::duration<double,std::milli>(middle-start).count();
        phases[1]+=std::chrono::duration<double,std::milli>(Clock::now()-middle).count();
    }
}
static std::size_t canonical(unsigned x) {return (x/1024)*1024+4*(x%256)+(x/256)%4;}
// Independent forward outer: direct BCH matrix, scalar GL32 transpose and
// inverse routing. No optimized BCH or GFNI map participates in this oracle.
static void outerForward(const k::block* message,k::block* canonicalOut,k::block* out,std::size_t n,
    const pd::Gl32& gl,const Packets<4>& packets) {
    alignas(64) k::block encoded[1024];
    for(std::size_t tile=0;tile<n/1024;++tile) {
        for(unsigned c=0;c<256;++c)for(unsigned lane=0;lane<4;++lane) {
            auto value=_mm_setzero_si128();
            for(unsigned j=0;j<128;++j)if((k::BchRows[j][c/64]>>(c%64))&1)
                value=_mm_xor_si128(value,message[tile*512+128*lane+j].mData);
            encoded[4*c+lane]=k::block(value);
        }
        auto* mixed=canonicalOut+1024*tile;std::fill(mixed,mixed+1024,k::block{});
        const auto* coeff=gl.coeff.data()+1024*tile;
        for(unsigned g=0;g<32;++g)for(unsigned row=0;row<32;++row) {
            const unsigned lane=row/8,j=row%8;unsigned mask=0;
            for(unsigned d=0;d<4;++d)mask|=unsigned((coeff[32*g+8*d+2*lane]>>(8*(7-j)))&255)<<(8*((lane+d)%4));
            for(unsigned c=0;c<32;++c)if((mask>>c)&1)
                mixed[4*(8*g+c%8)+c/8]^=encoded[4*(8*g+j)+lane];
        }
    }
    for(std::size_t i=0;i<n;++i)out[i]=canonicalOut[canonical(packets.inverse[i])];
}
static void adjoint(const k::block* x,const k::block* ey,std::size_t n,
    const k::block* y,const k::block* etx,std::size_t ksize,const char* why) {
    auto left=_mm_setzero_si128(),right=left;
    for(std::size_t i=0;i<n;++i)left=_mm_xor_si128(left,_mm_and_si128(x[i].mData,ey[i].mData));
    for(std::size_t i=0;i<ksize;++i)right=_mm_xor_si128(right,_mm_and_si128(y[i].mData,etx[i].mData));
    equal(&left,&right,16,why);
}
static void checkMaps() {
    for(unsigned p=0;p<64;++p) {
        alignas(32) __m128i values[64]{},syndrome[16],retained[16];
        values[p]=k::block(0x593a456ad635120eULL,0x8c7629f3564321a8ULL).mData;
        const auto unit=values[p];k::zeta<64>(values);
        Fixed::finish(values,syndrome);k::Map64S16::finish(values,retained);
        equal(syndrome,retained,sizeof(syndrome),"t64 copied/retained finisher mismatch");
        for(unsigned j=0;j<16;++j) {
            const auto expected=((Fixed::expansionRows[j]>>p)&1)?unit:_mm_setzero_si128();
            equal(syndrome+j,&expected,16,"t64 A/zeta basis mismatch");
            if(((Fixed::feedbackColumns[p]>>j)&1)!=((Fixed::expansionRows[j]>>p)&1))
                throw std::runtime_error("t64 C=A^T declaration mismatch");
        }
    }
    alignas(32) __m128i state[16],raw[64],table[4][16];
    alignas(32) k::block input[64],actual[64];
    for(unsigned j=0;j<16;++j)state[j]=k::block(1ULL<<j,1ULL<<(63-j)).mData;
    for(unsigned p=0;p<64;++p)input[p]=k::block(0xc09385391df38a76ULL^p,0x129cadbec00fee39ULL+p);
    Fixed::emissionTable(state,table);
    unsigned next=64;
    auto emit=[&](std::size_t p,k::block value) {
        if(!next || p!=--next)throw std::runtime_error("t64 emitter order mismatch");
        actual[p]=value;
    };
    Fixed::emitReverse(input,raw,table,0,emit,std::make_index_sequence<64>{});
    if(next)throw std::runtime_error("t64 emitter omitted coordinates");
    for(unsigned p=0;p<64;++p) {
        auto expected=input[p].mData;
        for(unsigned j=0;j<16;++j)if((Fixed::expansionRows[j]>>p)&1)expected=_mm_xor_si128(expected,state[j]);
        equal(actual+p,&expected,16,"t64 grouped-emission/state-basis mismatch");
        equal(raw+p,input+p,16,"t64 emission did not retain raw feedback input");
    }
}
static void checkSetup(const Updates& setup) {
    for(std::size_t epoch=0;epoch<setup.packed.size();++epoch) {
        if(rank(setup.forward[epoch])!=16)throw std::runtime_error("singular t64 state update");
        alignas(32) __m128i basis[16],expected[16],actual[16];
        for(unsigned j=0;j<16;++j)basis[j]=k::block(1ULL<<j,1ULL<<(63-j)).mData;
        for(unsigned t=0;t<2;++t) {
            std::memcpy(expected,basis,sizeof(basis));
            if(t)updateOracle<true>(expected,setup,epoch);else updateOracle<false>(expected,setup,epoch);
            gf::Packed<16> state;gf::pack<16>(basis,state);
            gf::denseStep16(state,t?setup.packed[epoch]:gf::makeDense16(setup.forward[epoch].data()));
            gf::unpack(state,actual);
            equal(actual,expected,sizeof(actual),"t64 every-epoch forward/transpose update basis mismatch");
        }
    }
}
static void experiment(unsigned exponent,std::uint64_t seed,unsigned calls,bool profile) {
    const unsigned n=2U<<exponent;
    const Packets<4> packets(n/256,seed,true,true);
    const Updates setup(n/64,seed);
    // Match the retained S16 driver's principal allocations and outer setup.
    std::vector<k::block> input(n),inner(n),actual(n),other(n),forward(n),expected(n),materialized(n);
    const auto scratchBlocks=(n/1024)*tileStride;std::vector<k::block> storage(scratchBlocks+4);
    auto* scratch=reinterpret_cast<k::block*>((reinterpret_cast<std::uintptr_t>(storage.data())+63)&~std::uintptr_t(63));
    k::workspace_routing::adviseOwned(scratch,scratchBlocks*16);
    const pd::Gl32 gl(n/1024,seed);const Coefficients coeff(gl);
    std::vector<k::block> outerInput(n),fullForward(n);
    checkMaps();checkSetup(setup);
    fill(other,971);innerOracle<false>(other.data(),forward.data(),n,setup);
    outerForward(other.data(),materialized.data(),outerInput.data(),n,gl,packets);
    innerOracle<false>(outerInput.data(),fullForward.data(),n,setup);
    auto selected=[&](double* phases){encode(input.data(),scratch,n,packets,coeff,setup.packed.data(),phases);};
    for(unsigned pattern=0;pattern<(calls?1U:4U);++pattern) {
        fill(input,pattern==2?997:913);
        if(pattern==1){std::fill(input.begin(),input.end(),k::block{});input[n-3]=k::block(1,3);}
        if(pattern==3) {
            std::fill(input.begin(),input.end(),k::block{});
            input[63]=k::block(1,2);input[64]=k::block(4,8);input[n-65]=k::block(16,32);
        }
        expected=input;innerOracle<true>(input.data(),inner.data(),n,setup);
        reverseInner(input.data(),n,setup.packed.data(),[&](std::size_t i,k::block value){actual[i]=value;});
        equal(actual.data(),inner.data(),n*16,"t64 inner/reference mismatch");
        adjoint(input.data(),forward.data(),n,other.data(),actual.data(),n,"t64 inner adjoint mismatch");
        for(std::size_t i=0;i<n;++i)materialized[canonical(packets.inverse[i])]=inner[i];
        alignas(64) k::block mixed[1024];
        for(std::size_t tile=0;tile<n/1024;++tile) {
            pd::scalarMix(materialized.data()+1024*tile,mixed,gl.coeff.data()+1024*tile);
            k::bchTranspose4(mixed,expected.data()+512*tile);
        }
        adjoint(input.data(),fullForward.data(),n,other.data(),expected.data(),n/2,"t64 complete encoder adjoint mismatch");
        selected(nullptr);equal(input.data(),expected.data(),n*16,"t64 complete encoder/reference/suffix mismatch");
    }
    std::cout<<"checks PASS: explicit A64 zeta basis, C=A^T, grouped emission/raw feedback, "
        "every64-step update/transpose basis, independent original updates, inner forward/reverse, "
        "scalar outer/direct forward BCH, inner+FULL adjoint, complete output/suffix; "
        <<"step=64; state=16; refresh=uniform-gl16; map_sha256="<<Fixed::mapSha256
        <<"; physical_epochs="<<n/64<<"; hot_update_bytes="<<setup.packed.size()*sizeof(setup.packed[0])
        <<"; sampled_matrix_draws="<<setup.draws<<"; tile_mode="<<k::packedCoeffTileMode()<<'\n';
    if(!calls)return;
    fill(input,913);double phases[2]{};auto run=[&]{selected(profile?phases:nullptr);};
    for(unsigned i=0;i<3;++i)run();std::vector<double> times;times.reserve(calls);
    for(unsigned i=0;i<calls;++i) {
        const auto start=Clock::now();run();times.push_back(std::chrono::duration<double,std::milli>(Clock::now()-start).count());
    }
    std::sort(times.begin(),times.end());
    std::cout<<"exponent,step,state,refresh,seed,calls,median_ms,p10_ms,p90_ms,checksum,tile_mode\n"
        <<exponent<<",64,16,uniform-gl16,"<<seed<<','<<calls<<','
        <<std::fixed<<std::setprecision(6)<<times[calls/2]<<','<<times[calls/10]<<','<<times[9*calls/10]<<','
        <<std::hex<<hash(input)<<std::dec<<','<<k::packedCoeffTileMode()<<'\n';
    if(profile)std::cout<<"phase_means_including_warmups,inner_route,"<<phases[0]/(calls+3)
        <<",packed_mixer_bch,"<<phases[1]/(calls+3)<<'\n';
}
}
int main(int argc,char** argv) {try {
    if(argc<4 || argc>5)throw std::invalid_argument(
        "usage: gfni_t64_probe exponent seed calls [profile]; calls=0 checks only; selected A64/C=A^T, fresh uniformGL16 every64-symbol step");
    if(!__builtin_cpu_supports("avx512f") || !__builtin_cpu_supports("avx512vl") ||
       !__builtin_cpu_supports("avx512bw") || !__builtin_cpu_supports("gfni"))throw std::runtime_error("AVX512/GFNI required");
    const unsigned exponent=std::stoul(argv[1]),calls=std::stoul(argv[3]);
    const auto seed=std::stoull(argv[2]);const bool profile=argc==5 && std::stoul(argv[4]);
    if(exponent<14 || exponent>20)throw std::invalid_argument("unsupported t64 exponent");
    t64probe::experiment(exponent,seed,calls,profile);
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
