// Diagnostic only: exact retained PackedCoeff1.o, varying working-set size.
// No certified kernel, setup distribution, or production source is changed.
#define main cost_unused_retained_driver_main
#include "packed_driver.cpp"
#undef main

namespace spin::detail::kernel {
void bchPackedCoeffCompact(const block*,block*,const std::uint64_t*);
unsigned packedCoeffTileMode();
}
namespace costprobe {
struct Blocks {
    std::vector<k::block> storage;
    k::block* data;
    explicit Blocks(std::size_t n):storage(n+4),data(reinterpret_cast<k::block*>(
        (reinterpret_cast<std::uintptr_t>(storage.data())+63)&~std::uintptr_t(63))) {}
};
struct Coefficients {
    std::vector<std::uint64_t> storage;
    std::uint64_t* data;
    explicit Coefficients(const pd::Gl32& gl):storage(gl.coeff.size()/2+8),data(reinterpret_cast<std::uint64_t*>(
        (reinterpret_cast<std::uintptr_t>(storage.data())+63)&~std::uintptr_t(63))) {
        for(std::size_t i=0;i<gl.coeff.size()/2;++i) {
            if(gl.coeff[2*i]!=gl.coeff[2*i+1])throw std::runtime_error("coefficient pair differs");
            data[i]=gl.coeff[2*i];
        }
    }
};
static SPIN_NOINLINE void run(const k::block* input,k::block* output,
    const std::uint64_t* coefficients,unsigned mask) {
    // One actual K20 outer pass always makes 2048 calls. Hot tests change
    // only the number of distinct tiles, never the called kernel or its work.
    for(unsigned tile=0;tile<2048;++tile) {
        const auto slot=tile&mask;
        k::bchPackedCoeffCompact(input+slot*tileStride,output+slot*512,coefficients+slot*512);
    }
}
static std::uint64_t checksum(const k::block* output,unsigned tiles) {
    std::uint64_t result=0;
    for(unsigned tile=0;tile<tiles;++tile) {
        std::uint64_t words[2];std::memcpy(words,output+tile*512+(tile&511),16);
        result=(result*0x9e3779b97f4a7c15ULL)^words[0]^words[1];
    }
    return result;
}
static void experiment(std::uint64_t seed,unsigned calls,unsigned selected) {
    constexpr unsigned tiles=2048;
    Blocks input(tiles*tileStride),output(tiles*512);
    k::workspace_routing::adviseOwned(input.data,tiles*tileStride*16);
    pd::Gl32 gl(tiles,seed);Coefficients coefficients(gl);
    k::setup::Words rng(seed^0xf317ec5aULL);
    for(unsigned i=0;i<tiles*tileStride;++i)input.data[i]=k::block(rng(),rng());
    // Verify the linked optimized routine and compact coefficient layout
    // against scalar GL32 plus the retained independent BCH reference.
    alignas(64) k::block mixed[1024],expected[512];
    for(unsigned tile: {0U,1U,17U,2047U}) {
        pd::scalarMix(input.data+tile*tileStride,mixed,gl.coeff.data()+1024*tile);
        k::bchTranspose4(mixed,expected);
        k::bchPackedCoeffCompact(input.data+tile*tileStride,output.data+tile*512,coefficients.data+512*tile);
        if(std::memcmp(expected,output.data+512*tile,sizeof(expected)))
            throw std::runtime_error("exact-current outer kernel/reference mismatch");
    }
    std::cout<<"checks PASS: exact compact GL32/BCH versus scalar mix+production BCH; tile_mode="
             <<k::packedCoeffTileMode()<<'\n';
    std::cout<<"seed,order,footprint_tiles,logical_footprint_bytes,calls,median_ms,p10_ms,p90_ms,checksum\n";
    constexpr unsigned footprints[]={1,4,16,64,256,1024,2048};
    for(unsigned order=0;order<2;++order)for(unsigned at=0;at<7;++at) {
        const auto count=footprints[order?6-at:at];
        if(selected && count!=selected)continue;
        for(unsigned i=0;i<3;++i)run(input.data,output.data,coefficients.data,count-1);
        std::vector<double> times;times.reserve(calls);
        for(unsigned i=0;i<calls;++i) {
            const auto start=Clock::now();run(input.data,output.data,coefficients.data,count-1);
            times.push_back(std::chrono::duration<double,std::milli>(Clock::now()-start).count());
        }
        std::sort(times.begin(),times.end());
        std::cout<<seed<<','<<order<<','<<count<<','<<std::size_t(count)*(tileStride*16+512*16+512*8)
                 <<','<<calls<<','<<std::fixed<<std::setprecision(6)<<times[calls/2]<<','
                 <<times[calls/10]<<','<<times[9*calls/10]<<','<<std::hex<<checksum(output.data,count)
                 <<std::dec<<'\n';
    }
}
}
int main(int argc,char** argv) {try {
    if(argc<3 || argc>4)throw std::invalid_argument("usage: encoder_cost_probe seed calls [footprint_tiles]");
    const auto seed=std::stoull(argv[1]);const unsigned calls=std::stoul(argv[2]);
    const unsigned selected=argc==4?std::stoul(argv[3]):0;
    if(!calls || (selected && selected!=1 && selected!=4 && selected!=16 && selected!=64
        && selected!=256 && selected!=1024 && selected!=2048))throw std::invalid_argument("unsupported diagnostic");
    if(!__builtin_cpu_supports("avx512f") || !__builtin_cpu_supports("avx512vl") ||
       !__builtin_cpu_supports("avx512bw") || !__builtin_cpu_supports("gfni"))throw std::runtime_error("AVX512/GFNI required");
    costprobe::experiment(seed,calls,selected);
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
