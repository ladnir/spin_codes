// Outer-only correctness checker, not a benchmark or an inner-kernel harness.
// Link the outer object under test plus the retained PackedMixer/BchAvx512
// objects needed by packed_driver.cpp's renamed (never called) driver.
// Usage: creative-basis SEED [SEED ...]
#define main creative_basis_unused_packed_driver_main
#include "packed_driver.cpp"
#undef main
#include "generated/BchCircuit.h"
#include <array>
#include <cstdint>
#include <cstring>

namespace spin::detail::kernel {
void bchPackedCoeffOriginal(const block*,block*,const std::uint64_t*);
void bchPackedCoeffAligned(const block*,block*,const std::uint64_t*);
void bchPackedCoeffCompact(const block*,block*,const std::uint64_t*);
}

namespace creativebasis {

static void requireEqual(const void* actual,const void* expected,std::size_t bytes,const char* reason) {
    if(std::memcmp(actual,expected,bytes))throw std::runtime_error(reason);
}

static void scalarOuter(const k::block* input,k::block* output,const std::uint64_t* originalCoefficients) {
    alignas(64) std::array<k::block,1024> mixed;
    pd::scalarMix(input,mixed.data(),originalCoefficients);
    // Literal binary generator rows, independent of GFNI matrix packing,
    // byte transposes, output assembly, and the generated BCH circuit.
    for(unsigned lane=0;lane<4;++lane)
        for(unsigned row=0;row<128;++row) {
            k::block value{};
            for(unsigned column=0;column<256;++column)
                if((k::BchRows[row][column/64]>>(column%64))&1U)
                    value^=mixed[4*column+lane];
            output[128*lane+row]=value;
        }
}

template<std::size_t Size>
static void checkGuards(const std::array<k::block,Size>& storage,
                        std::size_t first,std::size_t count,const k::block& canary) {
    for(std::size_t i=0;i<first;++i)
        requireEqual(&storage[i],&canary,sizeof(canary),"leading block canary");
    for(std::size_t i=first+count;i<Size;++i)
        requireEqual(&storage[i],&canary,sizeof(canary),"trailing block canary");
}

template<unsigned Export>
static void checkExport(const k::block* input,k::block* output,const std::uint64_t* coefficients,
                        const std::array<k::block,1024>& originalInput,
                        const std::array<k::block,512>& expected,const k::block& canary) {
    static_assert(Export<3);
    // Reset every output, so a missing store cannot inherit a correct zero.
    std::fill(output,output+512,canary);
    if constexpr(Export==0)k::bchPackedCoeffOriginal(input,output,coefficients);
    else if constexpr(Export==1)k::bchPackedCoeffAligned(input,output,coefficients);
    else k::bchPackedCoeffCompact(input,output,coefficients);
    requireEqual(output,expected.data(),sizeof(expected),
        Export==0?"Original physical basis":Export==1?"Aligned physical basis":"Compact physical basis");
    requireEqual(input,originalInput.data(),sizeof(originalInput),"outer modified input");
}

static void check(std::uint64_t seed) {
    const pd::Gl32 gl(1,seed);
    alignas(64) std::array<std::uint64_t,1024> aligned;
    alignas(64) std::array<std::uint64_t,512> compact;
    // Original promises unaligned coefficient loads; exercise that too.
    alignas(64) std::array<std::uint64_t,1032> originalStorage;
    constexpr std::uint64_t coefficientCanary=0xa5d184391ca47b63ULL;
    originalStorage.fill(coefficientCanary);
    auto* original=originalStorage.data()+1;
    std::copy(gl.coeff.begin(),gl.coeff.end(),aligned.begin());
    std::copy(gl.coeff.begin(),gl.coeff.end(),original);
    for(unsigned i=0;i<512;++i) {
        if(gl.coeff[2*i]!=gl.coeff[2*i+1])throw std::runtime_error("coefficient pair differs");
        compact[i]=gl.coeff[2*i];
    }

    const k::block canary(0xe173542834875129ULL,0x6315872054816739ULL);
    alignas(64) std::array<k::block,1032> inputStorage;
    alignas(64) std::array<k::block,520> outputStorage;
    alignas(64) std::array<k::block,1024> originalInput{};
    alignas(64) std::array<k::block,512> expected,masks;
    inputStorage.fill(canary);outputStorage.fill(canary);
    // Both pointers are valid block-aligned addresses, deliberately 16 mod64.
    auto* input=inputStorage.data()+1;
    auto* output=outputStorage.data()+1;
    if((reinterpret_cast<std::uintptr_t>(input)&63)!=16 ||
       (reinterpret_cast<std::uintptr_t>(output)&63)!=16)
        throw std::runtime_error("unaligned test fixture lost its alignment");
    std::fill(input,input+1024,k::block{});
    const k::block ones(_mm_set1_epi64x(-1)),zero{};
    std::size_t checked=0;
    for(unsigned coordinate=0;coordinate<1024;++coordinate) {
        // The scalar map acts identically and independently on all128 payload
        // bits. Its all-ones image gives this coordinate's binary column mask;
        // AND with each unit below is the literal scalar single-bit answer.
        input[coordinate]=ones;
        scalarOuter(input,masks.data(),gl.coeff.data());
        for(const auto& mask:masks)
            if(std::memcmp(&mask,&zero,sizeof(mask)) && std::memcmp(&mask,&ones,sizeof(mask)))
                throw std::runtime_error("scalar oracle mixed payload positions");
        for(unsigned bit=0;bit<128;++bit) {
            const auto unit=_mm_set_epi64x(bit>=64?std::uint64_t(1)<<(bit-64):0,
                                           bit<64?std::uint64_t(1)<<bit:0);
            input[coordinate]=originalInput[coordinate]=k::block(unit);
            for(unsigned i=0;i<512;++i)expected[i]=k::block(_mm_and_si128(masks[i].mData,unit));
            checkExport<0>(input,output,original,originalInput,expected,canary);
            checkGuards(inputStorage,1,1024,canary);checkGuards(outputStorage,1,512,canary);
            checkExport<1>(input,output,aligned.data(),originalInput,expected,canary);
            checkGuards(inputStorage,1,1024,canary);checkGuards(outputStorage,1,512,canary);
            checkExport<2>(input,output,compact.data(),originalInput,expected,canary);
            checkGuards(inputStorage,1,1024,canary);checkGuards(outputStorage,1,512,canary);
            ++checked;
        }
        input[coordinate]=originalInput[coordinate]=k::block{};
    }
    requireEqual(aligned.data(),gl.coeff.data(),sizeof(aligned),"aligned coefficients modified");
    requireEqual(original,gl.coeff.data(),sizeof(aligned),"original coefficients modified");
    for(unsigned i=0;i<512;++i)
        if(compact[i]!=gl.coeff[2*i])throw std::runtime_error("compact coefficients modified");
    if(originalStorage.front()!=coefficientCanary)throw std::runtime_error("coefficient prefix canary");
    for(unsigned i=1025;i<originalStorage.size();++i)
        if(originalStorage[i]!=coefficientCanary)throw std::runtime_error("coefficient suffix canary");
    if(checked!=131072)throw std::runtime_error("incomplete physical basis");
    std::cout<<"checks PASS outer physical basis="<<checked<<" per export seed="<<seed
             <<" exports=Original,Aligned,Compact total_calls="<<3*checked
             <<"; scalar GL32 + literal BCH rows; unaligned data/output; guards; immutable input/coefficients\n";
}

static std::uint64_t seedArgument(const char* text) {
    if(!*text)throw std::invalid_argument("empty seed");
    for(const char* p=text;*p;++p)
        if(*p<'0' || *p>'9')throw std::invalid_argument("seed must be an unsigned decimal integer");
    return std::stoull(text);
}
} // namespace creativebasis

int main(int argc,char** argv) {try {
    if(argc<2)throw std::invalid_argument("usage: creative-basis SEED [SEED ...]");
    if(!__builtin_cpu_supports("avx512f") || !__builtin_cpu_supports("avx512vl") ||
       !__builtin_cpu_supports("avx512bw") || !__builtin_cpu_supports("gfni"))
        throw std::runtime_error("AVX512 F/VL/BW and GFNI required");
#if defined(__AVX512VBMI__)
    if(!__builtin_cpu_supports("avx512vbmi"))throw std::runtime_error("VBMI build requires VBMI CPU");
#endif
    for(int i=1;i<argc;++i)creativebasis::check(creativebasis::seedArgument(argv[i]));
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
