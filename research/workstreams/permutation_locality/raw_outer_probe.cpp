// Separate raw-polynomial message-basis experiment. Mode0 retains the exact
// systematic map; mode1 has the same code image but DIFFERENT output coordinates.
// Load the authenticated reference first: its nested main-renaming macro must
// not consume the rename around the complete size-probe source below.
#include "PacketReference.h"
#define main raw_outer_unused_size_main
#include "inner_size_probe.cpp"
#undef main
#include <limits>

namespace spin::detail::kernel {
void bchPackedCoeffRawPolynomial(const block*, block*, const std::uint64_t*);
void bchPackedCoeffRawPolynomialScalar(const block*, block*, const std::uint64_t*);
}

namespace rawprobe {
using namespace sizeprobe;

// Both modes call this same compiled inner, not two separately allocated
// copies selected by outer-template register pressure. It is the selected
// StreamVbmi/cached/no-prefetch recipe and allocates nothing.
static SPIN_NOINLINE void commonInner(const k::block* input, k::block* scratch,
    std::size_t n, const Packets<4>& route, const ip::PreparedPacketUpdates16& updates) {
    sizeprobe::Route<true> emit{scratch, route.bases.data()};
    ip::transposeInner<ip::InnerRecipe::StreamVbmi>(input, n, updates.conjugatedRows(), emit);
}

template<bool Raw> static SPIN_NOINLINE void encode(k::block* input, k::block* scratch,
    std::size_t n, const Packets<4>& route, const Coefficients& coeff,
    const ip::PreparedPacketUpdates16& updates, double* phases) {
    const auto start = phases ? Clock::now() : Clock::time_point{};
    commonInner(input, scratch, n, route, updates);
    const auto middle = phases ? Clock::now() : Clock::time_point{};
    for(std::size_t tile = 0; tile < n / 1024; ++tile) {
        if constexpr(Raw)
            k::bchPackedCoeffRawPolynomial(scratch + tile * tileStride,
                input + tile * 512, coeff.compact + tile * 512);
        else
            k::bchPackedCoeffCompact(scratch + tile * tileStride,
                input + tile * 512, coeff.compact + tile * 512);
    }
    if(phases) {
        phases[0] += std::chrono::duration<double,std::milli>(middle - start).count();
        phases[1] += std::chrono::duration<double,std::milli>(Clock::now() - middle).count();
    }
}

// Untimed original-coordinate inner, explicit retained canonical route, and
// scalar outer. The raw oracle shares mathematical rows, not optimized GFNI
// instructions, polynomial graph, prepared packets, or state basis.
template<bool Raw> static void scalarFull(const std::vector<k::block>& input,
    std::vector<k::block>& expected, std::vector<k::block>& routed,
    const Packets<4>& route, const Updates& setup, const pd::Gl32& gl,
    const Coefficients& coeff) {
    const auto n = input.size();
    innerOracle<true>(input.data(), expected.data(), n, setup);
    for(std::size_t i = 0; i < n; ++i)
        routed[canonical(route.inverse[i])] = expected[i];
    alignas(64) k::block mixed[1024];
    for(std::size_t tile = 0; tile < n / 1024; ++tile) {
        if constexpr(Raw)
            k::bchPackedCoeffRawPolynomialScalar(routed.data() + tile * 1024,
                expected.data() + tile * 512, coeff.compact + tile * 512);
        else {
            pd::scalarMix(routed.data() + tile * 1024, mixed, gl.coeff.data() + tile * 1024);
            k::bchTranspose4(mixed, expected.data() + tile * 512);
        }
    }
    std::copy(input.begin() + n / 2, input.end(), expected.begin() + n / 2);
}

template<bool Raw> static void experiment(unsigned exponent, std::uint64_t seed,
                                          unsigned calls, bool profile) {
    // Keep construction and allocation order equal between both modes and
    // the retained size probe, including its unused expanded coefficient owner.
    const unsigned n = 2U << exponent;
    const Packets<4> route(n / 256, seed, true, true);
    const Updates setup(n / 64, seed); const packetprobe::Setup retained(setup);
    const ip::PreparedPacketUpdates16 transformed(setup.reverse);
    equal(transformed.conjugatedRows().data(), retained.packed.data(),
          retained.packed.size() * sizeof(gf::Dense16Row), "prepared basis tables");
    const pd::Gl32 gl(n / 1024, seed); const Coefficients coeff(gl);
    const sizeprobe::ExpandedCoefficients expanded(gl);
    std::vector<k::block> input(n), expected(n), initial(n), guarded((n / 1024) * tileStride + 12);
    auto* scratch = reinterpret_cast<k::block*>((reinterpret_cast<std::uintptr_t>(guarded.data() + 4) + 63) & ~std::uintptr_t(63));
    const auto scratchBlocks = (n / 1024) * tileStride;
    k::workspace_routing::adviseOwned(scratch, scratchBlocks * 16);
    sizeprobe::checkInnerSize<6>(setup, transformed);
    if(!calls && !ip::vbmiStateSelfCheck()) throw std::runtime_error("VBMI exhaustive basis");
    auto selected = [&](double* phases) {
        rawprobe::encode<Raw>(input.data(), scratch, n, route, coeff, transformed, phases);
    };
    const k::block canary(0x781aca29312357bfULL, 0x928fafcac352697bULL);
    for(unsigned pattern = 0; pattern < (calls ? 1U : 4U); ++pattern) {
        fill(input, pattern == 2 ? 997 : 913);
        if(pattern == 1) { std::fill(input.begin(), input.end(), k::block{}); input[n - 3] = k::block(1, 3); }
        if(pattern == 3) {
            std::fill(input.begin(), input.end(), k::block{});
            input[63] = k::block(1, 2); input[64] = k::block(4, 8); input[n - 65] = k::block(16, 32);
        }
        scalarFull<Raw>(input, expected, initial, route, setup, gl, coeff);
        scratch[-1] = canary; scratch[scratchBlocks] = canary;
        selected(nullptr);
        equal(input.data(), expected.data(), n * sizeof(k::block), "complete scalar output and untouched suffix");
        equal(scratch - 1, &canary, sizeof(canary), "leading scratch canary");
        equal(scratch + scratchBlocks, &canary, sizeof(canary), "trailing scratch canary");
    }
    std::cout << "checks PASS mode=" << unsigned(Raw) << " seed=" << seed << " exponent=" << exponent
              << "; original-coordinate inner; explicit canonical route; scalar GL32/BCH; full output/suffix/canaries"
              << (Raw ? "; RAW POLYNOMIAL MESSAGE BASIS, not the original output coordinates" : "; exact retained message basis") << '\n';
    if(!calls) return;
    fill(input, 913); double phases[2]{};
    for(unsigned i = 0; i < 5; ++i) selected(profile ? phases : nullptr);
    std::vector<double> times; times.reserve(calls);
    for(unsigned i = 0; i < calls; ++i) {
        const auto start = Clock::now(); selected(profile ? phases : nullptr);
        times.push_back(std::chrono::duration<double,std::milli>(Clock::now() - start).count());
    }
    std::sort(times.begin(), times.end());
    std::cout << "mode,exponent,seed,calls,median_ms,p10_ms,p90_ms,checksum\n"
              << unsigned(Raw) << ',' << exponent << ',' << seed << ',' << calls << ','
              << std::fixed << std::setprecision(6) << times[calls / 2] << ','
              << times[calls / 10] << ',' << times[9 * calls / 10] << ','
              << std::hex << hash(input) << std::dec << '\n';
    if(profile) std::cout << "phase_means_including_warmups,inner_route," << phases[0] / (calls + 5)
                         << ",packed_mixer_bch," << phases[1] / (calls + 5) << '\n';
}

// Complete physical basis, with the scalar all-ones image reused as a binary
// column mask for all128 payload bits. No timing or dynamic hot-path callbacks.
static SPIN_NOINLINE void basis(std::uint64_t seed) {
    const pd::Gl32 gl(1, seed); const Coefficients coeff(gl);
    const k::block canary(0xe173542834875129ULL, 0x6315872054816739ULL);
    std::vector<k::block> input(1032), actual(520), expected(512), masks(512);
    std::fill(input.begin(), input.end(), k::block{});
    std::fill(actual.begin(), actual.end(), canary);
    auto* source = input.data() + 3;
    auto* output = actual.data() + 3;
    input[2] = input[1027] = canary;
    std::size_t checked = 0;
    for(unsigned coordinate = 0; coordinate < 1024; ++coordinate) {
        source[coordinate] = k::block(_mm_set1_epi64x(-1));
        k::bchPackedCoeffRawPolynomialScalar(source, masks.data(), coeff.compact);
        for(unsigned bit = 0; bit < 128; ++bit) {
            const auto unit = _mm_set_epi64x(bit >= 64 ? std::uint64_t(1) << (bit - 64) : 0,
                                            bit < 64 ? std::uint64_t(1) << bit : 0);
            source[coordinate] = k::block(unit);
            for(unsigned i = 0; i < 512; ++i) expected[i] = k::block(_mm_and_si128(masks[i].mData, unit));
            k::bchPackedCoeffRawPolynomial(source, output, coeff.compact);
            equal(output, expected.data(), 512 * sizeof(k::block), "raw outer complete physical basis");
            equal(source - 1, &canary, sizeof(canary), "basis input prefix");
            equal(source + 1024, &canary, sizeof(canary), "basis input suffix");
            equal(output - 1, &canary, sizeof(canary), "basis output prefix");
            equal(output + 512, &canary, sizeof(canary), "basis output suffix");
            ++checked;
        }
        source[coordinate] = k::block{};
    }
    std::cout << "checks PASS raw outer physical basis=" << checked << " seed=" << seed
              << "; unaligned input/output; guards; different message basis\n";
}
}

int main(int argc, char** argv) { try {
    if(!__builtin_cpu_supports("avx512f") || !__builtin_cpu_supports("avx512vl") ||
       !__builtin_cpu_supports("avx512bw") || !__builtin_cpu_supports("gfni") ||
       !__builtin_cpu_supports("avx512vbmi")) throw std::runtime_error("AVX512/GFNI/VBMI required");
    if(argc == 3 && std::string(argv[1]) == "basis") {
        rawprobe::basis(std::stoull(argv[2])); return 0;
    }
    if(argc < 5 || argc > 6) throw std::invalid_argument("usage: raw-outer mode exponent seed calls [profile]; or basis seed");
    const auto mode = std::stoul(argv[1]), exponent = std::stoul(argv[2]), calls = std::stoul(argv[4]);
    if(mode > 1 || exponent < 14 || exponent > 20 || calls > std::numeric_limits<unsigned>::max())
        throw std::invalid_argument("mode0 exact/mode1 raw; exponent14..20; calls must fit unsigned");
    const auto seed = std::stoull(argv[3]); const bool profile = argc == 6 && std::stoul(argv[5]);
    if(!calls) sizeprobe::checkInnerBoundaries();
    if(mode == 0) rawprobe::experiment<false>(unsigned(exponent), seed, unsigned(calls), profile);
    else rawprobe::experiment<true>(unsigned(exponent), seed, unsigned(calls), profile);
    return 0;
} catch(const std::exception& e) { std::cerr << e.what() << '\n'; return 1; } }
