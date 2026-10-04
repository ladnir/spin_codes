#include "K18Certified22.h"
#include "kernel/K16NativeOuter.h"
#include "outer_variants/OuterVariants.h"
#include "K20Variants.h"
#include "../../src/packet/PacketPlan.h"
#include <spin/Code.h>
#include <spin/PacketCode.h>
#include <algorithm>
#include <chrono>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace rs = spin::research::rs;
namespace small = spin::research::k18cert22;
namespace outer = spin::research::k16codesign;
namespace tuned = spin::research::k20small;
using Clock = std::chrono::steady_clock;
using Block = rs::Block;

static Block* blocks(spin::Buffer& buffer) {
    return reinterpret_cast<Block*>(buffer.bytes().data());
}
static void fill(spin::Buffer& buffer) {
    std::uint64_t state = 913;
    for (std::size_t i = 0; i < buffer.bytes().size(); i += 8) {
        state += 0x9e3779b97f4a7c15ULL; auto v = state;
        v = (v ^ (v >> 30)) * 0xbf58476d1ce4e5b9ULL;
        v = (v ^ (v >> 27)) * 0x94d049bb133111ebULL; v ^= v >> 31;
        std::memcpy(buffer.bytes().data() + i, &v, 8);
    }
}

// Both constructions use the same in-place pointer, warmups, and clock loop.
// The callable is a concrete template parameter: no erased hot-path dispatch.
template<class Encode>
static void measure(const char* mode, std::size_t k, std::uint64_t seed,
                    unsigned calls, bool huge, std::size_t scratchBytes,
                    spin::Buffer& buffer, Encode&& encode) {
    fill(buffer);
    auto* x = blocks(buffer) + 4; // Same 64-byte offset for both kernels.
    for (unsigned i = 0; i < 5; ++i) encode(x);
    std::vector<double> times(calls);
    for (auto& time : times) {
        const auto begin = Clock::now();
        encode(x);
        time = std::chrono::duration<double, std::micro>(Clock::now()-begin).count();
    }
    std::sort(times.begin(), times.end());
    std::uint64_t hash = 0;
    for (std::size_t i = 0; i < k*16; i += 8) {
        std::uint64_t v; std::memcpy(&v, reinterpret_cast<const std::byte*>(x)+i, 8);
        hash = (hash ^ v) * 0x100000001b3ULL;
    }
    std::cout << "mode,K,seed,calls,median_us,p10_us,p90_us,checksum,huge_requested,input_allocation_bytes,scratch_allocation_bytes\n"
        << mode << ',' << k << ',' << seed << ',' << calls << ',' << std::fixed << std::setprecision(4)
        << times[calls/2] << ',' << times[calls/10] << ',' << times[9*calls/10] << ','
        << std::hex << hash << std::dec << ',' << huge << ','
        << buffer.allocation_bytes() << ',' << scratchBytes << '\n';
}

template<class Route, class Finish>
static void measurePhases(const char* mode, std::size_t k, std::uint64_t seed,
                          unsigned calls, spin::Buffer& input, Route&& route, Finish&& finish) {
    fill(input); auto* x = blocks(input) + 4;
    for (unsigned i = 0; i < 5; ++i) { route(x); finish(x); }
    std::vector<double> total(calls), routeTimes(calls), outerTimes(calls);
    for (unsigned i = 0; i < calls; ++i) {
        const auto begin = Clock::now();
        route(x);
        const auto middle = Clock::now();
        finish(x);
        const auto end = Clock::now();
        total[i] = std::chrono::duration<double, std::micro>(end-begin).count();
        routeTimes[i] = std::chrono::duration<double, std::micro>(middle-begin).count();
        outerTimes[i] = std::chrono::duration<double, std::micro>(end-middle).count();
    }
    std::sort(total.begin(), total.end());
    std::sort(routeTimes.begin(), routeTimes.end());
    std::sort(outerTimes.begin(), outerTimes.end());
    std::cout << "mode,K,seed,calls,total_us,route_us,outer_us\n" << mode << ','
        << k << ',' << seed << ',' << calls << ',' << std::fixed << std::setprecision(4)
        << total[calls/2] << ',' << routeTimes[calls/2] << ',' << outerTimes[calls/2] << '\n';
}

template<unsigned RouteNt, unsigned OuterNt>
static void runSmall(const char* mode, std::size_t k, std::uint64_t seed, unsigned calls, bool huge, bool phases) {
    rs::Plan plan(k, seed, seed, rs::Variant::Rs16Gf16, rs::InnerKernel::RetainedStreaming);
    small::Tables tables; small::prepare(plan, seed, tables);
    outer::NativeOuterTables native; outer::customizeNativeField16(plan, seed, native);
    const auto policy = huge ? spin::MemoryPolicy::PreferHugePages : spin::MemoryPolicy::Normal;
    spin::Buffer input((plan.n+12)*16, policy), scratch(plan.scratchBlocks()*16, policy);
    if (phases) {
        return measurePhases(mode, k, seed, calls, input,
            [&](Block* x) { tuned::route<RouteNt>(x, blocks(scratch), plan, tables); },
            [&](Block* x) { tuned::finish<OuterNt>(blocks(scratch), x, plan, native); });
    }
    measure(mode, k, seed, calls, huge, scratch.allocation_bytes(), input, [&](Block* x) {
        tuned::route<RouteNt>(x, blocks(scratch), plan, tables);
        tuned::finish<OuterNt>(blocks(scratch), x, plan, native);
    });
}

static void runPublic(std::size_t k, std::uint64_t seed, unsigned calls, bool huge) {
    const auto policy = huge ? spin::MemoryPolicy::PreferHugePages : spin::MemoryPolicy::Normal;
    spin::Code code({k, spin::Parameters::PacketRsT64S22, seed, seed}, {spin::Backend::Avx512});
    auto work = code.make_workspace(spin::Width::Bits128, policy);
    spin::Buffer input((code.code_size()+12)*16, policy);
    measure("public_s22", k, seed, calls, huge, work.bytes(), input, [&](Block* x) {
        code.transpose_inplace_bytes(std::span<std::byte>(reinterpret_cast<std::byte*>(x),code.code_size()*16), work);
    });
}

static void runLarge(std::size_t k, std::uint64_t seed, unsigned calls, bool huge, bool phases) {
    const auto policy = huge ? spin::MemoryPolicy::PreferHugePages : spin::MemoryPolicy::Normal;
    if (phases) {
        namespace packet = spin::detail::packet;
        packet::Plan plan(k, seed, seed);
        spin::Buffer input((plan.n+12)*16, policy), scratch(plan.scratchBlocks()*16, policy);
        // Diagnostic only: same private stages selected by the public API at K20.
        // Keep the public API call for all headline measurements.
        if (k <= 262144)
            return measurePhases("large_s20_phases", k, seed, calls, input,
                [&](Block* x) { packet::reverseRouteCached(x, blocks(scratch), plan); },
                [&](Block* x) { packet::outerFastCached(blocks(scratch), x, plan); });
        return measurePhases("large_s20_phases", k, seed, calls, input,
            [&](Block* x) { packet::reverseRoute(x, blocks(scratch), plan); },
            [&](Block* x) { packet::outerFast(blocks(scratch), x, plan); });
    }
    spin::Code code({k, spin::Parameters::PacketRsT64S20, seed, seed}, {spin::Backend::Avx512});
    auto work = code.make_workspace(spin::Width::Bits128, policy);
    spin::Buffer input((code.code_size()+12)*16, policy);
    measure("large_s20", k, seed, calls, huge, work.bytes(), input, [&](Block* x) {
        code.transpose_inplace_bytes(std::span<std::byte>(reinterpret_cast<std::byte*>(x),code.code_size()*16), work);
    });
}

int main(int argc, char** argv) { try {
    if (!spin::packet_fast_available() || !spin::capabilities().forward512) return 77;
    if (argc < 3 || argc > 7)
        throw std::invalid_argument("small|small-packed|small-packed-unroll|small-unroll|small-nt-route|small-nt-outer|small-nt|small-nt-split|small-nt-unroll|small-nt-packed|small-nt-packed-unroll|public|large K [seed=41] [calls=1001] [huge=1] [phases=0]");
    const std::string mode = argv[1];
    const auto k = std::stoull(argv[2]), seed = argc > 3 ? std::stoull(argv[3]) : 41;
    const auto calls = argc > 4 ? std::stoul(argv[4]) : 1001;
    const bool huge = argc <= 5 || std::stoul(argv[5]) != 0;
    const bool phases = argc > 6 && std::stoul(argv[6]) != 0;
    if (!calls || calls > 1000000) throw std::invalid_argument("calls must be in 1..1000000");
    if (mode == "small") runSmall<false,false>("small_s22", k, seed, unsigned(calls), huge, phases);
    else if (mode == "small-packed") runSmall<3,0>("small_packed_cached", k, seed, unsigned(calls), huge, phases);
    else if (mode == "small-packed-unroll") runSmall<3,4>("small_packed_cached_unroll", k, seed, unsigned(calls), huge, phases);
    else if (mode == "small-unroll") runSmall<0,4>("small_cached_unroll", k, seed, unsigned(calls), huge, phases);
    else if (mode == "small-nt-route") runSmall<true,false>("small_nt_route", k, seed, unsigned(calls), huge, phases);
    else if (mode == "small-nt-outer") runSmall<false,true>("small_nt_outer", k, seed, unsigned(calls), huge, phases);
    else if (mode == "small-nt") runSmall<true,true>("small_nt", k, seed, unsigned(calls), huge, phases);
    else if (mode == "small-nt-split") runSmall<1,2>("small_nt_split", k, seed, unsigned(calls), huge, phases);
    else if (mode == "small-nt-unroll") runSmall<1,3>("small_nt_unroll", k, seed, unsigned(calls), huge, phases);
    else if (mode == "small-nt-packed") runSmall<2,1>("small_nt_packed", k, seed, unsigned(calls), huge, phases);
    else if (mode == "small-nt-packed-unroll") runSmall<2,3>("small_nt_packed_unroll", k, seed, unsigned(calls), huge, phases);
    else if (mode == "public") {
        if (phases) throw std::invalid_argument("public mode measures full API calls only");
        runPublic(k, seed, unsigned(calls), huge);
    }
    else if (mode == "large") {
        runLarge(k, seed, unsigned(calls), huge, phases);
    }
    else throw std::invalid_argument("unknown mode");
    return 0;
} catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; } }
