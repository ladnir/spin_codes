#include "K18Certified22.h"
#include "kernel/K16NativeOuter.h"
#include "outer_variants/OuterVariants.h"
#include <spin/Code.h>
#include <spin/PacketCode.h>
#include <algorithm>
#include <bit>
#include <chrono>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace rs = spin::research::rs;
namespace outer = spin::research::k16codesign;
namespace inner = spin::research::k18cert22;
namespace spin::research::k16codesign {
void outerScalar(const rs::Block*, rs::Block*, const rs::Plan&);
void outerForwardScalar(const rs::Block*, rs::Block*, const rs::Plan&);
}
using Block = rs::Block;
using Clock = std::chrono::steady_clock;
static Block* blocks(spin::Buffer& b) { return reinterpret_cast<Block*>(b.bytes().data()); }
static void fill(spin::Buffer& b, std::uint64_t state = 913) {
    for (std::size_t i = 0; i < b.bytes().size(); i += 8) {
        state += 0x9e3779b97f4a7c15ULL; auto v = state;
        v = (v ^ (v >> 30)) * 0xbf58476d1ce4e5b9ULL;
        v = (v ^ (v >> 27)) * 0x94d049bb133111ebULL; v ^= v >> 31;
        std::memcpy(b.bytes().data() + i, &v, 8);
    }
}
static unsigned dot(const Block* a, const Block* b, std::size_t n) {
    unsigned p = 0;
    for (std::size_t i = 0; i < n; ++i) {
        std::uint64_t x[2], y[2]; std::memcpy(x, a + i, 16); std::memcpy(y, b + i, 16);
        p ^= std::popcount(x[0] & y[0]) ^ std::popcount(x[1] & y[1]);
    }
    return p & 1;
}
static void require(bool condition, const char* reason) {
    if (!condition) throw std::runtime_error(reason);
}
static void run(std::size_t k, std::uint64_t seed, unsigned calls, bool phases, bool huge) {
    rs::Plan plan(k, seed, seed, rs::Variant::Rs16Gf16, rs::InnerKernel::RetainedStreaming);
    inner::Tables tables; inner::prepare(plan, seed, tables);
    outer::NativeOuterTables native; outer::customizeNativeField16(plan, seed, native);
    const auto policy = huge ? spin::MemoryPolicy::PreferHugePages : spin::MemoryPolicy::Normal;
    spin::Buffer input((plan.n + 12) * 16, policy), output((plan.k + 12) * 16, policy);
    spin::Buffer scratch(plan.scratchBlocks() * 16, policy), scalarScratch(plan.scratchBlocks() * 16);
    spin::Buffer expected(plan.k * 16), message(plan.k * 16), encoded(plan.n * 16);
    fill(input); fill(message, 1949);
    const std::vector<std::byte> saved(input.bytes().begin(), input.bytes().end());
    const auto route = [&](const Block* x) { inner::reverseRouteFast(x, blocks(scratch), plan, tables); };
    const auto finish = [&](Block* y) {
        outer::outervariants::fieldLoopSharedParity(blocks(scratch), y, plan, native);
    };
    const auto encode = [&](const Block* x, Block* y) { route(x); finish(y); };

    // Compare the routed data as well as the final output: the outer has a
    // kernel and cannot by itself detect every possible inner implementation bug.
    std::fill(scratch.bytes().begin(), scratch.bytes().end(), std::byte{0x3c});
    std::fill(scalarScratch.bytes().begin(), scalarScratch.bytes().end(), std::byte{0x3c});
    inner::reverseRouteScalar(blocks(input) + 4, blocks(scalarScratch), plan, tables);
    route(blocks(input) + 4);
    require(!std::memcmp(scratch.bytes().data(), scalarScratch.bytes().data(), scratch.bytes().size()),
            "routed scalar mismatch (including padding guards)");
    outer::outerScalar(blocks(scalarScratch), blocks(expected), plan);
    for (unsigned offset = 0; offset < 4; ++offset) {
        std::fill(output.bytes().begin(), output.bytes().end(), std::byte{0xa5});
        std::memcpy(blocks(input) + 4 + offset, saved.data() + 64, plan.n * 16);
        auto* out = blocks(output) + 4 + offset;
        encode(blocks(input) + 4 + offset, out);
        require(!std::memcmp(out, blocks(expected), plan.k * 16), "scalar mismatch");
        const auto guard = [](std::byte x) { return x == std::byte{0xa5}; };
        require(std::all_of(output.bytes().begin(), output.bytes().begin() + (4 + offset) * 16, guard)
            && std::all_of(output.bytes().begin() + (4 + offset + plan.k) * 16, output.bytes().end(), guard),
            "output guard mismatch");
        require(!std::memcmp(blocks(input) + 4 + offset, saved.data() + 64, plan.n * 16), "input modified");
    }
    std::memcpy(input.bytes().data(), saved.data(), saved.size());
    outer::outerForwardScalar(blocks(message), blocks(scalarScratch), plan);
    inner::forwardInnerScalar(blocks(scalarScratch), blocks(encoded), plan, tables);
    require(dot(blocks(encoded), blocks(input) + 4, plan.n) == dot(blocks(message), blocks(expected), plan.k),
            "adjoint mismatch");
    encode(blocks(input) + 4, blocks(input) + 4);
    require(!std::memcmp(blocks(input) + 4, blocks(expected), plan.k * 16)
        && !std::memcmp(blocks(input) + 4 + plan.k, saved.data() + 64 + plan.k * 16, (plan.n - plan.k) * 16),
        "in-place mismatch");
    if (!calls) { std::cout << "PASS certified22 K=" << k << " seed=" << seed << '\n'; return; }

    // Same measurement protocol as the frozen s15 control: precomputed setup,
    // five warmups, repeated in-place transpose calls, no mid-stage clock by default.
    fill(input); auto* x = blocks(input) + 4;
    for (unsigned i = 0; i < 5; ++i) encode(x, x);
    std::vector<double> total(calls), routeTimes(calls), outerTimes(calls);
    for (unsigned i = 0; i < calls; ++i) {
        const auto a = Clock::now(); route(x);
        const auto b = phases ? Clock::now() : a; finish(x);
        const auto c = Clock::now();
        total[i] = std::chrono::duration<double, std::micro>(c - a).count();
        if (phases) {
            routeTimes[i] = std::chrono::duration<double, std::micro>(b - a).count();
            outerTimes[i] = std::chrono::duration<double, std::micro>(c - b).count();
        }
    }
    std::sort(total.begin(), total.end()); std::sort(routeTimes.begin(), routeTimes.end());
    std::sort(outerTimes.begin(), outerTimes.end());
    std::uint64_t hash = 0;
    for (std::size_t i = 0; i < plan.k * 16; i += 8) {
        std::uint64_t v; std::memcpy(&v, reinterpret_cast<std::byte*>(x) + i, 8);
        hash = (hash ^ v) * 0x100000001b3ULL;
    }
    std::cout << "mode,K,seed,calls,median_us,p10_us,p90_us,inner_us,outer_us,checksum,huge\n"
        << "certified22," << k << ',' << seed << ',' << calls << ',' << std::fixed << std::setprecision(4)
        << total[calls / 2] << ',' << total[calls / 10] << ',' << total[9 * calls / 10] << ','
        << routeTimes[calls / 2] << ',' << outerTimes[calls / 2] << ',' << std::hex << hash << std::dec << ',' << huge << '\n';
}
int main(int argc, char** argv) { try {
    if (!spin::packet_fast_available() || !spin::capabilities().forward512) return 77;
    if (argc < 2 || argc > 6) throw std::invalid_argument("K [seed=1] [calls=501;0=check] [phases=0] [huge=0]");
    const auto k = std::stoull(argv[1]), seed = argc > 2 ? std::stoull(argv[2]) : 1;
    const auto calls = argc > 3 ? std::stoul(argv[3]) : 501;
    if (calls > 1000000) throw std::invalid_argument("too many calls");
    run(k, seed, unsigned(calls), argc > 4 && std::stoul(argv[4]), argc > 5 && std::stoul(argv[5]));
    return 0;
} catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; } }
