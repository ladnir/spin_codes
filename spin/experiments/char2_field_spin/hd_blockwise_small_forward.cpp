#include "HdFastOuter.h"
#include "HdBlockwiseSmallInner.h"
#include "HdBlockwiseOuter256.h"
#include "BenchPayload.h"
#include "Punctured8Oracle.h"
#include <atomic>
#include <iomanip>
#include <iostream>

// Isolated research cost control. No distance certificate is claimed here.
// Default: GRS[256,128], t16/s4. Width64/s4 or s8 is retained only as a
// rejected control: the quiet-subcode bound rules out 40% distance there.
// Compile-time GRS outer and small state, one route, direct SPIN recurrence.
// Setup, allocation, input generation, and checksum are outside timed calls.
namespace c2 = spin::research::char2;
namespace bw = c2::hdbwsmall;
#ifndef BW_SMALL_STATE
#define BW_SMALL_STATE 4
#endif
#ifndef BW_OUTER_LENGTH
#define BW_OUTER_LENGTH 256
#endif
#ifndef BW_PACKET_WIDTH
#define BW_PACKET_WIDTH 16
#endif
static constexpr unsigned StateSize = BW_SMALL_STATE;
static constexpr unsigned PacketWidth = BW_PACKET_WIDTH;
static constexpr unsigned OuterLength = BW_OUTER_LENGTH;
static constexpr unsigned OuterDimension = OuterLength / 2;
static_assert(StateSize == 4 || StateSize == 8);
static_assert(PacketWidth == 16 || PacketWidth == 32 || PacketWidth == 64);
static_assert(PacketWidth >= 2 * StateSize);
static_assert(OuterLength == 128 || OuterLength == 256);
using State = bw::State<StateSize>;
using Mixer = bw::Mixer<StateSize>;
namespace out = c2::hdfast::outer_detail;
using Block = c2::Block;
using Clock = std::chrono::steady_clock;
#ifndef BW_PREFETCH_PACKETS
#define BW_PREFETCH_PACKETS 4
#endif
#ifndef BW_PAD_ROWS
#define BW_PAD_ROWS 4
#endif
#ifndef BW_OUTER_NT
#define BW_OUTER_NT 1
#endif
#ifndef BW_STAGE_REGION
#define BW_STAGE_REGION 1
#endif
static_assert(BW_PREFETCH_PACKETS >= 0 && BW_PREFETCH_PACKETS <= 16);
static_assert(BW_PAD_ROWS >= 0 && BW_PAD_ROWS % 4 == 0);
static_assert(BW_STAGE_REGION >= 0 && BW_STAGE_REGION <= 2);
#if defined(_MSC_VER)
#define BW_NOINLINE __declspec(noinline)
#else
#define BW_NOINLINE __attribute__((noinline))
#endif

static void require(bool yes, const char* message) {
    if (!yes) throw std::runtime_error(message);
}
static std::uint16_t mul(std::uint16_t a, std::uint16_t b) {
    return c2::medium32::quadraticMultiply(a, b);
}
struct RawMixer {
    std::array<std::uint16_t, StateSize> u{}, v{};
};
struct Plan {
    std::size_t k, n, rows, stride;
    std::vector<std::uint16_t> route;
    std::vector<RawMixer> raw;
    std::vector<Mixer> mixers;
    explicit Plan(std::size_t size, std::uint64_t seed)
        : k(size), n(2 * size), rows(size / OuterDimension), stride(rows + BW_PAD_ROWS) {
        if (size < OuterDimension * 64 || size % (OuterDimension * 64) || rows > 65536)
            throw std::invalid_argument("research geometry: K multiple of 64*outer_dimension, rows<=65536");
        std::mt19937_64 rr(seed ^ 0x91af823017ccb91aULL);
        std::mt19937_64 mr(seed ^ 0x14b2ec8399500213ULL);
        route.resize(n);
        for (unsigned c = 0; c != OuterLength; ++c) {
            auto* r = route.data() + c * rows;
            std::iota(r, r + rows, std::uint16_t(0));
            for (auto j = rows - 1; j; --j) {
                std::uniform_int_distribution<std::size_t> draw(0, j);
                std::swap(r[j], r[draw(rr)]);
            }
        }
        raw.resize(n / PacketWidth);
        mixers.resize(raw.size());
        for (std::size_t i = 0; i != raw.size(); ++i) {
            auto& m = raw[i];
            bool nonzero;
            do {
                nonzero = false;
                for (auto& x : m.u) { x = std::uint16_t(mr()); nonzero |= x != 0; }
            } while (!nonzero);
            std::uint16_t pairing;
            do {
                pairing = 0;
                for (unsigned j = 0; j != StateSize; ++j) {
                    m.v[j] = std::uint16_t(mr());
                    pairing ^= mul(m.u[j], m.v[j]);
                }
            } while (pairing == 1);
            mixers[i] = bw::prepareMixer<StateSize>(m.u.data(), m.v.data());
        }
    }
    std::size_t planeRecords() const { return OuterLength * stride; }
    std::size_t scratchRecords() const { return planeRecords() + (BW_STAGE_REGION ? rows : 0); }
};

BW_NOINLINE static void outer(const Plan& p, const Block* input, Block* planes) {
    if constexpr (OuterLength == 256) {
        c2::hdbwouter::outerForward256<64, true, BW_OUTER_NT != 0>(p.k, p.stride, input, planes);
    } else {
        for (std::size_t first = 0; first < p.rows; first += 64) {
            const auto last = std::min(first + 64, p.rows);
            for (auto row = first; row != last; row += 4)
                out::coset64<0, true>(input + row * 64, planes + row, p.stride);
            for (auto row = first; row != last; row += 4)
                out::coset64<64, true>(input + row * 64, planes + 64 * p.stride + row, p.stride);
        }
        if (p.stride != p.rows)
            for (unsigned c = 0; c != OuterLength; ++c)
                std::memset(planes + c * p.stride + p.rows, 0, (p.stride - p.rows) * sizeof(Block));
    }
}
BW_NOINLINE static void inner(const Plan& p, const Block* planes, Block* output, Block* regionScratch) {
    State state{};
    for (unsigned region = 0; region != OuterLength; ++region) {
        const auto begin = region * p.rows;
        const Block* source = planes + region * p.stride;
        if constexpr (BW_STAGE_REGION) {
            // Sequential staging is charged to every full call. The reusable
            // region workspace is allocated by the caller, outside the kernel.
            for (std::size_t h = 0; h != p.rows; h += 16) {
                const auto a = _mm512_loadu_si512(source + h);
                const auto b = _mm512_loadu_si512(source + h + 4);
                const auto c = _mm512_loadu_si512(source + h + 8);
                const auto d = _mm512_loadu_si512(source + h + 12);
                _mm512_storeu_si512(regionScratch + h, a);
                _mm512_storeu_si512(regionScratch + h + 4, b);
                _mm512_storeu_si512(regionScratch + h + 8, c);
                _mm512_storeu_si512(regionScratch + h + 12, d);
            }
            source = regionScratch;
        }
        if constexpr (BW_STAGE_REGION == 2) {
            // A separate gather into final output is also charged. The raw
            // packet kernel supports exact alias and shares the same state.
            for (std::size_t j = 0; j != p.rows; j += 4)
                _mm512_storeu_si512(output + begin + j,
                    c2::hdbw::detail::input4<true>(source, p.route.data() + begin + j, 0));
            for (std::size_t j = 0; j != p.rows; j += PacketWidth)
                bw::packetForwardRaw<StateSize, PacketWidth>(output + begin + j, output + begin + j,
                    state, p.mixers[(begin + j) / PacketWidth]);
            continue;
        }
        for (std::size_t j = 0; j != p.rows; j += PacketWidth) {
            if constexpr (BW_PREFETCH_PACKETS != 0) {
                const auto next = j + PacketWidth * BW_PREFETCH_PACKETS;
                if (next < p.rows)
                    for (unsigned h = 0; h != PacketWidth; h += 4)
                        c2::c8Unroll(std::make_index_sequence<4>{}, [&](auto v) {
                            _mm_prefetch(reinterpret_cast<const char*>(
                                source + p.route[begin + next + h + v]), _MM_HINT_T0);
                        });
            }
            bw::packetForwardRouted<StateSize, PacketWidth>(source, p.route.data() + begin + j,
                output + begin + j, state, p.mixers[(begin + j) / PacketWidth]);
        }
    }
}
// Pre-routed arithmetic-only control, with the same packet inputs and mixers.
// Materialization is outside this control's time, never outside a full time.
BW_NOINLINE static void plainInner(const Plan& p, const Block* planes, Block* output) {
    State state{};
    for (std::size_t i = 0; i != p.n; i += PacketWidth)
        bw::packetForwardRaw<StateSize, PacketWidth>(planes + i, output + i, state, p.mixers[i / PacketWidth]);
}
// Routing-only control: no coding operation or distance claim.
BW_NOINLINE static void routeOnly(const Plan& p, const Block* planes, Block* output) {
    for (unsigned c = 0; c != OuterLength; ++c) {
        const auto begin = c * p.rows;
        for (std::size_t j = 0; j != p.rows; j += 4)
            _mm512_storeu_si512(output + begin + j,
                c2::hdbw::detail::input4<true>(planes + c * p.stride, p.route.data() + begin + j, 0));
    }
}
struct Times { double outer = 0, inner = 0; };
template<bool Timed = false>
BW_NOINLINE static void forward(const Plan& p, const Block* input, Block* output,
    Block* scratch, Times* times = nullptr) {
    auto start = Clock::time_point{};
    if constexpr (Timed) start = Clock::now();
    outer(p, input, scratch);
    if constexpr (Timed) {
        const auto now = Clock::now();
        times->outer = std::chrono::duration<double, std::micro>(now - start).count();
        start = now;
    }
    inner(p, scratch, output, scratch + p.planeRecords());
    if constexpr (Timed)
        times->inner = std::chrono::duration<double, std::micro>(Clock::now() - start).count();
}

static void add(Block& a, Block b) { a.lo ^= b.lo; a.hi ^= b.hi; }
static Block scale(Block x, std::uint16_t c) {
    std::uint16_t lanes[8];
    std::memcpy(lanes, &x, sizeof(x));
    for (auto& v : lanes) v = mul(v, c);
    std::memcpy(&x, lanes, sizeof(x));
    return x;
}
static Block byteScale(Block x, std::uint8_t c) {
    const auto& o = c2::p8Oracle();
    auto* bytes = reinterpret_cast<std::uint8_t*>(&x);
    for (unsigned byte = 0; byte != 16; ++byte) bytes[byte] = o.mul[c][bytes[byte]];
    return x;
}
// Independent product definition of the normalized additive-polynomial basis.
static const auto& evaluations() {
    static const auto table = [] {
        std::array<std::array<std::uint8_t, OuterDimension>, OuterLength> e{};
        const auto& o = c2::p8Oracle();
        for (unsigned point = 0; point != OuterLength; ++point) {
            std::uint8_t v[8]{};
            for (unsigned bit = 0; (1U << bit) < OuterDimension; ++bit) {
                std::uint8_t n = 1, d = 1;
                for (unsigned j = 0; j != (1U << bit); ++j) {
                    n = o.mul[n][point ^ j]; d = o.mul[d][(1U << bit) ^ j];
                }
                v[bit] = o.mul[n][o.inverse[d]];
            }
            for (unsigned j = 0; j != OuterDimension; ++j) {
                e[point][j] = 1;
                for (unsigned bit = 0; (1U << bit) < OuterDimension; ++bit)
                    if (j & (1U << bit)) e[point][j] = o.mul[e[point][j]][v[bit]];
            }
        }
        return e;
    }();
    return table;
}
static const auto& emission() {
    static const auto table = [] {
        std::array<std::array<std::uint8_t, StateSize>, PacketWidth> a{};
        const auto& o = c2::p8Oracle();
        // Direct Lagrange products: no constants or FFT helpers from the kernel.
        for (unsigned i = 0; i != PacketWidth; ++i)
            for (unsigned j = 0; j != StateSize; ++j) {
                std::uint8_t n = 1, d = 1;
                for (unsigned h = 0; h != StateSize; ++h) if (h != j) {
                    n = o.mul[n][i ^ (128 + h)];
                    d = o.mul[d][j ^ h];
                }
                a[i][j] = o.mul[n][o.inverse[d]];
            }
        return a;
    }();
    return table;
}
static void scalarPacket(const Block* input, Block* output,
    std::array<Block, StateSize>& state, const RawMixer& m) {
    const auto& a = emission();
    std::array<Block, StateSize> next = state;
    Block dot{};
    for (unsigned j = 0; j != StateSize; ++j) add(dot, scale(state[j], m.v[j]));
    for (unsigned j = 0; j != StateSize; ++j) add(next[j], scale(dot, m.u[j]));
    for (unsigned i = 0; i != PacketWidth; ++i) {
        auto y = input[i];
        for (unsigned j = 0; j != StateSize; ++j) {
            add(y, byteScale(state[j], a[i][j]));
            add(next[j], byteScale(input[i], a[i][j]));
        }
        output[i] = y;
    }
    state = next;
}
static std::vector<Block> oracle(const Plan& p, const Block* input) {
    std::vector<Block> planes(p.planeRecords()), output(p.n);
    const auto& e = evaluations();
    for (std::size_t i = 0; i != p.k; ++i) {
        const auto row = (i / (4 * OuterDimension)) * 4 + i % 4, degree = (i % (4 * OuterDimension)) / 4;
        for (unsigned c = 0; c != OuterLength; ++c)
            add(planes[c * p.stride + row], byteScale(input[i], e[c][degree]));
    }
    std::array<Block, StateSize> state{};
    std::array<Block, PacketWidth> packet{};
    for (std::size_t i = 0; i != p.n; i += PacketWidth) {
        const auto begin = (i / p.rows) * p.stride;
        for (unsigned j = 0; j != PacketWidth; ++j) packet[j] = planes[begin + p.route[i + j]];
        scalarPacket(packet.data(), output.data() + i, state, p.raw[i / PacketWidth]);
    }
    return output;
}
static void checkPackets() {
    const auto& a = emission();
    std::array<Block, StateSize> q{}, injected{}, zero{};
    std::array<Block, PacketWidth> x{}, image{};
    const Block one{0x0101010101010101ULL, 0x0101010101010101ULL};
    for (unsigned j = 0; j != StateSize; ++j) {
        q.fill(Block{}); q[j] = one;
        bw::emit<StateSize, PacketWidth>(q.data(), image.data());
        for (unsigned i = 0; i != PacketWidth; ++i) {
            require(bw::aCoefficient<StateSize, PacketWidth>(i, j) == a[i][j], "independent Lagrange coefficient");
            require(image[i] == byteScale(one, a[i][j]), "emission basis");
        }
        bw::inject<StateSize, PacketWidth>(image.data(), injected.data());
        require(injected == zero, "C A = 0 basis");
    }
    for (unsigned i = 0; i != PacketWidth; ++i) {
        x.fill(Block{}); x[i] = one;
        bw::inject<StateSize, PacketWidth>(x.data(), injected.data());
        for (unsigned j = 0; j != StateSize; ++j)
            require(injected[j] == byteScale(one, a[i][j]), "injection basis");
    }
    Plan p(64 * OuterDimension, 919);
    std::mt19937_64 rng(919);
    for (unsigned trial = 0; trial != 32; ++trial) {
        for (auto& v : x) v = {rng(), rng()};
        for (auto& v : q) v = {rng(), rng()};
        auto scalarState = q;
        auto m = p.raw[trial];
        if (trial == 0) m.v.fill(0); // Valid identity mixer is part of the ensemble.
        const auto mixer = bw::prepareMixer<StateSize>(m.u.data(), m.v.data());
        std::array<Block, PacketWidth> want{};
        scalarPacket(x.data(), want.data(), scalarState, m);
        State state;
        bw::loadState(state, q.data());
        bw::packetForwardRaw<StateSize, PacketWidth>(x.data(), image.data(), state, mixer);
        bw::storeState(injected.data(), state);
        require(image == want && injected == scalarState, "nonzero-state packet");
        bw::loadState(state, q.data());
        auto aliased = x;
        bw::packetForwardRaw<StateSize, PacketWidth>(aliased.data(), aliased.data(), state, mixer);
        bw::storeState(injected.data(), state);
        require(aliased == want && injected == scalarState, "raw packet exact alias");
        std::array<std::uint16_t, PacketWidth> r16{};
        std::array<std::uint32_t, PacketWidth> r32{};
        std::iota(r16.begin(), r16.end(), std::uint16_t(0));
        std::shuffle(r16.begin(), r16.end(), rng);
        std::array<Block, PacketWidth> region{};
        for (unsigned i = 0; i != PacketWidth; ++i) { region[r16[i]] = x[i]; r32[i] = r16[i]; }
        bw::loadState(state, q.data());
        bw::packetForwardRouted<StateSize, PacketWidth>(region.data(), r16.data(), image.data(), state, mixer);
        bw::storeState(injected.data(), state);
        require(image == want && injected == scalarState, "uint16 routed packet");
        bw::loadState(state, q.data());
        bw::packetForwardRouted<StateSize, PacketWidth>(region.data(), r32.data(), image.data(), state, mixer);
        bw::storeState(injected.data(), state);
        require(image == want && injected == scalarState, "uint32 routed packet");
    }
    std::cout << "check,fixed-map-bases-CA0-nonzero-state-identity-mixer-raw-alias,pass\n";
}
static void check() {
    checkPackets();
    std::mt19937_64 rng(331);
    const Block guard{119, 811};
    for (auto k : {64 * OuterDimension, 128 * OuterDimension, 192 * OuterDimension}) {
        Plan p(k, 113);
        std::vector<Block> source(k);
        for (auto& x : source) x = {rng(), rng()};
        const auto want = oracle(p, source.data());
        {
            std::vector<Block> planes(p.planeRecords()), routed(p.n), rawResult(p.n);
            outer(p, source.data(), planes.data());
            routeOnly(p, planes.data(), routed.data());
            plainInner(p, routed.data(), rawResult.data());
            require(rawResult == want, "pre-routed arithmetic control is exact");
        }
        for (unsigned offset = 0; offset != 4; ++offset) {
            std::vector<Block> input(p.n + 8, guard), output(p.n + 8, guard), scratch(p.scratchRecords() + 8, guard);
            std::copy(source.begin(), source.end(), input.begin() + offset);
            forward(p, input.data() + offset, output.data() + offset, scratch.data() + offset);
            require(std::equal(want.begin(), want.end(), output.begin() + offset), "dense scalar forward equality");
            require(std::equal(source.begin(), source.end(), input.begin() + offset), "source preserved");
            Times times;
            forward<true>(p, input.data() + offset, input.data() + offset, scratch.data() + offset, &times);
            require(std::equal(want.begin(), want.end(), input.begin() + offset), "exact input/output alias");
            for (const auto* v : {&input, &output, &scratch}) {
                for (unsigned j = 0; j != offset; ++j) require((*v)[j] == guard, "prefix guard");
                const auto length = v == &scratch ? p.scratchRecords() : p.n;
                for (auto j = offset + length; j != v->size(); ++j) require((*v)[j] == guard, "suffix guard");
            }
        }
        std::cout << "check,K," << k << ",scalar-route-alias-offset-guards,pass\n";
    }
}
static void produce(Block* p, std::size_t n) {
    std::uint64_t x = 913, y = 113;
    for (std::size_t i = 0; i != n; ++i) {
        x ^= x << 13; x ^= x >> 7; x ^= x << 17; y += 0x9e3779b97f4a7c15ULL;
        _mm_stream_si128(reinterpret_cast<__m128i*>(p + i), _mm_set_epi64x(y, x));
    }
    _mm_sfence();
}
static std::uint64_t hash(const Block* p, std::size_t n) {
    std::uint64_t h = 0;
    for (std::size_t i = 0; i != n; ++i) {
        h = (h ^ p[i].lo) * 0x100000001b3ULL; h = (h ^ p[i].hi) * 0x100000001b3ULL;
    }
    return h;
}
static void stats(const char* label, std::vector<double> v) {
    std::sort(v.begin(), v.end());
    std::cout << label << ',' << v[v.size()/2] << ',' << v[v.size()/10] << ',' << v[9*v.size()/10] << '\n';
}
static void bench(unsigned calls, const std::string& mode, unsigned logK) {
    Plan p(std::size_t(1) << logK, 113);
    c2::BenchPayload input(p.k), output(p.planeRecords()), scratch(p.scratchRecords());
    produce(input.p, p.k);
    if (mode == "inner" || mode == "plain-inner" || mode == "route-only") outer(p, input.p, scratch.p);
    if (mode == "plain-inner") {
        routeOnly(p, scratch.p, output.p);
        std::memcpy(scratch.p, output.p, p.n * sizeof(Block));
    }
    auto run = [&](Times* t) {
        if (mode == "outer") outer(p, input.p, output.p);
        else if (mode == "inner") inner(p, scratch.p, output.p, scratch.p + p.planeRecords());
        else if (mode == "plain-inner") plainInner(p, scratch.p, output.p);
        else if (mode == "route-only") routeOnly(p, scratch.p, output.p);
        else if (mode == "phases") forward<true>(p, input.p, output.p, scratch.p, t);
        else forward(p, input.p, output.p, scratch.p);
    };
    for (unsigned j = 0; j != 3; ++j) { Times t; run(&t); }
    std::vector<double> total(calls), ot(calls), it(calls);
    std::uint64_t consume = 0;
    for (unsigned j = 0; j != calls; ++j) {
        produce(input.p, p.k);
        Times t;
        const auto start = Clock::now(); run(&t);
        total[j] = std::chrono::duration<double, std::micro>(Clock::now() - start).count();
        ot[j] = t.outer; it[j] = t.inner;
        std::atomic_signal_fence(std::memory_order_seq_cst);
        consume ^= output.p[(j * 104729U) % p.n].lo;
    }
    std::cout << "profile,blockwise-small,outer," << OuterLength << ",state," << StateSize << ",packet_width," << PacketWidth << ",K," << p.k << ",N," << p.n
        << ",calls," << calls << ",mode," << mode << ",prefetch_packets," << BW_PREFETCH_PACKETS
        << ",hardware_gather," << HDBW_HW_GATHER
        << ",padding_rows," << BW_PAD_ROWS << ",outer_nt," << BW_OUTER_NT
        << ",stage_region," << BW_STAGE_REGION
        << ",checksum," << hash(output.p, p.n)
        << ",consume," << consume << "\nname,median_us,p10_us,p90_us\n"
        << std::fixed << std::setprecision(3);
    stats(mode.c_str(), total);
    if (mode == "phases") { stats("outer", ot); stats("route_inner", it); }
}
int main(int argc, char** argv) { try {
    if (!c2::gfni512Available()) return 77;
    if (argc == 2 && std::string(argv[1]) == "--check") { check(); return 0; }
    if (argc < 2 || std::string(argv[1]) != "--bench")
        throw std::invalid_argument("--check | --bench [calls=31] [full|phases|outer|inner|plain-inner|route-only] [logK=20]");
    const unsigned calls = argc > 2 ? std::stoul(argv[2]) : 31;
    const std::string mode = argc > 3 ? argv[3] : "full";
    const unsigned logK = argc > 4 ? std::stoul(argv[4]) : 20;
    if (!calls || calls > 100000 || logK < 13 || logK > 22 ||
        (mode != "full" && mode != "phases" && mode != "outer" && mode != "inner" &&
         mode != "plain-inner" && mode != "route-only"))
        throw std::invalid_argument("benchmark arguments");
    bench(calls, mode, logK);
    return 0;
} catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; } }
