#include "HdFastOuter.h"
#include "HdBlockwiseRandomizedOuter.h"
#include "HdBlockwiseOuter256.h"
#include "BenchPayload.h"
#include "Punctured8Oracle.h"
#include <atomic>
#include <iomanip>
#include <iostream>

// Isolated randomized-outer beta64 cost prototype, not a distance certificate.
// Setup, allocation, generation, and checksums stay outside timed calls.
namespace c2 = spin::research::char2;
namespace bw = c2::hdbwshared;
namespace randomouter = c2::hdbwrandom;
static constexpr unsigned StateSize=4,PacketWidth=64,OuterLength=256,OuterDimension=128;
using State=bw::State;
using Mixer=bw::Mixer;
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
static_assert(BW_STAGE_REGION == 1, "randomized outer fuses into the sequential staging pass");
#if defined(_MSC_VER)
#define BW_NOINLINE __declspec(noinline)
#else
#define BW_NOINLINE __attribute__((noinline))
#endif

static void require(bool yes, const char* message) {
    if (!yes) throw std::runtime_error(message);
}
static std::uint16_t mul(std::uint16_t a, std::uint16_t b) {
    // Independent schoolbook product, reduced by u^2=u+0x20, using the
    // bit-polynomial AES oracle. No native tower multiplier enters checks.
    const auto& f = c2::p8Oracle().mul;
    const unsigned a0 = a & 255, a1 = a >> 8, b0 = b & 255, b1 = b >> 8;
    const auto high = f[a1][b1];
    return std::uint16_t((f[a0][b0] ^ f[0x20][high]) |
        (unsigned(f[a0][b1] ^ f[a1][b0] ^ high) << 8));
}
struct RawMixer {
    std::array<std::uint16_t, StateSize> u{}, v{};
};
struct Plan {
    std::size_t k, n, rows, stride;
    std::vector<std::uint16_t> route;
    // route32 is populated only by correctness checks, never benchmark setup.
    std::vector<std::uint32_t> route32;
    std::vector<std::uint16_t> diagonal;
    std::vector<c2::wide16::Coefficient> preparedDiagonal;
    std::vector<RawMixer> raw;
    std::vector<Mixer> mixers;
    explicit Plan(std::size_t size, std::uint64_t seed)
        : k(size), n(2 * size), rows(size / OuterDimension), stride(rows + BW_PAD_ROWS) {
        if (size < OuterDimension * 64 || size % (OuterDimension * 64) || rows > 65536)
            throw std::invalid_argument("research geometry: K multiple of 64*outer_dimension, rows<=65536");
        std::mt19937_64 rr(seed ^ 0x91af823017ccb91aULL);
        std::mt19937_64 mr(seed ^ 0x14b2ec8399500213ULL);
        std::mt19937_64 dr(seed ^ 0xe685639417f293b5ULL);
        // Reproducible research setup, with an independent generator stream.
        // One uniform draw from F^* per coordinate; no per-lane randomness.
        std::uniform_int_distribution<unsigned> nonzero(1, 65535);
        diagonal.resize(n);
        preparedDiagonal.resize(n);
        for (std::size_t i = 0; i != n; ++i)
            diagonal[i] = static_cast<std::uint16_t>(nonzero(dr));
        c2::wide16::prepare(diagonal.data(), preparedDiagonal.data(), n);
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
            mixers[i] = bw::prepareMixer(m.u.data(), m.v.data());
        }
    }
    template<bool WideRoute> const auto* routeData() const {
        if constexpr (WideRoute) return route32.data();
        else return route.data();
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
template<bool Fft=true, bool Diagonal=true, bool WideRoute=false>
BW_NOINLINE static void inner(const Plan& p, const Block* planes, Block* output, Block* regionScratch) {
    State state{};
    const auto* route = p.routeData<WideRoute>();
    for (unsigned region = 0; region != OuterLength; ++region) {
        const auto begin = region * p.rows;
        randomouter::stageRegion<Diagonal>(planes + region * p.stride,
            p.preparedDiagonal.data() + begin, regionScratch, p.rows);
        const Block* source = regionScratch;
        for (std::size_t j = 0; j != p.rows; j += PacketWidth) {
            if constexpr (BW_PREFETCH_PACKETS != 0) {
                const auto next = j + PacketWidth * BW_PREFETCH_PACKETS;
                if (next < p.rows)
                    for (unsigned h = 0; h != PacketWidth; h += 4)
                        c2::c8Unroll(std::make_index_sequence<4>{}, [&](auto v) {
                            _mm_prefetch(reinterpret_cast<const char*>(
                                source + route[begin + next + h + v]), _MM_HINT_T0);
                        });
            }
            bw::packetForwardRouted<Fft>(source, route + begin + j,
                output + begin + j, state, p.mixers[(begin + j) / PacketWidth]);
        }
    }
}
// Pre-routed arithmetic-only control, with the same packet inputs and mixers.
// Materialization is outside this control's time, never outside a full time.
template<bool Fft=true>
BW_NOINLINE static void plainInner(const Plan& p, const Block* planes, Block* output) {
    State state{};
    for (std::size_t i = 0; i != p.n; i += PacketWidth)
        bw::packetForwardRaw<Fft>(planes + i, output + i, state, p.mixers[i / PacketWidth]);
}
// Staging/diagonal/routing control: exactly the packet input materialization.
template<bool Diagonal=true, bool WideRoute=false>
BW_NOINLINE static void routeOnly(const Plan& p, const Block* planes, Block* output, Block* regionScratch) {
    const auto* route = p.routeData<WideRoute>();
    for (unsigned c = 0; c != OuterLength; ++c) {
        const auto begin = c * p.rows;
        randomouter::stageRegion<Diagonal>(planes + c * p.stride,
            p.preparedDiagonal.data() + begin, regionScratch, p.rows);
        for (std::size_t j = 0; j != p.rows; j += 4)
            _mm512_storeu_si512(output + begin + j,
                c2::hdbw::detail::input4<true>(regionScratch, route + begin + j, 0));
    }
}
struct Times { double outer = 0, inner = 0; };
template<bool Timed=false, bool Fft=true, bool Diagonal=true, bool WideRoute=false>
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
    inner<Fft,Diagonal,WideRoute>(p, scratch, output, scratch + p.planeRecords());
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
struct DenseMaps {
    std::array<std::array<std::uint8_t,64>,64> f{};
    std::array<std::array<std::uint8_t,4>,64> a{};
    std::array<std::array<std::uint8_t,64>,4> c{};
};
static const DenseMaps& maps() {
    static const auto data=[] {
        DenseMaps result;
        const auto& o=c2::p8Oracle();
        const auto& e=evaluations();
        // Independent scalar elimination of the product-defined evaluation
        // matrix. No native cube constants or packed FFT maps enter this oracle.
        for(unsigned alpha=0;alpha!=64;alpha+=8){
            std::uint8_t work[8][16]{};
            for(unsigned i=0;i!=8;++i){
                for(unsigned j=0;j!=8;++j)work[i][j]=e[alpha^i][j];
                work[i][8+i]=1;
            }
            for(unsigned col=0;col!=8;++col){
                unsigned pivot=col;
                while(pivot!=8&&!work[pivot][col])++pivot;
                require(pivot!=8,"scalar evaluation inverse pivot");
                for(unsigned j=0;j!=16;++j)std::swap(work[col][j],work[pivot][j]);
                const auto scale=o.inverse[work[col][col]];
                for(auto& x:work[col])x=o.mul[x][scale];
                for(unsigned i=0;i!=8;++i)if(i!=col){
                    const auto multiple=work[i][col];
                    for(unsigned j=0;j!=16;++j)work[i][j]^=o.mul[multiple][work[col][j]];
                }
            }
            for(unsigned i=0;i!=8;++i){
                for(unsigned j=0;j!=4;++j){
                    result.a[alpha+i][j]=e[alpha^64^i][j];
                    result.c[j][alpha+i]=work[j+4][8+i];
                }
                for(unsigned j=0;j!=8;++j)
                    for(unsigned degree=0;degree!=8;++degree)
                        result.f[alpha+i][alpha+j]^=o.mul[e[alpha^64^i][degree]][work[degree][8+j]];
            }
        }
        return result;
    }();
    return data;
}
static void scalarPacket(const Block* input,Block* output,
    std::array<Block,StateSize>& state,const RawMixer& mixer) {
    const auto& map=maps();
    auto next=state;
    Block dot{};
    for(unsigned j=0;j!=4;++j)add(dot,scale(state[j],mixer.v[j]));
    for(unsigned j=0;j!=4;++j)add(next[j],scale(dot,mixer.u[j]));
    for(unsigned i=0;i!=64;++i){
        Block y{};
        for(unsigned k=0;k!=64;++k)if(map.f[i][k])
            add(y,byteScale(input[k],map.f[i][k]));
        for(unsigned j=0;j!=4;++j){
            add(y,byteScale(state[j],map.a[i][j]));
            add(next[j],byteScale(input[i],map.c[j][i]));
        }
        output[i]=y;
    }
    state=next;
}
template<bool Diagonal>
static std::vector<Block> oracle(const Plan& p, const Block* input) {
    std::vector<Block> planes(p.planeRecords()), output(p.n);
    const auto& e = evaluations();
    for (std::size_t i = 0; i != p.k; ++i) {
        const auto row = (i / (4 * OuterDimension)) * 4 + i % 4, degree = (i % (4 * OuterDimension)) / 4;
        for (unsigned c = 0; c != OuterLength; ++c)
            add(planes[c * p.stride + row], byteScale(input[i], e[c][degree]));
    }
    if constexpr (Diagonal)
        for (unsigned c = 0; c != OuterLength; ++c)
            for (std::size_t row = 0; row != p.rows; ++row)
                planes[c * p.stride + row] = scale(planes[c * p.stride + row],
                    p.diagonal[c * p.rows + row]);
    std::array<Block, StateSize> state{};
    std::array<Block, PacketWidth> packet{};
    for (std::size_t i = 0; i != p.n; i += PacketWidth) {
        const auto begin = (i / p.rows) * p.stride;
        for (unsigned j = 0; j != PacketWidth; ++j) packet[j] = planes[begin + p.route[i + j]];
        scalarPacket(packet.data(), output.data() + i, state, p.raw[i / PacketWidth]);
    }
    return output;
}
template<bool Fft>
static void checkOnePacket(const std::array<Block,64>& x,const std::array<Block,4>& q,
                          const RawMixer& raw,std::mt19937_64& rng,bool wideRoute) {
    const auto mixer=bw::prepareMixer(raw.u.data(),raw.v.data());
    auto wantState=q;
    std::array<Block,64> want{};
    scalarPacket(x.data(),want.data(),wantState,raw);
    const Block guard{0x91e92c011775ULL,0x38291051aaULL};
    for(unsigned offset=0;offset!=4;++offset){
        std::vector<Block> input(72,guard),output(72,guard);
        std::copy(x.begin(),x.end(),input.begin()+offset);
        State state;
        bw::loadState(state,q.data());
        bw::packetForwardRaw<Fft>(input.data()+offset,output.data()+offset,state,mixer);
        std::array<Block,4> actualState{};
        bw::storeState(actualState.data(),state);
        require(std::equal(want.begin(),want.end(),output.begin()+offset)&&actualState==wantState,
                "independent dense nonzero-state packet");
        require(std::equal(x.begin(),x.end(),input.begin()+offset),"raw source preserved");
        bw::loadState(state,q.data());
        bw::packetForwardRaw<Fft>(input.data()+offset,input.data()+offset,state,mixer);
        bw::storeState(actualState.data(),state);
        require(std::equal(want.begin(),want.end(),input.begin()+offset)&&actualState==wantState,
                "raw exact alias");
        for(const auto* v:{&input,&output}){
            for(unsigned j=0;j!=offset;++j)require((*v)[j]==guard,"packet prefix guard");
            for(unsigned j=offset+64;j!=72;++j)require((*v)[j]==guard,"packet suffix guard");
        }
    }
    std::array<std::uint16_t,64> r16{};
    std::array<std::uint32_t,64> r32{};
    std::iota(r16.begin(),r16.end(),std::uint16_t(0));
    std::shuffle(r16.begin(),r16.end(),rng);
    const unsigned base=wideRoute?65536:0;
    std::vector<Block> region(base+68,guard);
    std::array<Block,64> out16{},out32{};
    for(unsigned i=0;i!=64;++i){region[base+r16[i]]=x[i];r32[i]=base+r16[i];}
    State state;
    std::array<Block,4> actualState{};
    bw::loadState(state,q.data());
    bw::packetForwardRouted<Fft>(region.data()+base,r16.data(),out16.data(),state,mixer);
    bw::storeState(actualState.data(),state);
    require(out16==want&&actualState==wantState,"uint16 routed packet");
    bw::loadState(state,q.data());
    bw::packetForwardRouted<Fft>(region.data(),r32.data(),out32.data(),state,mixer);
    bw::storeState(actualState.data(),state);
    require(out32==want&&actualState==wantState,"uint32 routed packet");
    for(unsigned i=base+64;i!=region.size();++i)require(region[i]==guard,"route read suffix guard");
}
static void checkPackets() {
    std::mt19937_64 rng(919);
    RawMixer identity{};
    identity.u[0]=1;
    // Every input/state coordinate, separately in both extension components.
    for(unsigned component=0;component!=2;++component)
        for(unsigned coordinate=0;coordinate!=68;++coordinate){
            std::array<Block,64> x{};
            std::array<Block,4> q{};
            const std::uint64_t one=0x0001000100010001ULL<<(8*component);
            if(coordinate<64)x[coordinate]={one,one};
            else q[coordinate-64]={one,one};
            checkOnePacket<false>(x,q,identity,rng,false);
            checkOnePacket<true>(x,q,identity,rng,false);
        }
    Plan p(8192,919);
    for(unsigned trial=0;trial!=24;++trial){
        std::array<Block,64> x{};
        std::array<Block,4> q{};
        for(auto& z:x)z={rng(),rng()};
        for(auto& z:q)z={rng(),rng()};
        auto raw=p.raw[trial];
        if(!trial)raw.v.fill(0);
        checkOnePacket<false>(x,q,raw,rng,trial==0);
        checkOnePacket<true>(x,q,raw,rng,trial==0);
    }
    // Direct independent C*A=0 and F^2=I matrix checks.
    const auto& map=maps();
    const auto& o=c2::p8Oracle();
    for(unsigned i=0;i!=4;++i)for(unsigned j=0;j!=4;++j){
        std::uint8_t sum=0;
        for(unsigned k=0;k!=64;++k)sum^=o.mul[map.c[i][k]][map.a[k][j]];
        require(!sum,"independent C A zero");
    }
    for(unsigned i=0;i!=64;++i)for(unsigned j=0;j!=64;++j){
        std::uint8_t sum=0;
        for(unsigned k=0;k!=64;++k)sum^=o.mul[map.f[i][k]][map.f[k][j]];
        require(sum==(i==j),"independent feedthrough involution");
    }
    std::cout<<"check,beta64-cube-FFT-all68bases-both-extension-components-nonzero-state-alias-offset-route16-route32,pass\n";
}
static void checkDiagonal() {
    // Every nonzero scalar, every lane, both tower basis components. Distinct
    // per-record scalars detect coefficient expansion and source-index errors.
    const Block guard{0x91e92c011775ULL, 0x38291051aaULL};
    std::array<c2::wide16::Coefficient, 16> coefficients{};
    std::array<std::uint16_t, 16> scalars{};
    std::array<Block, 16> source{}, want{};
    for (unsigned first = 1; first <= 65535; first += 16) {
        for (unsigned j = 0; j != 16; ++j) {
            scalars[j] = static_cast<std::uint16_t>(1 + (first - 1 + j) % 65535);
            coefficients[j] = c2::wide16::prepareCoefficient(scalars[j]);
        }
        for (unsigned component = 0; component != 2; ++component) {
            const std::uint64_t basis = 0x0001000100010001ULL << (8 * component);
            for (unsigned j = 0; j != 16; ++j) {
                source[j] = {basis, basis};
                want[j] = scale(source[j], scalars[j]);
            }
            std::array<Block, 16> output{};
            randomouter::stageRegion<true>(source.data(), coefficients.data(), output.data(), 16);
            require(output == want, "every nonzero diagonal scalar, both tower components");
        }
    }
    std::mt19937_64 rng(215);
    std::array<c2::wide16::Coefficient, 20> offsetCoefficients{};
    for (unsigned j = 0; j != 16; ++j) {
        source[j] = {rng(), rng()};
        scalars[j] = std::uint16_t(1 + rng() % 65535);
        want[j] = scale(source[j], scalars[j]);
    }
    for (unsigned offset = 0; offset != 4; ++offset) {
        std::vector<Block> input(24, guard), output(24, guard);
        std::copy(source.begin(), source.end(), input.begin() + offset);
        c2::wide16::prepare(scalars.data(), offsetCoefficients.data() + offset, 16);
        randomouter::stageRegion<true>(input.data() + offset,
            offsetCoefficients.data() + offset, output.data() + offset, 16);
        require(std::equal(want.begin(), want.end(), output.begin() + offset), "diagonal payload/coefficient offsets");
        require(std::equal(source.begin(), source.end(), input.begin() + offset), "diagonal source preserved");
        randomouter::stageRegion<false>(input.data() + offset, nullptr, output.data() + offset, 16);
        require(std::equal(source.begin(), source.end(), output.begin() + offset), "identity copy control");
        randomouter::stageRegion<true>(input.data() + offset,
            offsetCoefficients.data() + offset, input.data() + offset, 16);
        require(std::equal(want.begin(), want.end(), input.begin() + offset), "diagonal exact alias");
        for (const auto* values : {&input, &output}) {
            for (unsigned j = 0; j != offset; ++j) require((*values)[j] == guard, "diagonal prefix guard");
            for (unsigned j = offset + 16; j != values->size(); ++j) require((*values)[j] == guard, "diagonal suffix guard");
        }
    }
    std::cout << "check,all65535-nonzero-GF65536-scalars-eight-lanes-both-components-alias-offset-guards,pass\n";
}
template<bool Diagonal, bool WideRoute>
static void checkFull(const Plan& p, const std::vector<Block>& source, const std::vector<Block>& want) {
    const Block guard{119, 811};
    {
        std::vector<Block> planes(p.scratchRecords()), routed(p.n), rawResult(p.n);
        outer(p, source.data(), planes.data());
        routeOnly<Diagonal,WideRoute>(p, planes.data(), routed.data(), planes.data() + p.planeRecords());
        // Independently verify the route-conjugated diagonal at EVERY coordinate.
        const auto* route = p.routeData<WideRoute>();
        for (std::size_t i = 0; i != p.n; ++i) {
            const auto row = route[i];
            const auto c = i / p.rows;
            const auto x = planes[c * p.stride + row];
            const auto expected = Diagonal ? scale(x, p.diagonal[c * p.rows + row]) : x;
            require(routed[i] == expected, "route-conjugated outer diagonal index");
        }
        plainInner(p, routed.data(), rawResult.data());
        require(rawResult == want, "pre-routed arithmetic control is exact");
    }
    for (unsigned offset = 0; offset != 4; ++offset) {
        std::vector<Block> input(p.n + 8, guard), output(p.n + 8, guard), scratch(p.scratchRecords() + 8, guard);
        std::copy(source.begin(), source.end(), input.begin() + offset);
        forward<false,true,Diagonal,WideRoute>(p, input.data() + offset, output.data() + offset, scratch.data() + offset);
        require(std::equal(want.begin(), want.end(), output.begin() + offset), "independent FFT end-to-end scalar equality");
        require(std::equal(source.begin(), source.end(), input.begin() + offset), "source preserved");
        forward<false,false,Diagonal,WideRoute>(p, input.data() + offset, output.data() + offset, scratch.data() + offset);
        require(std::equal(want.begin(), want.end(), output.begin() + offset), "cube same-map control");
        forward<false,true,Diagonal,WideRoute>(p, input.data() + offset, input.data() + offset, scratch.data() + offset);
        require(std::equal(want.begin(), want.end(), input.begin() + offset), "FFT exact input/output alias");
        std::copy(source.begin(), source.end(), input.begin() + offset);
        forward<false,false,Diagonal,WideRoute>(p, input.data() + offset, input.data() + offset, scratch.data() + offset);
        require(std::equal(want.begin(), want.end(), input.begin() + offset), "cube exact input/output alias");
        for (const auto* v : {&input, &output, &scratch}) {
            for (unsigned j = 0; j != offset; ++j) require((*v)[j] == guard, "prefix guard");
            const auto length = v == &scratch ? p.scratchRecords() : p.n;
            for (auto j = offset + length; j != v->size(); ++j) require((*v)[j] == guard, "suffix guard");
        }
    }
}
static void check() {
    checkDiagonal();
    checkPackets();
    std::mt19937_64 rng(331);
    for (auto k : {64 * OuterDimension, 128 * OuterDimension, 192 * OuterDimension}) {
        Plan p(k, 113);
        p.route32.assign(p.route.begin(), p.route.end());
        require(std::all_of(p.diagonal.begin(), p.diagonal.end(), [](auto c) { return c != 0; }),
            "diagonal sampler excludes zero");
        std::vector<Block> source(k);
        for (auto& x : source) x = {rng(), rng()};
        const auto randomized = oracle<true>(p, source.data());
        const auto identity = oracle<false>(p, source.data());
        require(randomized != identity, "diagonal enabled changes the map");
        checkFull<true,false>(p, source, randomized);
        checkFull<true,true>(p, source, randomized);
        checkFull<false,false>(p, source, identity);
        checkFull<false,true>(p, source, identity);
        std::cout << "check,K," << k
            << ",randomized-and-identity-scalar-FFT-cube-route16-route32-alias-offset-guards,pass\n";
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
template<bool Fft, bool Diagonal>
static void bench(unsigned calls, const std::string& mode, unsigned logK) {
    Plan p(std::size_t(1) << logK, 113);
    c2::BenchPayload input(p.k), output(p.planeRecords()), scratch(p.scratchRecords());
    produce(input.p, p.k);
    if (mode == "inner" || mode == "plain-inner" || mode == "route-only") outer(p, input.p, scratch.p);
    if (mode == "plain-inner") {
        routeOnly<Diagonal>(p, scratch.p, output.p, scratch.p + p.planeRecords());
        std::memcpy(scratch.p, output.p, p.n * sizeof(Block));
    }
    auto run = [&](Times* t) {
        if (mode == "outer") outer(p, input.p, output.p);
        else if (mode == "inner") inner<Fft,Diagonal>(p, scratch.p, output.p, scratch.p + p.planeRecords());
        else if (mode == "plain-inner") plainInner<Fft>(p, scratch.p, output.p);
        else if (mode == "route-only") routeOnly<Diagonal>(p, scratch.p, output.p, scratch.p + p.planeRecords());
        else if (mode == "phases") forward<true,Fft,Diagonal>(p, input.p, output.p, scratch.p, t);
        else forward<false,Fft,Diagonal>(p, input.p, output.p, scratch.p);
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
    std::cout << "profile,blockwise-randomized-outer-beta64,outer," << OuterLength << ",state," << StateSize << ",packet_width," << PacketWidth << ",K," << p.k << ",N," << p.n
        << ",calls," << calls << ",mode," << mode << ",prefetch_packets," << BW_PREFETCH_PACKETS
        << ",hardware_gather," << HDBW_HW_GATHER
        << ",padding_rows," << BW_PAD_ROWS << ",outer_nt," << BW_OUTER_NT
        << ",stage_region," << BW_STAGE_REGION << ",algorithm," << (Fft?"fft":"cube")
        << ",diagonal," << (Diagonal ? "GF65536" : "identity")
        << ",diagonal_placement,outer-output-before-route,diagonal_phase,route_inner"
        << ",diagonal_scalars," << p.n << ",field_lanes_per_record,8"
        << ",checksum," << hash(output.p, p.n)
        << ",consume," << consume << "\nname,median_us,p10_us,p90_us\n"
        << std::fixed << std::setprecision(3);
    stats(mode.c_str(), total);
    if (mode == "phases") { stats("outer", ot); stats("route_inner", it); }
}
template<bool Fft>
static void selectDiagonal(unsigned calls, const std::string& mode, unsigned logK, const std::string& diagonal) {
    // Sequential calls by construction: never run two benchmarks at once.
    if (diagonal == "identity" || diagonal == "both") bench<Fft,false>(calls, mode, logK);
    if (diagonal == "diagonal" || diagonal == "both") bench<Fft,true>(calls, mode, logK);
}
int main(int argc, char** argv) { try {
    if (!c2::gfni512Available()) return 77;
    if (argc == 2 && std::string(argv[1]) == "--check") { check(); return 0; }
    if (argc < 2 || std::string(argv[1]) != "--bench")
        throw std::invalid_argument("--check | --bench [calls=31] [full|phases|outer|inner|plain-inner|route-only] [logK=20] [fft|cube] [diagonal|identity|both]");
    const unsigned calls = argc > 2 ? std::stoul(argv[2]) : 31;
    const std::string mode = argc > 3 ? argv[3] : "full";
    const unsigned logK = argc > 4 ? std::stoul(argv[4]) : 20;
    if (!calls || calls > 100000 || logK < 13 || logK > 22 ||
        (mode != "full" && mode != "phases" && mode != "outer" && mode != "inner" &&
         mode != "plain-inner" && mode != "route-only"))
        throw std::invalid_argument("benchmark arguments");
    const std::string algorithm=argc>5?argv[5]:"fft";
    const std::string diagonal=argc>6?argv[6]:"diagonal";
    if (argc > 7 || (diagonal != "diagonal" && diagonal != "identity" && diagonal != "both"))
        throw std::invalid_argument("diagonal must be diagonal, identity or both");
    if(algorithm=="cube")selectDiagonal<false>(calls,mode,logK,diagonal);
    else if(algorithm=="fft")selectDiagonal<true>(calls,mode,logK,diagonal);
    else throw std::invalid_argument("algorithm must be cube or fft");
    return 0;
} catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; } }
