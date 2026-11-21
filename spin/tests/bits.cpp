// Bit-vector transpose: equals every bit lane of the 128-bit record transpose,
// for every backend this CPU runs, natural sizes, Code and PreparedEncoder
// (Full and BankedHeuristic, across reseeds). The call allocates nothing.
#include <spin/Code.h>
#include <spin/PreparedEncoder.h>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <new>
#include <random>
#include <span>
#include <vector>
#ifdef _MSC_VER
#include <malloc.h>
#endif

static std::atomic<std::size_t> allocations = 0;
void* operator new(std::size_t n) { ++allocations; if (auto* p = std::malloc(n ? n : 1)) return p; throw std::bad_alloc(); }
void* operator new[](std::size_t n) { return ::operator new(n); }
void operator delete(void* p) noexcept { std::free(p); }
void operator delete[](void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }
void operator delete[](void* p, std::size_t) noexcept { std::free(p); }
void* operator new(std::size_t n, std::align_val_t a) {
    ++allocations;
#ifdef _MSC_VER
    auto* p = _aligned_malloc(n ? n : 1, std::size_t(a));
#else
    void* p = nullptr; if (posix_memalign(&p, std::size_t(a), n ? n : 1)) p = nullptr;
#endif
    if (p) return p; throw std::bad_alloc();
}
void* operator new[](std::size_t n, std::align_val_t a) { return ::operator new(n, a); }
#ifdef _MSC_VER
void operator delete(void* p, std::align_val_t) noexcept { _aligned_free(p); }
void operator delete[](void* p, std::align_val_t) noexcept { _aligned_free(p); }
void operator delete(void* p, std::size_t, std::align_val_t) noexcept { _aligned_free(p); }
void operator delete[](void* p, std::size_t, std::align_val_t) noexcept { _aligned_free(p); }
#else
void operator delete(void* p, std::align_val_t) noexcept { std::free(p); }
void operator delete[](void* p, std::align_val_t) noexcept { std::free(p); }
void operator delete(void* p, std::size_t, std::align_val_t) noexcept { std::free(p); }
void operator delete[](void* p, std::size_t, std::align_val_t) noexcept { std::free(p); }
#endif

struct alignas(16) Rec { std::uint64_t w[2]; };
static int failures = 0;
static void expect(bool ok, const char* what) { if (!ok) { std::cerr << "FAIL: " << what << '\n'; ++failures; } }

static constexpr unsigned lanes[] = {0, 1, 63, 64, 127};
static std::vector<std::uint64_t> lane(const std::vector<Rec>& r, unsigned j) {
    std::vector<std::uint64_t> v((r.size() + 63) / 64, 0);
    for (std::size_t i = 0; i < r.size(); ++i) v[i / 64] |= ((r[i].w[j / 64] >> (j % 64)) & 1) << (i % 64);
    return v;
}
static std::vector<Rec> randomRecords(std::size_t n, std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::vector<Rec> r(n);
    for (auto& x : r) { x.w[0] = rng(); x.w[1] = rng(); }
    return r;
}

// Encoder is Code or PreparedEncoder; Work its record workspace, BitWork its bit workspace.
template<class Encoder, class Work, class BitWork>
static void compare(const Encoder& code, Work& work, BitWork& bits, std::uint64_t seed, const char* what) {
    const auto n = code.code_size(), k = code.message_size();
    const auto in = randomRecords(n, seed);
    std::vector<Rec> out(k);
    code.template transpose<Rec>(std::span<const Rec>(in), std::span<Rec>(out), work);
    std::vector<std::uint64_t> got(k / 64);
    // Warm call: after a reseed the first call may rebind (Full builds the new plan's tables).
    { const auto x = lane(in, 0); code.transpose_bits(std::span<const std::uint64_t>(x), std::span<std::uint64_t>(got), bits); }
    for (unsigned j : lanes) {
        const auto x = lane(in, j), want = lane(out, j);
        const auto before = allocations.load();
        code.transpose_bits(std::span<const std::uint64_t>(x), std::span<std::uint64_t>(got), bits);
        expect(allocations.load() == before, "transpose_bits allocated");
        expect(got == want, what);
    }
    // Sparse regular-noise input (one bit per 64-bit word block of 4 words).
    std::vector<std::uint64_t> sparse(n / 64, 0);
    std::mt19937_64 rng(seed ^ 77);
    for (std::size_t w = 0; w < sparse.size(); w += 4) sparse[w + rng() % 4] = std::uint64_t(1) << (rng() % 64);
    std::vector<Rec> sin(n, Rec{{0, 0}}), sout(k);
    for (std::size_t i = 0; i < n; ++i) sin[i].w[0] = (sparse[i / 64] >> (i % 64)) & 1;
    code.template transpose<Rec>(std::span<const Rec>(sin), std::span<Rec>(sout), work);
    code.transpose_bits(std::span<const std::uint64_t>(sparse), std::span<std::uint64_t>(got), bits);
    expect(got == lane(sout, 0), what);
}

int main() {
    const auto caps = spin::capabilities();
    std::vector<spin::Backend> backends = {spin::Backend::Automatic, spin::Backend::Portable, spin::Backend::Reference};
    if (caps.avx2) backends.push_back(spin::Backend::Avx2);
    if (caps.avx2gfni) backends.push_back(spin::Backend::Avx2Gfni);
    if (caps.avx512) backends.push_back(spin::Backend::Avx512);
    for (std::size_t k : {128, 256, 384, 1152, 16384, 65536}) {
        for (auto b : backends) {
            if (b == spin::Backend::Reference && k > 16384) continue;
            spin::Code code({.message_size = k, .route_seed = 17, .inner_seed = 41}, {.backend = b});
            auto work = code.make_workspace();
            auto bits = code.make_bit_workspace();
            compare(code, work, bits, k + unsigned(b), "Code bit transpose equals record lanes");
        }
    }
    // Errors and workspace binding.
    {
        spin::Code a({.message_size = 1024}), b({.message_size = 1024, .route_seed = 9});
        auto wa = a.make_bit_workspace();
        std::vector<std::uint64_t> in(a.code_size() / 64), out(a.message_size() / 64), shortOut(1);
        bool threw = false;
        try { b.transpose_bits(in, out, wa); } catch (const std::invalid_argument&) { threw = true; }
        expect(threw, "foreign bit workspace rejected");
        threw = false;
        try { a.transpose_bits(in, shortOut, wa); } catch (const std::invalid_argument&) { threw = true; }
        expect(threw, "bit size mismatch rejected");
        b.prepare_bit_workspace(wa);
        b.transpose_bits(in, out, wa);
    }
    // PreparedEncoder: Full and BankedHeuristic, bit workspace made before reseeds.
    for (auto mode : {spin::SetupMode::Full, spin::SetupMode::BankedHeuristic})
        for (auto b : backends) {
            if (b == spin::Backend::Reference) continue;
            const std::size_t k = 16384;
            spin::PreparedEncoder enc({.message_size = k, .route_seed = 3, .inner_seed = 4}, {.mode = mode}, {.backend = b});
            auto work = enc.make_workspace();
            auto bits = enc.make_bit_workspace();
            compare(enc, work, bits, 11, "Prepared bit transpose equals record lanes");
            enc.setCodeSeed({5, 6});
            compare(enc, work, bits, 12, "Prepared bit transpose after reseed");
            enc.setCodeSeed({7, 8});
            compare(enc, work, bits, 13, "Prepared bit transpose after second reseed");
        }
    if (failures) { std::cerr << failures << " failures\n"; return 1; }
    std::cout << "spin_bits: all checks passed\n";
    return 0;
}
