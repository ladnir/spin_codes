#include <spin/Code.h>
#include <immintrin.h>
#include <algorithm>
#include <charconv>
#include <chrono>
#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <span>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>
#if defined(__linux__) && defined(__x86_64__)
#include <sys/mman.h>
#endif

// One process measures one case. Run benchmark processes serially.
// Each timed batch contains only forward_bytes calls over prebuilt spans.
// Allocation, setup, cache preparation, reference checks and checksums are
// outside timing. Plan/workspace contents are never flushed by this harness.
// The encode phase includes the encoder's own streaming-store fence. Cached
// writeback can remain deferred, as in the real encoding phase; no consumer
// is included. "Produced" models only NT-written input, not all of FLOCK.
namespace {
using Clock = std::chrono::steady_clock;
constexpr std::size_t k = 65536;
constexpr std::size_t cacheLine = 64;
constexpr unsigned warmupBatches = 3;
constexpr std::uint64_t inputSeed = 913;
constexpr std::byte guardValue{0xa5};
constexpr std::string_view usage =
    "usage: spin_forward_batch_bench imt|packet|paired15 128|256 hot|g24|g92 "
    "reused|cold|produced cached|streaming ROUTE_SEED INNER_SEED ITERATIONS normal|huge|automatic|flock\n"
    "       spin_forward_batch_bench --self-test";

enum class CacheMode { Reused, Cold, Produced };
struct Config {
    std::string_view profile, shape, cacheName, storesName, memoryName;
    spin::Parameters parameters;
    spin::Width width;
    CacheMode cache;
    spin::OutputStores stores;
    spin::MemoryPolicy memory;
    unsigned bits = 0;
    std::uint64_t routeSeed = 0, innerSeed = 0;
    std::size_t iterations = 0, calls = 0;
    bool flockAdvice = false;
};

std::uint64_t number(std::string_view text, const char* name) {
    std::uint64_t value = 0;
    const auto result = std::from_chars(text.data(), text.data()+text.size(), value, 10);
    if(text.empty() || result.ec != std::errc{} || result.ptr != text.data()+text.size())
        throw std::invalid_argument(std::string("invalid ")+name+": expected unsigned decimal integer");
    return value;
}
Config parse(int argc, char** argv) {
    if(argc != 10) throw std::invalid_argument(std::string(usage));
    Config c{};
    c.profile = argv[1];
    if(c.profile == "imt") c.parameters = spin::Parameters::T128S19;
    else if(c.profile == "packet") c.parameters = spin::Parameters::PacketRsT64S20;
    else if(c.profile == "paired15") c.parameters = spin::Parameters::PacketRsT64S15K16;
    else throw std::invalid_argument("profile must be imt, packet, or paired15");
    const auto bits = number(argv[2], "width");
    if(bits != 128 && bits != 256) throw std::invalid_argument("width must be 128 or 256");
    c.bits = unsigned(bits);
    c.width = static_cast<spin::Width>(c.bits/8);
    c.shape = argv[3];
    if(c.shape == "hot") c.calls = 1;
    else if(c.shape == "g24") c.calls = 24/(c.bits/128);
    else if(c.shape == "g92") c.calls = 92/(c.bits/128);
    else throw std::invalid_argument("shape must be hot, g24, or g92");
    c.cacheName = argv[4];
    if(c.cacheName == "reused") c.cache = CacheMode::Reused;
    else if(c.cacheName == "cold") c.cache = CacheMode::Cold;
    else if(c.cacheName == "produced") c.cache = CacheMode::Produced;
    else throw std::invalid_argument("cache mode must be reused, cold, or produced");
    c.storesName = argv[5];
    if(c.storesName == "cached") c.stores = spin::OutputStores::Cached;
    else if(c.storesName == "streaming") c.stores = spin::OutputStores::Streaming;
    else throw std::invalid_argument("stores must be cached or streaming");
    c.routeSeed = number(argv[6], "route seed");
    c.innerSeed = number(argv[7], "inner seed");
    const auto iterations = number(argv[8], "iterations");
    if(iterations == 0 || iterations > 10001)
        throw std::invalid_argument("iterations must be in [1,10001]");
    c.iterations = std::size_t(iterations);
    c.memoryName = argv[9];
    if(c.memoryName == "normal") c.memory = spin::MemoryPolicy::Normal;
    else if(c.memoryName == "huge") c.memory = spin::MemoryPolicy::PreferHugePages;
    else if(c.memoryName == "automatic") c.memory = spin::MemoryPolicy::Automatic;
    else if(c.memoryName == "flock") { c.memory = spin::MemoryPolicy::Normal; c.flockAdvice = true; }
    else throw std::invalid_argument("memory must be normal, huge, automatic, or flock");
    return c;
}

struct Advice {
    std::size_t bytes = 0;
    // 1 means not attempted, 0 means the request succeeded, negative is errno.
    // Successful advice does not prove every page has huge-page backing.
    int hugepageStatus = 1, collapseStatus = 1;
};
Advice adviseFlock(std::span<std::byte> bytes) {
    Advice result;
#if defined(__linux__) && defined(__x86_64__)
    // Match hypercat/src/pcs/brakedown.rs AlignedColumns::advise_huge_pages:
    // cache-line alignment, already first-touched ordinary allocation, advice
    // only on complete owned interior 2 MiB ranges, and both setup-time calls.
    constexpr std::uintptr_t hugePage = 2*1024*1024;
    constexpr int hugepageAdvice = 14, collapseAdvice = 25;
    const auto ptr = reinterpret_cast<std::uintptr_t>(bytes.data());
    const auto start = (ptr+hugePage-1)&~(hugePage-1);
    const auto end = (ptr+bytes.size())&~(hugePage-1);
    if(end>start) {
        result.bytes = std::size_t(end-start);
        result.hugepageStatus = ::madvise(reinterpret_cast<void*>(start),result.bytes,hugepageAdvice)==0 ? 0 : -errno;
        result.collapseStatus = ::madvise(reinterpret_cast<void*>(start),result.bytes,collapseAdvice)==0 ? 0 : -errno;
    }
#else
    // FLOCK's Linux advice is unavailable; preserve ordinary aligned buffers.
    (void)bytes;
#endif
    return result;
}

struct GuardedBuffer {
    spin::Buffer storage;
    std::size_t size;
    GuardedBuffer(std::size_t bytes, spin::MemoryPolicy memory)
        : storage(bytes+2*cacheLine, memory), size(bytes) {
        // Buffer construction first-touches and zeroes its entire allocation.
        // Both guards have separate cache lines; the payload remains aligned.
        auto all = storage.bytes();
        std::fill_n(all.data(), cacheLine, guardValue);
        std::fill_n(all.data()+cacheLine+size, cacheLine, guardValue);
    }
    std::span<std::byte> bytes() { return storage.bytes().subspan(cacheLine, size); }
    std::span<const std::byte> bytes() const { return storage.bytes().subspan(cacheLine, size); }
    void verify(const char* name) const {
        const auto all = storage.bytes();
        const auto intact = [](std::byte value) { return value == guardValue; };
        if(!std::all_of(all.begin(), all.begin()+cacheLine, intact) ||
           !std::all_of(all.begin()+cacheLine+size, all.end(), intact))
            throw std::runtime_error(std::string(name)+" guard was overwritten");
    }
};

inline std::uint64_t nextWord(std::uint64_t& state) {
    state += 0x9e3779b97f4a7c15ULL;
    auto value = state;
    value = (value^(value>>30))*0xbf58476d1ce4e5b9ULL;
    value = (value^(value>>27))*0x94d049bb133111ebULL;
    return value^(value>>31);
}
inline __m128i nextPair(std::uint64_t& state) {
    const auto lo = nextWord(state), hi = nextWord(state);
    return _mm_set_epi64x(static_cast<long long>(hi), static_cast<long long>(lo));
}
template<bool Stream>
inline void storePair(__m128i* output, __m128i value) {
    if constexpr(Stream) _mm_stream_si128(output, value);
    else _mm_store_si128(output, value);
}
template<bool Stream>
void fillInput(std::span<std::byte> bytes) {
    if((reinterpret_cast<std::uintptr_t>(bytes.data())&(cacheLine-1)) || bytes.size()%cacheLine)
        throw std::invalid_argument("input preparation requires complete aligned cache lines");
    auto* output = reinterpret_cast<__m128i*>(bytes.data());
    std::uint64_t state = inputSeed;
    for(std::size_t line=0; line<bytes.size()/cacheLine; ++line) {
        // Four explicit SSE2 stores fill each line; data comes from registers,
        // not a large source copy that would change the cache preparation.
        const auto a = nextPair(state), b = nextPair(state);
        const auto c = nextPair(state), d = nextPair(state);
        storePair<Stream>(output+4*line, a);
        storePair<Stream>(output+4*line+1, b);
        storePair<Stream>(output+4*line+2, c);
        storePair<Stream>(output+4*line+3, d);
    }
    if constexpr(Stream) _mm_sfence();
}
void flush(std::span<const std::byte> bytes) {
    for(std::size_t offset=0; offset<bytes.size(); offset+=cacheLine)
        _mm_clflush(bytes.data()+offset);
}
void prepare(CacheMode mode, std::span<std::byte> input, std::span<const std::byte> output) {
    switch(mode) {
    case CacheMode::Reused:
        // No sweep or checksum: previous calls determine cache residency.
        return;
    case CacheMode::Cold:
        flush(input);
        flush(output);
        _mm_mfence();
        return;
    case CacheMode::Produced:
        // Only input is rewritten. Output is deliberately left untouched.
        fillInput<true>(input);
        return;
    }
    throw std::invalid_argument("invalid cache mode");
}
std::uint64_t checksum(std::span<const std::byte> bytes) {
    std::uint64_t hash = 0xcbf29ce484222325ULL;
    for(std::size_t offset=0; offset<bytes.size(); offset+=8) {
        std::uint64_t word;
        std::memcpy(&word, bytes.data()+offset, 8);
        hash = (hash^word)*0x100000001b3ULL;
    }
    return hash;
}
struct Call {
    std::span<const std::byte> input;
    std::span<std::byte> output;
};
std::vector<Call> calls(std::span<const std::byte> input, std::span<std::byte> output,
                        std::size_t count, std::size_t inputBytes) {
    std::vector<Call> result(count);
    for(std::size_t i=0; i<count; ++i)
        result[i] = {input.subspan(i*inputBytes,inputBytes),
                     output.subspan(2*i*inputBytes,2*inputBytes)};
    return result;
}
void encode(const spin::Code& code, spin::Workspace& workspace,
            const std::vector<Call>& batch, spin::OutputStores stores) {
    for(const auto& call : batch) code.forward_bytes(call.input, call.output, workspace, stores);
}
const char* backend(spin::Backend value) {
    switch(value) {
    case spin::Backend::Automatic: return "automatic";
    case spin::Backend::Avx2: return "avx2";
    case spin::Backend::Avx512: return "avx512";
    case spin::Backend::Portable: return "portable";
    }
    throw std::runtime_error("unknown selected backend");
}
const char* preparation(CacheMode mode) {
    switch(mode) {
    case CacheMode::Reused: return "none_between_batches";
    case CacheMode::Cold: return "clflush_input_and_output_mfence";
    case CacheMode::Produced: return "nt_rewrite_input_sfence_output_untouched";
    }
    throw std::runtime_error("invalid cache mode");
}

void selfTest() {
    GuardedBuffer cached(512,spin::MemoryPolicy::Normal), streamed(512,spin::MemoryPolicy::Normal);
    fillInput<false>(cached.bytes());
    fillInput<true>(streamed.bytes());
    if(std::memcmp(cached.bytes().data(),streamed.bytes().data(),cached.size))
        throw std::runtime_error("cached and NT input generators differ");
    const auto expected = checksum(cached.bytes());
    for(auto mode : {CacheMode::Reused,CacheMode::Cold,CacheMode::Produced}) {
        prepare(mode,streamed.bytes(),cached.bytes());
        if(checksum(streamed.bytes()) != expected || checksum(cached.bytes()) != expected)
            throw std::runtime_error("cache preparation changed input data or output");
    }
    cached.verify("self-test cached buffer");
    streamed.verify("self-test streamed buffer");
    for(auto text : {"", "-1", "+1", "1x", "18446744073709551616"}) {
        bool rejected = false;
        try { (void)number(text,"self-test"); } catch(const std::invalid_argument&) { rejected = true; }
        if(!rejected) throw std::runtime_error("strict integer parser accepted invalid input");
    }
    if(number("18446744073709551615","self-test") != UINT64_MAX)
        throw std::runtime_error("integer parser rejected maximum seed");
    std::cout << "PASS: deterministic SSE2 cache preparation, guards, and strict integer parsing\n";
}

void run(const Config& c) {
    const auto bytesPerCall = k*(c.bits/8);
    const auto inputBytes = c.calls*bytesPerCall;
    const auto outputBytes = 2*inputBytes;
    spin::Code code({k,c.parameters,c.routeSeed,c.innerSeed});
    if(!code.supports_forward(c.width)) throw std::runtime_error("selected backend does not support this width");
    // Match the FFI: its native workspace always uses Automatic independently
    // of caller I/O page policy. Reference storage is outside the measured set.
    auto workspace = code.make_workspace(c.width,spin::MemoryPolicy::Automatic);
    GuardedBuffer input(inputBytes,c.memory), output(outputBytes,c.memory);
    GuardedBuffer reference(outputBytes,spin::MemoryPolicy::Normal);
    const auto inputAdvice = c.flockAdvice ? adviseFlock(input.bytes()) : Advice{};
    const auto outputAdvice = c.flockAdvice ? adviseFlock(output.bytes()) : Advice{};
    fillInput<false>(input.bytes());
    const auto expectedInput = checksum(input.bytes());
    const auto measuredCalls = calls(input.bytes(),output.bytes(),c.calls,bytesPerCall);
    const auto referenceCalls = calls(input.bytes(),reference.bytes(),c.calls,bytesPerCall);
    std::vector<double> times(c.iterations);

    // Establish exact cached-reference equality before warming the case.
    // Reference storage is not read again until every timed trial is over.
    encode(code,workspace,referenceCalls,spin::OutputStores::Cached);
    encode(code,workspace,measuredCalls,c.stores);
    if(std::memcmp(output.bytes().data(),reference.bytes().data(),outputBytes))
        throw std::runtime_error("selected policy differs from cached reference before timing");
    input.verify("input"); output.verify("output"); reference.verify("reference");
    for(unsigned i=0; i<warmupBatches; ++i) {
        prepare(c.cache,input.bytes(),output.bytes());
        encode(code,workspace,measuredCalls,c.stores);
    }
    for(std::size_t i=0; i<c.iterations; ++i) {
        prepare(c.cache,input.bytes(),output.bytes());
        const auto begin = Clock::now();
        encode(code,workspace,measuredCalls,c.stores);
        const auto end = Clock::now();
        times[i] = std::chrono::duration<double,std::milli>(end-begin).count();
    }

    // Deliberately no checksum, guard read, or reference read between trials.
    const auto outputHash = checksum(output.bytes());
    if(checksum(input.bytes()) != expectedInput)
        throw std::runtime_error("input changed during encoding or cache preparation");
    if(std::memcmp(output.bytes().data(),reference.bytes().data(),outputBytes))
        throw std::runtime_error("measured output differs from cached reference");
    input.verify("input"); output.verify("output"); reference.verify("reference");
    auto sorted = times;
    std::sort(sorted.begin(),sorted.end());
    const auto median = sorted[sorted.size()/2];
    std::cout << "profile,k,bits,shape,cache_mode,cache_preparation,stores,route_seed,inner_seed,input_seed,"
        "iterations,warmup_batches,calls_per_batch,input_bytes,output_bytes,total_bytes,backend,memory_request,workspace_memory,"
        "setup_bytes,scratch_bytes,input_allocation_bytes,output_allocation_bytes,reference_allocation_bytes,"
        "input_advice_bytes,output_advice_bytes,input_hugepage_status,input_collapse_status,output_hugepage_status,output_collapse_status,"
        "phase,median_batch_ms,p10_batch_ms,p90_batch_ms,median_call_ms,input_checksum,output_checksum,reference_match,samples_ms\n"
        << c.profile << ',' << k << ',' << c.bits << ',' << c.shape << ',' << c.cacheName << ','
        << preparation(c.cache) << ',' << c.storesName << ',' << c.routeSeed << ',' << c.innerSeed << ',' << inputSeed << ','
        << c.iterations << ',' << warmupBatches << ',' << c.calls << ',' << inputBytes << ',' << outputBytes << ','
        << inputBytes+outputBytes << ',' << backend(code.backend()) << ',' << c.memoryName << ",automatic,"
        << code.setup_bytes() << ',' << workspace.bytes() << ',' << input.storage.allocation_bytes() << ','
        << output.storage.allocation_bytes() << ',' << reference.storage.allocation_bytes() << ','
        << inputAdvice.bytes << ',' << outputAdvice.bytes << ',' << inputAdvice.hugepageStatus << ','
        << inputAdvice.collapseStatus << ',' << outputAdvice.hugepageStatus << ',' << outputAdvice.collapseStatus << ','
        << "encode_only_no_consumer," << std::fixed << std::setprecision(6) << median << ',' << sorted[sorted.size()/10] << ','
        << sorted[9*sorted.size()/10] << ',' << median/double(c.calls) << ','
        << std::hex << expectedInput << ',' << outputHash << std::dec << ",true,";
    for(std::size_t i=0; i<times.size(); ++i) {
        if(i) std::cout << ';';
        std::cout << times[i];
    }
    std::cout << '\n';
}
}

int main(int argc,char** argv) {
    try {
        if(argc==2 && std::string_view(argv[1])=="--self-test") { selfTest(); return 0; }
        run(parse(argc,argv));
        return 0;
    } catch(const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
