// Fixed BCH(256) endpoint experiment. No adaptive stopping or seed search.
// Build: cl /O2 /std:c++20 /EHsc /arch:AVX2 test_bch_endpoints.cpp bcrypt.lib
#define NOMINMAX
#include <windows.h>
#include <bcrypt.h>
#include <immintrin.h>
#include <array>
#include <bit>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

namespace {
using Byte = std::uint8_t;
using Poly = std::array<Byte, 19>;
constexpr std::size_t Chunk = 1 << 20;

Byte multiply(Byte a, Byte b) noexcept {
    unsigned x = a, y = b, v = 0;
    while (y) {
        if (y & 1) v ^= x;
        y >>= 1;
        x <<= 1;
        if (x & 256) x ^= 0x14d;
    }
    return static_cast<Byte>(v);
}

// Separate carryless-product/remainder implementation for field-table tests.
Byte multiply_reference(Byte a, Byte b) noexcept {
    unsigned v = 0;
    for (unsigned j = 0; j < 8; ++j) if ((b >> j) & 1) v ^= unsigned(a) << j;
    for (int j = 14; j >= 8; --j) if ((v >> j) & 1) v ^= 0x14dU << (j-8);
    return static_cast<Byte>(v);
}

struct alignas(32) Tables {
    Byte mul[256][256]{};
    Byte inverse[256]{};
    alignas(32) Byte target[256]{};
    alignas(32) Byte term[19][256][256]{};

    Byte power(Byte x, unsigned e) const noexcept {
        Byte v = 1;
        for (; e; e >>= 1, x = mul[x][x]) if (e & 1) v = mul[v][x];
        return v;
    }
    Tables() {
        for (unsigned a=0; a<256; ++a) for (unsigned b=0; b<256; ++b)
            mul[a][b] = multiply(static_cast<Byte>(a), static_cast<Byte>(b));
        for (unsigned x=1; x<256; ++x) {
            inverse[x] = power(static_cast<Byte>(x), 254);
            target[x] = power(static_cast<Byte>(x), 146);
        }
        for (unsigned x=0; x<256; ++x) {
            Byte p=1;
            for (unsigned d=0; d<=18; ++d, p=mul[p][x])
                for (unsigned c=0; c<256; ++c) term[d][c][x]=mul[c][p];
        }
    }
};

Byte scalar_value(const Poly& g, Byte x, const Tables& t, int degree=18) noexcept {
    Byte v=0;
    for (int j=degree; j>=0; --j) v=t.mul[v][x]^g[j];
    return v;
}

Poly interpolate_newton(const Poly& nodes, const Tables& t) noexcept {
    Poly divided{};
    divided[0]=1;
    for (unsigned i=1; i<=18; ++i) divided[i]=t.target[nodes[i]];
    for (unsigned d=1; d<=18; ++d)
        for (int i=18; i>=static_cast<int>(d); --i)
            divided[i]=t.mul[divided[i]^divided[i-1]][t.inverse[nodes[i]^nodes[i-d]]];
    Poly g{};
    g[0]=divided[18];
    unsigned degree=0;
    for (int i=17; i>=0; --i) {
        Poly next{};
        for (unsigned j=0; j<=degree; ++j) {
            next[j]^=t.mul[g[j]][nodes[i]];
            next[j+1]^=g[j];
        }
        next[0]^=divided[i]; g=next; ++degree;
    }
    return g;
}

// Independently build g by adding multiples of the vanishing polynomial.
Poly interpolate_incremental(const Poly& nodes, const Tables& t) noexcept {
    Poly g{}, basis{};
    g[0]=1; basis[1]=1;
    for (unsigned i=1; i<=18; ++i) {
        const Byte x=nodes[i];
        const Byte denom=scalar_value(basis,x,t,static_cast<int>(i));
        const Byte residual=t.target[x]^scalar_value(g,x,t,static_cast<int>(i-1));
        const Byte scale=t.mul[residual][t.inverse[denom]];
        for (unsigned j=0; j<=i; ++j) g[j]^=t.mul[scale][basis[j]];
        if (i<18) {
            for (int j=static_cast<int>(i+1); j>=1; --j)
                basis[j]=basis[j-1]^t.mul[x][basis[j]];
            basis[0]=t.mul[x][basis[0]];
        }
    }
    return g;
}

unsigned agreements_vector(const Poly& g, const Tables& t) noexcept {
    // Eight explicit lanes cover the field. All tables and loads are aligned.
    __m256i a0=_mm256_set1_epi8(g[0]), a1=a0, a2=a0, a3=a0;
    __m256i a4=a0, a5=a0, a6=a0, a7=a0;
    for (unsigned d=1; d<=18; ++d) {
        const auto* p=reinterpret_cast<const __m256i*>(t.term[d][g[d]]);
        a0=_mm256_xor_si256(a0,_mm256_load_si256(p+0));
        a1=_mm256_xor_si256(a1,_mm256_load_si256(p+1));
        a2=_mm256_xor_si256(a2,_mm256_load_si256(p+2));
        a3=_mm256_xor_si256(a3,_mm256_load_si256(p+3));
        a4=_mm256_xor_si256(a4,_mm256_load_si256(p+4));
        a5=_mm256_xor_si256(a5,_mm256_load_si256(p+5));
        a6=_mm256_xor_si256(a6,_mm256_load_si256(p+6));
        a7=_mm256_xor_si256(a7,_mm256_load_si256(p+7));
    }
    const auto* p=reinterpret_cast<const __m256i*>(t.target);
    const auto count=[](__m256i a, const __m256i* b) noexcept {
        return std::popcount(static_cast<unsigned>(_mm256_movemask_epi8(
            _mm256_cmpeq_epi8(a,_mm256_load_si256(b)))));
    };
    return count(a0,p)+count(a1,p+1)+count(a2,p+2)+count(a3,p+3)
         + count(a4,p+4)+count(a5,p+5)+count(a6,p+6)+count(a7,p+7);
}

unsigned agreements_horner(const Poly& g, const Tables& t) noexcept {
    unsigned count=0;
    for (unsigned x=0; x<256; x+=8) {
        Byte v0=0,v1=0,v2=0,v3=0,v4=0,v5=0,v6=0,v7=0;
        // Fixed-width independent chains expose instruction-level parallelism.
        for (int j=18; j>=0; --j) {
            v0=t.mul[x+0][v0]^g[j]; v1=t.mul[x+1][v1]^g[j];
            v2=t.mul[x+2][v2]^g[j]; v3=t.mul[x+3][v3]^g[j];
            v4=t.mul[x+4][v4]^g[j]; v5=t.mul[x+5][v5]^g[j];
            v6=t.mul[x+6][v6]^g[j]; v7=t.mul[x+7][v7]^g[j];
        }
        count+=(v0==t.target[x])+ (v1==t.target[x+1])+ (v2==t.target[x+2])
             + (v3==t.target[x+3])+ (v4==t.target[x+4])+ (v5==t.target[x+5])
             + (v6==t.target[x+6])+ (v7==t.target[x+7]);
    }
    return count;
}

void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

void selftest(const Tables& t) {
    for (unsigned a=0;a<256;++a) for (unsigned b=0;b<256;++b)
        require(t.mul[a][b]==multiply_reference(Byte(a),Byte(b)),"field table");
    for (unsigned a=1;a<256;++a) require(t.mul[a][t.inverse[a]]==1,"inverse");
    const Poly witness={1,0x8e,0x6b,0x77,0x78,0xe3,0xc4,0x93,0xad,0xb1,
                        0xea,0x39,0x29,0x22,0x69,0x2a,0x42,0x61,0xa8};
    require(agreements_vector(witness,t)==37,"known endpoint AVX2");
    require(agreements_horner(witness,t)==37,"known endpoint Horner");
    Poly nodes{};
    unsigned k=1;
    for (unsigned x=1;x<256 && k<=18;++x)
        if (scalar_value(witness,Byte(x),t)==t.target[x]) nodes[k++]=Byte(x);
    require(k==19,"witness base");
    require(interpolate_newton(nodes,t)==witness,"witness Newton");
    require(interpolate_incremental(nodes,t)==witness,"witness incremental");
    // Deterministic diagnostic vectors, not statistical experiment samples.
    std::uint64_t state=0x317a03efbe12719ULL;
    for (unsigned test=0;test<4096;++test) {
        std::array<bool,256> used{}; k=1; nodes={};
        while(k<=18) {
            state^=state<<13; state^=state>>7; state^=state<<17;
            const Byte x=static_cast<Byte>(state);
            if(x && !used[x]) {used[x]=true; nodes[k++]=x;}
        }
        auto a=interpolate_newton(nodes,t), b=interpolate_incremental(nodes,t);
        require(a==b && a[0]==1,"interpolation cross-check");
        for (unsigned i=1;i<=18;++i)
            require(scalar_value(a,nodes[i],t)==t.target[nodes[i]],"base agreement");
        const unsigned n=agreements_vector(a,t);
        require(n>=18 && n<=37 && n==agreements_horner(a,t),"root-count cross-check");
    }
    std::cout << "PASS: field, known endpoint, 4096 independent kernel cross-checks\n";
}

template<bool Replay> struct ByteStream {
    HANDLE file=INVALID_HANDLE_VALUE;
    std::unique_ptr<Byte[]> buffer=std::make_unique<Byte[]>(Chunk);
    std::size_t pos=0, limit=0;
    std::uint64_t consumed=0, loaded=0;
    explicit ByteStream(const std::string& path) {
        file=CreateFileA(path.c_str(),Replay?GENERIC_READ:GENERIC_WRITE,FILE_SHARE_READ,
                         nullptr,Replay?OPEN_EXISTING:CREATE_NEW,FILE_FLAG_SEQUENTIAL_SCAN,nullptr);
        require(file!=INVALID_HANDLE_VALUE,"cannot open random tape (refusing overwrite)");
    }
    ~ByteStream() { if(file!=INVALID_HANDLE_VALUE) CloseHandle(file); }
    Byte next() {
        if(pos==limit) {
            DWORD done=0;
            if constexpr(Replay) {
                require(ReadFile(file,buffer.get(),DWORD(Chunk),&done,nullptr)!=0 && done>0,
                        "random tape exhausted or unreadable");
            } else {
                require(BCryptGenRandom(nullptr,buffer.get(),ULONG(Chunk),
                        BCRYPT_USE_SYSTEM_PREFERRED_RNG)==0,"BCryptGenRandom failed");
                require(WriteFile(file,buffer.get(),DWORD(Chunk),&done,nullptr)!=0 && done==Chunk,
                        "random tape write failed");
            }
            pos=0; limit=done; loaded+=done;
        }
        ++consumed; return buffer[pos++];
    }
    void finish() {
        if constexpr(!Replay) require(FlushFileBuffers(file)!=0,"random tape flush failed");
    }
};

template<bool Replay> void experiment(std::uint64_t n,const std::string& tape,
                                     const std::string& output,const Tables& t) {
    require(n>0,"sample count must be positive");
    require(!std::filesystem::exists(output),"result already exists");
    require(!std::filesystem::exists(output+".audit.jsonl"),"audit file already exists");
    std::ofstream audit(output+".audit.jsonl"); require(bool(audit),"audit open failed");
    ByteStream<Replay> source(tape);
    std::array<std::uint64_t,20> histogram{};
    std::uint64_t zeros=0,duplicates=0;
    const auto start=std::chrono::steady_clock::now();
    for(std::uint64_t i=0;i<n;++i) {
        const auto offset=source.consumed;
        Poly nodes{};
        unsigned k=1;
        if constexpr(Replay) {
            std::array<bool,256> used{};
            while(k<=18) {
                const Byte x=source.next();
                if(!x) {++zeros; continue;}
                if(used[x]) {++duplicates; continue;}
                used[x]=true; nodes[k++]=x;
            }
        } else {
            std::array<std::uint64_t,4> seen{};
            while(k<=18) {
                const Byte x=source.next();
                if(!x) {++zeros; continue;}
                const auto bit=std::uint64_t(1)<<(x&63);
                if(seen[x>>6]&bit) {++duplicates; continue;}
                seen[x>>6]|=bit; nodes[k++]=x;
            }
        }
        Poly g;
        unsigned count;
        if constexpr(Replay) {
            g=interpolate_incremental(nodes,t); count=agreements_horner(g,t);
        } else {
            g=interpolate_newton(nodes,t); count=agreements_vector(g,t);
        }
        require(g[0]==1 && count>=18 && count<=37,"invalid polynomial or root count");
        ++histogram[count-18];
        if(i<256 || i%1000000==0 || count==37 || i+1==n) {
            audit << "{\"sample\":" << i << ",\"byte_offset\":" << offset
                  << ",\"next_byte_offset\":" << source.consumed
                  << ",\"agreements\":" << count << ",\"nodes\":[";
            for(unsigned j=0;j<19;++j) audit << (j?",":"") << unsigned(nodes[j]);
            audit << "],\"coefficients\":[";
            for(unsigned j=0;j<19;++j) audit << (j?",":"") << unsigned(g[j]);
            audit << "]}\n";
            require(bool(audit),"audit write failed");
        }
        if((i+1)%5000000==0) {
            const auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
            std::cout << (Replay?"replay":"sample") << " samples=" << i+1
                      << " endpoints=" << histogram[19] << " seconds=" << seconds << std::endl;
        }
    }
    source.finish(); audit.close(); require(!audit.fail(),"audit close failed");
    const auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    std::ofstream out(output); require(bool(out),"result open failed");
    out << std::setprecision(17) << "{\n\"mode\":\"" << (Replay?"independent_replay":"BCryptGenRandom")
        << "\",\n\"samples\":" << n << ",\n\"endpoint_hits\":" << histogram[19]
        << ",\n\"bytes_consumed\":" << source.consumed
        << ",\n\"bytes_loaded\":" << source.loaded << ",\n\"zero_bytes_rejected\":" << zeros
        << ",\n\"duplicate_bytes_rejected\":" << duplicates
        << ",\n\"seconds\":" << seconds << ",\n\"extra_agreement_histogram\":[";
    for(unsigned j=0;j<20;++j) out << (j?",":"") << histogram[j];
    out << "]\n}\n"; out.close(); require(!out.fail(),"result close failed");
    std::cout << "COMPLETE " << (Replay?"replay":"sample") << " n=" << n
              << " endpoints=" << histogram[19] << " seconds=" << seconds << std::endl;
}
} // namespace

int main(int argc,char** argv) try {
    const auto tables=std::make_unique<Tables>();
    if(argc==2 && std::string(argv[1])=="selftest") {selftest(*tables); return 0;}
    if(argc!=6 || (std::string(argv[1])!="run" && std::string(argv[1])!="replay") ||
       std::string(argv[5])!="fixed-bch256-endpoint-v1") {
        std::cerr << "usage: test_bch_endpoints run|replay N TAPE OUTPUT fixed-bch256-endpoint-v1\n";
        return 2;
    }
    const auto n=std::stoull(argv[2]);
    if(std::string(argv[1])=="run") experiment<false>(n,argv[3],argv[4],*tables);
    else experiment<true>(n,argv[3],argv[4],*tables);
    return 0;
} catch(const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << '\n'; return 1;
}
