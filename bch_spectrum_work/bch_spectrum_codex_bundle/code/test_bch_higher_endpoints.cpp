// Fixed BCH weight-40/42 endpoint sampling. Preflight data proves no shell cap.
// cl /O2 /std:c++20 /EHsc /arch:AVX2 test_bch_higher_endpoints.cpp bcrypt.lib
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
using Byte=std::uint8_t;
constexpr std::size_t Chunk=1<<20;
void require(bool ok,const char* why) {if(!ok) throw std::runtime_error(why);}
Byte multiply(Byte a,Byte b) noexcept {
    unsigned x=a,y=b,v=0;
    while(y) {if(y&1) v^=x; y>>=1; x<<=1; if(x&256) x^=0x14d;}
    return Byte(v);
}
Byte reference_multiply(Byte a,Byte b) noexcept {
    unsigned v=0;
    for(unsigned j=0;j<8;++j) if((b>>j)&1) v^=unsigned(a)<<j;
    for(int j=14;j>=8;--j) if((v>>j)&1) v^=0x14dU<<(j-8);
    return Byte(v);
}
struct alignas(32) Tables {
    Byte mul[256][256]{},inverse[256]{},target[3][256]{};
    alignas(32) Byte term[24][256][256]{};
    Byte power(Byte x,unsigned e) const noexcept {
        Byte v=1;
        for(;e;e>>=1,x=mul[x][x]) if(e&1) v=mul[v][x];
        return v;
    }
    Tables() {
        for(unsigned a=0;a<256;++a) for(unsigned b=0;b<256;++b)
            mul[a][b]=multiply(Byte(a),Byte(b));
        for(unsigned x=0;x<256;++x) {
            if(x) inverse[x]=power(Byte(x),254);
            for(unsigned j=0;j<3;++j) target[j][x]=power(Byte(x),146+j);
            for(unsigned d=0;d<24;++d) {
                const Byte p=power(Byte(x),d<21?d:146+d-21);
                for(unsigned c=0;c<256;++c) term[d][c][x]=mul[c][p];
            }
        }
    }
};
template<std::size_t N> Byte evaluate(const std::array<Byte,N>& g,Byte x,
                                     const Tables& t,unsigned degree=N-1) noexcept {
    Byte v=0;
    for(int j=int(degree);j>=0;--j) v=t.mul[x][v]^g[j];
    return v;
}
template<int R> struct Candidate {
    static constexpr unsigned N=18+2*R,M=18+R,D=37+2*R;
    std::array<Byte,M+1> g{};
    std::array<Byte,R+1> h{};
    bool valid=false;
    unsigned count=0;
    bool operator==(const Candidate&) const = default;
};

template<unsigned N> std::array<Byte,N> newton(const std::array<Byte,N>& nodes,
                                               const Tables& t) noexcept {
    std::array<Byte,N> divided{},g{};
    for(unsigned i=0;i<N;++i) divided[i]=t.target[0][nodes[i]];
    for(unsigned d=1;d<N;++d) for(int i=N-1;i>=int(d);--i)
        divided[i]=t.mul[divided[i]^divided[i-1]][t.inverse[nodes[i]^nodes[i-d]]];
    g[0]=divided[N-1];
    for(int i=N-2;i>=0;--i) {
        const unsigned degree=N-1-i;
        for(int j=int(degree);j>=1;--j) g[j]=g[j-1]^t.mul[nodes[i]][g[j]];
        g[0]=t.mul[nodes[i]][g[0]]^divided[i];
    }
    return g;
}
// Replay reconstructs each remainder independently by incremental interpolation.
template<unsigned N> std::array<Byte,N> incremental(const std::array<Byte,N>& nodes,
                                                   unsigned shift,const Tables& t) noexcept {
    std::array<Byte,N> g{},basis{};
    basis[0]=1;
    for(unsigned i=0;i<N;++i) {
        const Byte x=nodes[i];
        const Byte scale=t.mul[t.target[shift][x]^evaluate(g,x,t)][
            t.inverse[evaluate(basis,x,t,i)]];
        for(unsigned j=0;j<=i;++j) g[j]^=t.mul[scale][basis[j]];
        if(i+1<N) {
            for(int j=int(i+1);j>=1;--j) basis[j]=basis[j-1]^t.mul[x][basis[j]];
            basis[0]=t.mul[x][basis[0]];
        }
    }
    return g;
}
template<int R,bool Replay> Candidate<R> reconstruct(
    const std::array<Byte,Candidate<R>::N>& nodes,Byte sigma,const Tables& t) noexcept {
    constexpr auto N=Candidate<R>::N,M=Candidate<R>::M;
    using Poly=std::array<Byte,N>;
    std::array<Poly,R+1> rem{};
    if constexpr(Replay) {
        for(unsigned j=0;j<=R;++j) rem[j]=incremental<N>(nodes,j,t);
    } else {
        rem[0]=newton<N>(nodes,t);
        std::array<Byte,N+1> basis{};
        basis[0]=1;
        for(unsigned i=0;i<N;++i) {
            for(int j=int(i+1);j>=1;--j) basis[j]=basis[j-1]^t.mul[nodes[i]][basis[j]];
            basis[0]=t.mul[nodes[i]][basis[0]];
        }
        for(unsigned j=1;j<=R;++j) {
            const Byte leading=rem[j-1][N-1];
            rem[j][0]=t.mul[leading][basis[0]];
            for(unsigned i=1;i<N;++i) rem[j][i]=rem[j-1][i-1]^t.mul[leading][basis[i]];
        }
    }
    Candidate<R> result;
    std::array<Byte,R> solution{};
    const Byte rhs0=1^t.mul[sigma][rem[0][0]];
    if constexpr(Replay) {
        Byte matrix[R][R+1]{};
        for(unsigned i=0;i<R;++i) {
            const unsigned degree=i?M+i:0;
            for(unsigned j=0;j<R;++j) matrix[i][j]=rem[j+1][degree];
            matrix[i][R]=Byte(i==0)^t.mul[sigma][rem[0][degree]];
        }
        for(unsigned j=0;j<R;++j) {
            unsigned pivot=j;
            while(pivot<R && !matrix[pivot][j]) ++pivot;
            if(pivot==R) return result;
            for(unsigned k=0;k<=R;++k) std::swap(matrix[j][k],matrix[pivot][k]);
            const Byte inv=t.inverse[matrix[j][j]];
            for(unsigned k=0;k<=R;++k) matrix[j][k]=t.mul[inv][matrix[j][k]];
            for(unsigned i=0;i<R;++i) if(i!=j) {
                const Byte scale=matrix[i][j];
                for(unsigned k=0;k<=R;++k) matrix[i][k]^=t.mul[scale][matrix[j][k]];
            }
        }
        for(unsigned j=0;j<R;++j) solution[j]=matrix[j][R];
    } else if constexpr(R==1) {
        if(!rem[1][0]) return result;
        solution[0]=t.mul[rhs0][t.inverse[rem[1][0]]];
    } else {
        const Byte a=rem[1][0],b=rem[2][0],c=rem[1][N-1],d=rem[2][N-1];
        const Byte determinant=t.mul[a][d]^t.mul[b][c];
        if(!determinant) return result;
        const Byte rhs1=t.mul[sigma][rem[0][N-1]],inv=t.inverse[determinant];
        solution[0]=t.mul[t.mul[rhs0][d]^t.mul[b][rhs1]][inv];
        solution[1]=t.mul[t.mul[a][rhs1]^t.mul[rhs0][c]][inv];
    }
    result.h[0]=sigma;
    for(unsigned j=0;j<R;++j) result.h[j+1]=solution[j];
    for(unsigned j=0;j<=R;++j) for(unsigned i=0;i<=M;++i)
        result.g[i]^=t.mul[result.h[j]][rem[j][i]];
    result.valid=true;
    return result;
}

template<int R> unsigned count_vector(const Candidate<R>& c,const Tables& t) noexcept {
    __m256i a0=_mm256_set1_epi8(c.g[0]),a1=a0,a2=a0,a3=a0,a4=a0,a5=a0,a6=a0,a7=a0;
    const auto add=[&](const Byte* row) {
        const auto* p=reinterpret_cast<const __m256i*>(row);
        a0=_mm256_xor_si256(a0,_mm256_load_si256(p+0));
        a1=_mm256_xor_si256(a1,_mm256_load_si256(p+1));
        a2=_mm256_xor_si256(a2,_mm256_load_si256(p+2));
        a3=_mm256_xor_si256(a3,_mm256_load_si256(p+3));
        a4=_mm256_xor_si256(a4,_mm256_load_si256(p+4));
        a5=_mm256_xor_si256(a5,_mm256_load_si256(p+5));
        a6=_mm256_xor_si256(a6,_mm256_load_si256(p+6));
        a7=_mm256_xor_si256(a7,_mm256_load_si256(p+7));
    };
    for(unsigned j=1;j<=c.M;++j) add(t.term[j][c.g[j]]);
    for(unsigned j=0;j<=R;++j) add(t.term[21+j][c.h[j]]);
    const auto count=[](__m256i x) {
        return std::popcount(unsigned(_mm256_movemask_epi8(_mm256_cmpeq_epi8(x,_mm256_setzero_si256()))));
    };
    // u=0 never contributes because g(0)=1 and u^146 h(u)=0.
    return count(a0)+count(a1)+count(a2)+count(a3)+count(a4)+count(a5)+count(a6)+count(a7);
}
template<int R> unsigned count_horner(const Candidate<R>& c,const Tables& t) noexcept {
    unsigned count=0;
    for(unsigned x=0;x<256;x+=8) {
        Byte v0=0,v1=0,v2=0,v3=0,v4=0,v5=0,v6=0,v7=0;
        for(int j=c.M;j>=0;--j) {
            v0=t.mul[x+0][v0]^c.g[j]; v1=t.mul[x+1][v1]^c.g[j];
            v2=t.mul[x+2][v2]^c.g[j]; v3=t.mul[x+3][v3]^c.g[j];
            v4=t.mul[x+4][v4]^c.g[j]; v5=t.mul[x+5][v5]^c.g[j];
            v6=t.mul[x+6][v6]^c.g[j]; v7=t.mul[x+7][v7]^c.g[j];
        }
        count+=(v0==t.mul[t.target[0][x+0]][evaluate(c.h,Byte(x+0),t)]);
        count+=(v1==t.mul[t.target[0][x+1]][evaluate(c.h,Byte(x+1),t)]);
        count+=(v2==t.mul[t.target[0][x+2]][evaluate(c.h,Byte(x+2),t)]);
        count+=(v3==t.mul[t.target[0][x+3]][evaluate(c.h,Byte(x+3),t)]);
        count+=(v4==t.mul[t.target[0][x+4]][evaluate(c.h,Byte(x+4),t)]);
        count+=(v5==t.mul[t.target[0][x+5]][evaluate(c.h,Byte(x+5),t)]);
        count+=(v6==t.mul[t.target[0][x+6]][evaluate(c.h,Byte(x+6),t)]);
        count+=(v7==t.mul[t.target[0][x+7]][evaluate(c.h,Byte(x+7),t)]);
    }
    return count;
}
template<int R,bool Replay> Candidate<R> trial(
    const std::array<Byte,Candidate<R>::N>& nodes,Byte sigma,const Tables& t) {
    auto c=reconstruct<R,Replay>(nodes,sigma,t);
    if(c.valid) {
        require(c.g[0]==1,"bad constant");
        if constexpr(Replay) c.count=count_horner(c,t); else c.count=count_vector(c,t);
        require(c.count>=c.N && c.count<=c.D,"bad agreement count");
    }
    return c;
}
template<bool Replay> struct ByteStream {
    HANDLE file=INVALID_HANDLE_VALUE;
    std::unique_ptr<Byte[]> buffer=std::make_unique<Byte[]>(Chunk);
    std::size_t pos=0,limit=0;
    std::uint64_t consumed=0,loaded=0;
    explicit ByteStream(const std::string& path) {
        file=CreateFileA(path.c_str(),Replay?GENERIC_READ:GENERIC_WRITE,FILE_SHARE_READ,
            nullptr,Replay?OPEN_EXISTING:CREATE_NEW,FILE_FLAG_SEQUENTIAL_SCAN,nullptr);
        require(file!=INVALID_HANDLE_VALUE,"cannot open tape; refusing overwrite");
    }
    ~ByteStream() {if(file!=INVALID_HANDLE_VALUE) CloseHandle(file);}
    Byte next() {
        if(pos==limit) {
            DWORD done=0;
            if constexpr(Replay) require(ReadFile(file,buffer.get(),DWORD(Chunk),&done,nullptr) && done,"tape exhausted");
            else {
                require(BCryptGenRandom(nullptr,buffer.get(),ULONG(Chunk),BCRYPT_USE_SYSTEM_PREFERRED_RNG)==0,"CNG failed");
                require(WriteFile(file,buffer.get(),DWORD(Chunk),&done,nullptr) && done==Chunk,"tape write failed");
            }
            pos=0; limit=done; loaded+=done;
        }
        ++consumed; return buffer[pos++];
    }
    void finish() {if constexpr(!Replay) require(FlushFileBuffers(file)!=0,"tape flush failed");}
};
template<class Values> void json_array(std::ostream& out,const Values& values) {
    out<<'['; bool first=true;
    for(auto v:values) {if(!first) out<<','; first=false; out<<std::uint64_t(v);}
    out<<']';
}
template<int R,bool Replay> void experiment(std::uint64_t samples,const std::string& tape,
                                           const std::string& output,const Tables& t) {
    constexpr unsigned N=Candidate<R>::N,D=Candidate<R>::D;
    require(samples>0,"sample count must be positive");
    require(!std::filesystem::exists(output) && !std::filesystem::exists(output+".audit.jsonl"),"output exists");
    std::ofstream audit(output+".audit.jsonl"); require(bool(audit),"audit open failed");
    ByteStream<Replay> source(tape);
    std::array<std::uint64_t,D+1> histogram{};
    std::uint64_t singular=0,zeros=0,duplicates=0;
    const auto start=std::chrono::steady_clock::now();
    for(std::uint64_t i=0;i<samples;++i) {
        const auto offset=source.consumed;
        // One byte maps uniformly onto S={0,...,31}, before any subset bytes.
        const Byte sigma=source.next()&31;
        std::array<Byte,N> nodes{}; unsigned k=0;
        if constexpr(Replay) {
            std::array<bool,256> seen{};
            while(k<N) {
                const Byte x=source.next();
                if(!x) {++zeros; continue;} if(seen[x]) {++duplicates; continue;}
                seen[x]=true; nodes[k++]=x;
            }
        } else {
            std::array<std::uint64_t,4> seen{};
            while(k<N) {
                const Byte x=source.next();
                if(!x) {++zeros; continue;}
                const auto bit=std::uint64_t(1)<<(x&63);
                if(seen[x>>6]&bit) {++duplicates; continue;}
                seen[x>>6]|=bit; nodes[k++]=x;
            }
        }
        const auto c=trial<R,Replay>(nodes,sigma,t);
        if(c.valid) ++histogram[c.count]; else ++singular;
        // Bounded singular audits avoid millions of full Gaussian trace checks.
        // All endpoints are retained. Neither branch alters the sampling law.
        if(i<256 || i%1000000==0 || (!c.valid && singular<=32) || c.count==D || i+1==samples) {
            audit<<"{\"sample\":"<<i<<",\"byte_offset\":"<<offset<<",\"next_byte_offset\":"<<source.consumed
                 <<",\"sigma\":"<<unsigned(sigma)<<",\"valid\":"<<(c.valid?"true":"false")
                 <<",\"agreements\":"<<c.count<<",\"nodes\":";
            json_array(audit,nodes); audit<<",\"g\":"; json_array(audit,c.g);
            audit<<",\"h\":"; json_array(audit,c.h); audit<<"}\n";
            require(bool(audit),"audit write failed");
        }
    }
    source.finish(); audit.close(); require(!audit.fail(),"audit close failed");
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    std::ofstream out(output); require(bool(out),"result open failed");
    out<<std::setprecision(17)<<"{\"mode\":\""<<(Replay?"independent_replay":"BCryptGenRandom")
       <<"\",\"weight\":"<<D+1<<",\"samples\":"<<samples<<",\"endpoint_hits\":"<<histogram[D]
       <<",\"singular_queries\":"<<singular<<",\"bytes_consumed\":"<<source.consumed
       <<",\"bytes_loaded\":"<<source.loaded<<",\"zero_bytes_rejected\":"<<zeros
       <<",\"duplicate_bytes_rejected\":"<<duplicates<<",\"seconds\":"<<seconds<<",\"agreement_histogram\":";
    json_array(out,histogram); out<<"}\n"; out.close(); require(!out.fail(),"result close failed");
    std::cout<<"COMPLETE "<<(Replay?"replay":"sample")<<" weight="<<D+1<<" samples="<<samples
             <<" endpoints="<<histogram[D]<<" singular="<<singular<<" seconds="<<seconds<<std::endl;
}
template<int R> void check_vector(std::istream& in,const Tables& t) {
    unsigned sigma,valid,count,value;
    require(bool(in>>sigma>>valid>>count) && sigma<32 && valid<2,"bad vector header");
    std::array<Byte,Candidate<R>::N> nodes{};
    std::array<bool,256> seen{};
    for(auto& x:nodes) {require(bool(in>>value) && value>0 && value<256 && !seen[value],"bad vector nodes"); seen[value]=true; x=Byte(value);}
    Candidate<R> expected; expected.valid=valid!=0; expected.count=count;
    for(auto& x:expected.g) {require(bool(in>>value) && value<256,"bad g"); x=Byte(value);}
    for(auto& x:expected.h) {require(bool(in>>value) && value<256,"bad h"); x=Byte(value);}
    require(trial<R,false>(nodes,Byte(sigma),t)==expected,"primary vector mismatch");
    require(trial<R,true>(nodes,Byte(sigma),t)==expected,"replay vector mismatch");
}
std::uint64_t positive_integer(const std::string& text) {
    require(!text.empty() && text.find_first_not_of("0123456789")==std::string::npos,"bad positive integer");
    auto n=std::stoull(text); require(n>0,"zero sample count"); return n;
}
} // namespace
int main(int argc,char** argv) try {
    auto t=std::make_unique<Tables>();
    if(argc==3 && std::string(argv[1])=="check") {
        for(unsigned a=0;a<256;++a) for(unsigned b=0;b<256;++b)
            require(t->mul[a][b]==reference_multiply(Byte(a),Byte(b)),"field cross-check");
        std::ifstream in(argv[2]); require(bool(in),"vectors missing");
        unsigned weight=0,count=0;
        while(in>>weight) {
            if(weight==40) check_vector<1>(in,*t); else if(weight==42) check_vector<2>(in,*t);
            else throw std::runtime_error("bad vector weight");
            ++count;
        }
        require(in.eof() && count>0,"bad or empty vectors");
        std::cout<<"PASS field table and "<<count<<" Python vectors through both kernels\n"; return 0;
    }
    require(argc==7,"usage: run|replay WEIGHT N TAPE OUTPUT fixed-bch256-higher-endpoint-v1; or check VECTORS");
    require(std::string(argv[6])=="fixed-bch256-higher-endpoint-v1","protocol mismatch");
    const bool replay=std::string(argv[1])=="replay";
    require(replay || std::string(argv[1])=="run","bad mode");
    const auto weight=positive_integer(argv[2]),n=positive_integer(argv[3]);
    require(weight==40 || weight==42,"unsupported weight");
    if(weight==40) {if(replay) experiment<1,true>(n,argv[4],argv[5],*t); else experiment<1,false>(n,argv[4],argv[5],*t);}
    else {if(replay) experiment<2,true>(n,argv[4],argv[5],*t); else experiment<2,false>(n,argv[4],argv[5],*t);}
    return 0;
} catch(const std::exception& e) {std::cerr<<"ERROR: "<<e.what()<<'\n'; return 1;}
