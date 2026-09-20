#include "Spin.h"
#include "BinaryRaa.h"
#include "libOTe/Tools/Ppcg/ChosenBlockBaa.h"
#include "libOTe/Tools/ExConvCodeOld/ExConvCodeOld.h"
#include <algorithm>
#include <chrono>
#include <cctype>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>
#ifdef __linux__
#include <fcntl.h>
#include <sched.h>
#include <sys/file.h>
#include <unistd.h>
#endif
using namespace osuCrypto;
using Clock = std::chrono::steady_clock;
static volatile u64 sink;
static void guard() {
#ifdef __linux__
    for (const auto& e: std::filesystem::directory_iterator("/proc")) {
        auto name=e.path().filename().string();
        if (name.empty() || !std::all_of(name.begin(),name.end(),[](unsigned char c){return std::isdigit(c);}) ||
            std::stoul(name)==static_cast<unsigned long>(getpid())) continue;
        std::error_code ec;
        auto file=std::filesystem::read_symlink(e.path()/"exe",ec).filename().string();
        if (ec) continue;
        std::transform(file.begin(),file.end(),file.begin(),[](unsigned char c){return char(std::tolower(c));});
        if (file.find("bench")!=std::string::npos || file.find("frontend_libote")!=std::string::npos)
            throw std::runtime_error("another benchmark is active: "+file);
    }
#endif
}
static void lockAndPin() {
#ifdef __linux__
    static int fd=open("/tmp/bare-spin-benchmark.lock",O_CREAT|O_RDWR,0600);
    if(fd<0 || flock(fd,LOCK_EX|LOCK_NB)) throw std::runtime_error("benchmark lock unavailable");
    cpu_set_t cpus; CPU_ZERO(&cpus); CPU_SET(15,&cpus);
    if(sched_setaffinity(0,sizeof(cpus),&cpus)) throw std::runtime_error("cannot pin CPU 15");
#endif
    guard();
}
static void fill(std::vector<block>& v) {
    PRNG p(block(231234,321312)); p.get(v.data(), v.size());
}
static u64 hash(span<const block> v) {
    u64 h=0xcbf29ce484222325ULL;
    for(auto b:v) for(unsigned j=0;j<2;++j) h=(h^b.get<u64>(j))*0x100000001b3ULL;
    return h;
}
static void equal(span<const block> a, span<const block> b) {
    if(a.size()!=b.size() || !std::equal(a.begin(),a.end(),b.begin()))
        throw std::runtime_error("correctness mismatch");
}
template<unsigned R> static void checkRaa(u64 k) {
    comparison::BinaryRaa<R> code(k);
    std::vector<block> x(k*R), y(k), gx(k), gty(k*R); fill(x); fill(y);
    code.encode(x.data(),gx.data()); code.forward(y.data(),gty.data());
    block a=ZeroBlock,b=ZeroBlock;
    for(u64 i=0;i<k;++i) a^=gx[i]&y[i];
    for(u64 i=0;i<k*R;++i) b^=x[i]&gty[i];
    if(a!=b) throw std::runtime_error("RAA transpose identity failed");
    std::fill(x.begin(),x.end(),ZeroBlock); code.encode(x.data(),gx.data());
    if(std::any_of(gx.begin(),gx.end(),[](block z){return z!=ZeroBlock;})) throw std::runtime_error("RAA zero failed");
}
static void checks() {
    for(u64 k:{1,2,3,7,16,31,64,257}) { checkRaa<3>(k); checkRaa<4>(k); }
    for(auto c:{ChosenBlockCode::Golay24,ChosenBlockCode::Rm32}) {
        const u64 n=chosenBlockCodeParams(c).m*17;
        ChosenBlockBaa code(c,n,block(34123421,2134123),4,true);
        std::vector<block> x(n),a(code.outputSize()),b(a.size()); fill(x);
        code.encode(x,a); code.encodeReference(x,b); equal(a,b);
    }
    // Legacy EC: 24 random taps and one fixed tap give paper memory m=25.
    // Independent scalar tap traversal, then a scalar expansion in the
    // upstream random-index stream order. Small sizes do not refill the PRNG.
    for(u64 n:{64,128,256}) {
        ExConvCodeOld ec; ec.config(n/2,n,33,24,false);
        std::vector<block> input(n), ref, out(n/2), expected(n/2,ZeroBlock); fill(input); ref=input;
        PRNG taps(ec.mSeed ^ OneBlock);
        const auto bytes=reinterpret_cast<const u8*>(taps.mBuffer.data());
        for(u64 i=0;i<n;++i) {
            for(unsigned j=0;j<24 && i+j+1<n;++j)
                if((bytes[3*i+j/8]>>(7-j%8))&1) ref[i+j+1]^=ref[i];
            if(i+25<n) ref[i+25]^=ref[i];
        }
        detail::ExpanderModd indices(ec.mExpander.mSeed,n);
        for(u64 i=0;i<n/2;++i)
            for(unsigned j=0;j<33;++j)
                expected[i]^=ref[indices.get()];
        ec.dualEncode<block>(input,out); equal(out,expected); equal(input,ref);
    }
    bare_spin::Spin spin(bare_spin::Configuration::T128S19,16);
    std::vector<block> x(spin.codeBlocks()),a(spin.messageBlocks()),b(a.size()); fill(x);
    bare_spin::Spin::Workspace w(spin); spin.reference(x.data(),b.data());
    spin.compact(); spin.encodeUnchecked(x.data(),a.data(),w); equal(a,b);
    std::cout << "All adapter, binary-lane, reference, and transpose checks passed.\n";
}
template<class Encode> static void measure(const std::string& mode,unsigned m,unsigned trials,
    std::vector<block>& input,span<block> output, Encode&& encode) {
    std::vector<double> samples; samples.reserve(trials);
    const auto run=[&] {
        // Native buffers are reused without restoring or copying input.
        auto start=Clock::now(); encode();
        auto ms=std::chrono::duration<double,std::milli>(Clock::now()-start).count();
        sink=output[0].get<u64>(0); return ms;
    };
    guard();
    for(unsigned i=0;i<3;++i) run();
    for(unsigned i=0;i<trials;++i) samples.push_back(run());
    guard();
    std::cout<<std::setprecision(12)<<"{\"mode\":\""<<mode<<"\",\"m\":"<<m<<",\"n\":"<<input.size()
        <<",\"k\":"<<output.size()<<",\"warmups\":3,\"cpu\":15,\"input_policy\":\"no_reset\",\"samples_ms\":[";
    for(unsigned i=0;i<trials;++i) std::cout<<(i?",":"")<<samples[i];
    std::cout<<"],\"output_hash\":\""<<std::hex<<hash(output)<<std::dec<<"\"}\n";
}
template<unsigned R> static void raa(const std::string& mode,unsigned m,unsigned trials) {
    const u64 k=1ull<<m; comparison::BinaryRaa<R> code(k);
    std::vector<block> input(k*R),out(k); fill(input);
    measure(mode,m,trials,input,out,[&]{code.encode(input.data(),out.data());});
}
int main(int argc,char** argv) { try {
    if(argc==2 && std::string(argv[1])=="check") { checks(); return 0; }
    if(argc!=4) throw std::invalid_argument("usage: transposed_benchmark check | MODE M ODD_TRIALS");
    std::string mode=argv[1]; unsigned m=std::stoul(argv[2]),trials=std::stoul(argv[3]);
    if((m!=16 && m!=18 && m!=20) || trials<3 || !(trials&1)) throw std::invalid_argument("invalid size or trials");
    lockAndPin(); const u64 k=1ull<<m;
    if(mode=="spin") {
        bare_spin::Spin code(bare_spin::Configuration::T128S19,m); code.compact();
        bare_spin::Spin::Workspace w(code);
        std::vector<block> input(2*k),out(k); fill(input);
        measure(mode,m,trials,input,out,[&]{code.encodeUnchecked(input.data(),out.data(),w);});
    } else if(mode=="golay" || mode=="rm") {
        const u64 n=(2*k/192)*192;
        ChosenBlockBaa code(mode=="golay"?ChosenBlockCode::Golay24:ChosenBlockCode::Rm32,n,block(34123421,2134123),4,true);
        std::vector<block> input(n),out(code.outputSize()); fill(input);
        measure(mode,m,trials,input,out,[&]{code.encode(input,out,1);});
    } else if(mode=="exconv") {
        ExConvCodeOld code; code.config(k,2*k,33,24,false);
        std::vector<block> input(2*k),out(k); fill(input);
        measure(mode,m,trials,input,out,[&]{code.dualEncode<block>(input,out);});
    } else if(mode=="raa3") raa<3>(mode,m,trials);
    else if(mode=="raa4") raa<4>(mode,m,trials);
    else throw std::invalid_argument("unknown mode");
    return 0;
} catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; } }
