"""Extract the production packed encoder into a small compiler reproducer.

Usage: python spin/tools/probe_packed.py OUTPUT.cpp
Compile with /std:c++20 /EHsc /arch:AVX2 /Ispin/src (MSVC), or
-std=c++20 -mavx2 -Ispin/src (GCC/Clang). No library link is needed.
"""
from pathlib import Path
import sys

package = Path(__file__).resolve().parents[1]
source = (package / 'src/kernels/Spin.cpp').read_text()
start = source.index('template<class Map> void Spin::forwardBitsMap(')
body = source[start:source.rindex('\n}')].replace('Spin::forwardBitsMap', 'Probe::forwardBitsMap')
if '--constexpr' in sys.argv:
    body = body.replace('static const auto', 'static constexpr auto')
prefix = r'''
#include "kernels/Inner.h"
#include "kernels/generated/BchCircuit.h"
#include <cstdio>
#include <cstdlib>
#include <numeric>
#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
static LONG WINAPI crash(EXCEPTION_POINTERS* p) {
    const auto base=reinterpret_cast<std::uintptr_t>(GetModuleHandleA(nullptr));
    std::fprintf(stderr,"exception=%lx pc=%llx module_offset=%llx address=%llx\n",
        p->ExceptionRecord->ExceptionCode,p->ContextRecord->Rip,p->ContextRecord->Rip-base,
        p->ExceptionRecord->ExceptionInformation[1]);
    return EXCEPTION_EXECUTE_HANDLER;
}
#endif
using namespace spin::detail::kernel;
struct Probe {
    std::size_t mK=16384;
    std::vector<u32> mRoute,mForwardFieldRows;
    std::size_t codeBlocks() const {return 2*mK;}
    template<class Map> void forwardBitsMap(const u64*,u64*,u64*) const;
};
'''
suffix = r'''
static u64 rng(u64& s) {s^=s<<13;s^=s>>7;s^=s<<17;return s;}
template<class Map> static int test() {
    Probe p;
    const auto n=p.codeBlocks();
    constexpr unsigned rounds=isTwoRoundMap<Map>?2:1;
    u64 seed=18;
    p.mRoute.resize(n);std::iota(p.mRoute.begin(),p.mRoute.end(),0);
    for(std::size_t i=n-1;i>0;--i) std::swap(p.mRoute[i],p.mRoute[rng(seed)%(i+1)]);
    p.mForwardFieldRows.resize(2*rounds*n/Map::T);
    for(std::size_t i=0;i<p.mForwardFieldRows.size();i+=2) {
        auto u=u32(rng(seed))&((1U<<Map::S)-1),v=u32(rng(seed))&((1U<<Map::S)-1);
        if(!u) u=1;
        if(std::popcount(u&v)&1) v^=u&-u;
        p.mForwardFieldRows[i]=u;p.mForwardFieldRows[i+1]=v;
    }
    std::vector<u64> in(p.mK/64),out(n/64),scratch(n/64),outer(n/64);
    for(auto& x:in)x=rng(seed);
    std::printf("start T=%u S=%u rounds=%u\n",Map::T,Map::S,rounds);
    p.forwardBitsMap<Map>(in.data(),out.data(),scratch.data());
    for(std::size_t row=0;row<p.mK/128;++row)
        for(unsigned j=0;j<128;++j) if((in[2*row+j/64]>>(j%64))&1)
            for(unsigned w=0;w<4;++w) outer[4*row+w]^=BchRows[j][w];
    if(outer!=scratch) {std::puts("BCH mismatch");return 1;}
    u32 state=0;
    for(std::size_t base=0;base<n;base+=Map::T) {
        u32 syndrome=0;
        for(unsigned j=0;j<Map::T;++j) {
            const auto a=p.mRoute[base+j];
            const unsigned x=unsigned((outer[a/64]>>(a%64))&1);
            const unsigned y=x^(std::popcount(state&Map::columns[j])&1);
            const unsigned actual=unsigned((out[(base+j)/64]>>((base+j)%64))&1);
            if(x) syndrome^=Map::feedbackColumns[j];
            if(y!=actual) {
                std::printf("mismatch bit=%zu state=%x expected=%u actual=%u\n",base+j,state,y,actual);
                return 2;
            }
        }
        for(unsigned r=0;r<rounds;++r) {
            const auto i=2*(rounds*(base/Map::T)+r);
            if(std::popcount(state&p.mForwardFieldRows[i+1])&1) state^=p.mForwardFieldRows[i];
        }
        state^=syndrome;
    }
    std::puts("PASS");return 0;
}
int main() {
    std::setvbuf(stdout,nullptr,_IONBF,0);
#ifdef _WIN32
    SetUnhandledExceptionFilter(crash);
#endif
    int result=test<Map128S19>();
    result|=test<Map64S12>();result|=test<Map64S12R2>();
    return result;
}
'''
Path(sys.argv[1]).write_text(prefix + body + suffix)
