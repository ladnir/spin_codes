#include <spin/Code.h>
#include <array>
#include <iostream>
#include <vector>
struct alignas(16) Word {std::array<std::uint64_t,2> lanes;};
static std::uint64_t digest(const std::vector<Word>& words) {
    std::uint64_t h=0xcbf29ce484222325ULL;
    for(auto w:words)for(auto x:w.lanes)for(unsigned b=0;b<8;++b) {
        h^=(x>>(8*b))&255;h*=0x100000001b3ULL;
    }
    return h;
}
int main(int argc,char**) {
    if(!spin::capabilities().avx2)return 77;
    // Descriptor v1, seeds (17,29), minimum natural K, xorshift input from 139.
    // FNV-1a over little-endian output bytes is a regression checksum, not a hash API.
    constexpr std::uint64_t expectedForward[3]={0x4015ca039758b7c5ULL,0x7eed51ec660021d1ULL,0xb05d1190e0c95cbdULL};
    constexpr std::uint64_t expectedTranspose[3]={0x67b806197faf2bd5ULL,0xb6343401871868e1ULL,0xd99ec72636d79ee3ULL};
    unsigned test=0;
    for(auto p:{spin::Parameters::T128S19,spin::Parameters::T64S12,spin::Parameters::T64S12R2}) {
        spin::Code c({spin::message_alignment(p),p,17,29});auto w=c.make_workspace();
        std::vector<Word> in(c.code_size()),out(c.code_size()),transposed(c.message_size());
        std::uint64_t s=139;
        for(auto& v:in)for(auto& lane:v.lanes) {s^=s<<13;s^=s>>7;s^=s<<17;lane=s;}
        c.forward<Word>(std::span<const Word>(in).first(c.message_size()),out,w);
        c.transpose<Word>(in,transposed,w);
        const auto f=digest(out),t=digest(transposed);
        if(argc>1) std::cout<<std::hex<<f<<" "<<t<<'\n';
        else if(f!=expectedForward[test] || t!=expectedTranspose[test]) {
            std::cerr<<"descriptor v1 known-answer mismatch: "<<unsigned(p)<<'\n';return 1;
        }
        ++test;
    }
}
