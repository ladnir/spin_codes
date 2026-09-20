#include "AsymmetricMap.h"
#include "Weight5Inner.h"
#include <array>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace bare_spin;
static u64 randomWord(u64& s) {
    s+=0x9e3779b97f4a7c15ULL;
    auto x=s;x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;
    x=(x^(x>>27))*0x94d049bb133111ebULL;return x^(x>>31);
}
template<bool Masked> static void check(unsigned epochs,u64 seed,bool impulses) {
    const unsigned n=128*epochs;
    std::vector<block> x(n),u(n),y(n),v(n);
    std::vector<u32> masks(2*epochs);
    for(auto& value:x) {auto lo=randomWord(seed),hi=randomWord(seed);value=block(hi,lo);}
    for(auto& value:u) {auto lo=randomWord(seed),hi=randomWord(seed);value=block(hi,lo);}
    if(impulses) {
        std::fill(x.begin(),x.end(),block(0,0));
        for(unsigned p:{0U,127U,128U,n-1}) if(p<n) x[p]=block(1,3);
    }
    for(unsigned e=0;e<epochs;++e) {
        u32 a;do a=randomWord(seed)&((1U<<19)-1);while(!a);
        u32 b=randomWord(seed)&((1U<<19)-1);
        if(std::popcount(a&b)&1)b^=a&-a;
        masks[2*e]=a;masks[2*e+1]=b;
    }
    std::array<block,19> state{},next{};
    for(auto& z:state)z=block(0,0);
    // Deliberately scalar forward oracle: raw A/B columns, no zeta or XOR circuit.
    for(unsigned e=0;e<epochs;++e) {
        block dot(0,0);
        for(unsigned j=0;j<19;++j)if((masks[2*e+1]>>j)&1)dot^=state[j];
        for(unsigned j=0;j<19;++j)next[j]=state[j]^(((masks[2*e]>>j)&1)?dot:block(0,0));
        for(unsigned p=0;p<128;++p) {
            auto value=x[128*e+p];
            for(unsigned j=0;j<19;++j) {
                if((AsymmetricMap::feedbackColumns[p]>>j)&1)value^=state[j];
                if((AsymmetricMap::columns[p]>>j)&1)next[j]^=x[128*e+p];
            }
            y[128*e+p]=value;
        }
        state=next;
    }
    asymmetricReverse<AsymmetricMap,true,Masked>(u.data(),n,masks.data(),[&](std::size_t p,block value){v[p]=value;});
    block lhs(0,0),rhs(0,0);
    for(unsigned p=0;p<n;++p) {
        lhs^=block(_mm_and_si128(y[p].mData,u[p].mData));
        rhs^=block(_mm_and_si128(x[p].mData,v[p].mData));
    }
    if(lhs!=rhs)throw std::runtime_error("forward/transpose identity failed");
}
int main() {
    try {
        for(unsigned e:{1U,2U,3U,8U})for(unsigned seed=1;seed<=12;++seed)for(bool impulses:{false,true}) {
            check<false>(e,seed,impulses);check<true>(e,seed,impulses);
        }
        std::cout<<"inner_forward_transpose_identity=PASS\n";
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
