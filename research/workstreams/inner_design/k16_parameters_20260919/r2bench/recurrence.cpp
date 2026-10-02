#include "Spin.h"
#include "Inner.h"
#include "K16Inner.h"
#include <array>
#include <cstring>
#include <iostream>
using namespace bare_spin;
int main() {
    using Map=Map64S12;
    constexpr unsigned epochs=5,n=epochs*Map::T;
    for(u64 seed0=1;seed0<=20;++seed0) {
        u64 seed=seed0;
        std::array<block,n> x,y,q,z;
        std::array<u32,4*epochs> masks;
        for(auto& a:x) {auto lo=splitmix(seed);a=block(splitmix(seed),lo);}
        for(auto& a:q) {auto lo=splitmix(seed);a=block(splitmix(seed),lo);}
        for(unsigned i=0;i<masks.size();i+=2) {
            u32 u;do u=splitmix(seed)&4095;while(!u);
            u32 v=splitmix(seed)&4095;
            if(std::popcount(u&v)&1) v^=u&-u;
            masks[i]=u;masks[i+1]=v;
        }
        std::array<block,12> state,feedback;
        for(auto& a:state) a=block(0,0);
        for(unsigned e=0;e<epochs;++e) {
            for(auto& a:feedback) a=block(0,0);
            for(unsigned p=0;p<64;++p) {
                auto a=x[64*e+p];
                for(unsigned j=0;j<12;++j) {
                    if((Map::columns[p]>>j)&1) a^=state[j];
                    if((Map::feedbackColumns[p]>>j)&1) feedback[j]^=x[64*e+p];
                }
                y[64*e+p]=a;
            }
            for(unsigned r=0;r<2;++r) {
                auto dot=block(0,0);
                for(unsigned j=0;j<12;++j) if((masks[4*e+2*r+1]>>j)&1) dot^=state[j];
                for(unsigned j=0;j<12;++j) if((masks[4*e+2*r]>>j)&1) state[j]^=dot;
            }
            for(unsigned j=0;j<12;++j) state[j]^=feedback[j];
        }
        inner64Reverse(q.data(),n,masks.data(),[&](std::size_t i,block a){z[i]=a;});
        auto lhs=_mm_setzero_si128(),rhs=lhs;
        for(unsigned i=0;i<n;++i) {
            lhs=_mm_xor_si128(lhs,_mm_and_si128(x[i].mData,z[i].mData));
            rhs=_mm_xor_si128(rhs,_mm_and_si128(y[i].mData,q[i].mData));
        }
        if(_mm_movemask_epi8(_mm_cmpeq_epi8(lhs,rhs))!=65535) return 1;
    }
    std::cout<<"Two-round independent forward/optimized transpose adjoint: 20 setups PASS\n";
}
