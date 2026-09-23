#include "Masks.h"
#include <iostream>
using namespace spin::experimental::bank;
int main(){
    unsigned rejected=0;
    for(unsigned seed=0;seed<2048;++seed) {
        bool fallback=false;const auto masks=fastMasks(seed,&fallback);
        if(masks!=spin::experimental::feistel::makeMasksK18(seed))return 1;
        rejected+=fallback;
    }
#if defined(__AVX512DQ__) && defined(__AVX512VPOPCNTDQ__) && defined(__AVX512VL__)
    if(!rejected)return 2;
#endif
    std::cout<<"2048 exact mask streams PASS; rejection fallbacks="<<rejected<<'\n';
}
