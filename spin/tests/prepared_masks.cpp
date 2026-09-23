#include <spin/Code.h>
#include <spin/detail/BankState.h>
#include "../src/Cpu.h"
#include "../src/kernels/BankKernel.h"
#include "../src/kernels/SetupRandom.h"
#include <bit>
#include <iostream>
#include <stdexcept>
#include <vector>

int main(){try {
    if(!spin::capabilities().avx2)return 77;
    unsigned rejected=0,vectorized=0;
    for(auto p:{spin::Parameters::T128S19,spin::Parameters::T64S12,spin::Parameters::T64S12R2}) {
        const auto k=p==spin::Parameters::T128S19?1U<<18:1U<<16;
        spin::detail::BankState bank({k,p,17,29},913);
        const unsigned bits=p==spin::Parameters::T128S19?19:12,mask=(1U<<bits)-1;
        std::vector<std::uint32_t> expected(bank.masks.size()),fast(expected.size());
        for(unsigned seed=0;seed<2048;++seed) {
            spin::detail::kernel::setup::Words words(seed);
            for(std::size_t i=0;i<expected.size();i+=2) {
                std::uint32_t u;do{u=std::uint32_t(words())&mask;}while(!u);
                auto v=std::uint32_t(words())&mask;
                if(std::popcount(u&v)&1)v^=u&-u;
                expected[i]=u;expected[i+1]=v;
            }
            bank.refresh(seed+17,seed);
            if(bank.masks!=expected)throw std::runtime_error("prepared mask stream changed");
#if SPIN_BCH_AVX512
            if(spin::detail::cpu_mask512()) {
                if(spin::detail::kernel::bankMasks512(fast.data(),fast.size(),bits,seed)) {
                    ++vectorized;
                    if(fast!=expected)throw std::runtime_error("SIMD mask stream changed");
                } else ++rejected;
            }
#endif
        }
    }
    if(spin::detail::cpu_mask512() && (!rejected || !vectorized))
        throw std::runtime_error("missing SIMD or rejection coverage");
    std::cout<<"6144 exact prepared mask streams PASS; vectorized="<<vectorized<<" rejected="<<rejected<<'\n';
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
