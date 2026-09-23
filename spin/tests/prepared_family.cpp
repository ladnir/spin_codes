#include <spin/PreparedEncoder.h>
#include "../experiments/permutation_bank/Composed.h"
#include <iostream>

int main() {try {
    using namespace spin::experimental::bank;
    const ComposedBank<1> prototype(913);
    spin::detail::BankState current({1U<<18,spin::Parameters::T128S19,17,29},913);
    for(unsigned seed=17;seed<25;++seed) {
        ComposedRouting<1,true,true,true,true> route(prototype,seed);
        current.refresh(seed,seed+12);
        for(unsigned i=0;i<CodeSize;++i)
            if(current.destination(i/2048,i%2048)!=route.outer(i))
                throw std::runtime_error("prepared routing differs from the measured row-rotate1 family");
        if(current.masks!=spin::experimental::feistel::makeMasksK18(seed+12))
            throw std::runtime_error("prepared masks differ from the measured family");
    }
    std::cout<<"Prepared family matches K18 prototype: eight seeds, complete routes and masks PASS\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
