#include "Feistel.h"
#include <spin/Code.h>
#include "../../src/kernels/Spin.h"
#include <cstring>
#include <iostream>
#include <stdexcept>
using namespace spin::experimental::feistel;
namespace k=spin::detail::kernel;
static void check(bool b,const char* why){if(!b)throw std::runtime_error(why);}
template<unsigned Bits,unsigned Rounds> void permutation() {
    for(unsigned seed=0;seed<16;++seed) {
        kernelWords words(seed);Permutation<Bits,Rounds> p;p.fill(words);
        std::vector<bool> seen(1U<<Bits);
        for(unsigned x=0;x<seen.size();++x) {
            const auto y=p.forward(x);
            check(y<seen.size() && !seen[y],"Feistel bijection");seen[y]=true;
            check(p.inverse(y)==x && p.forward(p.inverse(x))==x,"Feistel inverse");
        }
    }
}
template<unsigned Bits,unsigned Rounds> void routing() {
    Routing<Bits,Rounds> routing(17);auto route=routing.materialize();
    std::array<unsigned,16> batch;
    for(unsigned i=0;i<route.size();++i) {
        if((i&15)==0)routing.outerBatch(i,batch.data());
        check(route[i]==batch[i&15],"batched route differs");
        check(route[i]==routing.outer(i),"materialized route differs from direct evaluation");
        check(routing.inner(route[i])==i,"route inverse");
    }
    const auto size=route.size(),message=size/2;
    const auto config=message==65536?k::Configuration::T64S12R2:k::Configuration::T128S19;
    k::Spin code(config,k::MessageLength{message},17,29,256,k::BchBackend::Auto,false,std::move(route));
    code.validateSetup();
    k::Spin::Workspace work(code);
    std::vector<k::block> input(size),actual(message),expected(message);
    kernelWords words(183);for(auto& x:input){auto a=words(),b=words();x=k::block(a,b);}
    code.reference(input.data(),expected.data());
    code.compact();code.encode(input.data(),size,actual.data(),message,work);
    check(!std::memcmp(actual.data(),expected.data(),message*sizeof(k::block)),"Feistel transpose dense reference");
}
int main(){try {
    if(!spin::capabilities().avx2)return 77;
    permutation<6,4>();permutation<7,4>();permutation<8,4>();permutation<9,4>();
    permutation<10,6>();permutation<11,6>();permutation<12,8>();permutation<13,8>();
    permutation<8,6>();permutation<8,8>();permutation<11,4>();permutation<11,8>();
    permutation<8,2>();permutation<11,2>();permutation<16,6>();permutation<17,6>();
    routing<9,4>();routing<9,6>();routing<9,8>();
    routing<11,4>();routing<11,6>();routing<11,8>();
    std::cout<<"Feistel bijections, inverses, route structure, and dense encoding references PASS\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
