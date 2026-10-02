#include "Schedule.h"
#include <iostream>
#include <stdexcept>
using namespace bare_spin;

template<bool Packed,unsigned U,unsigned PF,unsigned Write>
void check(std::size_t n) {
    constexpr unsigned base=3;
    std::vector<block> input(n+base), actual(n+2), expected(n+2);
    std::vector<u8> p24(3*(n+base)+4);
    std::vector<u32> p32(n+base);
    const block guard(0x88776655,0x44332211);
    actual.front()=actual.back()=guard;
    expected.front()=expected.back()=guard;
    for(std::size_t j=0;j<n;++j) {
        input[base+j]=block(j+19,j+7);
        const auto offset=u32(n-1-j);
        p32[base+j]=offset;
        auto p=p24.data()+3*(base+j);
        p[0]=u8(offset);p[1]=u8(offset>>8);p[2]=u8(offset>>16);
        expected[1+offset]=input[base+j];
    }
    schedule::scatter<Packed,U,PF,Write>(input.data(),actual.data()+1,p24.data(),p32.data(),base,n);
    if(actual!=expected) throw std::runtime_error("scatter or canary mismatch");
}

int main() {
    for(unsigned n=0;n<=513;++n) {
        check<true,4,128,0>(n);check<false,4,128,0>(n);
        check<true,8,128,0>(n);check<false,8,128,0>(n);
        check<true,4,32,0>(n);check<false,4,32,0>(n);
        check<true,8,0,1>(n);check<false,8,0,1>(n);
    }
    std::cout<<"schedule tails 0..513, packed/32-bit, nonzero base, canaries PASS\n";
}
