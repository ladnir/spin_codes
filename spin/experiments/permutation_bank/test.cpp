#include "Bank.h"
#include "Routing.h"
#include <iostream>
using namespace spin::experimental::bank;
static void require(bool value){if(!value)throw std::runtime_error("bank test failed");}
template<unsigned Count> void routingTest() {
    Tables<Count> tables(913);Routing<Count> routing(tables,17);Routing<Count,true> vector(tables,17);
    Bank rows(8,Count,913^0x726f7773ULL),regions(11,Count,913^0x72656773ULL);
    Words rw(17^0x726f7773ULL),gw(17^0x72656773ULL);
    std::vector<Parameters> rp(2048),gp(256);
    for(auto& p:rp)p=rows.sample(rw);for(auto& p:gp)p=regions.sample(gw);
    std::array<unsigned,16> batch;
    std::vector<bool> seen(Routing<Count>::Size);
    for(unsigned i=0;i<Routing<Count>::Size;++i) {
        if(!(i&15))routing.outerBatch(i,batch.data());
        const unsigned g=i>>11,row=regions.backward<Mode::RotateAdd>(i&2047,gp[g]);
        const unsigned expected=row*256+rows.backward<Mode::RotateAdd>(g,rp[row]);
        require(expected==routing.outer(i) && expected==batch[i&15] && !seen[expected]);seen[expected]=true;
    }
    std::array<unsigned,2048> scalar,simd;
    for(bool packed:{false,true})for(unsigned region=0;region<256;++region) {
        routing.outerRegion(region,scalar.data(),packed);vector.outerRegion(region,simd.data(),packed);
        for(unsigned i=0;i<2048;++i) {
            const unsigned x=routing.outer(region*2048+i);
            const unsigned expected=packed?((x&~1023U)|((x&255U)<<2)|((x>>8)&3U)):x;
            require(scalar[i]==expected && simd[i]==expected);
        }
    }
}
template<Mode M> void check(unsigned bits) {
    Bank bank(bits,4,12345);Words words(6789);
    for(unsigned trial=0;trial<16;++trial) {
        const auto k=bank.sample(words);std::vector<bool> seen(bank.n);
        for(unsigned x=0;x<bank.n;++x) {
            const auto y=bank.forward<M>(x,k);
            require(y<bank.n && !seen[y]);seen[y]=true;
            require(bank.backward<M>(y,k)==x);
        }
    }
}
int main(){try {
    routingTest<16>();routingTest<64>();
    for(unsigned bits:{8U,11U,13U}){check<Mode::Xor>(bits);check<Mode::XorAdd>(bits);check<Mode::RotateAdd>(bits);check<Mode::InputRotateAdd>(bits);}
    Bank bank(8,1,123);Words words(456);
    for(unsigned trial=0;trial<8;++trial) {
        const auto k=bank.sample(words);
        std::vector<unsigned> expected(8),observed(8);
        for(unsigned x=0;x<256;++x) {
            ++expected[std::countr_zero(bank.p[x]^bank.p[x^128])];
            ++observed[std::countr_zero(bank.forward<Mode::XorAdd>(x,k)^bank.forward<Mode::XorAdd>(x^128,k))];
            for(unsigned y=0;y<256;++y)if(x!=y) {
                const unsigned v=std::countr_zero(x^y);
                require(unsigned(std::countr_zero(bank.in<Mode::XorAdd>(x,k)^bank.in<Mode::XorAdd>(y,k)))==v);
                require(unsigned(std::countr_zero(bank.out<Mode::XorAdd>(x,k)^bank.out<Mode::XorAdd>(y,k)))==v);
            }
        }
        require(expected==observed);
        // Input/output XOR translations only relabel the per-permutation XOR
        // difference histogram. Test every difference, not just sampled ones.
        for(unsigned dx=1;dx<256;++dx) {
            std::vector<unsigned> base(256),wrapped(256);
            for(unsigned x=0;x<256;++x) {
                ++base[bank.p[x]^bank.p[x^dx]];
                ++wrapped[bank.forward<Mode::Xor>(x,k)^bank.forward<Mode::Xor>(x^dx,k)];
            }
            require(base==wrapped);
        }
    }
    std::cout<<"bank bijection/inverse, XOR histogram, XOR+add valuation identities PASS\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
