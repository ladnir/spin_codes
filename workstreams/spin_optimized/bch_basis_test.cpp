#include "Spin.h"
#include "generated/BchCircuit.h"
#include <array>
#include <cstring>
#include <iostream>
#include <stdexcept>
namespace bare_spin {void bchTranspose4(const block*,block*);}
using namespace bare_spin;
static bool same(const block& a,const block& b) {return std::memcmp(&a,&b,sizeof(a))==0;}
int main() {
 if(!bchAvx512Available()) {std::cout<<"AVX-512 unavailable: skip\n";return 77;}
 try {
  alignas(64) std::array<block,1024> input;
  alignas(64) std::array<block,514> output;
  const block zero(0,0),probe(0xabcdef0123456789ULL,0x123456789abcdef0ULL),guard(17,29);
  auto checkGuards=[&] {if(!same(output.front(),guard)||!same(output.back(),guard)) throw std::runtime_error("BCH output canary");};
  for(unsigned lane=0;lane<4;++lane) for(unsigned bit=0;bit<256;++bit) {
   input.fill(zero);output.fill(guard);input[4*bit+lane]=probe;
   bchTranspose4(input.data(),output.data()+1);checkGuards();
   for(unsigned l=0;l<4;++l) for(unsigned j=0;j<128;++j) {
    const auto expected=l==lane && ((BchRows[j][bit/64]>>(bit%64))&1)?probe:zero;
    if(!same(output[1+128*l+j],expected)) throw std::runtime_error("BCH basis/lane mismatch");
   }
  }
  u64 seed=19;
  for(unsigned repeat=0;repeat<12;++repeat) {
   for(auto& v:input) {auto lo=splitmix(seed),hi=splitmix(seed);v=block(hi,lo);}
   const auto original=input;output.fill(guard);
   bchTranspose4(input.data(),output.data()+1);checkGuards();
   if(std::memcmp(input.data(),original.data(),sizeof(input))) throw std::runtime_error("BCH changed input");
   for(unsigned l=0;l<4;++l) for(unsigned j=0;j<128;++j) {
    block expected=zero;
    for(unsigned bit=0;bit<256;++bit) if((BchRows[j][bit/64]>>(bit%64))&1) expected^=input[4*bit+l];
    if(!same(output[1+128*l+j],expected)) throw std::runtime_error("BCH dense mismatch");
   }
  }
  std::cout<<"1024 basis coordinates, 12 dense inputs, input preservation, output canaries PASS\n";
 } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
