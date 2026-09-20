#include "Spin.h"
#include <algorithm>
#include <cstring>
#include <iostream>
#include <stdexcept>
using namespace bare_spin;
static void check(bool v,const char* s) {if(!v) throw std::runtime_error(s);}
static bool equal(const std::vector<block>& a,const std::vector<block>& b) {return !std::memcmp(a.data(),b.data(),a.size()*16);}
static __m128i dot(const std::vector<block>& a,const std::vector<block>& b) {
 auto x=_mm_setzero_si128();for(std::size_t i=0;i<a.size();++i) x=_mm_xor_si128(x,_mm_and_si128(a[i].mData,b[i].mData));return x;
}
int main() {
 try {
  for(unsigned m:{14U,16U,18U,20U}) {
   Spin oracle(Configuration::T128S19,m,17,29,0,BchBackend::Avx2);
   std::vector<block> x(oracle.messageBlocks()),q(oracle.codeBlocks()),y(q.size()),z(x.size()),ey(q.size()),ez(x.size());
   u64 seed=123;for(auto* a:{&x,&q}) for(auto& v:*a) {const auto lo=splitmix(seed),hi=splitmix(seed);v=block(hi,lo);}
   oracle.forwardReference(x.data(),ey.data());oracle.reference(q.data(),ez.data());
   for(auto layout:{Layout::Packed24,Layout::Indices32}) for(unsigned tile:m>=16?std::vector<unsigned>{0,2,4,16}:std::vector<unsigned>{0})
    for(auto backend:{BchBackend::Auto,BchBackend::Avx2}) {
     Spin code(Configuration::T128S19,m,17,29,tile,backend);code.validateSetup();
     if(layout==Layout::Packed24) {
      const auto view=code.wideView();
      Spin legacy(Configuration::T128S19,m,17,29,tile,BchBackend::Avx2);
      const auto old=legacy.wideView();
      check(view.n==old.n && view.tile==old.tile &&
            !std::memcmp(view.offsets,old.offsets,3*view.n) &&
            !std::memcmp(view.slots,old.slots,3*view.n) &&
            !std::memcmp(view.fieldRows,old.fieldRows,2*(view.n/128)*sizeof(u32)),"WideView changed");
      for(std::size_t i=0;i<view.n;++i) {
       const u32 offset=u32(view.offsets[3*i])|(u32(view.offsets[3*i+1])<<8)|(u32(view.offsets[3*i+2])<<16);
       check(offset<view.tile,"WideView offset out of range");
      }
     }
     code.compact(layout);Spin::Workspace w(code);
     code.forward(x.data(),x.size(),y.data(),y.size(),w,layout);code.encode(q.data(),q.size(),z.data(),z.size(),w,layout);
     check(equal(y,ey)&&equal(z,ez),"forward/transpose oracle mismatch");
     check(_mm_movemask_epi8(_mm_cmpeq_epi8(dot(y,q),dot(x,z)))==0xffff,"adjoint failed");
     auto inplace=q;code.encodeInplace(inplace.data(),inplace.size(),w,layout);
     check(!std::memcmp(inplace.data(),ez.data(),ez.size()*16),"inplace output failed");
     check(!std::memcmp(inplace.data()+ez.size(),q.data()+ez.size(),ez.size()*16),"inplace suffix failed");
     if(m>=16 && tile==0) {
      std::vector<u64> bits(x.size()/64),out(q.size()/64),scratch(out.size());
      std::vector<block> bx(x.size()),by(q.size());
      for(auto& v:bits) v=splitmix(seed);
      for(std::size_t i=0;i<bx.size();++i) bx[i]=block(0,(bits[i/64]>>(i%64))&1);
      code.forwardBits(bits.data(),bits.size(),out.data(),out.size(),scratch.data(),scratch.size());
      code.forward(bx.data(),bx.size(),by.data(),by.size(),w,layout);
      for(std::size_t i=0;i<by.size();++i) check((u64(_mm_cvtsi128_si64(by[i].mData))&1)==((out[i/64]>>(i%64))&1),"forwardBits failed");
     }
    }
   std::cout<<"m="<<m<<" forward, transpose, compact, adjoint PASS\n";
  }
  for(auto cfg:{Configuration::T64S16,Configuration::T64S20,Configuration::T256S14}) {
   Spin code(cfg,16);check(code.bchBackend()==BchBackend::Avx2,"legacy config selected new backend");
   code.validateSetup();Spin::Workspace w(code);
   std::vector<block> in(code.codeBlocks()),out(code.messageBlocks()),expected(out.size());u64 seed=7;
   for(auto& v:in) {auto lo=splitmix(seed),hi=splitmix(seed);v=block(hi,lo);}
   code.reference(in.data(),expected.data());code.encode(in.data(),in.size(),out.data(),out.size(),w);
   check(equal(out,expected),"legacy config changed");
  }
  bool rejected=false;
  try {Spin bad(Configuration::T128S19,16,1,2,1);} catch(const std::invalid_argument&) {rejected=true;}
  check(rejected,"one-row tile must be rejected");
  rejected=false;
  try {Spin bad(Configuration::T128S19,16,1,2,2,BchBackend::Avx512);} catch(const std::invalid_argument&) {rejected=true;}
  check(rejected,"four-row kernel accepted a two-row tile");
 } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
