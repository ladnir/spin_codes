#include "Spin.h"
#include <algorithm>
#include <iostream>
#include <stdexcept>
using namespace bare_spin;
static void require(bool ok,const char* what) {if(!ok) throw std::runtime_error(what);}
int main() {
 try {
  const bool available=bchAvx512Available();
  for(unsigned m:{16U,18U,20U}) {
   Spin ref(Configuration::T128S19,m,17,29,0,Outer::Bch256x128,BchBackend::Avx2);
   Spin::Workspace rw(ref);
   std::vector<block> input(ref.codeBlocks()),expected(ref.messageBlocks()),out(expected.size());
   u64 seed=3;
   for(auto& v:input) {const auto lo=splitmix(seed),hi=splitmix(seed);v=block(hi,lo);}
   ref.encode(input.data(),input.size(),expected.data(),expected.size(),rw);
   // Cover direct K18 and tiled K20 with tiny tiles and alternate schedules.
   for(unsigned tile: m<=18?std::vector<unsigned>{0,2,4,8,16,512}:std::vector<unsigned>{0,2,4,512,1024,4096})
    for(auto layout:{Layout::Packed24,Layout::Indices32})
     for(auto backend:{BchBackend::Auto,BchBackend::Avx2,BchBackend::Avx512}) {
      const bool canFour=available && tile!=2;
      if(backend==BchBackend::Avx512 && !canFour) {
       bool rejected=false;
       try {Spin bad(Configuration::T128S19,m,17,29,tile,Outer::Bch256x128,backend);}
       catch(const std::invalid_argument&) {rejected=true;}
       require(rejected,"unavailable forced backend accepted");continue;
      }
      Spin code(Configuration::T128S19,m,17,29,tile,Outer::Bch256x128,backend);
      require(code.bchBackend()==((canFour && backend!=BchBackend::Avx2)?BchBackend::Avx512:BchBackend::Avx2),"wrong backend");
      require(code.routeHash()==ref.routeHash(),"backend changed route");
      code.validateSetup();code.compact(layout);Spin::Workspace work(code);
      code.encode(input.data(),input.size(),out.data(),out.size(),work,layout);
      require(out==expected,"backend or tile changed output");
      auto inplace=input;code.encodeInplace(inplace.data(),inplace.size(),work,layout);
      require(std::equal(expected.begin(),expected.end(),inplace.begin()),"inplace output mismatch");
      require(std::equal(input.begin()+expected.size(),input.end(),inplace.begin()+expected.size()),"suffix changed");
     }
  }
  bool rejected=false;
  try {Spin bad(Configuration::T128S19,16,1,2,0,Outer::Bch256x128,static_cast<BchBackend>(99));}
  catch(const std::invalid_argument&) {rejected=true;}
  require(rejected,"unknown backend accepted");
  std::cout<<"backend equivalence PASS; AVX512 available="<<available<<'\n';
 } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
