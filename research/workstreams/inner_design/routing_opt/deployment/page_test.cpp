#include "WorkspaceRouting.h"
#include <array>
#include <cerrno>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace bare_spin::workspace_routing;
static void require(bool ok) {if(!ok) throw std::runtime_error("page helper test failed");}
int main() {
    static_assert(wholePages(4096,4096,4096).bytes==4096);
    static_assert(wholePages(4097,8192,4096).begin==8192);
    static_assert(wholePages(4097,8192,4096).bytes==4096);
    static_assert(wholePages(4097,4096,4096).bytes==0);
    static_assert(wholePages(4096,0,4096).bytes==0);
    static_assert(wholePages(4096,4096,0).bytes==0);
    static_assert(wholePages(4096,4096,3).bytes==0);
    static_assert(wholePages(std::numeric_limits<std::uintptr_t>::max()-10,100,4096).bytes==0);
    require(!eligible(true,1U<<16));require(!eligible(false,1U<<20));
    require(eligible(true,1U<<18)==compiled);require(eligible(true,1U<<20)==compiled);
    std::vector<unsigned char> data(8U<<20,0x5a);
    errno=E2BIG;
    const auto advice=adviseOwned(data.data(),data.size());
    require(errno==E2BIG);
    for(auto value:data) require(value==0x5a);
    if(!compiled) require(advice.bytes==0 && !advice.hinted && !advice.collapsed);
#ifdef SPIN_WORKSPACE_ROUTING_ADVICE_FAILURE_TEST
    require(advice.bytes>0 && !advice.hinted && !advice.collapsed);
#endif
    // Tiny unaligned owned span must not cause advice on the surrounding stack.
    std::array<unsigned char,17> small{};
    require(adviseOwned(small.data()+1,15).bytes==0);
    std::cout<<"page helper PASS compiled="<<compiled<<" hint="<<advice.hinted<<" collapse="<<advice.collapsed<<'\n';
}
