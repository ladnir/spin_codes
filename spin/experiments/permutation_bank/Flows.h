#pragma once
#include "Routing.h"
#include "../feistel/Direct.h"

namespace spin::experimental::bank {
namespace kernel=detail::kernel;
inline constexpr unsigned MessageSize=1U<<18,CodeSize=2*MessageSize;

// Flow 1: fixed code, complete routing schedules and masks, persistent scratch.
// An empty route selects the original shuffle family. A supplied route permits
// comparison with exactly the same map as flow 2.
class PrecomputedCode {
    kernel::Spin code_;
    kernel::Spin::Workspace scratch_;
public:
    PrecomputedCode(std::uint64_t routeSeed,std::uint64_t maskSeed,std::vector<std::uint32_t> route={})
        :code_(kernel::Configuration::T128S19,kernel::MessageLength{MessageSize},routeSeed,maskSeed,256,kernel::BchBackend::Auto,true,std::move(route)) {
        code_.compact();scratch_=kernel::Spin::Workspace(code_);
    }
    void encode(kernel::block* input){code_.encodeInplace(input,CodeSize,scratch_);}
    std::size_t planBytes()const{return code_.setupBytes();}
    std::size_t workspaceBytes()const{return scratch_.bytes();}
};

// Flow 2: reusable bank, fresh instance parameters and masks, no full route.
// A returned FreshCode borrows the bank. The BankFlow must outlive the code;
// construction from a temporary BankFlow is deliberately disabled.
template<unsigned Count> class BankFlow {
    Tables<Count> tables_;
public:
    using FreshCode=feistel::Direct<0,false,Routing<Count>>;
    explicit BankFlow(std::uint64_t bankSeed):tables_(bankSeed) {}
    BankFlow(const BankFlow&)=delete;
    BankFlow& operator=(const BankFlow&)=delete;
    BankFlow(BankFlow&&)=delete;
    BankFlow& operator=(BankFlow&&)=delete;
    FreshCode fresh(std::uint64_t routeSeed,std::uint64_t maskSeed)const & {
        return FreshCode(Routing<Count>(tables_,routeSeed),maskSeed);
    }
    FreshCode fresh(std::uint64_t,std::uint64_t)const &&=delete;
    std::vector<std::uint32_t> materialize(std::uint64_t seed)const {
        Routing<Count> routing(tables_,seed);std::vector<std::uint32_t> out(CodeSize);
        for(unsigned i=0;i<CodeSize;i+=16)routing.outerBatch(i,out.data()+i);
        return out;
    }
    std::size_t bytes()const{return tables_.bytes();}
};
}
