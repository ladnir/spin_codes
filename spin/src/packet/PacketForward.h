#pragma once
#include "PacketPlan.h"
namespace spin::detail::packet {
// Prepared once. These are binary adjoints of the same immutable RS packet map.
struct ForwardPlan {
    std::vector<BorderRow> updates;
    std::vector<std::array<std::uint64_t,9>> outer;
    std::vector<std::uint32_t> inverseRoute;
    explicit ForwardPlan(const Plan&);
    std::size_t setupBytes() const noexcept {
        return updates.capacity()*sizeof(updates[0])+outer.capacity()*sizeof(outer[0])
            +inverseRoute.capacity()*sizeof(inverseRoute[0]);
    }
};
void forwardOuter(const Block*,Block*,const Plan&,const ForwardPlan&);
void forwardInner(const Block*,Block*,const Plan&,const ForwardPlan&);
void forwardOuterGather(const Block*,Block*,const Plan&,const ForwardPlan&);
void forwardInnerGather(const Block*,Block*,const Plan&,const ForwardPlan&);
inline void forwardPacket(const Block* in,Block* out,Block* scratch,
                          const Plan& p,const ForwardPlan& f) {
    forwardOuter(in,scratch,p,f);
    forwardInner(scratch,out,p,f);
}
// Measured size-aware path. Caller supplies disjoint message/output/scratch;
// scratch has 64-byte alignment; input/output need only 16-byte alignment.
inline void forwardSelected(const Block* in,Block* out,Block* scratch,
                            const Plan& p,const ForwardPlan& f) {
    if(p.k<=262144) {
        forwardOuterGather(in,scratch,p,f);forwardInnerGather(scratch,out,p,f);
    } else if((reinterpret_cast<std::uintptr_t>(out)&63)==0) {
        forwardOuter(in,out,p,f);forwardInner(out,out,p,f);
    } else {
        // Scatter requires complete aligned cache lines. Foreign outputs retain
        // the public 16-byte contract by using the already-owned workspace.
        forwardOuter(in,scratch,p,f);forwardInner(scratch,out,p,f);
    }
}
}
