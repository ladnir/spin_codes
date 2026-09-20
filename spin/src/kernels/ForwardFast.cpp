#include <algorithm>
#include "Spin.h"
#include "Inner.h"
#include <cstring>
namespace spin::detail::kernel {
void bchForward4(const block*,block*);
static SPIN_FORCEINLINE u32 unpack(const u8* p) {u32 v;std::memcpy(&v,p,4);return v&0xffffff;}
template<class Map,bool Packed> void Spin::runForwardFour(const block* in,block* out,Workspace& w) const {
    block* values=w.buckets.data(); block* tile=w.tile.data();
    const auto tileSize=tileBlocks(),n=codeBlocks();
    for(std::size_t base=0;base<n;base+=tileSize) {
        for(std::size_t j=0;j<tileSize;j+=1024)
            bchForward4(in+base/2+j/2,tile+j);
        for(std::size_t j=0;j<tileSize;++j) {
            const auto offset=Packed?unpack(mBchOffsets24.data()+3*(base+j)):mBchOffsets32[base+j];
            values[base+j]=tile[offset];
        }
    }
    innerForward<Map>(n,mForwardFieldRows.data(),[&](std::size_t i) {
        if constexpr(Packed) return values[unpack(mSlots24.data()+3*i)];
        else return values[mSlots32[i]];
    },out);
}
template<class Map,bool Packed> void Spin::runForwardFourTail(const block* in,block* out,Workspace& w) const {
    block* values=w.buckets.data(); block* tile=w.tile.data();
    const auto tileSize=tileBlocks(),n=codeBlocks();
    for(std::size_t base=0;base<n;base+=tileSize) {
        const auto activeSize=std::min(tileSize,n-base);
        for(std::size_t j=0;j<activeSize;j+=1024)
            bchForward4(in+base/2+j/2,tile+j);
        for(std::size_t j=0;j<activeSize;++j) {
            const auto offset=Packed?unpack(mBchOffsets24.data()+3*(base+j)):mBchOffsets32[base+j];
            values[base+j]=tile[offset];
        }
    }
    innerForward<Map>(n,mForwardFieldRows.data(),[&](std::size_t i) {
        if constexpr(Packed) return values[unpack(mSlots24.data()+3*i)];
        else return values[mSlots32[i]];
    },out);
}
template void Spin::runForwardFour<Map128S19,true>(const block*,block*,Workspace&) const;
template void Spin::runForwardFour<Map128S19,false>(const block*,block*,Workspace&) const;
template void Spin::runForwardFour<Map64S12,true>(const block*,block*,Workspace&) const;
template void Spin::runForwardFour<Map64S12,false>(const block*,block*,Workspace&) const;
template void Spin::runForwardFour<Map64S12R2,true>(const block*,block*,Workspace&) const;
template void Spin::runForwardFour<Map64S12R2,false>(const block*,block*,Workspace&) const;
template<class Map> void Spin::runForwardDirect(const block* in,block* out,Workspace& w) const {
    auto* values=w.buckets.data();const auto n=codeBlocks();
    for(std::size_t j=0;j<n;j+=1024) bchForward4(in+j/2,values+j);
    innerForward<Map>(n,mForwardFieldRows.data(),[&](std::size_t i) {
        return values[mForwardDirect[i]];
    },out);
}
template void Spin::runForwardDirect<Map128S19>(const block*,block*,Workspace&) const;
template void Spin::runForwardDirect<Map64S12>(const block*,block*,Workspace&) const;
template void Spin::runForwardDirect<Map64S12R2>(const block*,block*,Workspace&) const;
template void Spin::runForwardFourTail<Map128S19,true>(const block*,block*,Workspace&) const;
template void Spin::runForwardFourTail<Map128S19,false>(const block*,block*,Workspace&) const;
template void Spin::runForwardFourTail<Map64S12,true>(const block*,block*,Workspace&) const;
template void Spin::runForwardFourTail<Map64S12,false>(const block*,block*,Workspace&) const;
template void Spin::runForwardFourTail<Map64S12R2,true>(const block*,block*,Workspace&) const;
template void Spin::runForwardFourTail<Map64S12R2,false>(const block*,block*,Workspace&) const;
}
