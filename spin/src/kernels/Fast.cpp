#include <algorithm>
#include "Spin.h"
#include "Inner.h"
#include "K16Inner.h"
#include "K16R2Inner.h"
#include <cstring>
namespace spin::detail::kernel {
void bchTranspose4(const block*,block*);
static SPIN_FORCEINLINE u32 unpack(const u8* p) {u32 v;std::memcpy(&v,p,4);return v&0xffffff;}
template<bool Packed> void Spin::runFour(const block* in,block* out,Workspace& w) const {
    using Map=Map128S19;
    block* values=w.buckets.data();
    innerReverse<Map>(in,codeBlocks(),mFieldRows.data(),[&](std::size_t i,block v) {
        if constexpr(Packed) values[unpack(mSlots24.data()+3*i)]=v;
        else values[mSlots32[i]]=v;
    });
    block* tile=w.tile.data();
    const auto tileSize=tileBlocks(),n=codeBlocks();
    for(std::size_t base=0;base<n;base+=tileSize) {
        for(std::size_t j=0;j<tileSize;++j) {
            if(j+32<tileSize) {
                const auto future=Packed?unpack(mBchOffsets24.data()+3*(base+j+32)):mBchOffsets32[base+j+32];
                _mm_prefetch(reinterpret_cast<const char*>(tile+future),_MM_HINT_T0);
            }
            const auto offset=Packed?unpack(mBchOffsets24.data()+3*(base+j)):mBchOffsets32[base+j];
            tile[offset]=values[base+j];
        }
        for(std::size_t j=0;j<tileSize;j+=1024)
            bchTranspose4(tile+j,out+base/2+j/2);
    }
}
template<bool Packed> void Spin::runFourTail(const block* in,block* out,Workspace& w) const {
    using Map=Map128S19;
    block* values=w.buckets.data();
    innerReverse<Map>(in,codeBlocks(),mFieldRows.data(),[&](std::size_t i,block v) {
        if constexpr(Packed) values[unpack(mSlots24.data()+3*i)]=v;
        else values[mSlots32[i]]=v;
    });
    block* tile=w.tile.data();
    const auto tileSize=tileBlocks(),n=codeBlocks();
    for(std::size_t base=0;base<n;base+=tileSize) {
        const auto activeSize=std::min(tileSize,n-base);
        for(std::size_t j=0;j<activeSize;++j) {
            if(j+32<activeSize) {
                const auto future=Packed?unpack(mBchOffsets24.data()+3*(base+j+32)):mBchOffsets32[base+j+32];
                _mm_prefetch(reinterpret_cast<const char*>(tile+future),_MM_HINT_T0);
            }
            const auto offset=Packed?unpack(mBchOffsets24.data()+3*(base+j)):mBchOffsets32[base+j];
            tile[offset]=values[base+j];
        }
        for(std::size_t j=0;j<activeSize;j+=1024)
            bchTranspose4(tile+j,out+base/2+j/2);
    }
}

template void Spin::runFour<true>(const block*,block*,Workspace&) const;
template void Spin::runFour<false>(const block*,block*,Workspace&) const;
void Spin::runSmall(const block* in,block* out,Workspace& w) const {
        using Map=Map128S19;
        block* direct=w.buckets.data();
        innerReverse<Map>(in,131072,mFieldRows.data(),[&](std::size_t i,block v) {
            direct[mSmallRoute32[i]]=v;
        });
        for(std::size_t j=0;j<131072;j+=1024) bchTranspose4(direct+j,out+j/2);
        return;
    }
template<bool Packed> void Spin::runK16Four(const block* in,block* out,Workspace& w) const {
    using Map=Map64S12;
    block* values=w.buckets.data();
    inner64Reverse(in,codeBlocks(),mFieldRows.data(),[&](std::size_t i,block v) {
        if constexpr(Packed) values[unpack(mSlots24.data()+3*i)]=v;
        else values[mSlots32[i]]=v;
    });
    block* tile=w.tile.data();
    const auto tileSize=tileBlocks(),n=codeBlocks();
    for(std::size_t base=0;base<n;base+=tileSize) {
        for(std::size_t j=0;j<tileSize;++j) {
            if(j+32<tileSize) {
                const auto future=Packed?unpack(mBchOffsets24.data()+3*(base+j+32)):mBchOffsets32[base+j+32];
                _mm_prefetch(reinterpret_cast<const char*>(tile+future),_MM_HINT_T0);
            }
            const auto offset=Packed?unpack(mBchOffsets24.data()+3*(base+j)):mBchOffsets32[base+j];
            tile[offset]=values[base+j];
        }
        for(std::size_t j=0;j<tileSize;j+=1024)
            bchTranspose4(tile+j,out+base/2+j/2);
    }
}
template<bool Packed> void Spin::runK16FourTail(const block* in,block* out,Workspace& w) const {
    using Map=Map64S12;
    block* values=w.buckets.data();
    inner64Reverse(in,codeBlocks(),mFieldRows.data(),[&](std::size_t i,block v) {
        if constexpr(Packed) values[unpack(mSlots24.data()+3*i)]=v;
        else values[mSlots32[i]]=v;
    });
    block* tile=w.tile.data();
    const auto tileSize=tileBlocks(),n=codeBlocks();
    for(std::size_t base=0;base<n;base+=tileSize) {
        const auto activeSize=std::min(tileSize,n-base);
        for(std::size_t j=0;j<activeSize;++j) {
            if(j+32<activeSize) {
                const auto future=Packed?unpack(mBchOffsets24.data()+3*(base+j+32)):mBchOffsets32[base+j+32];
                _mm_prefetch(reinterpret_cast<const char*>(tile+future),_MM_HINT_T0);
            }
            const auto offset=Packed?unpack(mBchOffsets24.data()+3*(base+j)):mBchOffsets32[base+j];
            tile[offset]=values[base+j];
        }
        for(std::size_t j=0;j<activeSize;j+=1024)
            bchTranspose4(tile+j,out+base/2+j/2);
    }
}


template void Spin::runK16Four<true>(const block*,block*,Workspace&) const;
template void Spin::runK16Four<false>(const block*,block*,Workspace&) const;
void Spin::runSmallK16(const block* in,block* out,Workspace& w) const {
        using Map=Map64S12;
        block* direct=w.buckets.data();
        inner64Reverse(in,131072,mFieldRows.data(),[&](std::size_t i,block v) {
            direct[mSmallRoute32[i]]=v;
        });
        for(std::size_t j=0;j<131072;j+=1024) bchTranspose4(direct+j,out+j/2);
        return;
    }
template<bool Packed> void Spin::runK16R2Four(const block* in,block* out,Workspace& w) const {
    using Map=Map64S12;
    block* values=w.buckets.data();
    inner64R2Reverse(in,codeBlocks(),mFieldRows.data(),[&](std::size_t i,block v) {
        if constexpr(Packed) values[unpack(mSlots24.data()+3*i)]=v;
        else values[mSlots32[i]]=v;
    });
    block* tile=w.tile.data();
    const auto tileSize=tileBlocks(),n=codeBlocks();
    for(std::size_t base=0;base<n;base+=tileSize) {
        for(std::size_t j=0;j<tileSize;++j) {
            if(j+32<tileSize) {
                const auto future=Packed?unpack(mBchOffsets24.data()+3*(base+j+32)):mBchOffsets32[base+j+32];
                _mm_prefetch(reinterpret_cast<const char*>(tile+future),_MM_HINT_T0);
            }
            const auto offset=Packed?unpack(mBchOffsets24.data()+3*(base+j)):mBchOffsets32[base+j];
            tile[offset]=values[base+j];
        }
        for(std::size_t j=0;j<tileSize;j+=1024)
            bchTranspose4(tile+j,out+base/2+j/2);
    }
}
template<bool Packed> void Spin::runK16R2FourTail(const block* in,block* out,Workspace& w) const {
    using Map=Map64S12;
    block* values=w.buckets.data();
    inner64R2Reverse(in,codeBlocks(),mFieldRows.data(),[&](std::size_t i,block v) {
        if constexpr(Packed) values[unpack(mSlots24.data()+3*i)]=v;
        else values[mSlots32[i]]=v;
    });
    block* tile=w.tile.data();
    const auto tileSize=tileBlocks(),n=codeBlocks();
    for(std::size_t base=0;base<n;base+=tileSize) {
        const auto activeSize=std::min(tileSize,n-base);
        for(std::size_t j=0;j<activeSize;++j) {
            if(j+32<activeSize) {
                const auto future=Packed?unpack(mBchOffsets24.data()+3*(base+j+32)):mBchOffsets32[base+j+32];
                _mm_prefetch(reinterpret_cast<const char*>(tile+future),_MM_HINT_T0);
            }
            const auto offset=Packed?unpack(mBchOffsets24.data()+3*(base+j)):mBchOffsets32[base+j];
            tile[offset]=values[base+j];
        }
        for(std::size_t j=0;j<activeSize;j+=1024)
            bchTranspose4(tile+j,out+base/2+j/2);
    }
}


template void Spin::runK16R2Four<true>(const block*,block*,Workspace&) const;
template void Spin::runK16R2Four<false>(const block*,block*,Workspace&) const;
void Spin::runSmallK16R2(const block* in,block* out,Workspace& w) const {
        using Map=Map64S12;
        block* direct=w.buckets.data();
        inner64R2Reverse(in,131072,mFieldRows.data(),[&](std::size_t i,block v) {
            direct[mSmallRoute32[i]]=v;
        });
        for(std::size_t j=0;j<131072;j+=1024) bchTranspose4(direct+j,out+j/2);
        return;
    }
void Spin::runK18(const block* in,block* out,Workspace& w) const {
    using Map=Map128S19;
    block* direct=w.buckets.data();
    innerReverse<Map>(in,524288,mFieldRows.data(),[&](std::size_t i,block v) {
        direct[mK18Route32[i]]=v;
    });
    for(std::size_t j=0;j<524288;j+=1024) bchTranspose4(direct+j,out+j/2);
}
template void Spin::runFourTail<true>(const block*,block*,Workspace&) const;
template void Spin::runFourTail<false>(const block*,block*,Workspace&) const;
template void Spin::runK16FourTail<true>(const block*,block*,Workspace&) const;
template void Spin::runK16FourTail<false>(const block*,block*,Workspace&) const;
template void Spin::runK16R2FourTail<true>(const block*,block*,Workspace&) const;
template void Spin::runK16R2FourTail<false>(const block*,block*,Workspace&) const;
// Range-direct transpose from the standalone target; compile-time map selection.
template<class Map> void Spin::runRange(const block* in,block* out,Workspace& w) const {
    block* direct=w.buckets.data();
    const auto n=codeBlocks();
    const auto* route=mForwardDirect.empty()?mRangeRoute32.data():mForwardDirect.data();
    auto emit=[&](std::size_t i,block v) {direct[route[i]]=v;};
    if constexpr(std::is_same_v<Map,Map64S12>) inner64Reverse(in,n,mFieldRows.data(),emit);
    else if constexpr(std::is_same_v<Map,Map64S12R2>) inner64R2Reverse(in,n,mFieldRows.data(),emit);
    else innerReverse<Map>(in,n,mFieldRows.data(),emit);
    for(std::size_t j=0;j<n;j+=1024) bchTranspose4(direct+j,out+j/2);
}
template void Spin::runRange<Map128S19>(const block*,block*,Workspace&) const;
template void Spin::runRange<Map64S12>(const block*,block*,Workspace&) const;
template void Spin::runRange<Map64S12R2>(const block*,block*,Workspace&) const;
}
