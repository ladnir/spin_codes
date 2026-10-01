#pragma once
#include <algorithm>
#include <bit>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace spin::research {
// Execution schedule only: inverse is the already sampled, checked route.
// Pairing never changes which inner position reaches which BCH coordinate.
struct TwoBitTiled {
    unsigned tileBlocks;
    unsigned bucketStride;
    std::vector<unsigned> pairSlots;
    std::vector<unsigned> offsets;
    std::vector<unsigned> gatherSources;

    explicit TwoBitTiled(std::span<const unsigned> inverse, unsigned tileRows, bool padded=false, bool gather=false)
        :tileBlocks(0),bucketStride(0) {
        if(inverse.empty() || inverse.size()%1024 ||
           inverse.size()>std::numeric_limits<unsigned>::max() ||
           tileRows<4 || !std::has_single_bit(tileRows))
            throw std::invalid_argument("invalid two-bit tile geometry");
        tileBlocks=256*std::min(tileRows,unsigned(inverse.size()/256));
        bucketStride=tileBlocks+(padded?4:0);
        if(inverse.size()%tileBlocks)
            throw std::invalid_argument("two-bit tiles must divide the route");
        if((inverse.size()/tileBlocks)*bucketStride>std::numeric_limits<unsigned>::max())
            throw std::invalid_argument("two-bit bucket offsets exceed 32 bits");
        pairSlots.resize(inverse.size()/2);
        offsets.resize(inverse.size());
        if(gather)gatherSources.resize(inverse.size());
        std::vector<unsigned> remaining(inverse.size()/tileBlocks,tileBlocks);
        std::vector<bool> seen(inverse.size());
        for(std::size_t i=inverse.size();i;) {
            i-=2;
            const auto x=inverse[i],y=inverse[i+1];
            if(x>=inverse.size() || y>=inverse.size() || x==y || seen[x] || seen[y] ||
               x/512!=y/512 || x/256==y/256)
                throw std::invalid_argument("invalid two-row packet");
            seen[x]=seen[y]=true;
            const auto bucket=x/tileBlocks;
            if(remaining[bucket]<2)throw std::invalid_argument("tile overflow");
            remaining[bucket]-=2;
            const auto slot=bucket*tileBlocks+remaining[bucket];
            pairSlots[i/2]=bucket*bucketStride+remaining[bucket];
            offsets[slot]=packed(x%tileBlocks);
            offsets[slot+1]=packed(y%tileBlocks);
            if(gather) {
                gatherSources[bucket*tileBlocks+offsets[slot]]=remaining[bucket];
                gatherSources[bucket*tileBlocks+offsets[slot+1]]=remaining[bucket]+1;
            }
        }
        for(auto count:remaining)if(count)throw std::invalid_argument("incomplete tile");
        if(gather)std::vector<unsigned>().swap(offsets);
    }

    static unsigned packed(unsigned x) noexcept {
        return (x&~1023U)|((x&255U)<<2)|((x>>8)&3U);
    }
};
}
