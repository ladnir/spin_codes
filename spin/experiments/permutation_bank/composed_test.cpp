#include "Composed.h"
#include <bitset>
#include <iostream>
using namespace spin::experimental::bank;
static void require(bool b){if(!b)throw std::runtime_error("composed test failed");}
template<unsigned Count,bool ShiftOnly,bool RowShift,bool RowAffine=false,bool RowRotate=false> void test() {
    ComposedBank<Count> bank(913);
    std::array<std::bitset<256>,2> columns;
    // Bank-adaptive supports, fixed before the instance is sampled. Choose
    // the columns in rows 0 and 1 assigned to the same 32 stored regions.
    for(unsigned g=0;g<32;++g)for(unsigned p=0;p<2048;++p) {
        const unsigned x=unpackDestination(bank.routes[g*Count*2048+p]);
        if((x>>8)<2)columns[x>>8].set(x&255);
    }
    require(columns[0].count()==32 && columns[1].count()==32);
    unsigned minOverlap=256,maxOverlap=0,total=0;
    for(unsigned seed=17;seed<49;++seed) {
        ComposedRouting<Count,ShiftOnly,RowShift,RowAffine,RowRotate> route(bank,seed);
        std::array<std::bitset<256>,2> regions;
        std::vector<bool> seen(CodeSize);
        std::array<unsigned,2048> plain,packed;
        for(unsigned g=0;g<256;++g) {
            route.outerRegion(g,plain.data());route.outerRegion(g,packed.data(),true);
            if constexpr(ShiftOnly && RowShift) {
                std::array<unsigned,2048> bytes;
                route.outerRegionBytes(g,bytes.data(),true);
                for(unsigned i=0;i<2048;++i)require(bytes[i]==16*packed[i]);
                route.outerRegionBytes(g,bytes.data(),false);
                for(unsigned i=0;i<2048;++i)require(bytes[i]==16*plain[i]);
            }
            if constexpr(ShiftOnly) {
                std::array<unsigned,128> epoch;
                for(unsigned offset=0;offset<2048;offset+=128) {
                    route.outerEpoch(g,offset,epoch.data(),true);
                    require(std::equal(epoch.begin(),epoch.end(),packed.begin()+offset));
                    route.outerEpoch(g,offset,epoch.data(),false);
                    require(std::equal(epoch.begin(),epoch.end(),plain.begin()+offset));
                }
                auto checkChunk=[&]<unsigned Size>() {
                    std::array<unsigned,Size> chunk;
                    for(unsigned offset=0;offset<2048;offset+=Size) {
                        route.template outerChunk<Size>(g,offset,chunk.data(),true);
                        require(std::equal(chunk.begin(),chunk.end(),packed.begin()+offset));
                        route.template outerChunk<Size>(g,offset,chunk.data(),false);
                        require(std::equal(chunk.begin(),chunk.end(),plain.begin()+offset));
                    }
                };
                checkChunk.template operator()<512>();checkChunk.template operator()<1024>();
            }
            std::bitset<2048> rows;
            for(unsigned p=0;p<2048;++p) {
                const unsigned x=route.outer(g*2048+p);
                require(x<CodeSize && !seen[x] && !rows[x>>8]);seen[x]=true;rows.set(x>>8);
                require(plain[p]==x && packed[p]==packDestination(x));
                if((x>>8)<2 && columns[x>>8][x&255])regions[x>>8].set(g);
            }
        }
        require(regions[0].count()==32 && regions[1].count()==32);
        const unsigned overlap=unsigned((regions[0]&regions[1]).count());
        if constexpr(!RowShift)require(overlap==32);
        minOverlap=std::min(minOverlap,overlap);maxOverlap=std::max(maxOverlap,overlap);total+=overlap;
    }
    if constexpr(RowShift)require(minOverlap<32); // Removes the exact invariant, not a security test.
    std::cout<<"count="<<Count<<" shift_only="<<ShiftOnly<<" row_shift="<<RowShift<<" row_affine="<<RowAffine<<" row_rotate="<<RowRotate
      <<" correctness PASS; adaptive overlap min="<<minOverlap<<" mean="<<double(total)/32<<" max="<<maxOverlap<<'\n';
}
int main(){try {
    test<1,false,false>();test<4,true,false>();test<4,true,true>();test<1,true,true>();test<1,true,true,true>();
    test<1,true,true,true,true>();
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
