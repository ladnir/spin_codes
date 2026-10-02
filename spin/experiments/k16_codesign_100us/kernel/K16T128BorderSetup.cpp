#include "K16T128Border.h"
#include "../../../src/kernels/SetupRandom.h"
#include <algorithm>
#include <stdexcept>
namespace spin::research::k16codesign {
namespace {
bool fullRank(std::array<std::uint32_t,20> rows,unsigned s) {
    for(unsigned column=0;column<s;++column) {
        unsigned pivot=column;
        while(pivot<s && !(rows[pivot]&(1U<<column)))++pivot;
        if(pivot==s)return false;
        std::swap(rows[pivot],rows[column]);
        for(unsigned row=column+1;row<s;++row)
            if(rows[row]&(1U<<column))rows[row]^=rows[column];
    }
    return true;
}
}
void prepareT128Border(const rs::Plan& plan,std::uint64_t seed,T128BorderTables& tables,unsigned s) {
    if((s!=19 && s!=20) || plan.n%128)
        throw std::invalid_argument("t128border requires s19/s20 and whole 128-output steps");
    tables.stateBits=s;
    tables.reverseMatrices.resize(plan.n/128);
    tables.packedUpdates.resize(plan.n/128);
    detail::kernel::setup::Words words(seed^0x3f625a92ULL);
    const auto mask=(1U<<s)-1;
    for(std::size_t epoch=0;epoch<plan.n/128;++epoch) {
        std::array<std::uint32_t,20> forward{};
        do {for(unsigned row=0;row<s;++row)forward[row]=static_cast<std::uint32_t>(words())&mask;}
        while(!fullRank(forward,s));
        auto& reverse=tables.reverseMatrices[epoch];reverse.fill(0);
        for(unsigned row=0;row<s;++row)
            for(unsigned col=0;col<s;++col)reverse[col]|=((forward[row]>>col)&1U)<<row;
        auto& compact=tables.packedUpdates[epoch];compact.fill(0);
        for(unsigned out=0;out<3;++out)
            for(unsigned in=0;in<3;++in)
                for(unsigned j=0;j<8 && 8*out+j<s;++j)
                    compact[3*out+in]|=std::uint64_t((reverse[8*out+j]>>(8*in))&255U)<<(8*(7-j));
    }
}
}
