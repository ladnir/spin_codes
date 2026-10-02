#include "K16Paired15.h"
#include "../../../src/kernels/SetupRandom.h"
#include <array>
#include <stdexcept>
namespace spin::research::k16codesign {
namespace {
bool fullRank15(std::array<std::uint16_t,16> rows) {
    unsigned pivot=0;
    for(unsigned column=0;column<15;++column) {
        unsigned selected=pivot;
        while(selected<15 && !((rows[selected]>>column)&1U))++selected;
        if(selected==15)continue;
        std::swap(rows[pivot],rows[selected]);
        for(unsigned row=pivot+1;row<15;++row)
            if((rows[row]>>column)&1U)rows[row]^=rows[pivot];
        ++pivot;
    }
    return pivot==15;
}
}
void customizePaired15(rs::Plan& plan,std::uint64_t seed,PairedTables& tables,PairedOptimizedTables& optimized) {
    if(plan.variant!=rs::Variant::Rs16Gf16)throw std::invalid_argument("Paired15 requires RS16x8");
    detail::kernel::setup::Words words(seed^0x3f625a92ULL);
    for(auto& rows:plan.reverseMatrices) {
        do {
            for(unsigned row=0;row<15;++row)rows[row]=std::uint16_t(words()&0x7fffU);
            rows[15]=0x8000U;
        } while(!fullRank15(rows));
    }
    preparePairedBasis(plan,tables,optimized);
}
}
