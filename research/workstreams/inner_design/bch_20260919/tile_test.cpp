#include "Spin.h"
#include <iostream>
#include <stdexcept>
#include <algorithm>
using namespace bare_spin;
int main() {
    try {
        Spin ref(Configuration::T128S19,16,17,29,4);
        std::vector<block> in(ref.codeBlocks()),expected(ref.messageBlocks()),out(expected.size());
        u64 seed=3;
        for(auto& x:in) {const auto lo=splitmix(seed),hi=splitmix(seed);x=block(hi,lo);}
        ref.reference(in.data(),expected.data());
        for(unsigned tile:{4U,8U,16U,512U}) for(auto layout:{Layout::Packed24,Layout::Indices32}) {
            Spin code(Configuration::T128S19,16,17,29,tile);
            code.validateSetup();code.compact(layout);
            Spin::Workspace work(code);
            code.encode(in.data(),in.size(),out.data(),out.size(),work,layout);
            if(out!=expected) throw std::runtime_error("tile changed map");
            auto alias=in;code.encodeInplace(alias.data(),alias.size(),work,layout);
            if(!std::equal(expected.begin(),expected.end(),alias.begin()) ||
               !std::equal(in.begin()+expected.size(),in.end(),alias.begin()+expected.size()))
                throw std::runtime_error("inplace mismatch");
        }
        bool rejected=false;
        try {Spin bad(Configuration::T128S19,16,1,2,2);} catch(const std::invalid_argument&) {rejected=true;}
        if(!rejected) throw std::runtime_error("four-row minimum not checked");
        std::cout<<"Four-row tiles 4/8/16/512, both schedules, compaction, inplace PASS\n";
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
