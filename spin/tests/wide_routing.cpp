#include <spin/Code.h>
#include "../src/kernels/Spin.h"
#include <array>
#include <cstring>
#include <iostream>
#include <vector>

extern "C" void* spin_internal_wide_workspace256(const void*) noexcept;
extern "C" void spin_internal_wide_destroy256(void*) noexcept;
extern "C" int spin_internal_wide_encode256(const void*,void*,const void*,std::size_t,void*,std::size_t,bool) noexcept;

int main() {
    if(!spin::capabilities().avx2)return 77;
    using namespace spin::detail::kernel;
    struct alignas(32) Row {std::array<std::uint64_t,4> words;};
    unsigned cases=0;
    for(auto configuration:{Configuration::T128S19,Configuration::T64S12R2})
    for(auto layout:{Layout::Packed24,Layout::Indices32})
    for(auto backend:{BchBackend::Avx2,BchBackend::Auto}) {
        // Three outer-region units leave a partial final tile at this tile size.
        const std::size_t k=configuration==Configuration::T128S19?49152:24576;
        Spin code(configuration,MessageLength{k},17,29,1024,backend);
        code.compact(layout);
        void* work=spin_internal_wide_workspace256(&code);
        if(!work)return 1;
        std::vector<Row> input(k),output(2*k);
        std::vector<block> lane(k),expected(2*k);Spin::Workspace reference(code);
        std::uint64_t rng=1934771;
        for(unsigned kind=0;kind<3;++kind) {
            for(auto& row:input)for(auto& x:row.words) {
                rng^=rng<<13;rng^=rng>>7;rng^=rng<<17;x=kind==0?rng:0;
            }
            if(kind==2)input.back().words[3]=1;
            if(spin_internal_wide_encode256(&code,work,input.data(),2*k,output.data(),4*k,false))return 2;
            for(unsigned j=0;j<2;++j) {
                for(std::size_t i=0;i<k;++i)std::memcpy(&lane[i],input[i].words.data()+2*j,16);
                code.forward(lane.data(),k,expected.data(),2*k,reference,layout);
                for(std::size_t i=0;i<2*k;++i)
                    if(std::memcmp(&expected[i],output[i].words.data()+2*j,16))return 3;
            }
            ++cases;
        }
        spin_internal_wide_destroy256(work);
    }
    std::cout<<"PASS: "<<cases<<" wide routing cases, including 32-bit indices, partial tiles, and scratch reuse.\n";
}
