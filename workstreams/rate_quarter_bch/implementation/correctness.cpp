#include "Spin.h"
#include "generated/QuarterCircuit.h"
#include <algorithm>
#include <iostream>
#include <stdexcept>
using namespace bare_spin;
static void require(bool ok,const char* message) {if(!ok) throw std::runtime_error(message);}
template<class F> static void rejects(F&& f) {
    bool ok=false;try {f();} catch(const std::invalid_argument&) {ok=true;}
    require(ok,"invalid argument accepted");
}
static void fill(std::vector<block>& v,u64 seed) {
    for(auto& x:v) {const auto lo=splitmix(seed),hi=splitmix(seed);x=block(hi,lo);}
}
int main() {
    try {
        // All 128 coordinate impulses, in independent lanes of the paired kernel.
        for(unsigned c=0;c<128;++c) {
            block a[128]{},b[128]{},x[32],y[32];
            for(auto& v:a) v=block(0,0);
            for(auto& v:b) v=block(0,0);
            a[c]=block(1,2); b[127-c]=block(3,4);
            quarterTranspose2(a,b,x,y);
            for(unsigned j=0;j<32;++j) {
                require(x[j]==(((QuarterRows[j][c/64]>>(c%64))&1)?block(1,2):block(0,0)),"outer impulse mismatch");
                unsigned d=127-c;
                require(y[j]==(((QuarterRows[j][d/64]>>(d%64))&1)?block(3,4):block(0,0)),"paired outer lane mismatch");
            }
        }
        for(unsigned m:{16U,18U,20U}) {
            Spin code(Configuration::T128S19,m,1,2,0,Outer::Bch128x32);
            require(code.codeBlocks()==4*code.messageBlocks(),"wrong rate");
            code.validateSetup(); Spin::Workspace work(code);
            std::vector<block> input(code.codeBlocks()),actual(code.messageBlocks()),expected(actual.size()),other(actual.size());
            fill(input,123);
            code.reference(input.data(),expected.data());
            for(auto layout:{Layout::Packed24,Layout::Indices32}) {
                code.encode(input.data(),input.size(),actual.data(),actual.size(),work,layout);
                require(actual==expected,"full dense oracle mismatch");
                auto inplace=input;
                code.encodeInplace(inplace.data(),inplace.size(),work,layout);
                require(std::equal(expected.begin(),expected.end(),inplace.begin()),"inplace output mismatch");
                require(std::equal(input.begin()+expected.size(),input.end(),inplace.begin()+expected.size()),"inplace changed suffix");
            }
            if(m==16) {
                // Different tile sizes and layouts must leave the same linear map.
                for(unsigned tile:{2U,128U,512U,1024U,2048U}) {
                    Spin alt(Configuration::T128S19,m,1,2,tile,Outer::Bch128x32);
                    Spin::Workspace w(alt);alt.validateSetup();alt.compact(Layout::Indices32);
                    alt.encode(input.data(),input.size(),other.data(),other.size(),w,Layout::Indices32);
                    require(other==expected,"tile changed map");
                    rejects([&]{alt.encode(input.data(),input.size(),other.data(),other.size(),w);});
                }
                rejects([&]{code.encode(input.data(),input.size()-1,actual.data(),actual.size(),work);});
                rejects([&]{code.encode(input.data(),input.size(),input.data(),actual.size(),work);});
                Spin::Workspace empty;
                rejects([&]{code.encode(input.data(),input.size(),actual.data(),actual.size(),empty);});
                std::fill(input.begin(),input.end(),block(0,0));
                code.encode(input.data(),input.size(),actual.data(),actual.size(),work);
                require(std::all_of(actual.begin(),actual.end(),[](block v){return v==block(0,0);}),"nonzero image of zero");
                Spin second(Configuration::T128S19,m,17,29,0,Outer::Bch128x32);
                Spin::Workspace w(second); second.validateSetup();
                const auto region=code.messageBlocks()/32;
                for(auto p:std::vector<std::size_t>{0,127,128,region-1,region,input.size()-1}) input[p]=block(1,2);
                second.reference(input.data(),actual.data());
                second.encode(input.data(),input.size(),other.data(),other.size(),w);
                require(actual==other,"boundary impulse mismatch");
                // Linearity on an independent dense input and setup.
                auto secondInput=input; fill(secondInput,456);
                second.encode(secondInput.data(),secondInput.size(),actual.data(),actual.size(),w);
                for(std::size_t i=0;i<input.size();++i) input[i]^=secondInput[i];
                std::vector<block> sum(other.size());
                second.encode(input.data(),input.size(),sum.data(),sum.size(),w);
                for(std::size_t i=0;i<sum.size();++i) require(sum[i]==(actual[i]^other[i]),"linearity mismatch");
            }
            fill(input,123); code.compact();
            code.encode(input.data(),input.size(),actual.data(),actual.size(),work);
            require(actual==expected,"compaction changed map");
            rejects([&]{code.compact(Layout::Indices32);});
            rejects([&]{code.encodeInplace(input.data(),input.size()-1,work);});
            rejects([&]{code.encodeInplace(input.data(),input.size(),work,Layout::Indices32);});
            std::cout<<"quarter t128_s19 m="<<m<<" PASS\n"<<std::flush;
        }
        rejects([]{Spin c(Configuration::T64S20,20,1,2,0,Outer::Bch128x32);});
        rejects([]{Spin c(Configuration::T128S19,20,1,2,0,static_cast<Outer>(99));});
        rejects([]{Spin c(Configuration::T128S19,20,1,2,1,Outer::Bch128x32);});
        std::cout<<"quarter correctness=PASS\n";
    } catch(const std::exception& e) {std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}
}
