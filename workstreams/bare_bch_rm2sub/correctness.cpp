#include "Spin.h"
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace bare_spin;

static void require(bool value,const char* message) { if(!value) throw std::runtime_error(message); }
template<class F> static void rejects(F&& f) {
    bool rejected=false;
    try { f(); } catch(const std::invalid_argument&) { rejected=true; }
    require(rejected,"invalid argument accepted");
}
static u64 hash(const std::vector<block>& a) {
    u64 h=0xcbf29ce484222325ULL;
    for(const auto& v:a) for(auto w:v.get<u64>()) {h^=w; h*=0x100000001b3ULL;}
    return h;
}
int main() {
    try {
        for(unsigned c=0;c<4;++c) for(unsigned m:{16U,18U,20U}) {
            Spin code(static_cast<Configuration>(c),m,1,2);
            code.validateSetup(); Spin::Workspace work(code);
            std::vector<block> input(code.codeBlocks()),actual(code.messageBlocks()),expected(actual.size()),other(actual.size());
            u64 seed=123;
            for(auto& v:input) {auto lo=splitmix(seed);auto hi=splitmix(seed);v=block(hi,lo);}
            code.reference(input.data(),expected.data());
            code.encode(input.data(),input.size(),actual.data(),actual.size(),work);
            require(actual==expected,"packed implementation differs from dense oracle");
            code.encode(input.data(),input.size(),other.data(),other.size(),work,Layout::Indices32);
            require(other==expected,"32-bit route differs from dense oracle");
            const auto routeHash=code.routeHash();
            if(m==16) {
                Spin alternate(static_cast<Configuration>(c),m,1,2,512);
                Spin::Workspace small(alternate);
                alternate.encode(input.data(),input.size(),other.data(),other.size(),small);
                require(other==expected,"tile size changed map");
                rejects([&]{code.encode(input.data(),input.size()-1,actual.data(),actual.size(),work);});
                rejects([&]{code.encode(input.data(),input.size(),input.data(),actual.size(),work);});
                Spin::Workspace empty;
                rejects([&]{code.encode(input.data(),input.size(),actual.data(),actual.size(),empty);});
                rejects([&]{code.encode(input.data(),input.size(),actual.data(),actual.size(),small);});
                std::fill(input.begin(),input.end(),block(0,0));
                code.encode(input.data(),input.size(),actual.data(),actual.size(),work);
                require(std::all_of(actual.begin(),actual.end(),[](block v){return v==block(0,0);}),"zero is not mapped to zero");
                // Boundary impulses cross epoch and region transitions; use a second setup.
                Spin second(static_cast<Configuration>(c),m,17,29);
                Spin::Workspace secondWork(second);
                const auto t=second.step();
                const auto rows=second.messageBlocks()/128;
                for(auto p:std::vector<std::size_t>{0,t-1,t,rows-1,rows,input.size()-1}) input[p]=block(1,2);
                second.reference(input.data(),actual.data());
                second.encode(input.data(),input.size(),other.data(),other.size(),secondWork);
                require(actual==other,"boundary impulse oracle mismatch");
            }
            const auto bytes=code.setupBytes();
            // Restore the original random input after boundary tests.
            seed=123;
            for(auto& v:input) {auto lo=splitmix(seed);auto hi=splitmix(seed);v=block(hi,lo);}
            code.compact();
            require(code.setupBytes()<bytes,"compaction did not release storage");
            code.encode(input.data(),input.size(),actual.data(),actual.size(),work);
            require(actual==expected,"compaction changed the map");
            rejects([&]{code.encode(input.data(),input.size(),actual.data(),actual.size(),work,Layout::Indices32);});
            rejects([&]{code.compact(Layout::Indices32);});
            std::cout<<code.name()<<" m="<<m<<" hash="<<std::hex<<hash(expected)<<" route="<<routeHash<<std::dec<<" PASS\n";
        }
        rejects([]{Spin code(Configuration::T64S20,15);});
        rejects([]{Spin code(Configuration::T64S20,20,1,2,3);});
        rejects([]{Spin code(static_cast<Configuration>(99),20);});
        std::cout<<"correctness=PASS\n";
    } catch(const std::exception& e) {std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}
}
