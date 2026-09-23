#include <spin/Code.h>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <optional>
#include <string_view>
#include <vector>
struct alignas(16) Word {std::uint64_t a,b;};
using Clock=std::chrono::steady_clock;
static double elapsed(Clock::time_point a,Clock::time_point b) {
    return std::chrono::duration<double,std::milli>(b-a).count();
}
static double median(std::vector<double> x) {std::sort(x.begin(),x.end());return x[x.size()/2];}
int main(int argc,char** argv) {
    if(!spin::capabilities().avx2)return 77;
    const bool reuse=argc==2 && std::string_view(argv[1])=="--reuse";
#ifndef SPIN_BENCH_REUSE
    if(reuse){std::cerr<<"Scratch reuse unavailable in this baseline build\n";return 1;}
#endif
    std::cout<<"log_k,plan_ms,workspace_ms,first_encode_ms,setup_plus_first_ms,warm_encode_ms,setup_bytes,workspace_bytes,checksum\n";
    for(unsigned lg:{16,18,20}) {
        const auto k=std::size_t{1}<<lg;
        std::vector<Word> input(2*k);
        std::optional<spin::Workspace> workspace;
        std::vector<double> plans,workspaces,first,total,warm;
        std::uint64_t digest=0;std::size_t setupBytes=0,workspaceBytes=0;
        for(unsigned rep=0;rep<9;++rep) {
            for(std::size_t i=0;i<input.size();++i)input[i]={i*0x9e3779b97f4a7c15ULL+rep,i^0xa531b78420ULL};
            const auto a=Clock::now();
            spin::Code c({k,lg==16?spin::Parameters::T64S12R2:spin::Parameters::T128S19,17+rep,29+rep});
            const auto b=Clock::now();
#ifdef SPIN_BENCH_REUSE
            if(reuse && workspace)c.prepare_workspace(*workspace);
            else
#endif
                workspace.emplace(c.make_workspace());
            auto& w=*workspace;const auto d=Clock::now();
            c.transpose_inplace<Word>(input,w);const auto e=Clock::now();
            for(unsigned j=0;j<3;++j)c.transpose_inplace<Word>(input,w);
            const auto f=Clock::now();
            if(rep>=2) {
                plans.push_back(elapsed(a,b));workspaces.push_back(elapsed(b,d));first.push_back(elapsed(d,e));
                total.push_back(elapsed(a,e));warm.push_back(elapsed(e,f)/3);
            }
            for(const auto x:input) {digest=(digest^x.a)*0x100000001b3ULL;digest=(digest^x.b)*0x100000001b3ULL;}
            setupBytes=c.setup_bytes();workspaceBytes=w.bytes();
            if(!reuse)workspace.reset();
        }
        std::cout<<std::fixed<<std::setprecision(6)<<lg<<','<<median(plans)<<','<<median(workspaces)<<','<<median(first)
                 <<','<<median(total)<<','<<median(warm)<<','<<setupBytes<<','<<workspaceBytes<<','<<std::hex<<digest<<std::dec<<'\n';
    }
}
