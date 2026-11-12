// Single-use cost: fresh code seed, then one 128-bit transposed encode.
// "turbo"/"imt": full setup from a new seed (Code construction) into a reused
// workspace. "imtbank"/"turbobank": BankedHeuristic PreparedEncoder refresh. One process
// per invocation; run serially. Workspace creation and input fill are excluded.
#include <spin/Code.h>
#include <spin/PreparedEncoder.h>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
struct alignas(16) Word { std::uint64_t a,b; };
using Clock=std::chrono::steady_clock;
static double ms(Clock::time_point a,Clock::time_point b){return std::chrono::duration<double,std::milli>(b-a).count();}
static double median(std::vector<double> x){std::sort(x.begin(),x.end());return x[x.size()/2];}
int main(int argc,char** argv){try{
    if(argc!=4 && argc!=5) throw std::invalid_argument("usage: spin_setup_once_bench turbo|imt|imtbank|turbobank K reps [auto|avx2|avx512|portable|reference]");
    const std::string b=argc==5?argv[4]:"auto";
    const spin::ExecutionOptions exec{b=="avx2"?spin::Backend::Avx2:b=="avx512"?spin::Backend::Avx512:
        b=="portable"?spin::Backend::Portable:b=="reference"?spin::Backend::Reference:spin::Backend::Automatic};
    const std::string mode=argv[1];const auto k=std::stoull(argv[2]);const auto reps=std::stoul(argv[3]);
    std::vector<Word> input(2*k);
    std::vector<double> setup,first,warm;
    const auto fill=[&](unsigned r){for(std::size_t i=0;i<input.size();++i)input[i]={i*0x9e3779b97f4a7c15ULL+r,i^0xa531b78420ULL};};
    if(mode=="imtbank" || mode=="turbobank") {
        const auto p=mode=="turbobank"?spin::Parameters::PacketRs8T64:spin::Parameters::T128S19;
        spin::PreparedEncoder enc({k,p,17,29},{spin::SetupMode::BankedHeuristic,913},exec);
        auto w=enc.make_workspace();
        for(unsigned r=0;r<reps+3;++r){
            fill(r);const auto a=Clock::now();enc.setCodeSeed({100+r,200+r});enc.prepare_workspace(w);
            const auto b=Clock::now();enc.transpose_inplace<Word>(input,w);const auto c=Clock::now();
            enc.transpose_inplace<Word>(input,w);const auto d=Clock::now();
            if(r>=3){setup.push_back(ms(a,b));first.push_back(ms(b,c));warm.push_back(ms(c,d));}
        }
    } else {
        const auto p=mode=="turbo"?spin::Parameters::PacketRs8T64:spin::Parameters::T128S19;
        spin::Code base({k,p,1,2},exec);auto w=base.make_workspace();
        for(unsigned r=0;r<reps+3;++r){
            fill(r);const auto a=Clock::now();spin::Code c({k,p,100+r,200+r},exec);c.prepare_workspace(w);
            const auto b=Clock::now();c.transpose_inplace<Word>(input,w);const auto e=Clock::now();
            c.transpose_inplace<Word>(input,w);const auto f=Clock::now();
            if(r>=3){setup.push_back(ms(a,b));first.push_back(ms(b,e));warm.push_back(ms(e,f));}
        }
    }
    std::cout<<mode<<','<<b<<','<<k<<",setup_ms="<<std::fixed<<std::setprecision(4)<<median(setup)
        <<",first_encode_ms="<<median(first)<<",single_use_ms="<<median(setup)+median(first)
        <<",warm_encode_ms="<<median(warm)<<'\n';
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
