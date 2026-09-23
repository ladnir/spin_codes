#include <spin/PreparedEncoder.h>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

struct alignas(16) Word {std::uint64_t a,b;};
using Clock=std::chrono::steady_clock;
static double ms(Clock::time_point a,Clock::time_point b) {
    return std::chrono::duration<double,std::milli>(b-a).count();
}
static double median(std::vector<double> v){std::sort(v.begin(),v.end());return v[v.size()/2];}
int main(int argc,char** argv) {try {
    if(!spin::capabilities().avx2)return 77;
    if(argc!=3 || (std::string(argv[1])!="full" && std::string(argv[1])!="bank"))
        throw std::invalid_argument("usage: spin_prepared_bench full|bank K");
    const auto k=std::stoull(argv[2]);
    const bool bank=std::string(argv[1])=="bank";
    const auto parameter=k<=65536?spin::Parameters::T64S12R2:spin::Parameters::T128S19;
    const auto start=Clock::now();
    spin::PreparedEncoder c({k,parameter,17,29},
        {bank?spin::SetupMode::BankedHeuristic:spin::SetupMode::Full,913});
    const auto planned=Clock::now();auto w=c.make_workspace();const auto prepared=Clock::now();
    std::vector<Word> input(2*k);std::vector<double> refresh,encode,combined;
    refresh.reserve(101);encode.reserve(101);combined.reserve(101);
    std::uint64_t digest=0;
    for(unsigned rep=0;rep<121;++rep) {
        for(std::size_t i=0;i<input.size();++i)input[i]={i*0x9e3779b97f4a7c15ULL+rep,i^0xa531b78420ULL};
        const auto a=Clock::now();c.setCodeSeed({17+(bank?rep:0),29+(bank?rep:0)});
        const auto b=Clock::now();c.transpose_inplace<Word>(input,w);const auto d=Clock::now();
        if(rep>=20){refresh.push_back(ms(a,b));encode.push_back(ms(b,d));combined.push_back(ms(a,d));}
        digest^=input[rep%k].a;
    }
    std::cout<<"mode,k,setup_ms,workspace_ms,refresh_ms,encode_ms,combined_ms,setup_bytes,workspace_bytes,checksum\n"
        <<argv[1]<<','<<k<<','<<std::fixed<<std::setprecision(6)<<ms(start,planned)<<','<<ms(planned,prepared)
        <<','<<median(refresh)<<','<<median(encode)<<','<<median(combined)<<','<<c.setup_bytes()<<','<<w.bytes()<<','<<digest<<'\n';
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
