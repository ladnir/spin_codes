// Fixed code, fixed input, separate output: setup is outside the timed calls.
#include <spin/Code.h>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

struct alignas(16) Word {
    std::uint64_t a,b;
    bool operator==(const Word&) const = default;
};
static std::uint64_t word(std::uint64_t& seed) {
    auto x=(seed+=0x9e3779b97f4a7c15ULL);
    x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;
    x=(x^(x>>27))*0x94d049bb133111ebULL;
    return x^(x>>31);
}
int main(int argc,char** argv){try {
    if(!spin::capabilities().avx2)return 77;
    if(argc!=3 && argc!=4)throw std::invalid_argument("usage: spin_paper_bench forward|transpose exponent [r2]");
    const std::string direction=argv[1];const auto m=std::stoul(argv[2]);
    if((direction!="forward" && direction!="transpose") || (m!=16 && m!=18 && m!=20))
        throw std::invalid_argument("unsupported paper cell");
    const std::size_t k=std::size_t{1}<<m;
    const bool forward=direction=="forward";
    const bool r2=argc==4;
    if(r2 && (m!=16 || std::string(argv[3])!="r2"))throw std::invalid_argument("r2 requires exponent 16");
    const spin::CodeSpec spec{k,r2?spin::Parameters::T64S12R2:spin::Parameters::T128S19,1,2};
    // Preserve the measured transpose tile choice; forward prefers 256 rows.
    const unsigned tile=(!forward && m==20)?2048:256;
    spin::Code code(spec,{spin::Backend::Automatic,tile});auto workspace=code.make_workspace();
    std::vector<Word> input(forward?k:2*k),output(forward?2*k:k);
    std::uint64_t seed=123;
    for(auto& x:input)x={word(seed),word(seed)};
    auto run=[&]{if(forward)code.forward<Word>(input,output,workspace);
                else code.transpose<Word>(input,output,workspace);};
    {
        spin::Code reference(spec,{spin::Backend::Avx2});auto scratch=reference.make_workspace();
        std::vector<Word> expected(output.size());
        if(forward)reference.forward<Word>(input,expected,scratch);
        else reference.transpose<Word>(input,expected,scratch);
        run();if(output!=expected)throw std::runtime_error("backend mismatch");
    }
    for(unsigned i=0;i<3;++i)run();
    std::vector<double> samples(101);
    for(auto& sample:samples) {
        const auto start=std::chrono::steady_clock::now();run();
        sample=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
    }
    auto sorted=samples;std::sort(sorted.begin(),sorted.end());
    std::uint64_t checksum=0;
    for(const auto& x:output){checksum=(checksum^x.a)*0x100000001b3ULL;checksum=(checksum^x.b)*0x100000001b3ULL;}
    std::cout<<std::setprecision(12)<<"{\"direction\":\""<<direction<<"\",\"m\":"<<m
        <<",\"parameters\":\""<<(r2?"T64S12R2":"T128S19")<<"\",\"route_seed\":1,\"inner_seed\":2,\"warmups\":3,\"trials\":101"
        <<",\"backend\":"<<int(code.backend())<<",\"record_bytes\":16,\"buffer_policy\":\"fixed_input_separate_output\""
        <<",\"tile_rows\":"<<tile
        <<",\"median_ms\":"<<sorted[50]<<",\"setup_bytes\":"<<code.setup_bytes()
        <<",\"workspace_bytes\":"<<workspace.bytes()<<",\"checksum\":"<<checksum<<",\"samples_ms\":[";
    for(unsigned i=0;i<samples.size();++i)std::cout<<(i?",":"")<<samples[i];
    std::cout<<"]}\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
