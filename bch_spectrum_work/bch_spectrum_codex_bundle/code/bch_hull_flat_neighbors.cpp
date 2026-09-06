// Scan affine four-flats around a small set of orbit seeds, single-threaded.
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
struct Word {
    uint64_t a,b,c,d;
    bool operator==(const Word&) const = default;
    bool operator<(const Word& x) const {
        if(d!=x.d) return d<x.d;
        if(c!=x.c) return c<x.c;
        if(b!=x.b) return b<x.b;
        return a<x.a;
    }
};
template<class T> void read_exact(std::ifstream& f,T* p,size_t n) {
    f.read(reinterpret_cast<char*>(p),sizeof(T)*n);
    if(!f) throw std::runtime_error("truncated input");
}
int main(int argc,char** argv) {
    if(argc!=3) return 2;
    static_assert(sizeof(Word)==32 && std::endian::native==std::endian::little);
    std::ifstream input(argv[1],std::ios::binary);
    uint32_t header[2]; read_exact(input,header,2);
    if(header[0]!=200787 || header[1]>20000) return 3;
    std::vector<Word> directions(header[0]),seeds(header[1]);
    read_exact(input,directions.data(),directions.size());
    read_exact(input,seeds.data(),seeds.size());
    char extra; if(input.get(extra)) return 4;
    std::array<uint8_t,256> field{},position{};
    unsigned alpha=1;
    for(unsigned i=0;i<255;++i) {
        field[i]=alpha; position[alpha]=i; alpha<<=1; if(alpha&256) alpha^=0x14d;
    }
    if(alpha!=1) return 5;
    field[255]=0; position[0]=255;
    std::vector<Word> flats; flats.reserve(directions.size()*16);
    for(const Word v:directions) {
        std::array<uint8_t,16> members{}; unsigned used=0;
        const uint64_t limbs[4]={v.a,v.b,v.c,v.d};
        for(unsigned l=0;l<4;++l) {
            uint64_t x=limbs[l];
            while(x) {
                if(used>=16) return 6;
                members[used++]=field[64*l+std::countr_zero(x)]; x&=x-1;
            }
        }
        if(used!=16) return 6;
        std::array<uint8_t,256> visited{};
        unsigned count=0;
        for(unsigned rep=0;rep<256;++rep) if(!visited[rep]) {
            uint64_t out[4]={0,0,0,0};
            for(uint8_t x:members) {
                unsigned y=x^rep,p=position[y]; visited[y]=1;
                out[p>>6]|=uint64_t(1)<<(p&63);
            }
            flats.push_back({out[0],out[1],out[2],out[3]}); ++count;
        }
        if(count!=16) return 7;
    }
    std::vector<Word> candidates; candidates.reserve(1u<<20);
    uint64_t queries=0;
    for(const Word seed:seeds) {
        for(const Word v:flats) {
            const Word x={seed.a^v.a,seed.b^v.b,seed.c^v.c,seed.d^v.d}; ++queries;
            if(std::popcount(x.a)+std::popcount(x.b)+std::popcount(x.c)+std::popcount(x.d)==20)
                candidates.push_back(x);
        }
    }
    const auto raw=candidates.size();
    std::sort(candidates.begin(),candidates.end());
    candidates.erase(std::unique(candidates.begin(),candidates.end()),candidates.end());
    std::ofstream output(argv[2],std::ios::binary|std::ios::trunc);
    output.write(reinterpret_cast<const char*>(candidates.data()),32*candidates.size());
    if(!output) return 8;
    std::cout<<"{\"affine_four_flats\":"<<flats.size()<<",\"queries\":"<<queries
        <<",\"raw_weight20_hits\":"<<raw<<",\"distinct_weight20_seeds\":"<<candidates.size()<<"}\n";
}
