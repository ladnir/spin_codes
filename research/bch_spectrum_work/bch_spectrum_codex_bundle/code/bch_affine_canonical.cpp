// Affine orbit representatives via ordered-pair normalization, single-threaded.
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
struct Canon { Word word; uint32_t stabilizer; };
struct Field {
    std::array<uint8_t,256> coord{},position{},inverse{};
    std::array<std::array<uint8_t,256>,256> mult{};
    std::array<Word,256> bits{};
    static uint8_t mul(unsigned a,unsigned b) {
        unsigned out=0;
        while(b) { if(b&1) out^=a; b>>=1; a<<=1; if(a&256) a^=0x14d; }
        return out;
    }
    Field() {
        unsigned alpha=1;
        for(unsigned i=0;i<255;++i) { coord[i]=alpha; position[alpha]=i; alpha=mul(alpha,2); }
        if(alpha!=1) throw std::runtime_error("field generator");
        coord[255]=0; position[0]=255;
        for(unsigned a=0;a<256;++a) {
            for(unsigned b=0;b<256;++b) {
                const auto c=mul(a,b); mult[a][b]=c; if(c==1) inverse[a]=b;
            }
            const unsigned p=position[a]; uint64_t limbs[4]={0,0,0,0};
            limbs[p>>6]=uint64_t(1)<<(p&63); bits[a]={limbs[0],limbs[1],limbs[2],limbs[3]};
        }
    }
};
template<unsigned W> static Canon canonical(Word word,const Field& f) {
    std::array<uint8_t,W> support{};
    const uint64_t limbs[4]={word.a,word.b,word.c,word.d}; unsigned used=0;
    for(unsigned l=0;l<4;++l) {
        uint64_t x=limbs[l];
        while(x) {
            if(used>=W) throw std::runtime_error("excess weight");
            support[used++]=f.coord[64*l+std::countr_zero(x)]; x&=x-1;
        }
    }
    if(used!=W) throw std::runtime_error("wrong weight");
    Word best={UINT64_MAX,UINT64_MAX,UINT64_MAX,UINT64_MAX}; unsigned count=0;
    for(unsigned i=0;i<W;++i) for(unsigned j=0;j<W;++j) if(i!=j) {
        const unsigned x=support[i],a=f.inverse[x^support[j]];
        const auto& row=f.mult[a]; Word image={0,0,0,0};
        for(unsigned k=0;k<W;++k) {
            const Word bit=f.bits[row[support[k]^x]];
            image.a|=bit.a; image.b|=bit.b; image.c|=bit.c; image.d|=bit.d;
        }
        if(image<best) { best=image; count=1; }
        else if(image==best) ++count;
    }
    if(!count || 65280%count) throw std::runtime_error("stabilizer");
    return {best,count};
}
int main(int argc,char** argv) {
    if(argc!=4) return 2;
    const int weight=std::stoi(argv[1]);
    if(weight!=14 && weight!=16 && weight!=18 && weight!=20 && weight!=22) return 3;
    static_assert(sizeof(Word)==32 && std::endian::native==std::endian::little);
    std::ifstream input(argv[2],std::ios::binary|std::ios::ate);
    if(!input) return 4;
    const auto size=input.tellg(); if(size<=0 || size%32) return 5;
    std::vector<Word> seeds(static_cast<size_t>(size)/32);
    input.seekg(0); input.read(reinterpret_cast<char*>(seeds.data()),size); if(!input) return 6;
    const Field field;
    std::vector<Canon> values; values.reserve(seeds.size());
    for(const Word seed:seeds) {
        switch(weight) {
            case 14: values.push_back(canonical<14>(seed,field)); break;
            case 16: values.push_back(canonical<16>(seed,field)); break;
            case 18: values.push_back(canonical<18>(seed,field)); break;
            case 20: values.push_back(canonical<20>(seed,field)); break;
            case 22: values.push_back(canonical<22>(seed,field)); break;
        }
    }
    std::sort(values.begin(),values.end(),[](const Canon& x,const Canon& y){return x.word<y.word;});
    std::ofstream output(argv[3],std::ios::binary|std::ios::trunc);
    uint64_t orbits=0,words=0; const Canon* previous=nullptr;
    for(const auto& c:values) {
        if(previous && c.word==previous->word) {
            if(c.stabilizer!=previous->stabilizer) throw std::runtime_error("inconsistent stabilizer");
            continue;
        }
        output.write(reinterpret_cast<const char*>(&c.word),32);
        output.write(reinterpret_cast<const char*>(&c.stabilizer),4);
        ++orbits; words+=65280/c.stabilizer; previous=&c;
    }
    if(!output) return 7;
    std::cout<<"{\"input_seeds\":"<<seeds.size()<<",\"affine_orbits\":"<<orbits
        <<",\"distinct_orbit_words\":"<<words<<"}\n";
}
