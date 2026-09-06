// Exact, single-threaded syndrome lookup on triples of 256-bit indicators.
#include <array>
#include <bit>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>
#include <algorithm>

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
struct Record { Word w; uint64_t lo,hi; };
struct Slot { uint64_t lo=0,hi=0; uint32_t head=UINT32_MAX; };
static inline uint64_t hash_key(uint64_t lo,uint64_t hi) {
    uint64_t x=lo^std::rotl(hi,29);
    x^=x>>30; x*=0xbf58476d1ce4e5b9ULL;
    x^=x>>27; x*=0x94d049bb133111ebULL;
    return x^(x>>31);
}
template<class T> void read_exact(std::ifstream& f,T* p,size_t count) {
    f.read(reinterpret_cast<char*>(p),sizeof(T)*count);
    if(!f) throw std::runtime_error("truncated input");
}
int main(int argc,char** argv) {
    if(argc!=3) return 2;
    static_assert(sizeof(Record)==48 && sizeof(Word)==32);
    static_assert(std::endian::native==std::endian::little);
    std::ifstream input(argv[1],std::ios::binary);
    uint32_t header[2]; read_exact(input,header,2);
    const uint32_t n=header[0],r=header[1];
    if(n!=97155 || r!=381) return 3;
    std::vector<Record> records(n);
    std::vector<uint32_t> reps(r),next(n,UINT32_MAX);
    read_exact(input,records.data(),n); read_exact(input,reps.data(),r);
    char extra; if(input.get(extra)) return 4;
    constexpr uint32_t capacity=1u<<18,mask=capacity-1;
    std::vector<Slot> table(capacity);
    for(uint32_t i=0;i<n;++i) {
        const auto& x=records[i];
        uint32_t p=hash_key(x.lo,x.hi)&mask;
        while(table[p].head!=UINT32_MAX && (table[p].lo!=x.lo || table[p].hi!=x.hi)) p=(p+1)&mask;
        next[i]=table[p].head;
        table[p]={x.lo,x.hi,i};
    }
    std::vector<Word> seeds;
    seeds.reserve(1u<<18);
    uint64_t queries=0,hits=0;
    std::array<uint64_t,25> histogram{};
    for(uint32_t i:reps) {
        if(i>=n) return 5;
        const auto x=records[i];
        for(uint32_t j=0;j<n;++j) {
            const auto& y=records[j];
            const uint64_t lo=x.lo^y.lo,hi=x.hi^y.hi;
            uint32_t p=hash_key(lo,hi)&mask;
            while(table[p].head!=UINT32_MAX && (table[p].lo!=lo || table[p].hi!=hi)) p=(p+1)&mask;
            ++queries;
            if(table[p].head==UINT32_MAX) continue;
            // Fixed-width XOR/popcount; no callbacks or allocation in lookup.
            const Word xy{x.w.a^y.w.a,x.w.b^y.w.b,x.w.c^y.w.c,x.w.d^y.w.d};
            for(uint32_t k=table[p].head;k!=UINT32_MAX;k=next[k]) {
                if(j>=k) continue;
                const auto& z=records[k];
                Word word{xy.a^z.w.a,xy.b^z.w.b,xy.c^z.w.c,xy.d^z.w.d};
                unsigned weight=std::popcount(word.a)+std::popcount(word.b)+std::popcount(word.c)+std::popcount(word.d);
                if(weight>24) return 6;
                ++hits; ++histogram[weight];
                if(weight==18) seeds.push_back(word);
            }
        }
    }
    std::sort(seeds.begin(),seeds.end());
    seeds.erase(std::unique(seeds.begin(),seeds.end()),seeds.end());
    std::ofstream output(argv[2],std::ios::binary|std::ios::trunc);
    output.write(reinterpret_cast<const char*>(seeds.data()),seeds.size()*sizeof(Word));
    if(!output) return 7;
    std::cout<<"{\"queries\":"<<queries<<",\"matching_triples\":"<<hits<<",\"distinct_weight18_seeds\":"<<seeds.size()<<",\"weight_histogram\":{";
    bool first=true;
    for(unsigned w=0;w<=24;++w) if(histogram[w]) {
        if(!first) std::cout<<",";
        first=false;
        std::cout<<"\""<<w<<"\":"<<histogram[w];
    }
    std::cout<<"}}\n";
}
