// Candidate search on translated 3-flats; fixed-width, single-threaded kernels.
#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
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
struct Triple { uint32_t i,j,k; };
static inline Word wxor(Word x,Word y) { return {x.a^y.a,x.b^y.b,x.c^y.c,x.d^y.d}; }
static inline Word wand(Word x,Word y) { return {x.a&y.a,x.b&y.b,x.c&y.c,x.d&y.d}; }
static inline unsigned weight(Word x) { return std::popcount(x.a)+std::popcount(x.b)+std::popcount(x.c)+std::popcount(x.d); }
static inline uint64_t hash_key(uint64_t lo,uint64_t hi) {
    uint64_t x=lo^std::rotl(hi,29);
    x^=x>>30; x*=0xbf58476d1ce4e5b9ULL;
    x^=x>>27; x*=0x94d049bb133111ebULL;
    return x^(x>>31);
}
template<class T> void read_exact(std::ifstream& f,T* p,size_t n) {
    f.read(reinterpret_cast<char*>(p),sizeof(T)*n);
    if(!f) throw std::runtime_error("truncated input");
}
static bool possible18(Word x,Word y,Word z) {
    const unsigned a=weight(wand(x,y)),b=weight(wand(x,z)),c=weight(wand(y,z));
    const unsigned t=weight(wand(wand(x,y),z));
    // A translated intersection is empty or has the direction-intersection size.
    for(unsigned mask=0;mask<8;++mask) {
        const unsigned u=(mask&1)?a:0,v=(mask&2)?b:0,w=(mask&4)?c:0;
        if(24-2*int(u+v+w)==18) return true;
        if(u>=t && v>=t && w>=t && 24-2*int(u+v+w)+4*int(t)==18) return true;
    }
    return false;
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
        const auto& x=records[i]; uint32_t p=hash_key(x.lo,x.hi)&mask;
        while(table[p].head!=UINT32_MAX && (table[p].lo!=x.lo || table[p].hi!=x.hi)) p=(p+1)&mask;
        next[i]=table[p].head; table[p]={x.lo,x.hi,i};
    }
    std::vector<Triple> triples; triples.reserve(4096);
    std::vector<uint8_t> needed(n,0);
    uint64_t queries=0,hits=0;
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
            for(uint32_t k=table[p].head;k!=UINT32_MAX;k=next[k]) {
                if(j>=k) continue;
                ++hits;
                if(possible18(x.w,y.w,records[k].w)) {
                    triples.push_back({i,j,k}); needed[j]=needed[k]=1;
                }
            }
        }
    }
    std::array<uint8_t,256> field{},position{};
    uint16_t alpha=1;
    for(unsigned i=0;i<255;++i) {
        field[i]=alpha; position[alpha]=i;
        alpha<<=1; if(alpha&256) alpha^=0x14d;
    }
    if(alpha!=1) return 6;
    field[255]=0; position[0]=255;
    std::vector<uint32_t> dense(n,UINT32_MAX);
    uint32_t count=0;
    for(uint32_t i=0;i<n;++i) if(needed[i]) dense[i]=count++;
    std::vector<std::array<Word,32>> cosets(count);
    for(uint32_t i=0;i<n;++i) if(needed[i]) {
        std::array<uint8_t,8> members{};
        const Word v=records[i].w;
        const uint64_t limbs[4]={v.a,v.b,v.c,v.d};
        unsigned used=0;
        for(unsigned limb=0;limb<4;++limb) {
            uint64_t x=limbs[limb];
            while(x) { unsigned bit=std::countr_zero(x); members[used++]=field[64*limb+bit]; x&=x-1; }
        }
        if(used!=8) return 7;
        std::array<uint8_t,256> visited{};
        unsigned c=0;
        for(unsigned rep=0;rep<256;++rep) if(!visited[rep]) {
            uint64_t out[4]={0,0,0,0};
            for(uint8_t x:members) {
                const unsigned y=x^rep,p=position[y]; visited[y]=1;
                out[p>>6]|=uint64_t(1)<<(p&63);
            }
            cosets[dense[i]][c++]={out[0],out[1],out[2],out[3]};
        }
        if(c!=32) return 8;
    }
    std::vector<Word> seeds; seeds.reserve(1u<<20);
    uint64_t affine_checks=0,raw18=0;
    for(const auto t:triples) {
        const Word x=records[t.i].w;
        const auto& ys=cosets[dense[t.j]];
        const auto& zs=cosets[dense[t.k]];
        for(const auto y:ys) {
            const Word xy=wxor(x,y);
            for(const auto z:zs) {
                const Word word=wxor(xy,z); ++affine_checks;
                if(weight(word)==18) { ++raw18; seeds.push_back(word); }
            }
        }
    }
    std::sort(seeds.begin(),seeds.end());
    seeds.erase(std::unique(seeds.begin(),seeds.end()),seeds.end());
    std::ofstream output(argv[2],std::ios::binary|std::ios::trunc);
    output.write(reinterpret_cast<const char*>(seeds.data()),seeds.size()*sizeof(Word));
    if(!output) return 9;
    std::cout<<"{\"queries\":"<<queries<<",\"matching_direction_triples\":"<<hits
        <<",\"eligible_direction_triples\":"<<triples.size()<<",\"coset_families\":"<<count
        <<",\"affine_triples_checked\":"<<affine_checks<<",\"raw_weight18_hits\":"<<raw18
        <<",\"distinct_weight18_seeds\":"<<seeds.size()<<"}\n";
}
