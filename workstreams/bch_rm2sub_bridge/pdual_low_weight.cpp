// Bounded information-set search. Fixed four-limb XOR/popcount kernels.
// Each trial enumerates one- and two-row sums in a randomized systematic basis.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>

struct Word { uint64_t a,b,c,d; };
static inline Word add(Word x, Word y) { return {x.a^y.a,x.b^y.b,x.c^y.c,x.d^y.d}; }
static inline Word complement(Word x) { return {~x.a,~x.b,~x.c,~x.d}; }
static inline unsigned weight(Word x) {
    return __builtin_popcountll(x.a)+__builtin_popcountll(x.b)+
           __builtin_popcountll(x.c)+__builtin_popcountll(x.d);
}
static inline bool bit(Word x,unsigned col) {
    const unsigned limb=col>>6,shift=col&63;
    return ((limb==0?x.a:limb==1?x.b:limb==2?x.c:x.d)>>shift)&1;
}
static Word parse(std::string s) {
    if(s.substr(0,2)=="0x")s.erase(0,2);
    s=std::string(64-s.size(),'0')+s;
    return {std::stoull(s.substr(48,16),nullptr,16),std::stoull(s.substr(32,16),nullptr,16),
            std::stoull(s.substr(16,16),nullptr,16),std::stoull(s.substr(0,16),nullptr,16)};
}
static std::string hex(Word w) {
    std::ostringstream out;out<<"0x"<<std::hex<<std::setfill('0')
        <<std::setw(16)<<w.d<<std::setw(16)<<w.c<<std::setw(16)<<w.b<<std::setw(16)<<w.a;
    return out.str();
}
int main(int argc,char**argv) {
    if(argc!=3)return 2;
    const double seconds=std::stod(argv[1]);const uint64_t seed=std::stoull(argv[2]);
    std::array<Word,125> original,rows;std::array<unsigned,125> tags,original_tags;
    for(unsigned i=0;i<125;++i){std::string s;if(!(std::cin>>s>>original_tags[i]))return 3;original[i]=parse(s);}
    std::array<unsigned,256> columns;for(unsigned i=0;i<256;++i)columns[i]=i;
    std::mt19937_64 random(seed);
    Word best_word=original[0],best_nonzero_word=original[0];
    unsigned best=256,best_nonzero=256;uint64_t trials=0,candidates=0;
    const auto start=std::chrono::steady_clock::now();
    auto elapsed=[&](){return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();};
    auto consider=[&](Word w,unsigned tag){
        ++candidates;unsigned count=weight(w);
        if(count>128){w=complement(w);count=256-count;}
        if(count==0)return;
        if(count<best){best=count;best_word=w;}
        if(tag && count<best_nonzero){best_nonzero=count;best_nonzero_word=w;}
    };
    do {
        rows=original;tags=original_tags;std::shuffle(columns.begin(),columns.end(),random);
        unsigned rank=0;
        for(unsigned col:columns){
            unsigned j=rank;while(j<125 && !bit(rows[j],col))++j;
            if(j==125)continue;
            std::swap(rows[rank],rows[j]);std::swap(tags[rank],tags[j]);
            const Word pivot=rows[rank];const unsigned tag=tags[rank];
            for(unsigned i=0;i<125;++i)if(i!=rank && bit(rows[i],col)){
                rows[i]=add(rows[i],pivot);tags[i]^=tag;
            }
            if(++rank==125)break;
        }
        if(rank!=125)return 4;
        for(unsigned i=0;i<125;++i){
            consider(rows[i],tags[i]);
            for(unsigned j=i+1;j<125;++j)consider(add(rows[i],rows[j]),tags[i]^tags[j]);
        }
        ++trials;
    }while(best>30 && elapsed()<seconds);
    std::cout<<"{\"trials\":"<<trials<<",\"candidates\":"<<candidates
        <<",\"elapsed_seconds\":"<<elapsed()<<",\"best_weight\":"<<best
        <<",\"best_word\":\""<<hex(best_word)<<"\",\"best_nonzero_F7_weight\":"<<best_nonzero
        <<",\"best_nonzero_F7_word\":\""<<hex(best_nonzero_word)<<"\"}\n";
}
