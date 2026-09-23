#include "Bank.h"
#include <array>
#include <iostream>
#include <iomanip>
using namespace spin::experimental::bank;
using Set=std::array<std::uint64_t,4>;
static unsigned apply(unsigned x,unsigned a,unsigned b,unsigned e){return ((a*x+b)&255)^e;}
static std::vector<Set> images(const Bank& rows,unsigned row,const std::array<unsigned,32>& target,unsigned rotation=0) {
    // (a,b,e) and (a,b+128,e XOR 128) define the same map. So do
    // (a,b,e) and (256-a,255-b,e XOR 255). Restrict a,b<128 and
    // retain all e: each representative accounts for exactly four keys.
    std::vector<Set> out;out.reserve(1U<<21);
    for(unsigned a=1;a<128;a+=2) {
        unsigned inverse=1;while(((a*inverse)&255)!=1)inverse+=2;
        for(unsigned e=0;e<256;++e) {
            std::array<unsigned,32> pre;
            for(unsigned i=0;i<32;++i)pre[i]=((rows.inverse[row*256+target[i]]^e)*inverse)&255;
            for(unsigned b=0;b<128;++b) {
                Set bits{};const unsigned shift=b*inverse;
                for(auto x:pre){const unsigned c=std::rotl(std::uint8_t(x-shift),int(rotation));const unsigned g=rows.p[row*256+c];bits[g>>6]|=1ULL<<(g&63);}
                out.push_back(bits);
            }
        }
    }
    std::sort(out.begin(),out.end());return out;
}
int main(int argc,char** argv){try {
    const unsigned seed=argc>1?unsigned(std::stoul(argv[1])):913;
    // Exhaustive identity checks justify the weighted representative space.
    for(unsigned a=1;a<256;a+=2)for(unsigned b=0;b<256;++b)for(unsigned x=0;x<256;++x) {
        if(apply(x,a,b,37)!=apply(x,a,(b+128)&255,37^128) ||
           apply(x,a,b,37)!=apply(x,256-a,255-b,37^255))throw std::runtime_error("key equivalence");
    }
    Bank rows(8,2048,std::uint64_t(seed)^0x726f7773ULL);
    std::array<unsigned,32> target;std::iota(target.begin(),target.end(),0);unsigned rightRow=1;
    const bool cycle=argc>2 && std::string(argv[2])=="cycle";
    if(cycle) {
        bool found=false;
        for(unsigned row=1;row<2048 && !found;++row) {
            std::array<bool,256> seen{};
            for(unsigned start=0;start<256 && !found;++start)if(!seen[start]) {
                std::vector<unsigned> component{start};seen[start]=true;
                for(unsigned i=0;i<component.size();++i)for(unsigned r:{0U,row}) {
                    const auto g=rows.p[r*256+((rows.inverse[r*256+component[i]]+128)&255)];
                    if(!seen[g]){seen[g]=true;component.push_back(g);}
                }
                if(component.size()==32){std::copy(component.begin(),component.end(),target.begin());rightRow=row;found=true;}
            }
        }
        if(!found)throw std::runtime_error("no 32-cycle support found");
    }
    const unsigned rotations=argc>3 && std::string(argv[3])=="rotate"?8:1;
    std::vector<std::vector<Set>> rights;
    for(unsigned r=0;r<rotations;++r) {
        rights.push_back(images(rows,rightRow,target,r));
        if(rights.back().size()!=(1U<<21))throw std::runtime_error("key count");
        std::cerr<<"right rotation "<<r<<" ready\n";
    }
    std::uint64_t equalPairs=0,matchedSets=0;
    for(unsigned r=0;r<rotations;++r) {
        const auto left=images(rows,0,target,r);
        if(left.size()!=(1U<<21))throw std::runtime_error("key count");
        for(const auto& right:rights) {
            std::size_t i=0,j=0;
            while(i<left.size() && j<right.size()) {
                if(left[i]<right[j]){++i;continue;}if(right[j]<left[i]){++j;continue;}
                const auto value=left[i];const auto a=i,b=j;
                while(i<left.size() && left[i]==value)++i;while(j<right.size() && right[j]==value)++j;
                equalPairs+=std::uint64_t(i-a)*(j-b);++matchedSets;
            }
        }
        std::cerr<<"left rotation "<<r<<" counted\n";
    }
    const std::uint64_t denominator=(1ULL<<42)*rotations*rotations;
    std::cout<<std::setprecision(15)<<"seed,right_row,cycle,rotations,matched_rotation_groups,equal_pairs,denominator,probability\n"<<seed<<','<<rightRow<<','<<cycle<<','<<rotations<<','<<matchedSets<<','<<equalPairs<<','<<denominator<<','<<double(equalPairs)/denominator<<'\n';
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
