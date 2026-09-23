#include "Composed.h"
#include <bitset>
#include <cmath>
#include <iomanip>
#include <iostream>
using namespace spin::experimental::bank;
static long double choose(unsigned n,unsigned k) {
    long double x=1;for(unsigned i=1;i<=k;++i)x=x*(n-k+i)/i;return x;
}
static void rowTails(unsigned seed) {
    Bank rows(8,2048,std::uint64_t(seed)^0x726f7773ULL);
    for(unsigned shape=0;shape<3;++shape) {
        std::array<std::array<std::bitset<256>,256>,2> supports;
        for(unsigned row=0;row<2;++row)for(unsigned d=0;d<256;++d)for(unsigned i=0;i<32;++i) {
            const unsigned c=shape==0?rows.inverse[row*256+i]:shape==1?i:8*i;
            supports[row][d].set(rows.p[row*256+((c-d)&255)]);
        }
        std::array<unsigned,33> histogram{};
        for(unsigned d=0;d<256;++d)for(unsigned e=0;e<256;++e)++histogram[(supports[0][d]&supports[1][e]).count()];
        double mean=0;for(unsigned k=0;k<=32;++k)mean+=k*double(histogram[k])/65536;
        for(unsigned threshold:{8U,12U,16U,24U,32U}) {
            unsigned count=0;long double uniform=0;
            for(unsigned k=threshold;k<=32;++k){count+=histogram[k];uniform+=choose(32,k)*choose(224,32-k)/choose(256,32);}
            std::cout<<"row_tail,"<<seed<<','<<shape<<','<<threshold<<','<<count<<",65536,"<<double(count)/65536<<','<<double(uniform)<<','<<mean<<'\n';
        }
    }
}
template<unsigned Count> static void regionScreen(unsigned seed) {
    ComposedBank<Count> bank(seed);constexpr unsigned n=2048;
    std::vector<unsigned> labels(bank.routes.size());
    for(unsigned i=0;i<labels.size();++i)labels[i]=unpackDestination(bank.routes[i])>>8;
    const std::uint64_t total=labels.size();
    // Exact conditional marginal for one actual region: sigma makes the stored
    // region uniform, j selects a table, and b enumerates its cyclic shifts.
    for(unsigned dx:{1U,2U,4U,8U,16U,32U,64U,128U,256U,512U,1024U,2047U}) {
        std::array<std::uint64_t,n> hist{};
        for(unsigned base=0;base<labels.size();base+=n)for(unsigned b=0;b<n;++b)++hist[labels[base+b]^labels[base+((b+dx)&(n-1))]];
        long double sq=0;for(auto count:hist)sq+=static_cast<long double>(count)*count;
        const auto excess=(n-1)*sq/(static_cast<long double>(total)*total)-1;
        std::cout<<"region_pair,"<<seed<<','<<Count<<','<<dx<<",0,"<<total<<','<<double(excess)<<",0,0\n";
    }
    const std::array<std::array<unsigned,4>,4> shapes{{{0,1,2,3},{0,1,1024,1025},{0,256,512,768},{0,17,743,1801}}};
    for(unsigned shape=0;shape<shapes.size();++shape) {
        std::uint64_t xors=0,quarters=0;
        for(unsigned base=0;base<labels.size();base+=n)for(unsigned b=0;b<n;++b) {
            std::array<unsigned,4> y;for(unsigned i=0;i<4;++i)y[i]=labels[base+((b+shapes[shape][i])&(n-1))];
            xors+=(y[0]^y[1]^y[2]^y[3])==0;
            quarters+=(y[0]/512==y[1]/512 && y[0]/512==y[2]/512 && y[0]/512==y[3]/512);
        }
        const double quarter=double(511)*510*509/(double(2047)*2046*2045);
        std::cout<<"region_four_xor,"<<seed<<','<<Count<<','<<shape<<','<<xors<<','<<total<<','<<double(xors)/total<<','<<1.0/2045<<",0\n";
        std::cout<<"region_four_quarter,"<<seed<<','<<Count<<','<<shape<<','<<quarters<<','<<total<<','<<double(quarters)/total<<','<<quarter<<",0\n";
    }
}
int main(){try {
    std::cout<<std::setprecision(12)<<"test,seed,shape_or_count,threshold_or_shape,count,denominator,value,uniform,mean\n";
    for(unsigned seed=913;seed<921;++seed){rowTails(seed);regionScreen<1>(seed);regionScreen<4>(seed);}
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
