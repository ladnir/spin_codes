#include "Feistel.h"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <numeric>
using namespace spin::experimental::feistel;
template<unsigned Bits,unsigned Rounds> void screen(unsigned samples) {
    constexpr unsigned N=1U<<Bits,W=32;
    std::vector<unsigned> a(N),b(N),hist(N),counts;
    std::uint64_t totalMaxDiff=0,totalPeak=0,totalOverlap=0;
    unsigned worstDiff=0,worstPeak=0,worstOverlap=0,tail=0;
    std::vector<spin::detail::kernel::setup::Divisor> divisors(N+1);
    for(unsigned d=2;d<=N;++d)divisors[d]=spin::detail::kernel::setup::Divisor(d);
    for(unsigned seed=0;seed<samples;++seed) {
        kernelWords words(seed^0x73637265656eULL);
        if constexpr(Rounds) {
            Permutation<Bits,Rounds> p,q;p.fill(words);q.fill(words);
            for(unsigned x=0;x<N;++x){a[x]=p.forward(x);b[x]=q.forward(x);}
        }else {
            std::iota(a.begin(),a.end(),0);std::iota(b.begin(),b.end(),0);
            for(unsigned i=N;i>1;--i)std::swap(a[i-1],a[divisors[i].sample(words,i)]);
            for(unsigned i=N;i>1;--i)std::swap(b[i-1],b[divisors[i].sample(words,i)]);
        }
        unsigned maxDiff=0;
        for(unsigned dx:{1U,(1U<<(Bits/2))-1,1U<<(Bits/2),N/2,N-1}) {
            std::fill(hist.begin(),hist.end(),0);
            for(unsigned x=0;x<N;++x)++hist[a[x]^a[x^dx]];
            maxDiff=std::max(maxDiff,*std::max_element(hist.begin(),hist.end()));
        }
        totalMaxDiff+=maxDiff;worstDiff=std::max(worstDiff,maxDiff);
        // Two fixed 32-coordinate supports: contiguous and evenly spaced.
        // These are structural probes, not assertions that they are BCH words.
        for(unsigned stride:{1U,N/W}) {
            std::array<unsigned,4> bins{};
            std::fill(hist.begin(),hist.end(),0);
            for(unsigned j=0;j<W;++j){++bins[a[j*stride]/(N/4)];hist[b[j*stride]]=1;}
            unsigned overlap=0;for(unsigned j=0;j<W;++j)overlap+=hist[a[j*stride]];
            const auto peak=*std::max_element(bins.begin(),bins.end());
            totalPeak+=peak;worstPeak=std::max(worstPeak,peak);tail+=peak>=16;
            totalOverlap+=overlap;worstOverlap=std::max(worstOverlap,overlap);
        }
    }
    std::cout<<Bits<<','<<Rounds<<','<<samples<<','<<double(totalMaxDiff)/samples<<','<<worstDiff<<','
             <<double(totalPeak)/(2*samples)<<','<<worstPeak<<','<<tail<<','<<double(totalOverlap)/(2*samples)<<','<<worstOverlap<<'\n';
}
template<unsigned Bits> void all(unsigned n){screen<Bits,0>(n);screen<Bits,2>(n);screen<Bits,4>(n);screen<Bits,6>(n);screen<Bits,8>(n);}
int main(){std::cout<<std::fixed<<std::setprecision(4)
    <<"bits,rounds,samples,mean_max_xor_count,worst_xor_count,mean_peak_quarter,worst_peak_quarter,quarter_ge16,mean_overlap,worst_overlap\n";
    all<8>(1024);all<11>(256);all<13>(128);
}
