#include "Bank.h"
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
using namespace spin::experimental::bank;
static double excess(const std::vector<std::uint64_t>& hist,std::uint64_t samples) {
    // Unbiased estimate of (N-1)*sum_y p_y^2 - 1, excluding impossible zero.
    long double pairs=0;for(unsigned i=1;i<hist.size();++i)pairs+=static_cast<long double>(hist[i])*(hist[i]-double(hist[i]!=0));
    return double((hist.size()-1)*pairs/(static_cast<long double>(samples)*(samples-1))-1);
}
template<Mode M> void check(const Bank& bank) {
    Words w(736);
    for(unsigned trial=0;trial<32;++trial) {
        const auto k=bank.sample(w);std::vector<bool> seen(bank.n);
        for(unsigned x=0;x<bank.n;++x) {
            const auto y=bank.forward<M>(x,k);
            if(y>=bank.n || seen[y] || bank.backward<M>(y,k)!=x)throw std::runtime_error("bijection/inverse");
            seen[y]=true;
            if constexpr(M==Mode::Xor || M==Mode::XorAdd) {
                const auto a=bank.in<M>(x,k),b=bank.in<M>(x^(bank.n/2),k);
                if((a^b)!=bank.n/2)throw std::runtime_error("top-pair invariant");
            }
        }
    }
}
template<Mode M,bool Uniform=false> void screen(const Bank& bank,unsigned bankSeed,unsigned samples,bool allPairs) {
    constexpr unsigned W=32;
    const unsigned n=bank.n;
    std::vector<unsigned> differences;
    for(unsigned d=1;d<n;d<<=1)differences.push_back(d);
    differences.push_back(n-1);differences.push_back((n/2)-1);
    if(allPairs){differences.resize(n-1);std::iota(differences.begin(),differences.end(),1);}
    std::vector<unsigned> queries{0};
    queries.insert(queries.end(),differences.begin(),differences.end());
    std::array<std::vector<unsigned>,3> support;
    for(unsigned i=0;i<W;++i){support[0].push_back(i);support[1].push_back(i*n/W);support[2].push_back(bank.inverse[i]);}
    for(auto& s:support)queries.insert(queries.end(),s.begin(),s.end());
    std::sort(queries.begin(),queries.end());queries.erase(std::unique(queries.begin(),queries.end()),queries.end());
    std::vector<unsigned> y(n),z(n),seen(n),member(n);
    std::vector<std::vector<std::uint64_t>> xhist(differences.size(),std::vector<std::uint64_t>(n)),ahist=xhist;
    std::array<std::uint64_t,3> peak{},peakLow{},tail{},overlap{},overlapTail{};
    std::vector<std::uint64_t> valuation(bank.bits);
    Words words(0x7374756479000000ULL+bankSeed*773+bank.count*33+unsigned(M)*19+Uniform*12345);
    for(unsigned sample=0;sample<samples;++sample) {
        if constexpr(Uniform) {
            // Uniform injective images for all queried inputs; exactly the
            // restriction of a uniform full permutation, without building it.
            for(unsigned pass=0;pass<2;++pass) {
                const unsigned stamp=2*sample+pass+1;
                for(auto x:queries){unsigned v;do{v=unsigned(words())&bank.mask;}while(seen[v]==stamp);seen[v]=stamp;(pass?z:y)[x]=v;}
            }
        }else {
            const auto a=bank.sample(words);auto b=bank.sample(words);
            b.entry=a.entry; // Stress reuse of the same table, not different tables.
            for(auto x:queries){y[x]=bank.forward<M>(x,a);z[x]=bank.forward<M>(x,b);}
        }
        for(unsigned d=0;d<differences.size();++d) {
            ++xhist[d][y[0]^y[differences[d]]];
            ++ahist[d][(y[differences[d]]-y[0])&bank.mask];
        }
        ++valuation[std::countr_zero(y[0]^y[n/2])];
        for(unsigned s=0;s<support.size();++s) {
            std::array<unsigned,4> bins{},low{};
            for(auto x:support[s]){++bins[y[x]/(n/4)];++low[y[x]&3];member[z[x]]=sample*3+s+1;}
            unsigned ov=0;for(auto x:support[s])ov+=member[y[x]]==sample*3+s+1;
            const auto pk=*std::max_element(bins.begin(),bins.end());
            peak[s]+=pk;peakLow[s]+=*std::max_element(low.begin(),low.end());tail[s]+=pk>=16;
            overlap[s]+=ov;overlapTail[s]+=ov>=(n==256?12:5);
        }
    }
    double mx=-1,ma=-1;unsigned worst=0;
    for(unsigned d=0;d<differences.size();++d) {
        const double e=excess(xhist[d],samples);if(e>mx){mx=e;worst=differences[d];}
        ma=std::max(ma,excess(ahist[d],samples));
    }
    // Exact conditional valuation profile for top-bit input pairs under XOR
    // and XOR+add. It is fixed by the bank; enumerate every top-bit pair.
    std::vector<unsigned> exact(bank.bits);
    for(unsigned j=0;j<bank.count;++j)for(unsigned x=0;x<n/2;++x)++exact[std::countr_zero(bank.p[j*n+x]^bank.p[j*n+(x^(n/2))])];
    double chi=0,exactChi=0,maxError=0;
    for(unsigned v=0;v<bank.bits;++v) {
        const double p=double(n>>(v+1))/(n-1),observed=double(valuation[v])/samples,e=double(exact[v])/(bank.count*(n/2));
        chi+=(observed-p)*(observed-p)/p;exactChi+=(e-p)*(e-p)/p;
        maxError=std::max(maxError,std::abs(observed-e));
    }
    const char* name=Uniform?"uniform":M==Mode::Xor?"xor":M==Mode::XorAdd?"xor_add":M==Mode::InputRotateAdd?"input_rotate_add":"rotate_add";
    std::cout<<bank.bits<<','<<bank.count<<','<<bankSeed<<','<<name<<','<<samples<<','<<mx<<','<<worst<<','<<ma<<','<<chi-double(bank.bits-1)/samples<<','<<exactChi<<','<<maxError;
    for(unsigned s=0;s<3;++s)std::cout<<','<<double(peak[s])/samples<<','<<double(peakLow[s])/samples<<','<<double(tail[s])/samples<<','<<double(overlap[s])/samples<<','<<double(overlapTail[s])/samples;
    std::cout<<'\n';
}
int main(int argc,char** argv){try {
    const unsigned samples=argc>1?unsigned(std::stoul(argv[1])):32768;
    const unsigned banks=argc>2?unsigned(std::stoul(argv[2])):8;
    const unsigned onlyBits=argc>3?unsigned(std::stoul(argv[3])):0;
    const unsigned onlySize=argc>4?unsigned(std::stoul(argv[4])):0;
    const bool allPairs=argc>5 && std::string(argv[5])=="all-pairs";
    const bool inputOnly=argc>6 && std::string(argv[6])=="input-only";
    if(samples<2 || !banks)throw std::invalid_argument("sample count");
    std::cout<<std::fixed<<std::setprecision(7)<<"bits,bank_size,bank_seed,family,samples,max_xor_excess,worst_dx,max_add_excess,top_valuation_excess,exact_bank_valuation_excess,valuation_prediction_error";
    for(auto name:{"interval","subspace","bank_preimage"})for(auto metric:{"peak_high","peak_low","peak_ge16","overlap","overlap_tail"})std::cout<<','<<name<<'_'<<metric;
    std::cout<<'\n';
    const std::vector<unsigned> sizes=onlySize?std::vector<unsigned>{onlySize}:std::vector<unsigned>{1,4,16,64};
    for(unsigned bits:{8U,11U})for(unsigned seed=0;seed<banks;++seed)for(unsigned size:sizes) {
        if((onlyBits && onlyBits!=bits)||(onlySize && onlySize!=size))continue;
        Bank bank(bits,size,0x62616e6b00000000ULL+seed);
        if(inputOnly){check<Mode::InputRotateAdd>(bank);screen<Mode::Xor,true>(bank,seed,samples,allPairs);screen<Mode::InputRotateAdd>(bank,seed,samples,allPairs);continue;}
        if(!seed){check<Mode::Xor>(bank);check<Mode::XorAdd>(bank);check<Mode::RotateAdd>(bank);}
        screen<Mode::Xor,true>(bank,seed,samples,allPairs);screen<Mode::Xor>(bank,seed,samples,allPairs);
        screen<Mode::XorAdd>(bank,seed,samples,allPairs);screen<Mode::RotateAdd>(bank,seed,samples,allPairs);
        std::cerr<<"done bits="<<bits<<" bank="<<size<<" seed="<<seed<<'\n';
    }
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
