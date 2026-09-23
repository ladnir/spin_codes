#include "Bank.h"
#include <array>
#include <iomanip>
#include <iostream>
#include <string>
using namespace spin::experimental::bank;
using Quad=std::array<unsigned,4>;
struct Uniform {
    std::vector<unsigned> pool;
    std::vector<spin::detail::kernel::setup::Divisor> div;
    explicit Uniform(unsigned n):pool(n),div(n+1){std::iota(pool.begin(),pool.end(),0);for(unsigned i=2;i<=n;++i)div[i]=spin::detail::kernel::setup::Divisor(i);}
    void sample(Words& words,const std::vector<unsigned>& queries,std::vector<unsigned>& out) {
        // Partial Fisher--Yates gives an independent uniform injection even
        // when the initial pool is the permutation left by the preceding draw.
        for(unsigned i=0;i<queries.size();++i) {
            const unsigned left=unsigned(pool.size())-i;
            const unsigned j=i+(left>1?unsigned(div[left].sample(words,left)):0);
            std::swap(pool[i],pool[j]);out[queries[i]]=pool[i];
        }
    }
};
template<Mode M,bool U> void draw(const Bank& bank,Words& words,Uniform& uniform,const std::vector<unsigned>& q,std::vector<unsigned>& y) {
    if constexpr(U)uniform.sample(words,q,y);
    else {const auto k=bank.sample(words);for(auto x:q)y[x]=bank.forward<M>(x,k);}
}
template<Mode M,bool U> const char* name(){return U?"uniform":M==Mode::XorAdd?"xor_add":M==Mode::InputRotateAdd?"input_rotate_add":"rotate_add";}
std::vector<Quad> quads(const Bank& bank) {
    std::vector<Quad> out;
    auto add=[&](Quad q){auto s=q;std::sort(s.begin(),s.end());if(std::adjacent_find(s.begin(),s.end())==s.end())out.push_back(q);};
    for(unsigned i=0;i<bank.bits;++i)for(unsigned j=i+1;j<bank.bits;++j)add({0,1U<<i,1U<<j,(1U<<i)|(1U<<j)});
    for(unsigned a=1;a<bank.n/2;a<<=1)add({0,a,2*a,3*a});
    for(unsigned j=0;j<std::min(bank.count,16U);++j)for(unsigned stride:{1U,bank.n/4}) {
        Quad q;for(unsigned i=0;i<4;++i)q[i]=bank.inverse[j*bank.n+i*stride];add(q);
    }
    Words words(78651+bank.bits);
    while(out.size()<128) {Quad q;for(auto& x:q)x=unsigned(words())&bank.mask;add(q);}
    return out;
}
template<Mode M,bool U=false> void four(const Bank& bank,unsigned seed,unsigned samples) {
    const auto shapes=quads(bank);
    std::vector<unsigned> queries,y(bank.n);
    for(const auto& q:shapes)queries.insert(queries.end(),q.begin(),q.end());
    std::sort(queries.begin(),queries.end());queries.erase(std::unique(queries.begin(),queries.end()),queries.end());
    std::vector<std::uint64_t> xr(shapes.size()),ar(shapes.size()),qr(shapes.size());
    Words words(0x666f757200000000ULL+seed*777+bank.count*9+unsigned(M)+U*31);Uniform uniform(bank.n);
    for(unsigned trial=0;trial<samples;++trial) {
        draw<M,U>(bank,words,uniform,queries,y);
        for(unsigned i=0;i<shapes.size();++i) {
            const auto& q=shapes[i];const auto a=y[q[0]],b=y[q[1]],c=y[q[2]],d=y[q[3]];
            xr[i]+=(a^b^c^d)==0;ar[i]+=((a+d-b-c)&bank.mask)==0;
            qr[i]+=((a^b)|(a^c)|(a^d))<(bank.n/4);
        }
    }
    for(unsigned i=0;i<shapes.size();++i) {
        std::cout<<"four,"<<bank.bits<<','<<bank.count<<','<<seed<<','<<name<M,U>()<<','<<i<<','<<samples<<','<<double(xr[i])/samples<<','<<double(ar[i])/samples<<','<<double(qr[i])/samples;
        for(auto x:shapes[i])std::cout<<','<<x;
        std::cout<<'\n';
    }
}
template<Mode M,bool U=false> void adaptive(const Bank& bank,unsigned seed,unsigned samples,unsigned epoch) {
    const unsigned n=bank.n,W=32,train=4096;
    std::vector<unsigned> queries(n),y(n),affinity(n*n);
    std::iota(queries.begin(),queries.end(),0);
    Words words(0x747261696e000000ULL+seed*777+bank.count*9+unsigned(M)+U*31);Uniform uniform(n);
    for(unsigned trial=0;trial<train;++trial) {
        draw<M,U>(bank,words,uniform,queries,y);
        std::array<std::vector<unsigned>,4> bins;
        for(unsigned x=0;x<n;++x)bins[y[x]/(n/4)].push_back(x);
        for(auto& bin:bins)for(unsigned i=0;i<bin.size();++i)for(unsigned j=0;j<i;++j) {
            ++affinity[bin[i]*n+bin[j]];++affinity[bin[j]*n+bin[i]];
        }
    }
    // Maximize the number of same-quarter pairs in a 32-position support.
    // Sixteen greedy-swap restarts; select once, before drawing holdout data.
    std::vector<unsigned> best;std::uint64_t bestScore=0;
    for(unsigned restart=0;restart<16;++restart) {
        std::vector<unsigned> support,scores(n);std::vector<bool> member(n);
        for(unsigned i=0;i<W;++i) {
            unsigned x;
            if(restart<bank.count/4)x=bank.inverse[(restart%bank.count)*n+i];
            else {do{x=unsigned(words())&bank.mask;}while(member[x]);}
            support.push_back(x);member[x]=true;
        }
        for(unsigned v=0;v<n;++v)for(auto x:support)scores[v]+=affinity[v*n+x];
        for(unsigned iteration=0;iteration<256;++iteration) {
            int gain=0;unsigned remove=0,insert=0;
            for(unsigned i=0;i<W;++i)for(unsigned v=0;v<n;++v)if(!member[v]) {
                const unsigned u=support[i];const int g=int(scores[v])-int(scores[u])-int(affinity[u*n+v]);
                if(g>gain){gain=g;remove=i;insert=v;}
            }
            if(!gain)break;
            const auto old=support[remove];member[old]=false;member[insert]=true;support[remove]=insert;
            for(unsigned v=0;v<n;++v)scores[v]=scores[v]-affinity[v*n+old]+affinity[v*n+insert];
        }
        std::uint64_t score=0;for(auto x:support)score+=scores[x];
        if(score>bestScore){bestScore=score;best=support;}
    }
    if(best.size()!=W)throw std::runtime_error("adaptive support");
    Words holdout(0x686f6c646f757400ULL+seed*1777+bank.count*19+unsigned(M)+U*31+epoch*0x9e3779b97f4a7c15ULL);
    std::uint64_t pairs=0,tail=0,peak=0;
    for(unsigned trial=0;trial<samples;++trial) {
        draw<M,U>(bank,holdout,uniform,best,y);std::array<unsigned,4> bins{};
        for(auto x:best)++bins[y[x]/(n/4)];
        const auto pk=*std::max_element(bins.begin(),bins.end());tail+=pk>=16;peak+=pk;
        for(auto c:bins)pairs+=c*(c-1)/2;
    }
    std::cout<<"adaptive,"<<bank.bits<<','<<bank.count<<','<<seed<<','<<name<M,U>()<<",0,"<<samples<<','<<double(bestScore)/(train*W*(W-1))<<','<<double(pairs)/(samples*(W*(W-1)/2ULL))<<','<<double(tail)/samples<<','<<double(peak)/samples;
    for(auto x:best)std::cout<<','<<x;
    std::cout<<'\n';
}
int main(int argc,char** argv){try {
    const unsigned samples=argc>1?unsigned(std::stoul(argv[1])):262144,banks=argc>2?unsigned(std::stoul(argv[2])):8;
    const bool adaptiveOnly=argc>3 && std::string(argv[3])=="adaptive";
    const unsigned epoch=argc>4?unsigned(std::stoul(argv[4])):0;
    const bool inputOnly=argc>5 && std::string(argv[5])=="input-only";
    if(!samples||!banks)throw std::invalid_argument("samples");
    std::cout<<std::fixed<<std::setprecision(8)<<"test,bits,bank_size,seed,family,shape,samples,metric1,metric2,metric3,extra...\n";
    for(unsigned bits:{8U,11U})for(unsigned size:{16U,64U})for(unsigned seed=0;seed<banks;++seed) {
        if(adaptiveOnly && bits!=8)continue;
        Bank bank(bits,size,0x62616e6b00000000ULL+seed);
        if(inputOnly){
            if(!adaptiveOnly){four<Mode::RotateAdd,true>(bank,seed,samples);four<Mode::InputRotateAdd>(bank,seed,samples);}
            if(bits==8){adaptive<Mode::RotateAdd,true>(bank,seed,samples,epoch);adaptive<Mode::InputRotateAdd>(bank,seed,samples,epoch);}
            continue;
        }
        if(!adaptiveOnly){four<Mode::RotateAdd,true>(bank,seed,samples);four<Mode::XorAdd>(bank,seed,samples);four<Mode::RotateAdd>(bank,seed,samples);}
        if(bits==8){adaptive<Mode::RotateAdd,true>(bank,seed,samples,epoch);adaptive<Mode::XorAdd>(bank,seed,samples,epoch);adaptive<Mode::RotateAdd>(bank,seed,samples,epoch);}
        std::cerr<<"higher done "<<bits<<' '<<size<<' '<<seed<<'\n';
    }
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
