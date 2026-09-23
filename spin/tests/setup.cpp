#include "../src/kernels/SetupRandom.h"
#include <algorithm>
#include <bit>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>
using namespace spin::detail::kernel::setup;
using U=std::uint64_t;
static void check(bool b){if(!b)throw std::runtime_error("setup sampler regression");}
static U oldWord(U& s) {
    auto v=(s+=0x9e3779b97f4a7c15ULL);
    v=(v^(v>>30))*0xbf58476d1ce4e5b9ULL;
    v=(v^(v>>27))*0x94d049bb133111ebULL;
    return v^(v>>31);
}
static U oldBelow(U& s,U d) {const auto threshold=(U(0)-d)%d;for(;;){auto x=oldWord(s);if(x>=threshold)return x%d;}}
int main(){try {
    for(U seed:{U(0),U(1),U(17),~U(0)}) {
        Words words(seed);auto old=seed;
        for(unsigned i=0;i<10001;++i)check(words()==oldWord(old));
    }
    U seed=912;
    for(U d=2;d<=65536;++d) {
        Divisor divisor(d);check(divisor.threshold==(U(0)-d)%d);
        for(auto x:{U(0),U(1),d-1,d,d+1,~U(0),~U(0)-d,oldWord(seed),oldWord(seed)})
            check(divisor.remainder(x,d)==x%d);
    }
    for(U d:{U(1)<<24,(U(1)<<24)+1,(U(1)<<31)-1,(U(1)<<32)-1}) {
        Divisor divisor(d);
        for(unsigned i=0;i<10000;++i){auto x=oldWord(seed);check(divisor.remainder(x,d)==x%d);}
    }
    struct Source {unsigned calls=0;U operator()(){return calls++?5:0;}} source;
    check(Divisor(3).sample(source,3)==2 && source.calls==2);
    for(unsigned n:{64,128,256,384,2048,3072,8192}) {
        Words words(17);U old=17;
        std::vector<Divisor> divisors(n+1);for(unsigned d=2;d<=n;++d)divisors[d]=Divisor(d);
        std::vector<unsigned> a(n),b(n);
        for(unsigned rep=0;rep<9;++rep) {
            std::iota(a.begin(),a.end(),0);b=a;
            for(unsigned i=n;i>1;--i)std::swap(a[i-1],a[oldBelow(old,i)]);
            for(unsigned i=n;i>1;--i)std::swap(b[i-1],b[divisors[i].sample(words,i)]);
            check(a==b);
        }
    }
    for(unsigned bits:{12,19}) {
        Words words(29);U old=29,mask=(U(1)<<bits)-1;
        for(unsigned i=0;i<100000;++i) {
            U u,v,a,b;do{u=oldWord(old)&mask;}while(!u);v=oldWord(old)&mask;
            do{a=words()&mask;}while(!a);b=words()&mask;
            if(std::popcount(u&v)&1)v^=u&-u;
            if(std::popcount(a&b)&1)b^=a&-a;
            check(a==u && b==v && !(std::popcount(a&b)&1));
        }
    }
    std::cout<<"Exact stream, modulo, rejection, shuffles, and transvection pairs PASS\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
