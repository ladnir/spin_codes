// Research inner-only gate. Not a full encoder or a sampled-code certificate.
// Accepted key pools are reused to bound untimed setup; expanded coefficient
// streams and independently shuffled partitions cover every timed update.
#include "WideDiagonal32.h"
#include <algorithm>
#include <bit>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>

namespace m=spin::research::char2::medium32;
namespace w=spin::research::char2::wide32;
using E=m::Record128;
using U=std::uint32_t;
using V=std::array<U,4>;
using P=std::array<U,2>;
constexpr U mu=0x20000000U, delta=2;
constexpr unsigned S=512,L=4096,N=S*L;
#define FI inline __attribute__((always_inline))
void check(bool b){if(!b)throw std::runtime_error("quartet native check failed");}
U mul(U a,U b){return m::towerMultiplyScalar(a,b);}
U inv(U a){U r=1;for(U p=0xfffffffeU;p;p>>=1){if(p&1)r=mul(r,a);a=mul(a,a);}return r;}
P qm(P a,P b){U p=mul(a[0],b[0]),q=mul(a[1],b[1]);return {p^mul(mu,q),mul(a[0]^a[1],b[0]^b[1])^p};}
V fm(V a,V b){
    P p=qm({a[0],a[1]},{b[0],b[1]}),q=qm({a[2],a[3]},{b[2],b[3]});
    P r=qm({a[0]^a[2],a[1]^a[3]},{b[0]^b[2],b[1]^b[3]}),d=qm({0,mu},q);
    return {p[0]^d[0],p[1]^d[1],r[0]^p[0],r[1]^p[1]};
}
std::array<U,9> prepare(V c){
    P d=qm({0,mu},{c[2],c[3]});
    auto q=[](P a){return std::array<U,3>{a[0],mul(mu,a[1]),a[0]^a[1]};};
    auto a=q({c[0],c[1]}),b=q(d),e=q({c[0]^c[2],c[1]^c[3]});
    return {a[0],a[1],a[2],b[0],b[1],b[2],e[0],e[1],e[2]};
}
bool superregular(std::array<V,4> a,unsigned g){
    for(unsigned rows=1;rows<(1U<<g);++rows)for(unsigned cols=1;cols<(1U<<g);++cols){
        unsigned n=std::popcount(rows);if(n!=unsigned(std::popcount(cols)))continue;
        U z[4][4]{};unsigned i=0;
        for(unsigned r=0;r<g;++r)if(rows>>r&1){unsigned j=0;for(unsigned c=0;c<g;++c)if(cols>>c&1)z[i][j++]=a[r][c];++i;}
        for(unsigned c=0;c<n;++c){unsigned pivot=c;while(pivot<n&&!z[pivot][c])++pivot;if(pivot==n)return false;
            for(unsigned j=c;j<n;++j)std::swap(z[pivot][j],z[c][j]);
            U reciprocal=inv(z[c][c]);
            for(unsigned r=c+1;r<n;++r){U scale=mul(z[r][c],reciprocal);for(unsigned j=c+1;j<n;++j)z[r][j]^=mul(scale,z[c][j]);}
        }
    }return true;
}
bool accepted(std::array<V,4> a,unsigned g){if(!superregular(a,g))return false;for(unsigned j=0;j<g;++j)a[j][j]^=delta;return superregular(a,g);}
struct alignas(64) Key {std::array<w::Coefficient,4> c;};
static_assert(sizeof(Key)==64);
FI __m512i prod(__m512i x,const Key& k){return w::multiply512(x,k.c.data());}
FI __m512i add(__m512i a,__m512i b){return _mm512_xor_si512(a,b);}
FI void quadratic(__m512i x,__m512i y,const Key* k,__m512i& a,__m512i& b){
    auto p=prod(x,k[0]);a=add(p,prod(y,k[1]));b=add(p,prod(add(x,y),k[2]));
}
FI void quartet(__m512i& a,__m512i& b,__m512i& c,__m512i& d,const Key* k){
    __m512i p0,p1,q0,q1,r0,r1;
    quadratic(a,b,k,p0,p1);quadratic(c,d,k+3,q0,q1);quadratic(add(a,c),add(b,d),k+6,r0,r1);
    a=add(p0,q0);b=add(p1,q1);c=add(p0,r0);d=add(p1,r1);
}
FI void pair(__m512i& a,__m512i& b,const Key* k){b=add(b,prod(a,k[0]));a=add(a,prod(b,k[1]));b=add(b,prod(a,k[2]));}
template<unsigned G,unsigned J> FI __m512i gather(const E* q,const std::uint16_t* p){
    auto x=_mm512_castsi128_si512(_mm_loadu_si128((const __m128i*)(q+p[J])));
    x=_mm512_inserti32x4(x,_mm_loadu_si128((const __m128i*)(q+p[G+J])),1);
    x=_mm512_inserti32x4(x,_mm_loadu_si128((const __m128i*)(q+p[2*G+J])),2);
    return _mm512_inserti32x4(x,_mm_loadu_si128((const __m128i*)(q+p[3*G+J])),3);
}
template<unsigned G,unsigned J> FI void scatter(E* q,const std::uint16_t* p,__m512i x){
    _mm_storeu_si128((__m128i*)(q+p[J]),_mm512_castsi512_si128(x));
    _mm_storeu_si128((__m128i*)(q+p[G+J]),_mm512_extracti32x4_epi32(x,1));
    _mm_storeu_si128((__m128i*)(q+p[2*G+J]),_mm512_extracti32x4_epi32(x,2));
    _mm_storeu_si128((__m128i*)(q+p[3*G+J]),_mm512_extracti32x4_epi32(x,3));
}
template<unsigned G,unsigned J> FI void finish(E* q,const E* x,E* y,const std::uint16_t* p,__m512i old,__m512i next){
    next=add(next,gather<G,J>(x,p));scatter<G,J>(q,p,next);
    scatter<G,J>(y,p,add(next,_mm512_gf2p8mul_epi8(old,_mm512_set1_epi8(delta))));
}
template<unsigned G,bool Recurrence> __attribute__((noinline)) void run(const Key* keys,const std::uint16_t* perm,const E* input,E* output){
    alignas(64) E state[S]{};
    constexpr unsigned Products=G==4?9:3;
    for(unsigned t=0;t<L;++t){
        for(unsigned off=0;off<S;off+=G*4,keys+=Products){
            if constexpr(Recurrence){
                auto a=gather<G,0>(state,perm+off),b=gather<G,1>(state,perm+off),oa=a,ob=b;
                if constexpr(G==4){auto c=gather<4,2>(state,perm+off),d=gather<4,3>(state,perm+off),oc=c,od=d;
                    quartet(a,b,c,d,keys);finish<4,2>(state,input,output,perm+off,oc,c);finish<4,3>(state,input,output,perm+off,od,d);
                }else pair(a,b,keys);
                finish<G,0>(state,input,output,perm+off,oa,a);finish<G,1>(state,input,output,perm+off,ob,b);
            }else{
                auto a=_mm512_loadu_si512(input+off),b=_mm512_loadu_si512(input+off+4);
                if constexpr(G==4){auto c=_mm512_loadu_si512(input+off+8),d=_mm512_loadu_si512(input+off+12);quartet(a,b,c,d,keys);
                    _mm512_storeu_si512(output+off+8,c);_mm512_storeu_si512(output+off+12,d);
                }else pair(a,b,keys);
                _mm512_storeu_si512(output+off,a);_mm512_storeu_si512(output+off+4,b);
            }
        }input+=S;output+=S;perm+=S;
    }
}
FI __m512i gather4(const E* q,const std::uint16_t* p){return gather<1,0>(q,p);}
FI void finishPacked(E* q,const E* x,E* y,__m512i old,__m512i next){
    next=add(next,_mm512_loadu_si512(x));_mm512_storeu_si512(q,next);
    _mm512_storeu_si512(y,add(next,_mm512_gf2p8mul_epi8(old,_mm512_set1_epi8(delta))));
}
// Store each step in its partition's SIMD order. Input routing must compose
// this fixed relabeling; output coordinates are likewise relabeled. Only the
// previous state is gathered. Double buffering prevents overwrite hazards.
template<unsigned G> __attribute__((noinline)) void packedRun(const Key* keys,const std::uint16_t* perm,const E* input,E* output){
    alignas(64) E states[2][S]{};E* previous=states[0];E* next=states[1];
    constexpr unsigned Products=G==4?9:3;
    for(unsigned t=0;t<L;++t){
        for(unsigned off=0;off<S;off+=G*4,keys+=Products){
            auto a=gather4(previous,perm+off),b=gather4(previous,perm+off+4),oa=a,ob=b;
            if constexpr(G==4){auto c=gather4(previous,perm+off+8),d=gather4(previous,perm+off+12),oc=c,od=d;
                quartet(a,b,c,d,keys);finishPacked(next+off+8,input+off+8,output+off+8,oc,c);finishPacked(next+off+12,input+off+12,output+off+12,od,d);
            }else pair(a,b,keys);
            finishPacked(next+off,input+off,output+off,oa,a);finishPacked(next+off+4,input+off+4,output+off+4,ob,b);
        }
        std::swap(previous,next);input+=S;output+=S;perm+=S;
    }
}
template<unsigned G> std::vector<Key> makeKeys(const std::vector<std::array<U,9>>& pool,std::mt19937& rng){
    constexpr unsigned Products=G==4?9:3;
    std::vector<Key> keys(N/G/4*Products);
    for(unsigned i=0;i<keys.size();i+=Products)for(unsigned lane=0;lane<4;++lane){
        const auto& c=pool[rng()%pool.size()];for(unsigned j=0;j<Products;++j)keys[i+j].c[lane]=w::prepareCoefficient(c[j]);
    }return keys;
}
template<unsigned G> void verify(const std::vector<Key>& keys,const std::vector<std::uint16_t>& perm,const std::vector<E>& input,std::vector<E>& output){
    run<G,true>(keys.data(),perm.data(),input.data(),output.data());
    // Independent scalar replay of the first two entire updates, all payload lanes.
    E state[S]{};constexpr unsigned Products=G==4?9:3;
    for(unsigned t=0;t<2;++t)for(unsigned off=0;off<S;off+=G*4)for(unsigned lane=0;lane<4;++lane){
        auto p=perm.data()+t*S+off+lane*G;std::array<U,9> k{};
        for(unsigned j=0;j<Products;++j)k[j]=w::recoverScalar(keys[(t*S+off)/G/4*Products+j].c[lane]);
        for(unsigned l=0;l<4;++l){V a{},b{};for(unsigned j=0;j<G;++j)a[j]=state[p[j]].lane[l];b=a;
            if constexpr(G==4){V c{k[0],k[0]^k[2],0,0};c[2]=k[6]^c[0];c[3]=(k[6]^k[8])^c[1];b=fm(c,a);}
            else{b[1]^=mul(b[0],k[0]);b[0]^=mul(b[1],k[1]);b[1]^=mul(b[0],k[2]);}
            for(unsigned j=0;j<G;++j){b[j]^=input[t*S+p[j]].lane[l];check(output[t*S+p[j]].lane[l]==(b[j]^mul(delta,a[j])));state[p[j]].lane[l]=b[j];}
        }
    }
}
template<unsigned G> void verifyPacked(const std::vector<Key>& keys,const std::vector<std::uint16_t>& perm,const std::vector<E>& input,std::vector<E>& output){
    packedRun<G>(keys.data(),perm.data(),input.data(),output.data());
    E states[2][S]{};E* previous=states[0];E* next=states[1];constexpr unsigned Products=G==4?9:3;
    for(unsigned t=0;t<3;++t){for(unsigned off=0;off<S;off+=G*4)for(unsigned lane=0;lane<4;++lane){
        std::array<U,9> k{};for(unsigned j=0;j<Products;++j)k[j]=w::recoverScalar(keys[(t*S+off)/G/4*Products+j].c[lane]);
        for(unsigned l=0;l<4;++l){V a{},b{};for(unsigned j=0;j<G;++j)a[j]=previous[perm[t*S+off+4*j+lane]].lane[l];b=a;
            if constexpr(G==4){V c{k[0],k[0]^k[2],0,0};c[2]=k[6]^c[0];c[3]=(k[6]^k[8])^c[1];b=fm(c,a);}
            else{b[1]^=mul(b[0],k[0]);b[0]^=mul(b[1],k[1]);b[1]^=mul(b[0],k[2]);}
            for(unsigned j=0;j<G;++j){unsigned index=off+4*j+lane;b[j]^=input[t*S+index].lane[l];
                check(output[t*S+index].lane[l]==(b[j]^mul(delta,a[j])));next[index].lane[l]=b[j];}
        }
    }std::swap(previous,next);}
}
int main(){
    std::mt19937 rng(0x91abc);std::vector<std::array<U,9>> pairs,quartets;
    U trace=mu,x=mu;for(unsigned i=1;i<32;++i){x=mul(x,x);trace^=x;}check(trace==1);
    while(quartets.size()<16){V c{U(rng()),U(rng()),U(rng()),U(rng())};std::array<V,4> a{};
        for(unsigned j=0;j<4;++j){V e{};e[j]=1;auto y=fm(c,e);for(unsigned i=0;i<4;++i)a[i][j]=y[i];}
        if(accepted(a,4)){
            if(quartets.empty()){
                // Extra witness for the optional stronger setup diagnostic.
                // This does not change the benchmark's accepted-pool rule.
                for(unsigned mask=0;mask<16;++mask){auto shifted=a;
                    for(unsigned j=0;j<4;++j)if(mask>>j&1)shifted[j][j]^=delta;
                    check(superregular(shifted,4));
                }
                std::cout<<"native_all_diagonal_minor_witness="<<std::hex<<c[0]<<','<<c[1]<<','<<c[2]<<','<<c[3]<<std::dec<<'\n';
            }
            quartets.push_back(prepare(c));
        }
    }
    while(pairs.size()<16){U a=rng(),b=rng(),c=rng();std::array<V,4> matrix{};
        matrix[0]={1^mul(a,b),b,0,0};matrix[1]={a^mul(c,matrix[0][0]),1^mul(b,c),0,0};
        if(accepted(matrix,2))pairs.push_back({a,b,c,0,0,0,0,0,0});
    }
    auto pk=makeKeys<2>(pairs,rng),qk=makeKeys<4>(quartets,rng);
    std::vector<std::uint16_t> perm(N);for(unsigned t=0;t<L;++t){auto first=perm.begin()+t*S;std::iota(first,first+S,0);std::shuffle(first,first+S,rng);}
    std::vector<E> input(N),output(N);for(auto& e:input)for(auto& a:e.lane)a=rng();
    verify<2>(pk,perm,input,output);verify<4>(qk,perm,input,output);
    verifyPacked<2>(pk,perm,input,output);verifyPacked<4>(qk,perm,input,output);
    std::cout<<"checks=pass accepted_keys_per_family=16 field_mu=0x20000000 delta=2 inner_only=1 pool_reused=1\n";
    using Clock=std::chrono::steady_clock;
    auto time=[&](unsigned which){auto start=Clock::now();
        switch(which){case 0:run<2,false>(pk.data(),perm.data(),input.data(),output.data());break;
        case 1:run<4,false>(qk.data(),perm.data(),input.data(),output.data());break;
        case 2:run<2,true>(pk.data(),perm.data(),input.data(),output.data());break;
        case 3:run<4,true>(qk.data(),perm.data(),input.data(),output.data());break;
        case 4:packedRun<2>(pk.data(),perm.data(),input.data(),output.data());break;
        case 5:packedRun<4>(qk.data(),perm.data(),input.data(),output.data());break;}
        auto end=Clock::now();U checksum=0;for(unsigned i=0;i<N;i+=4093)checksum^=output[i].lane[0];
        asm volatile("" : : "r"(checksum) : "memory");
        return std::chrono::duration<double,std::milli>(end-start).count();};
    for(unsigned i=0;i<6;++i)time(i);
    std::array<std::vector<double>,6> times;
    for(unsigned round=0;round<10;++round)for(unsigned j=0;j<6;++j){unsigned i=round%2?5-j:j;times[i].push_back(time(i));}
    const char* names[]={"pair_arithmetic","quartet_arithmetic","pair_recurrence","quartet_recurrence","pair_packed_recurrence","quartet_packed_recurrence"};
    for(unsigned i=0;i<6;++i){auto& t=times[i];std::sort(t.begin(),t.end());std::cout<<names[i]<<" min_ms="<<t.front()<<" median_ms="<<(t[4]+t[5])/2<<" max_ms="<<t.back()<<'\n';}
}
