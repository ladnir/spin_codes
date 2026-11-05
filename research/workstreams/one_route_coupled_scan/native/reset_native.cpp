// Research reset-only timing gate. No feedback, outer, route, or full encoder timing.
#define SPIN_OR_FULL_CANTOR 1
#include "AdditivePersistent32.h"
#include <algorithm>
#include <bit>
#include <chrono>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

namespace reset {
namespace m=spin::research::char2::medium32;
namespace w=spin::research::char2::wide32;
using E=m::Record128;
using Basis=spin::research::one_route::AdditivePersistent32<32>;
using Clock=std::chrono::steady_clock;
#define FI inline __attribute__((always_inline))
#define NI __attribute__((noinline))
void require(bool ok){if(!ok)throw std::runtime_error("reset algebra mismatch");}
std::uint32_t mul(std::uint32_t a,std::uint32_t b){return m::towerMultiplyScalar(a,b);}
std::uint32_t delta(std::uint32_t a){return mul(a,a)^a;}
std::uint32_t inverse(std::uint32_t a){
    std::uint32_t out=1;
    for(std::uint32_t p=0xfffffffeU;p;p>>=1){if(p&1)out=mul(out,a);a=mul(a,a);}
    return out;
}
std::uint32_t point(unsigned bits){
    std::uint32_t out=0;
    for(unsigned i=0;bits;bits>>=1,++i)if(bits&1)out^=Basis::shiftDelta(i);
    return out;
}
struct Packed {std::array<w::Coefficient,4> c;unsigned kind;};
static_assert(sizeof(Packed)==68);
Packed prepare(std::array<std::uint32_t,4> a){
    Packed p;unsigned bits=0;bool ones=true;
    for(unsigned i=0;i<4;++i){p.c[i]=w::prepareCoefficient(a[i]);bits|=a[i];ones&=a[i]==1;}
    p.kind=bits==0?0:ones?1:bits<256?8:bits<65536?16:32;
    return p;
}
Packed broadcast(std::uint32_t a){return prepare({a,a,a,a});}
FI __m512i product(__m512i x,const Packed& p){
    if(p.kind==0)return _mm512_setzero_si512();
    if(p.kind==1)return x;
    const auto c=_mm512_loadu_si512(p.c.data());
    const auto p0=_mm512_gf2p8mul_epi8(x,w::expandDiagonal<0>(c));
    if(p.kind==8)return p0;
    if(p.kind==16){
        const auto swap=_mm512_broadcast_i32x4(_mm_set_epi32(0x0e0f0c0d,0x0a0b0809,0x06070405,0x02030001));
        const auto off=_mm512_mask_blend_epi8(0xaaaaaaaaaaaaaaaaULL,
                         w::expandDiagonal<1>(c),w::expandDiagonal<3>(c));
        return _mm512_xor_si512(p0,_mm512_gf2p8mul_epi8(_mm512_shuffle_epi8(x,swap),off));
    }
    const auto p1=_mm512_gf2p8mul_epi8(_mm512_ror_epi32(x,8),w::expandDiagonal<1>(c));
    const auto p2=_mm512_gf2p8mul_epi8(_mm512_ror_epi32(x,16),w::expandDiagonal<2>(c));
    const auto p3=_mm512_gf2p8mul_epi8(_mm512_ror_epi32(x,24),w::expandDiagonal<3>(c));
    return _mm512_xor_si512(_mm512_ternarylogic_epi32(p0,p1,p2,0x96),p3);
}
template<unsigned I> FI __m512i smallTranspose(__m512i x,const Packed& p){
    constexpr int lo=I==0?0xa0:0x44,hi=I==0?0xf5:0xee;
    constexpr __mmask8 mask=I==0?0xcc:0xf0;
    const auto a=_mm512_shuffle_i32x4(x,x,lo),b=_mm512_shuffle_i32x4(x,x,hi);
    const auto s=_mm512_xor_si512(a,b);
    return _mm512_mask_blend_epi64(mask,s,_mm512_xor_si512(product(s,p),b));
}
template<unsigned H> struct Kernel {
    static constexpr unsigned B=1U<<H,Keys=3*B/4-1;
    using Buffer=std::array<E,B>;
    using Fft=std::array<Packed,Keys>;
    using Translation=std::array<Packed,H>;
    Fft fixed;
    std::array<Translation,22-H> translations;
    std::array<Packed,B/4> scaling;
    std::array<std::uint32_t,B> scale;
    unsigned normalization;
    explicit Kernel(unsigned normal=0):fixed(prepareFft(0)),normalization(normal){
        std::array<std::uint32_t,H> inv{};
        for(unsigned i=0;i<H;++i)inv[i]=(normal==1||(normal==2&&H-i<=7))?
                                          inverse(Basis::shiftDelta(H-i)):1;
        scale[0]=1;
        for(unsigned i=0;i<H;++i)for(unsigned j=0;j<(1U<<i);++j)
            scale[j+(1U<<i)]=mul(scale[j],inv[i]);
        for(unsigned j=0;j<B;j+=4)scaling[j/4]=prepare({scale[j],scale[j+1],scale[j+2],scale[j+3]});
        for(unsigned k=0;k<translations.size();++k)for(unsigned i=0;i<H;++i)
            translations[k][i]=broadcast(mul(Basis::shiftDelta(H+k-i),inv[i]));
    }
    static Fft prepareFft(std::uint32_t a){
        Fft result;std::array<std::uint32_t,H> s{};s[0]=a;
        for(unsigned i=1;i<H;++i)s[i]=delta(s[i-1]);
        for(unsigned i=1;i<H;++i)for(unsigned offset=0;offset<B;offset+=2U<<i){
            // Delta^i(point(offset)) shifts its Cantor indices down by i.
            auto c=s[i]^point(offset>>i);
            result[(B>>(i+1))-1+(offset>>(i+1))]=broadcast(c);
        }
        for(unsigned offset=0;offset<B;offset+=4){
            const auto a0=a^point(offset),a1=a^point(offset+2);
            result[B/2-1+offset/4]=prepare({a0,a0,a1,a1});
        }
        return result;
    }
    template<unsigned I> static FI void upper(E* x,const Fft& plan){
        if constexpr(I<H){
            constexpr unsigned half=1U<<I;
            for(unsigned offset=0;offset<B;offset+=2*half){
                const auto& key=plan[(B>>(I+1))-1+(offset>>(I+1))];
                for(unsigned j=0;j<half;j+=4){
                    const auto a=_mm512_loadu_si512(x+offset+j);
                    const auto b=_mm512_loadu_si512(x+offset+half+j);
                    const auto s=_mm512_xor_si512(a,b);
                    _mm512_storeu_si512(x+offset+j,s);
                    _mm512_storeu_si512(x+offset+half+j,_mm512_xor_si512(product(s,key),b));
                }
            }
            upper<I+1>(x,plan);
        }
    }
    static FI void fft(E* x,const Fft& plan){
        for(unsigned offset=0;offset<B;offset+=4){
            auto v=_mm512_loadu_si512(x+offset);
            v=smallTranspose<0>(v,plan[B/2-1+offset/4]);
            v=smallTranspose<1>(v,plan[B/4-1+offset/4]);
            _mm512_storeu_si512(x+offset,v);
        }
        upper<2>(x,plan);
    }
    template<unsigned I=0> static FI void translate(E* state,const Translation& t){
        if constexpr(I<H){
            constexpr unsigned half=1U<<I;
            if constexpr(I<2){
                for(unsigned off=0;off<B;off+=4){
                    auto v=_mm512_loadu_si512(state+off);
                    auto low=_mm512_shuffle_i32x4(v,v,I==0?0xa0:0x44);
                    v=_mm512_mask_xor_epi64(v,I==0?0xcc:0xf0,v,product(low,t[I]));
                    _mm512_storeu_si512(state+off,v);
                }
            }else{
                for(unsigned off=0;off<B;off+=2*half)for(unsigned j=0;j<half;j+=4){
                    auto a=_mm512_loadu_si512(state+off+j);
                    auto b=_mm512_loadu_si512(state+off+half+j);
                    _mm512_storeu_si512(state+off+half+j,_mm512_xor_si512(b,product(a,t[I])));
                }
            }
            translate<I+1>(state,t);
        }
    }
    static void scalarTranslate(Buffer& z,std::uint32_t a){
        for(unsigned i=0;i<H;++i){
            for(unsigned off=0;off<B;off+=2U<<i)for(unsigned j=0;j<(1U<<i);++j)
                for(unsigned lane=0;lane<4;++lane)
                    z[off+(1U<<i)+j].lane[lane]^=mul(a,z[off+j].lane[lane]);
            a=delta(a);
        }
    }
    template<unsigned Mode> NI Buffer stream(const E* input,unsigned n,const Fft* plans=nullptr)const{
        alignas(64) Buffer state{},work{};
        for(unsigned t=0;t<n;++t){
            std::memcpy(work.data(),input+t*B,sizeof(work));
            if constexpr(Mode==0)fft(work.data(),plans[t]);
            else {
                fft(work.data(),fixed);
                translate(state.data(),translations[std::countr_zero(t+1)]);
            }
            for(unsigned j=0;j<B;j+=4){
                auto v=_mm512_loadu_si512(work.data()+j);
                if constexpr(Mode==2)v=product(v,scaling[j/4]);
                _mm512_storeu_si512(state.data()+j,_mm512_xor_si512(v,_mm512_loadu_si512(state.data()+j)));
            }
        }
        return state;
    }
    template<unsigned D> static NI Buffer striped(const E* input,unsigned n,const Kernel<D>& local){
        static_assert(D<H);
        constexpr unsigned size=1U<<D;
        alignas(64) Buffer state{};
        alignas(64) typename Kernel<D>::Buffer work{};
        for(unsigned t=0;t<n;++t){
            const auto& shift=local.translations[std::countr_zero(t+1)];
            for(unsigned group=0;group<B;group+=size){
                std::memcpy(work.data(),input+t*B+group,sizeof(work));
                Kernel<D>::fft(work.data(),local.fixed);
                Kernel<D>::translate(state.data()+group,shift);
                for(unsigned j=0;j<size;j+=4)
                    _mm512_storeu_si512(state.data()+group+j,_mm512_xor_si512(
                        _mm512_loadu_si512(work.data()+j),_mm512_loadu_si512(state.data()+group+j)));
            }
        }
        return state;
    }
};
std::uint64_t digest(const E* input,unsigned n){
    std::uint64_t hash=1;
    for(unsigned j=0;j<n;++j)for(auto lane:input[j].lane)hash=(hash^lane)*1099511628211ULL;
    return hash;
}
template<unsigned H> void check(std::mt19937_64& rng){
    using K=Kernel<H>;constexpr auto B=K::B;
    K plain(0),normal(1),partial(2);
    std::vector<E> input(B*8);for(auto& e:input)for(auto& x:e.lane)x=rng();
    std::vector<typename K::Fft> plans(8);
    for(unsigned t=0;t<8;++t)plans[t]=K::prepareFft(point(((t+1)^((t+1)>>1))<<H));
    auto direct=plain.template stream<0>(input.data(),8,plans.data());
    auto gray=plain.template stream<1>(input.data(),8);
    auto norm=normal.template stream<2>(input.data(),8);
    auto part=partial.template stream<2>(input.data(),8);
    for(unsigned j=0;j<B;++j)for(unsigned lane=0;lane<4;++lane){
        norm[j].lane[lane]=mul(norm[j].lane[lane],inverse(normal.scale[j]));
        part[j].lane[lane]=mul(part[j].lane[lane],inverse(partial.scale[j]));
    }
    const auto last=point((8^(8>>1))<<H);
    K::scalarTranslate(gray,last);K::scalarTranslate(norm,last);K::scalarTranslate(part,last);
    require(!std::memcmp(direct.data(),gray.data(),sizeof(direct)));
    require(!std::memcmp(direct.data(),norm.data(),sizeof(direct)));
    require(!std::memcmp(direct.data(),part.data(),sizeof(direct)));
    // Independent polynomial evaluation of selected columns at actual dimensions.
    for(unsigned column:{0U,1U,B/2,B-1}){
        alignas(64) typename K::Buffer x{};for(auto& a:x[column].lane)a=1;
        K::fft(x.data(),plans[3]);
        std::array<std::uint32_t,B> expected{};expected[0]=1;
        auto a=point(column)^point((4^(4>>1))<<H);
        for(unsigned i=0;i<H;++i){
            for(unsigned j=0;j<(1U<<i);++j)expected[j+(1U<<i)]=mul(a,expected[j]);
            a=delta(a);
        }
        for(unsigned j=0;j<B;++j)for(auto lane:x[j].lane)require(lane==expected[j]);
    }
    std::cout<<"check B="<<B<<" streams/columns PASS\n";
    if constexpr(H>=8){
        Kernel<3> local;
        const auto striped=K::template striped<3>(input.data(),8,local);
        for(unsigned g=0;g<B;g+=8){
            std::array<E,64> gather{};
            for(unsigned t=0;t<8;++t)std::memcpy(gather.data()+8*t,input.data()+B*t+g,8*sizeof(E));
            auto expected=local.template stream<1>(gather.data(),8);
            require(!std::memcmp(expected.data(),striped.data()+g,sizeof(expected)));
        }
        std::cout<<"check B="<<B<<" striped reset stream PASS\n";
    }
}
template<unsigned H> void bench(unsigned mode,unsigned reps){
    using K=Kernel<H>;constexpr unsigned N=1U<<21,n=N/K::B;
    K kernel(mode==2?1:mode==3?2:0);
    Kernel<3> small8;Kernel<5> small32;
    std::vector<E> input(N);std::mt19937_64 rng(173);
    for(auto& e:input)for(auto& a:e.lane)a=rng();
    std::vector<typename K::Fft> plans;
    if(mode==0){
        plans.resize(n);
        for(unsigned t=0;t<n;++t)plans[t]=K::prepareFft(point(((t+1)^((t+1)>>1))<<H));
    }
    std::vector<double> samples;std::uint64_t hash=0;
    for(unsigned r=0;r<reps+3;++r){
        const auto start=Clock::now();
        const auto state=mode==4?K::template striped<3>(input.data(),n,small8):
                         mode==5?K::template striped<5>(input.data(),n,small32):
                         mode==0?kernel.template stream<0>(input.data(),n,plans.data()):
                         mode==1?kernel.template stream<1>(input.data(),n):
                                 kernel.template stream<2>(input.data(),n);
        const auto stop=Clock::now();
        hash^=digest(state.data(),K::B)+r;
        if(r>=3)samples.push_back(std::chrono::duration<double,std::milli>(stop-start).count());
    }
    std::sort(samples.begin(),samples.end());
    std::cout<<std::fixed<<std::setprecision(3)<<"RESET_ONLY B="<<K::B<<" mode="<<mode
             <<" median_ms="<<samples[reps/2]<<" min_ms="<<samples.front()<<" max_ms="<<samples.back()
             <<" checksum="<<hash<<" plan_bytes="<<plans.size()*sizeof(typename K::Fft)<<"\n";
}
}
int main(int argc,char** argv){
    try{
        if(argc==2&&std::string(argv[1])=="--check"){
            std::mt19937_64 rng(71);
            for(unsigned kind:{0U,1U,8U,16U,32U})for(unsigned r=0;r<32;++r){
                std::array<std::uint32_t,4> c{};
                for(auto& x:c)x=kind==0?0:kind==1?1:std::uint32_t(rng())&(kind==8?255:kind==16?65535:0xffffffffU);
                auto p=reset::prepare(c);std::array<reset::E,4> x{},y{};
                for(auto& e:x)for(auto& a:e.lane)a=rng();
                _mm512_storeu_si512(y.data(),reset::product(_mm512_loadu_si512(x.data()),p));
                for(unsigned j=0;j<4;++j)for(unsigned k=0;k<4;++k)
                    reset::require(y[j].lane[k]==reset::m::towerMultiplyOracle(x[j].lane[k],c[j]));
            }
            reset::check<3>(rng);reset::check<5>(rng);reset::check<8>(rng);reset::check<9>(rng);
            std::cout<<"160 SIMD coefficient fixtures PASS\n";return 0;
        }
        if(argc!=5||std::string(argv[1])!="--bench")return 2;
        unsigned b=std::stoul(argv[2]),mode=std::stoul(argv[3]),reps=std::stoul(argv[4]);
        if(mode>5||reps<3||reps>31||reps%2==0)return 2;
        if(b==256)reset::bench<8>(mode,reps);else if(b==512)reset::bench<9>(mode,reps);else return 2;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
