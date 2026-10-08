#pragma once
#include "Cubic8Code.h"
#include "HdBlockwiseInner.h"
#include "RateSpeed8Inner.h"
#include <vector>

namespace spin::research::char2::hdsumrank {
using Scalar=std::uint16_t;
using State=hdbw::State16;
using Mixer=hdbw::Mixer16;
using Coefficient=wide16::Coefficient;

// Native s16 research kernels. The reference has uniform GL(16,F) feedback;
// the cheaper diagonal-Cauchy variant uses an explicit point-density bound.
// The tower-basis variant preserves the proof's ranks, not identical outputs.
// All maps remain F-linear. Setup constructs the fixed maps; hot kernels
// use four-record GFNI batches and no allocation. See HD_ENCODING_OPTIMIZATION.md.
using Matrix=std::vector<std::vector<Scalar>>;
inline Scalar mul(Scalar a,Scalar b) {return medium32::quadraticMultiply(a,b);}
inline Scalar power(Scalar a,unsigned n) {
    Scalar r=1;for(;n;n>>=1,a=mul(a,a))if(n&1)r=mul(r,a);return r;
}
inline Matrix inverse(Matrix a) {
    const auto n=a.size();Matrix b(n,std::vector<Scalar>(n));
    for(unsigned i=0;i<n;++i)b[i][i]=1;
    for(unsigned j=0;j<n;++j){
        unsigned p=j;while(p<n&&!a[p][j])++p;
        if(p==n)throw std::invalid_argument("singular sumrank matrix");
        std::swap(a[j],a[p]);std::swap(b[j],b[p]);
        const auto c=power(a[j][j],65534);
        for(unsigned k=0;k<n;++k){a[j][k]=mul(a[j][k],c);b[j][k]=mul(b[j][k],c);}
        for(unsigned i=0;i<n;++i)if(i!=j){
            const auto v=a[i][j];
            for(unsigned k=0;k<n;++k){a[i][k]^=mul(v,a[j][k]);b[i][k]^=mul(v,b[j][k]);}
        }
    }return b;
}
inline Matrix product(const Matrix& a,const Matrix& b) {
    Matrix c(a.size(),std::vector<Scalar>(b[0].size()));
    for(unsigned i=0;i<a.size();++i)for(unsigned j=0;j<b[0].size();++j)
        for(unsigned k=0;k<b.size();++k)c[i][j]^=mul(a[i][k],b[k][j]);
    return c;
}
inline Matrix transpose(const Matrix& a) {
    Matrix b(a[0].size(),std::vector<Scalar>(a.size()));
    for(unsigned i=0;i<a.size();++i)for(unsigned j=0;j<a[i].size();++j)b[j][i]=a[i][j];
    return b;
}
inline Matrix blockFeedback(unsigned localRank,bool globalNodes=false,unsigned localShift=0){
    if(localRank!=4&&localRank!=8)throw std::invalid_argument("block-feedback rank");
    Matrix c(16,std::vector<Scalar>(64));
    for(unsigned i=0;i!=64;++i)for(unsigned j=0;j!=localRank;++j){
        Scalar value=1;
        for(unsigned bit=0;(1U<<bit)<localRank;++bit)if(j>>bit&1)value=mul(value,hdbw::detail::normalized(bit,(globalNodes?i:i%(4*localRank))^localShift));
        c[(i/(4*localRank))*localRank+j][i]=value;
    }
    return c;
}

struct alignas(64) FourMap {
    // Four inputs -> four outputs, in output-packed form for each input.
    std::array<std::array<Coefficient,4>,4> columns{};
    alignas(64) std::uint8_t direct[4][64]{},cross[4][64]{};
    void prepare(const Matrix& a) {
        for(unsigned i=0;i!=4;++i)for(unsigned j=0;j!=4;++j){
            columns[j][i]=wide16::prepareCoefficient(a[i][j]);
            const auto c=columns[j][i];
            for(unsigned b=0;b!=16;b+=2){
                direct[j][16*i+b]=c.a;direct[j][16*i+b+1]=c.aPlusB;
                cross[j][16*i+b]=c.b;cross[j][16*i+b+1]=c.deltaB;
            }
        }
    }
};
struct alignas(64) Maps {
    Matrix f,a,c;
    FourMap emit,feedback;
    // tile, group of four consecutive state coefficients, packed GF256
    // powers of a_i. The Moore matrix has period four in its row index.
    alignas(64) std::uint8_t weights[8][4][64]{};
    alignas(64) std::uint64_t affineWeights[8][4][8]{};
    std::uint8_t nodes[16]{};
    Scalar normalized(unsigned level,unsigned offset) const {
        Scalar a=1,b=1;
        for(unsigned j=0;j<(1U<<level);++j){a=mul(a,nodes[offset^j]);b=mul(b,nodes[(1U<<level)^j]);}
        return mul(a,power(b,65534));
    }
    explicit Maps(bool towerBasis=false,bool fftWeights=false,bool systematic=false,bool globalMds=false) {
        if(globalMds){
            f=Matrix(64,std::vector<Scalar>(64));a=Matrix(64,std::vector<Scalar>(16));
            for(unsigned i=0;i!=64;++i){
                f[i][i]=1;
                for(unsigned j=0;j!=16;++j){
                    Scalar value=1;
                    for(unsigned bit=0;bit!=4;++bit)if(j>>bit&1)value=mul(value,hdbw::detail::normalized(bit,i));
                    a[i][j]=value;
                }
            }
            c=transpose(a);return;
        }
        const auto g=power(3,17);
        for(unsigned i=0;i!=16;++i)for(unsigned j=0;j!=4;++j)
            if(i>>j&1)nodes[i]^=static_cast<std::uint8_t>(power(g,j));
        Matrix e0(8,std::vector<Scalar>(8)),e1=e0,novel=e0,novelOut=e0;
        for(unsigned i=0;i!=8;++i)for(unsigned j=0;j!=8;++j){
            e0[i][j]=power(nodes[i],j);e1[i][j]=power(nodes[i+8],j);
            novel[i][j]=novelOut[i][j]=1;
            for(unsigned b=0;b!=3;++b)if(j>>b&1){novel[i][j]=mul(novel[i][j],normalized(b,i));novelOut[i][j]=mul(novelOut[i][j],normalized(b,i+8));}
        }
        const auto ei=inverse(e0), ni=inverse(novel);
        auto localF=product(e1,ei);
        Matrix localA(8,std::vector<Scalar>(4)),localC(4,std::vector<Scalar>(8));
        Matrix monoLow(8,std::vector<Scalar>(4));
        for(unsigned i=0;i!=8;++i)for(unsigned j=0;j!=4;++j){
            localA[i][j]=towerBasis?novelOut[i][j]:e1[i][j];localC[j][i]=towerBasis?ni[j+4][i]:ei[j+4][i];monoLow[i][j]=e0[i][j];
        }
        if(systematic){
            if(!towerBasis||!fftWeights)throw std::invalid_argument("systematic maps require FFT tower basis");
            for(unsigned i=0;i!=8;++i){
                for(unsigned j=0;j!=8;++j)localF[i][j]=Scalar(i==j);
                for(unsigned j=0;j!=4;++j)localA[i][j]=localC[j][i]=novel[i][j];
            }
        }
        const auto toNovel=product(ni,monoLow),fromNovel=product(localC,novel);
        Matrix low(4,std::vector<Scalar>(4)),high=low,moore=low;
        for(unsigned i=0;i!=4;++i)for(unsigned j=0;j!=4;++j){
            low[i][j]=towerBasis?Scalar(i==j):toNovel[i][j];high[i][j]=towerBasis?Scalar(i==j):fromNovel[i][j+4];
            const Scalar basis[4]={1,79,256,20224};
            moore[i][j]=power(towerBasis?basis[j]:power(256,j),1U<<(4*i));
        }
        emit.prepare(product(low,transpose(moore)));
        feedback.prepare(product(moore,high));
        f=Matrix(64,std::vector<Scalar>(64));a=Matrix(64,std::vector<Scalar>(16));c=Matrix(16,std::vector<Scalar>(64));
        for(unsigned tile=0;tile!=8;++tile){
            const auto representative=fftWeights?Scalar(nodes[tile+8]):power(g,tile);
            Matrix r(16,std::vector<Scalar>(4));
            for(unsigned j=0;j!=16;++j){
                Scalar w=power(representative,j);
                if(fftWeights){
                    Scalar n=1;for(unsigned b=0;b!=2;++b)if((j/4)>>b&1)n=mul(n,normalized(b,tile+8));
                    w=mul(power(representative,j%4),power(n,4));
                }
                for(unsigned b=0;b!=16;++b)weights[tile][j/4][16*(j%4)+b]=static_cast<std::uint8_t>(w);
                affineWeights[tile][j/4][2*(j%4)]=affineWeights[tile][j/4][2*(j%4)+1]=rs8AffineMatrix(static_cast<std::uint8_t>(w));
                for(unsigned l=0;l!=4;++l)r[j][l]=mul(w,moore[j%4][l]);
            }
            const auto aa=product(localA,transpose(r)),cc=product(r,localC);
            for(unsigned i=0;i!=8;++i){
                for(unsigned j=0;j!=8;++j)f[8*tile+i][8*tile+j]=localF[i][j];
                for(unsigned j=0;j!=16;++j){a[8*tile+i][j]=aa[i][j];c[j][8*tile+i]=cc[j][i];}
            }
        }
    }
};

struct alignas(64) Refresh {
    std::array<std::array<Coefficient,16>,16> columns{};
    void prepare(const Matrix& a) {
        for(unsigned i=0;i!=16;++i)for(unsigned j=0;j!=16;++j)
            columns[j][i]=wide16::prepareCoefficient(a[i][j]);
    }
};
struct DiagonalRefresh {
    std::array<Coefficient,16> before{},after{};
    void prepare(const Scalar* d1,const Scalar* d2){
        wide16::prepare(d1,before.data(),16);wide16::prepare(d2,after.data(),16);
    }
};
inline Matrix cauchyRefreshMatrix(){
    Matrix a(16,std::vector<Scalar>(16)),b=a;
    for(unsigned i=0;i!=16;++i)for(unsigned j=0;j!=16;++j){a[i][j]=power(i,j);b[i][j]=power(i+16,j);}
    return product(b,inverse(a));
}

#if defined(_MSC_VER)
#define HDS_INLINE __forceinline
#else
#define HDS_INLINE inline __attribute__((always_inline))
#endif
template<unsigned I> HDS_INLINE __m512i repeat(__m512i x) {
    return _mm512_shuffle_i32x4(x,x,_MM_SHUFFLE(I,I,I,I));
}
HDS_INLINE __m512i map4(__m512i x,const FourMap& m) {
    // Fixed maps reuse expanded coefficients. Sum direct and cross terms
    // separately so byte swapping occurs once after the linear reduction.
    const auto a=repeat<0>(x),b=repeat<1>(x),c=repeat<2>(x),d=repeat<3>(x);
    const auto direct=_mm512_xor_si512(
        _mm512_xor_si512(_mm512_gf2p8mul_epi8(a,_mm512_load_si512(m.direct[0])),_mm512_gf2p8mul_epi8(b,_mm512_load_si512(m.direct[1]))),
        _mm512_xor_si512(_mm512_gf2p8mul_epi8(c,_mm512_load_si512(m.direct[2])),_mm512_gf2p8mul_epi8(d,_mm512_load_si512(m.direct[3]))));
    const auto cross=_mm512_xor_si512(
        _mm512_xor_si512(_mm512_gf2p8mul_epi8(a,_mm512_load_si512(m.cross[0])),_mm512_gf2p8mul_epi8(b,_mm512_load_si512(m.cross[1]))),
        _mm512_xor_si512(_mm512_gf2p8mul_epi8(c,_mm512_load_si512(m.cross[2])),_mm512_gf2p8mul_epi8(d,_mm512_load_si512(m.cross[3]))));
    return _mm512_xor_si512(direct,medium32::selectBytes<1,0,3,2>(cross));
}
constexpr std::uint8_t node(unsigned i){
    constexpr auto g=medium32::bytePower(3,17);
    std::uint8_t x=0;for(unsigned j=0;j!=4;++j)if(i>>j&1)x^=medium32::bytePower(g,j);
    return x;
}
constexpr std::uint8_t normalized(unsigned level,unsigned offset){
    std::uint8_t a=1,b=1;
    for(unsigned j=0;j<(1U<<level);++j){a=sfMultiply(a,node(offset^j));b=sfMultiply(b,node((1U<<level)^j));}
    return sfMultiply(a,hdbw::gf8Inverse(b));
}
template<bool Affine,unsigned A,unsigned B=A,unsigned C=A,unsigned D=A>
HDS_INLINE __m512i byteProduct4(__m512i x){
    if constexpr(Affine){
        const auto m=_mm512_setr_epi64(rs8AffineMatrix(A),rs8AffineMatrix(A),rs8AffineMatrix(B),rs8AffineMatrix(B),
                                      rs8AffineMatrix(C),rs8AffineMatrix(C),rs8AffineMatrix(D),rs8AffineMatrix(D));
        return _mm512_gf2p8affine_epi64_epi8(x,m,0);
    }else return _mm512_gf2p8mul_epi8(x,c256_inner_fft_detail::coefficients<A,B,C,D>());
}
template<unsigned Half,unsigned Offset,bool Inverse,bool Affine=false>
HDS_INLINE __m512i fixedStage(__m512i x){
    constexpr auto c0=normalized(Half==1?0:1,Offset),c1=normalized(Half==1?0:1,Offset^(Half==1?2:0));
    constexpr auto hi=Half==1?0xcc:0xf0,lo=Half==1?0x33:0x0f;
    const auto a=_mm512_shuffle_i32x4(x,x,Half==1?0xa0:0x44);
    if constexpr(c0==0&&c1==0)return _mm512_mask_xor_epi64(x,hi,x,a);
    else{
        const auto b=_mm512_shuffle_i32x4(x,x,Half==1?0xf5:0xee);
        if constexpr(Inverse){const auto sum=_mm512_xor_si512(a,b);
            const auto product=[&](){if constexpr(Half==1)return byteProduct4<Affine,c0,c0,c1,c1>(sum);else return byteProduct4<Affine,c0>(sum);}();
            return _mm512_mask_xor_epi64(sum,lo,a,product);
        }
        else{
            const auto product=[&](){
                if constexpr(Half==1)return byteProduct4<Affine,c0,c0^1,c1,c1^1>(b);
                else return byteProduct4<Affine,c0,c0,c0^1,c0^1>(b);
            }();
            return _mm512_xor_si512(a,product);
        }
    }
}
template<unsigned Offset,bool Inverse,bool Affine=false>
HDS_INLINE void fixedTransform(__m512i& a,__m512i& b){
    constexpr auto c=normalized(2,Offset);
    if constexpr(Inverse){
        a=fixedStage<1,Offset,true,Affine>(a);b=fixedStage<1,Offset^4,true,Affine>(b);
        a=fixedStage<2,Offset,true,Affine>(a);b=fixedStage<2,Offset^4,true,Affine>(b);
        b=_mm512_xor_si512(a,b);
        if constexpr(c==1)a=_mm512_xor_si512(a,b);
        else if constexpr(c!=0)a=_mm512_xor_si512(a,byteProduct4<Affine,c>(b));
    }else{
        if constexpr(c==1)a=_mm512_xor_si512(a,b);
        else if constexpr(c!=0)a=_mm512_xor_si512(a,byteProduct4<Affine,c>(b));
        b=_mm512_xor_si512(a,b);
        a=fixedStage<2,Offset,false,Affine>(a);b=fixedStage<2,Offset^4,false,Affine>(b);
        a=fixedStage<1,Offset,false,Affine>(a);b=fixedStage<1,Offset^4,false,Affine>(b);
    }
}
struct Fixed {
    Maps maps;
    explicit Fixed(bool towerBasis=false,bool fftWeights=false,bool systematic=false,bool globalMds=false):maps(towerBasis,fftWeights,systematic,globalMds){}
};

// The same natural-basis FFT as hdbw::transform16, with fixed products
// expressed as GFNI affine maps. Keep each stage/unrolled width explicit.
template<unsigned Coefficient,hdbw::detail::Kind K>
HDS_INLINE void affinePair(__m512i& a,__m512i& b){
    using Kind=hdbw::detail::Kind;
    const auto product=[](__m512i x){return byteProduct4<true,Coefficient>(x);};
    if constexpr(K==Kind::Forward){
        if constexpr(Coefficient==1)a=_mm512_xor_si512(a,b);
        else if constexpr(Coefficient!=0)a=_mm512_xor_si512(a,product(b));
        b=_mm512_xor_si512(a,b);
    }else if constexpr(K==Kind::Inverse){
        b=_mm512_xor_si512(a,b);
        if constexpr(Coefficient==1)a=_mm512_xor_si512(a,b);
        else if constexpr(Coefficient!=0)a=_mm512_xor_si512(a,product(b));
    }else if constexpr(K==Kind::ForwardTranspose){
        a=_mm512_xor_si512(a,b);
        if constexpr(Coefficient==1)b=_mm512_xor_si512(b,a);
        else if constexpr(Coefficient!=0)b=_mm512_xor_si512(b,product(a));
    }else{
        if constexpr(Coefficient==1)b=_mm512_xor_si512(b,a);
        else if constexpr(Coefficient!=0)b=_mm512_xor_si512(b,product(a));
        a=_mm512_xor_si512(a,b);
    }
}
template<unsigned Half,unsigned Offset,hdbw::detail::Kind K>
HDS_INLINE __m512i affineSmall(__m512i x){
    using Kind=hdbw::detail::Kind;
    constexpr auto c0=hdbw::detail::normalized(Half==1?0:1,Offset);
    constexpr auto c1=hdbw::detail::normalized(Half==1?0:1,Offset^(Half==1?2:0));
    constexpr int lowShuffle=Half==1?0xa0:0x44,highShuffle=Half==1?0xf5:0xee;
    constexpr __mmask8 lowMask=Half==1?0x33:0x0f,highMask=static_cast<__mmask8>(~lowMask);
    if constexpr(c0==0&&c1==0){
        if constexpr(K==Kind::Forward||K==Kind::Inverse)return _mm512_mask_xor_epi64(x,highMask,x,_mm512_shuffle_i32x4(x,x,lowShuffle));
        else return _mm512_mask_xor_epi64(x,lowMask,x,_mm512_shuffle_i32x4(x,x,highShuffle));
    }else{
        const auto a=_mm512_shuffle_i32x4(x,x,lowShuffle),b=_mm512_shuffle_i32x4(x,x,highShuffle);
        if constexpr(K==Kind::Forward||K==Kind::InverseTranspose){
            constexpr unsigned flip=K==Kind::Forward?0:1;
            const auto product=[](__m512i value){
                if constexpr(Half==1)return byteProduct4<true,c0^flip,c0^flip^1,c1^flip,c1^flip^1>(value);
                else return byteProduct4<true,c0^flip,c0^flip,c0^flip^1,c0^flip^1>(value);
            };
            if constexpr(K==Kind::Forward)return _mm512_xor_si512(a,product(b));
            else return _mm512_xor_si512(b,product(a));
        }else{
            const auto sum=_mm512_xor_si512(a,b);
            const auto product=[&](){
                if constexpr(Half==1)return byteProduct4<true,c0,c0,c1,c1>(sum);
                else return byteProduct4<true,c0>(sum);
            }();
            if constexpr(K==Kind::Inverse)return _mm512_mask_xor_epi64(sum,lowMask,a,product);
            else return _mm512_mask_xor_epi64(sum,highMask,b,product);
        }
    }
}
template<unsigned Offset,hdbw::detail::Kind K,bool Affine>
HDS_INLINE void naturalFft16(__m512i& a,__m512i& b,__m512i& c,__m512i& d){
    using Kind=hdbw::detail::Kind;
    if constexpr(!Affine)hdbw::detail::transform16<Offset,K>(a,b,c,d);
    else if constexpr(K==Kind::Forward||K==Kind::InverseTranspose){
        affinePair<hdbw::detail::normalized(3,Offset),K>(a,c);affinePair<hdbw::detail::normalized(3,Offset),K>(b,d);
        affinePair<hdbw::detail::normalized(2,Offset),K>(a,b);affinePair<hdbw::detail::normalized(2,Offset^8),K>(c,d);
        a=affineSmall<2,Offset,K>(a);b=affineSmall<2,Offset^4,K>(b);c=affineSmall<2,Offset^8,K>(c);d=affineSmall<2,Offset^12,K>(d);
        a=affineSmall<1,Offset,K>(a);b=affineSmall<1,Offset^4,K>(b);c=affineSmall<1,Offset^8,K>(c);d=affineSmall<1,Offset^12,K>(d);
    }else{
        a=affineSmall<1,Offset,K>(a);b=affineSmall<1,Offset^4,K>(b);c=affineSmall<1,Offset^8,K>(c);d=affineSmall<1,Offset^12,K>(d);
        a=affineSmall<2,Offset,K>(a);b=affineSmall<2,Offset^4,K>(b);c=affineSmall<2,Offset^8,K>(c);d=affineSmall<2,Offset^12,K>(d);
        affinePair<hdbw::detail::normalized(2,Offset),K>(a,b);affinePair<hdbw::detail::normalized(2,Offset^8),K>(c,d);
        affinePair<hdbw::detail::normalized(3,Offset),K>(a,c);affinePair<hdbw::detail::normalized(3,Offset),K>(b,d);
    }
}

HDS_INLINE __m512i towerScale(__m512i x){
    // v=256, u=79, v^16=v+u, u^16=u+1; rows use v,v+u,v+1,v+u+1.
    static constexpr Coefficient c[4]={wide16::prepareCoefficient(256),wide16::prepareCoefficient(335),
                                      wide16::prepareCoefficient(257),wide16::prepareCoefficient(334)};
    return wide16::multiply512(x,c);
}
template<bool Transpose,bool Affine=false> HDS_INLINE __m512i towerMoore(__m512i x){
    if constexpr(!Transpose){
        const auto low=_mm512_shuffle_i32x4(x,x,0xa0),high=_mm512_shuffle_i32x4(x,x,0xf5);
        auto z=_mm512_xor_si512(low,byteProduct4<Affine,79>(high));
        z=_mm512_mask_xor_epi64(z,0xcc,z,high);
        return _mm512_xor_si512(_mm512_shuffle_i32x4(z,z,0x44),towerScale(_mm512_shuffle_i32x4(z,z,0xee)));
    }else{
        const auto low=_mm512_shuffle_i32x4(x,x,0x44),high=_mm512_shuffle_i32x4(x,x,0xee);
        const auto sum=_mm512_xor_si512(low,high);
        // towerScale's high lanes contain v+1 and v+u+1. XOR low lanes,
        // rather than high lanes, to realize v*sum+old_high.
        auto z=_mm512_mask_xor_epi64(sum,0xf0,towerScale(sum),low);
        const auto a=_mm512_shuffle_i32x4(z,z,0xa0),b=_mm512_shuffle_i32x4(z,z,0xf5);
        const auto ab=_mm512_xor_si512(a,b);
        return _mm512_mask_xor_epi64(ab,0xcc,byteProduct4<Affine,79>(ab),b);
    }
}

template<unsigned J> HDS_INLINE void refreshColumn(const State& x,State& y,const Refresh& t) {
    const auto v=repeat<J%4>(x.values[J/4]);
    y.values[0]=_mm512_xor_si512(y.values[0],wide16::multiply512(v,t.columns[J].data()));
    y.values[1]=_mm512_xor_si512(y.values[1],wide16::multiply512(v,t.columns[J].data()+4));
    y.values[2]=_mm512_xor_si512(y.values[2],wide16::multiply512(v,t.columns[J].data()+8));
    y.values[3]=_mm512_xor_si512(y.values[3],wide16::multiply512(v,t.columns[J].data()+12));
}
HDS_INLINE State refresh(const State& x,const Refresh& t) {
    State y{};
    c8Unroll(std::make_index_sequence<16>{},[&](auto j){refreshColumn<j>(x,y,t);});
    return y;
}
template<bool Affine=false> HDS_INLINE State refresh(const State& x,const DiagonalRefresh& t) {
    State y;
    c8Unroll(std::make_index_sequence<4>{},[&](auto i){y.values[i]=wide16::multiply512(x.values[i],t.before.data()+4*i);});
    using Kind=hdbw::detail::Kind;
    naturalFft16<0,Kind::Inverse,Affine>(y.values[0],y.values[1],y.values[2],y.values[3]);
    naturalFft16<16,Kind::Forward,Affine>(y.values[0],y.values[1],y.values[2],y.values[3]);
    c8Unroll(std::make_index_sequence<4>{},[&](auto i){y.values[i]=wide16::multiply512(y.values[i],t.after.data()+4*i);});
    return y;
}
template<unsigned I,bool Affine=false> HDS_INLINE __m512i weight(__m512i x,const Maps& m,unsigned tile) {
    if constexpr(Affine)return _mm512_gf2p8affine_epi64_epi8(x,_mm512_load_si512(m.affineWeights[tile][I]),0);
    else return _mm512_gf2p8mul_epi8(x,_mm512_load_si512(m.weights[tile][I]));
}
template<unsigned Tile,bool Routed,class Index,bool TowerBasis=false,bool Scale=false>
HDS_INLINE void tile(const Block* input,const Index* route,Block* output,
                    const State& q,State& feedback,const Fixed& f,const Coefficient* diagonal) {
    auto a=hdbw::detail::input4<Routed>(input,route,8*Tile);
    auto b=hdbw::detail::input4<Routed>(input,route,8*Tile+4);
    if constexpr(Scale){a=wide16::multiply512(a,diagonal+8*Tile);b=wide16::multiply512(b,diagonal+8*Tile+4);}
    fixedTransform<0,true>(a,b);
    const auto h=[&](){if constexpr(TowerBasis)return towerMoore<false>(b);else return map4(b,f.maps.feedback);}();
    feedback.values[0]=_mm512_xor_si512(feedback.values[0],weight<0>(h,f.maps,Tile));
    feedback.values[1]=_mm512_xor_si512(feedback.values[1],weight<1>(h,f.maps,Tile));
    feedback.values[2]=_mm512_xor_si512(feedback.values[2],weight<2>(h,f.maps,Tile));
    feedback.values[3]=_mm512_xor_si512(feedback.values[3],weight<3>(h,f.maps,Tile));
    const auto z=_mm512_xor_si512(_mm512_xor_si512(weight<0>(q.values[0],f.maps,Tile),weight<1>(q.values[1],f.maps,Tile)),
                                _mm512_xor_si512(weight<2>(q.values[2],f.maps,Tile),weight<3>(q.values[3],f.maps,Tile)));
    a=_mm512_xor_si512(a,[&](){if constexpr(TowerBasis)return towerMoore<true>(z);else return map4(z,f.maps.emit);}());
    fixedTransform<8,false>(a,b);
    _mm512_storeu_si512(output+8*Tile,a);_mm512_storeu_si512(output+8*Tile+4,b);
}
template<bool Routed,class Index,bool RefreshEnabled=true,bool TowerBasis=false,bool Scale=false,class Feedback>
HDS_INLINE void packet(const Block* input,const Index* route,Block* output,State& state,
                       const Mixer& mixer,const Feedback& t,const Fixed& f,const Coefficient* diagonal=nullptr) {
    const auto q=state;State feedback{};
    c8Unroll(std::make_index_sequence<8>{},[&](auto i){tile<i,Routed,Index,TowerBasis,Scale>(input,route,output,q,feedback,f,diagonal);});
    hdbw::detail::mix(state,mixer);
    if constexpr(RefreshEnabled)feedback=refresh(feedback,t);
    state.values[0]=_mm512_xor_si512(state.values[0],feedback.values[0]);
    state.values[1]=_mm512_xor_si512(state.values[1],feedback.values[1]);
    state.values[2]=_mm512_xor_si512(state.values[2],feedback.values[2]);
    state.values[3]=_mm512_xor_si512(state.values[3],feedback.values[3]);
}

template<unsigned Level,unsigned Offset,bool Transpose,bool Affine=false>
HDS_INLINE void weightPair(__m512i& a,__m512i& b){
    constexpr auto n=normalized(Level,Offset);
    constexpr auto c=sfMultiply(sfMultiply(n,n),sfMultiply(n,n));
    if constexpr(Transpose){
        a=_mm512_xor_si512(a,b);b=_mm512_xor_si512(b,byteProduct4<Affine,c>(a));
    }else{
        a=_mm512_xor_si512(a,byteProduct4<Affine,c>(b));b=_mm512_xor_si512(a,b);
    }
}
template<unsigned Offset,bool Transpose,bool Affine=false>
HDS_INLINE void weightFft(__m512i& a,__m512i& b,__m512i& c,__m512i& d){
    if constexpr(Transpose){
        weightPair<0,Offset,true,Affine>(a,b);weightPair<0,Offset+2,true,Affine>(c,d);
        weightPair<1,Offset,true,Affine>(a,c);weightPair<1,Offset,true,Affine>(b,d);
    }else{
        weightPair<1,Offset,false,Affine>(a,c);weightPair<1,Offset,false,Affine>(b,d);
        weightPair<0,Offset,false,Affine>(a,b);weightPair<0,Offset+2,false,Affine>(c,d);
    }
}
template<bool Routed,class Index,unsigned Load>
HDS_INLINE __m512i fftInput4(const Block* input,const Index* route,unsigned offset){
    if constexpr(!Routed||Load==0)return hdbw::detail::input4<Routed>(input,route,offset);
    else if constexpr(Load==1){
        const auto a=_mm_loadu_si128(reinterpret_cast<const __m128i*>(input+route[offset]));
        const auto b=_mm_loadu_si128(reinterpret_cast<const __m128i*>(input+route[offset+1]));
        const auto c=_mm_loadu_si128(reinterpret_cast<const __m128i*>(input+route[offset+2]));
        const auto d=_mm_loadu_si128(reinterpret_cast<const __m128i*>(input+route[offset+3]));
        const auto low=_mm256_inserti128_si256(_mm256_castsi128_si256(a),b,1);
        const auto high=_mm256_inserti128_si256(_mm256_castsi128_si256(c),d,1);
        return _mm512_inserti64x4(_mm512_castsi256_si512(low),high,1);
    }else if constexpr(std::is_same_v<Index,std::uint16_t>){
        const auto packed=_mm_loadl_epi64(reinterpret_cast<const __m128i*>(route+offset));
        const auto duplicate=_mm_unpacklo_epi16(packed,packed);
        const auto word=_mm256_slli_epi32(_mm256_cvtepu16_epi32(duplicate),1);
        const auto slots=_mm256_or_si256(word,_mm256_setr_epi32(0,1,0,1,0,1,0,1));
        return _mm512_i32gather_epi64(slots,input,8);
    }else{
        const auto packed=_mm_loadu_si128(reinterpret_cast<const __m128i*>(route+offset));
        const auto duplicate=_mm256_permutevar8x32_epi32(_mm256_broadcastsi128_si256(packed),_mm256_setr_epi32(0,0,1,1,2,2,3,3));
        const auto word=_mm512_slli_epi64(_mm512_cvtepu32_epi64(duplicate),1);
        return _mm512_i64gather_epi64(_mm512_or_si512(word,_mm512_setr_epi64(0,1,0,1,0,1,0,1)),input,8);
    }
}
template<unsigned Half,unsigned Offset,bool Affine=false>
HDS_INLINE __m512i fixedTransposeStage(__m512i x){
    constexpr auto c0=normalized(Half==1?0:1,Offset),c1=normalized(Half==1?0:1,Offset^(Half==1?2:0));
    constexpr auto hi=Half==1?0xcc:0xf0;
    const auto a=_mm512_shuffle_i32x4(x,x,Half==1?0xa0:0x44);
    const auto b=_mm512_shuffle_i32x4(x,x,Half==1?0xf5:0xee);
    const auto sum=_mm512_xor_si512(a,b);
    if constexpr(c0==0&&c1==0)return _mm512_mask_mov_epi64(sum,hi,b);
    else{
        const auto product=[&](){
            if constexpr(Half==1)return byteProduct4<Affine,c0,c0,c1,c1>(sum);
            else return byteProduct4<Affine,c0>(sum);
        }();
        return _mm512_mask_xor_epi64(sum,hi,b,product);
    }
}
template<unsigned Tile,bool Routed,class Index,unsigned Load=0,bool Stream=false,bool Systematic=false,bool Affine=false>
HDS_INLINE __m512i fftTile(const Block* input,const Index* route,Block* output,
                          __m512i emission,const Fixed& f){
    auto a=fftInput4<Routed,Index,Load>(input,route,8*Tile);
    auto b=fftInput4<Routed,Index,Load>(input,route,8*Tile+4);
    __m512i feedback;
    if constexpr(Systematic){
        const auto h=_mm512_xor_si512(fixedTransposeStage<2,0,Affine>(fixedTransposeStage<1,0,Affine>(a)),
                                      fixedTransposeStage<2,4,Affine>(fixedTransposeStage<1,4,Affine>(b)));
        feedback=weight<0,Affine>(towerMoore<false,Affine>(h),f.maps,Tile);
        auto e=towerMoore<true,Affine>(weight<0,Affine>(emission,f.maps,Tile)),other=_mm512_setzero_si512();
        fixedTransform<0,false,Affine>(e,other);
        a=_mm512_xor_si512(a,e);b=_mm512_xor_si512(b,other);
    }else{
        fixedTransform<0,true,Affine>(a,b);
        feedback=weight<0,Affine>(towerMoore<false,Affine>(b),f.maps,Tile);
        a=_mm512_xor_si512(a,towerMoore<true,Affine>(weight<0,Affine>(emission,f.maps,Tile)));
        fixedTransform<8,false,Affine>(a,b);
    }
    if constexpr(Stream){
        _mm512_stream_si512(reinterpret_cast<__m512i*>(output+8*Tile),a);
        _mm512_stream_si512(reinterpret_cast<__m512i*>(output+8*Tile+4),b);
    }else{
        _mm512_storeu_si512(output+8*Tile,a);_mm512_storeu_si512(output+8*Tile+4,b);
    }
    return feedback;
}
template<unsigned Base,bool Routed,class Index,unsigned Load=0,bool Stream=false,bool Systematic=false,bool Affine=false>
HDS_INLINE State fftHalf(const Block* input,const Index* route,Block* output,
                         const State& q,const Fixed& f){
    auto a=q.values[0],b=q.values[1],c=q.values[2],d=q.values[3];
    weightFft<8+Base,false,Affine>(a,b,c,d);
    a=fftTile<Base,Routed,Index,Load,Stream,Systematic,Affine>(input,route,output,a,f);
    b=fftTile<Base+1,Routed,Index,Load,Stream,Systematic,Affine>(input,route,output,b,f);
    c=fftTile<Base+2,Routed,Index,Load,Stream,Systematic,Affine>(input,route,output,c,f);
    d=fftTile<Base+3,Routed,Index,Load,Stream,Systematic,Affine>(input,route,output,d,f);
    weightFft<8+Base,true,Affine>(a,b,c,d);
    return {{a,b,c,d}};
}
template<bool Routed,class Index,unsigned Load=0,bool Stream=false,bool Systematic=false,bool Affine=false>
HDS_INLINE void fftPacket(const Block* input,const Index* route,Block* output,State& state,
                         const Mixer& mixer,const DiagonalRefresh& t,const Fixed& f){
    auto feedback=fftHalf<0,Routed,Index,Load,Stream,Systematic,Affine>(input,route,output,state,f);
    const auto high=fftHalf<4,Routed,Index,Load,Stream,Systematic,Affine>(input,route,output,state,f);
    c8Unroll(std::make_index_sequence<4>{},[&](auto i){feedback.values[i]=_mm512_xor_si512(feedback.values[i],high.values[i]);});
    hdbw::detail::mix(state,mixer);feedback=refresh(feedback,t);
    c8Unroll(std::make_index_sequence<4>{},[&](auto i){state.values[i]=_mm512_xor_si512(state.values[i],feedback.values[i]);});
}
template<unsigned Offset,bool Routed,class Index,bool Affine=false>
HDS_INLINE void mdsTile(const Block* input,const Index* route,Block* output,
                        const State& q,State& feedback){
    using Kind=hdbw::detail::Kind;
    auto x0=hdbw::detail::input4<Routed>(input,route,Offset);
    auto x1=hdbw::detail::input4<Routed>(input,route,Offset+4);
    auto x2=hdbw::detail::input4<Routed>(input,route,Offset+8);
    auto x3=hdbw::detail::input4<Routed>(input,route,Offset+12);
    auto a=q.values[0],b=q.values[1],c=q.values[2],d=q.values[3];
    naturalFft16<Offset,Kind::Forward,Affine>(a,b,c,d);
    _mm512_storeu_si512(output+Offset,_mm512_xor_si512(x0,a));
    _mm512_storeu_si512(output+Offset+4,_mm512_xor_si512(x1,b));
    _mm512_storeu_si512(output+Offset+8,_mm512_xor_si512(x2,c));
    _mm512_storeu_si512(output+Offset+12,_mm512_xor_si512(x3,d));
    naturalFft16<Offset,Kind::ForwardTranspose,Affine>(x0,x1,x2,x3);
    feedback.values[0]=_mm512_xor_si512(feedback.values[0],x0);
    feedback.values[1]=_mm512_xor_si512(feedback.values[1],x1);
    feedback.values[2]=_mm512_xor_si512(feedback.values[2],x2);
    feedback.values[3]=_mm512_xor_si512(feedback.values[3],x3);
}
HDS_INLINE void transposeRecords4(__m512i& a,__m512i& b,__m512i& c,__m512i& d){
    const auto ab0=_mm512_shuffle_i32x4(a,b,0x44),ab1=_mm512_shuffle_i32x4(a,b,0xee);
    const auto cd0=_mm512_shuffle_i32x4(c,d,0x44),cd1=_mm512_shuffle_i32x4(c,d,0xee);
    a=_mm512_shuffle_i32x4(ab0,cd0,0x88);b=_mm512_shuffle_i32x4(ab0,cd0,0xdd);
    c=_mm512_shuffle_i32x4(ab1,cd1,0x88);d=_mm512_shuffle_i32x4(ab1,cd1,0xdd);
}

template<unsigned Level,bool Transpose,unsigned LocalWidth=64,unsigned Shift=0>
HDS_INLINE void stripedStage(__m512i* x){
    c8Unroll(std::make_index_sequence<8>{},[&](auto index){
        constexpr unsigned j=index,half=1U<<Level,a=(j/half)*(2*half)+j%half,b=a+half;
        constexpr unsigned offset=a&~(2*half-1);
        constexpr auto c0=hdbw::detail::normalized(Level,(offset%LocalWidth)^Shift),c1=hdbw::detail::normalized(Level,((offset+16)%LocalWidth)^Shift);
        constexpr auto c2=hdbw::detail::normalized(Level,((offset+32)%LocalWidth)^Shift),c3=hdbw::detail::normalized(Level,((offset+48)%LocalWidth)^Shift);
        if constexpr(Transpose){
            x[a]=_mm512_xor_si512(x[a],x[b]);
            x[b]=_mm512_xor_si512(x[b],byteProduct4<true,c0,c1,c2,c3>(x[a]));
        }else{
            x[a]=_mm512_xor_si512(x[a],byteProduct4<true,c0,c1,c2,c3>(x[b]));
            x[b]=_mm512_xor_si512(x[b],x[a]);
        }
    });
}

// Accumulate the two tower-product terms separately. Byte-pair swapping
// commutes with XOR, so swap the cross-term sum once instead of each of
// the four state vectors. This preserves the sampled coefficients exactly.
HDS_INLINE __m512i dotProduct16(const State& state,const Coefficient* coefficients){
#if defined(HD_SUMRANK_DELAYED_DOT)
    const auto crossSelect=_mm512_setr_epi64(
        0x0203020302030203LL,0x0203020302030203LL,
        0x0607060706070607LL,0x0607060706070607LL,
        0x0a0b0a0b0a0b0a0bLL,0x0a0b0a0b0a0b0a0bLL,
        0x0e0f0e0f0e0f0e0fLL,0x0e0f0e0f0e0f0e0fLL);
    __m512i direct[4],cross[4];
    c8Unroll(std::make_index_sequence<4>{},[&](auto i){
        const auto packed=_mm512_broadcast_i32x4(_mm_loadu_si128(
            reinterpret_cast<const __m128i*>(coefficients+4*i)));
        direct[i]=_mm512_gf2p8mul_epi8(state.values[i],wide16::expandPairs<0>(packed));
        cross[i]=_mm512_gf2p8mul_epi8(state.values[i],_mm512_shuffle_epi8(packed,crossSelect));
    });
    const auto a=_mm512_xor_si512(_mm512_xor_si512(direct[0],direct[1]),_mm512_xor_si512(direct[2],direct[3]));
    const auto b=_mm512_xor_si512(_mm512_xor_si512(cross[0],cross[1]),_mm512_xor_si512(cross[2],cross[3]));
    return hdbw::detail::repeatRecordSum(_mm512_xor_si512(a,medium32::selectBytes<1,0,3,2>(b)));
#else
    const auto a=wide16::multiply512(state.values[0],coefficients);
    const auto b=wide16::multiply512(state.values[1],coefficients+4);
    const auto c=wide16::multiply512(state.values[2],coefficients+8);
    const auto d=wide16::multiply512(state.values[3],coefficients+12);
    return hdbw::detail::repeatRecordSum(_mm512_xor_si512(_mm512_xor_si512(a,b),_mm512_xor_si512(c,d)));
#endif
}

template<bool Affine,bool Shared,unsigned Carried=0,bool RankTwo=false,class MixerData=Mixer>
HDS_INLINE void mdsFinish(State& state,State feedback,const MixerData& mixer,const DiagonalRefresh& t){
    if constexpr(Carried){
        // New ensemble: carry the state through the same Cauchy circuit
        // used for feedback. This is not the certified rank-one update.
        c8Unroll(std::make_index_sequence<4>{},[&](auto i){
            if constexpr(Carried>=2)
                feedback.values[i]=wide16::multiply512(_mm512_xor_si512(state.values[i],feedback.values[i]),t.before.data()+4*i);
            else feedback.values[i]=_mm512_xor_si512(state.values[i],wide16::multiply512(feedback.values[i],t.before.data()+4*i));
        });
        using Kind=hdbw::detail::Kind;
        naturalFft16<0,Kind::Inverse,Affine>(feedback.values[0],feedback.values[1],feedback.values[2],feedback.values[3]);
        naturalFft16<16,Kind::Forward,Affine>(feedback.values[0],feedback.values[1],feedback.values[2],feedback.values[3]);
        c8Unroll(std::make_index_sequence<4>{},[&](auto i){state.values[i]=wide16::multiply512(feedback.values[i],t.after.data()+4*i);});
    }else if constexpr(Shared){
        const auto dot=dotProduct16(state,mixer.v.data());
        __m512i second{};
        if constexpr(RankTwo){
            // Shared refresh does not read mixer.u. For this explicit
            // variant it stores the second row of V, not a rank-one u.
            second=dotProduct16(state,mixer.u.data());
        }
        c8Unroll(std::make_index_sequence<4>{},[&](auto i){feedback.values[i]=wide16::multiply512(feedback.values[i],t.before.data()+4*i);});
        using Kind=hdbw::detail::Kind;
        naturalFft16<0,Kind::Inverse,Affine>(feedback.values[0],feedback.values[1],feedback.values[2],feedback.values[3]);
        naturalFft16<16,Kind::Forward,Affine>(feedback.values[0],feedback.values[1],feedback.values[2],feedback.values[3]);
        c8Unroll(std::make_index_sequence<4>{},[&](auto i){
            // Setup ties mixer.u to this final diagonal. This is a new
            // joint ensemble, not an algebraic rewrite of independent maps.
            auto mixed=_mm512_xor_si512(feedback.values[i],dot);
            if constexpr(RankTwo){
                constexpr unsigned j=16+4*decltype(i)::value;
                mixed=_mm512_xor_si512(mixed,byteProduct4<true,j,j+1,j+2,j+3>(second));
            }
            state.values[i]=_mm512_xor_si512(state.values[i],wide16::multiply512(mixed,t.after.data()+4*i));
        });
    }else{
        hdbw::detail::mix(state,mixer);feedback=refresh<Affine>(feedback,t);
        c8Unroll(std::make_index_sequence<4>{},[&](auto i){state.values[i]=_mm512_xor_si512(state.values[i],feedback.values[i]);});
    }
}

// Distinct type prevents interpreting a byte offset as a record index.
// Precomputing the scale can remove one address-generation shift per load.
struct ByteOffset { std::uint32_t value; };
struct WordOffset { std::uint16_t value; };
#if defined(HD_SUMRANK_ROUTE_UNITS8)
using PreparedOffset=WordOffset;
#else
using PreparedOffset=ByteOffset;
#endif
template<bool Routed,class Index>
HDS_INLINE __m512i stripedInput4(const Block* input,const Index* route,unsigned offset){
    if constexpr(Routed&&(std::is_same_v<Index,ByteOffset>||std::is_same_v<Index,WordOffset>)){
        const auto* bytes=reinterpret_cast<const char*>(input);
        constexpr unsigned scale=std::is_same_v<Index,WordOffset>?8:1;
        auto x=_mm512_castsi128_si512(_mm_loadu_si128(reinterpret_cast<const __m128i*>(bytes+scale*route[offset].value)));
        x=_mm512_inserti32x4(x,_mm_loadu_si128(reinterpret_cast<const __m128i*>(bytes+scale*route[offset+1].value)),1);
        x=_mm512_inserti32x4(x,_mm_loadu_si128(reinterpret_cast<const __m128i*>(bytes+scale*route[offset+2].value)),2);
        return _mm512_inserti32x4(x,_mm_loadu_si128(reinterpret_cast<const __m128i*>(bytes+scale*route[offset+3].value)),3);
    }else return hdbw::detail::input4<Routed>(input,route,offset);
}

// Same RS maps, SIMD-striped over the four sixteen-point cosets. FFT
// butterflies are now lane-local; transpose only at the packet interface.
template<bool Routed,class Index,bool Shared=false,unsigned LocalRank=16,unsigned EmissionRank=16,unsigned Carried=0,bool RankTwo=false,bool Stream=false,bool StripedInput=false,class MixerData=Mixer>
HDS_INLINE void mdsStriped(const Block* input,const Index* route,Block* output,State& state,
                           const MixerData& mixer,const DiagonalRefresh& t,bool eightChecks=false){
    static_assert(LocalRank==0||LocalRank==4||LocalRank==8||LocalRank==16);
    static_assert(EmissionRank==16||EmissionRank==8||EmissionRank==4);
    static_assert(!StripedInput||Routed);
    __m512i x[16];
    if constexpr(EmissionRank==4){
        auto a=state.values[0],b=state.values[1],c=state.values[2],d=state.values[3];
        transposeRecords4(a,b,c,d);
        x[0]=a;x[1]=b;x[2]=c;x[3]=d;
        c8Unroll(std::make_index_sequence<4>{},[&](auto i){x[i+4]=x[i];x[i+8]=x[i];x[i+12]=x[i];});
        stripedStage<1,false,16>(x);stripedStage<0,false,16>(x);
    }else if constexpr(EmissionRank==8){
        // Two independent RS[32,8] emission maps, one per packet half.
        // Four SIMD records carry two cosets of each map. The omitted
        // degree-eight stage only replicates the eight coefficient vectors.
        c8Unroll(std::make_index_sequence<8>{},[&](auto i){
            constexpr unsigned r=i%4;
            x[i]=_mm512_shuffle_i32x4(state.values[i/4],state.values[2+i/4],_MM_SHUFFLE(r,r,r,r));
            x[i+8]=x[i];
        });
        stripedStage<2,false,32>(x);stripedStage<1,false,32>(x);stripedStage<0,false,32>(x);
    }else{
        c8Unroll(std::make_index_sequence<16>{},[&](auto i){x[i]=repeat<i%4>(state.values[i/4]);});
        stripedStage<3,false>(x);stripedStage<2,false>(x);stripedStage<1,false>(x);stripedStage<0,false>(x);
    }
    c8Unroll(std::make_index_sequence<4>{},[&](auto g){
        constexpr unsigned j=4*g;
        auto a=x[j],b=x[j+1],c=x[j+2],d=x[j+3];
        if constexpr(StripedInput){
            // Setup stores route[4*j+lane] = original[16*lane+j].
            // The gathers therefore produce the feedback's striped layout
            // directly. Only the completed output needs a register transpose.
            x[j]=stripedInput4<true>(input,route,4*j);
            x[j+1]=stripedInput4<true>(input,route,4*j+4);
            x[j+2]=stripedInput4<true>(input,route,4*j+8);
            x[j+3]=stripedInput4<true>(input,route,4*j+12);
            a=_mm512_xor_si512(a,x[j]);b=_mm512_xor_si512(b,x[j+1]);
            c=_mm512_xor_si512(c,x[j+2]);d=_mm512_xor_si512(d,x[j+3]);
            transposeRecords4(a,b,c,d);
            if constexpr(Stream){
                _mm512_stream_si512(reinterpret_cast<__m512i*>(output+j),a);
                _mm512_stream_si512(reinterpret_cast<__m512i*>(output+16+j),b);
                _mm512_stream_si512(reinterpret_cast<__m512i*>(output+32+j),c);
                _mm512_stream_si512(reinterpret_cast<__m512i*>(output+48+j),d);
            }else{
                _mm512_storeu_si512(output+j,a);_mm512_storeu_si512(output+16+j,b);
                _mm512_storeu_si512(output+32+j,c);_mm512_storeu_si512(output+48+j,d);
            }
        }else{
            transposeRecords4(a,b,c,d);
            auto v0=stripedInput4<Routed>(input,route,j);
            auto v1=stripedInput4<Routed>(input,route,16+j);
            auto v2=stripedInput4<Routed>(input,route,32+j);
            auto v3=stripedInput4<Routed>(input,route,48+j);
            a=_mm512_xor_si512(a,v0);b=_mm512_xor_si512(b,v1);
            c=_mm512_xor_si512(c,v2);d=_mm512_xor_si512(d,v3);
            // Store before forming the feedback transpose. Keeping the
            // output registers live across it regressed the frozen path.
            if constexpr(Stream){
                _mm512_stream_si512(reinterpret_cast<__m512i*>(output+j),a);
                _mm512_stream_si512(reinterpret_cast<__m512i*>(output+16+j),b);
                _mm512_stream_si512(reinterpret_cast<__m512i*>(output+32+j),c);
                _mm512_stream_si512(reinterpret_cast<__m512i*>(output+48+j),d);
            }else{
                _mm512_storeu_si512(output+j,a);_mm512_storeu_si512(output+16+j,b);
                _mm512_storeu_si512(output+32+j,c);_mm512_storeu_si512(output+48+j,d);
            }
            transposeRecords4(v0,v1,v2,v3);
            x[j]=v0;x[j+1]=v1;x[j+2]=v2;x[j+3]=v3;
        }
    });
    State feedback;
    if constexpr(LocalRank==0){
        // Both block codes use the actual coordinate nodes. Translating
        // a block's nodes changes its polynomial basis, not its MDS ranks.
        // This shares emission, two feedback stages, and refresh; only the
        // short rank-selection tail branches.
        stripedStage<0,true>(x);stripedStage<1,true>(x);
        if(eightChecks){
            stripedStage<2,true>(x);
            c8Unroll(std::make_index_sequence<8>{},[&](auto j){x[j]=_mm512_xor_si512(x[j],x[j+8]);});
            c8Unroll(std::make_index_sequence<2>{},[&](auto g){
                auto a=x[4*g],b=x[4*g+1],c=x[4*g+2],d=x[4*g+3];transposeRecords4(a,b,c,d);
                feedback.values[g]=_mm512_xor_si512(a,b);feedback.values[2+g]=_mm512_xor_si512(c,d);
            });
        }else{
            c8Unroll(std::make_index_sequence<4>{},[&](auto j){feedback.values[j]=_mm512_xor_si512(_mm512_xor_si512(x[j],x[j+4]),_mm512_xor_si512(x[j+8],x[j+12]));});
            transposeRecords4(feedback.values[0],feedback.values[1],feedback.values[2],feedback.values[3]);
        }
    }else if constexpr(LocalRank==16){
        stripedStage<0,true>(x);stripedStage<1,true>(x);stripedStage<2,true>(x);stripedStage<3,true>(x);
        c8Unroll(std::make_index_sequence<4>{},[&](auto g){
            auto a=x[4*g],b=x[4*g+1],c=x[4*g+2],d=x[4*g+3];
            transposeRecords4(a,b,c,d);
            feedback.values[g]=_mm512_xor_si512(_mm512_xor_si512(a,b),_mm512_xor_si512(c,d));
        });
    }else if constexpr(LocalRank==8){
        stripedStage<0,true,32>(x);stripedStage<1,true,32>(x);stripedStage<2,true,32>(x);
        c8Unroll(std::make_index_sequence<8>{},[&](auto j){x[j]=_mm512_xor_si512(x[j],x[j+8]);});
        c8Unroll(std::make_index_sequence<2>{},[&](auto g){
            auto a=x[4*g],b=x[4*g+1],c=x[4*g+2],d=x[4*g+3];
            transposeRecords4(a,b,c,d);
            feedback.values[g]=_mm512_xor_si512(a,b);
            feedback.values[2+g]=_mm512_xor_si512(c,d);
        });
    }else{
        stripedStage<0,true,16,Carried==3?16:0>(x);stripedStage<1,true,16,Carried==3?16:0>(x);
        c8Unroll(std::make_index_sequence<4>{},[&](auto j){
            feedback.values[j]=_mm512_xor_si512(_mm512_xor_si512(x[j],x[j+4]),_mm512_xor_si512(x[j+8],x[j+12]));
        });
        if constexpr(Carried!=3)transposeRecords4(feedback.values[0],feedback.values[1],feedback.values[2],feedback.values[3]);
    }
    mdsFinish<true,Shared,Carried,RankTwo>(state,feedback,mixer,t);
}

template<unsigned Offset,bool Routed,class Index>
HDS_INLINE void blockFeedbackTile(const Block* input,const Index* route,Block* output,const State& state,State& feedback,bool eightChecks){
    using Kind=hdbw::detail::Kind;
    auto x0=hdbw::detail::input4<Routed>(input,route,Offset),x1=hdbw::detail::input4<Routed>(input,route,Offset+4);
    auto x2=hdbw::detail::input4<Routed>(input,route,Offset+8),x3=hdbw::detail::input4<Routed>(input,route,Offset+12);
    auto a=state.values[0],b=state.values[1],c=state.values[2],d=state.values[3];
    naturalFft16<Offset,Kind::Forward,true>(a,b,c,d);
    _mm512_storeu_si512(output+Offset,_mm512_xor_si512(x0,a));_mm512_storeu_si512(output+Offset+4,_mm512_xor_si512(x1,b));
    _mm512_storeu_si512(output+Offset+8,_mm512_xor_si512(x2,c));_mm512_storeu_si512(output+Offset+12,_mm512_xor_si512(x3,d));
    x0=affineSmall<1,Offset,Kind::ForwardTranspose>(x0);x1=affineSmall<1,Offset+4,Kind::ForwardTranspose>(x1);
    x2=affineSmall<1,Offset+8,Kind::ForwardTranspose>(x2);x3=affineSmall<1,Offset+12,Kind::ForwardTranspose>(x3);
    x0=affineSmall<2,Offset,Kind::ForwardTranspose>(x0);x1=affineSmall<2,Offset+4,Kind::ForwardTranspose>(x1);
    x2=affineSmall<2,Offset+8,Kind::ForwardTranspose>(x2);x3=affineSmall<2,Offset+12,Kind::ForwardTranspose>(x3);
    if(eightChecks){
        affinePair<hdbw::detail::normalized(2,Offset),Kind::ForwardTranspose>(x0,x1);
        affinePair<hdbw::detail::normalized(2,Offset+8),Kind::ForwardTranspose>(x2,x3);
        constexpr unsigned group=2*(Offset/32);
        feedback.values[group]=_mm512_xor_si512(feedback.values[group],_mm512_xor_si512(x0,x2));
        feedback.values[group+1]=_mm512_xor_si512(feedback.values[group+1],_mm512_xor_si512(x1,x3));
    }else feedback.values[Offset/16]=_mm512_xor_si512(_mm512_xor_si512(x0,x1),_mm512_xor_si512(x2,x3));
}

template<bool Routed,class Index>
HDS_INLINE void blockFeedbackPacket(const Block* input,const Index* route,Block* output,State& state,
                                    const Mixer& mixer,const DiagonalRefresh& t,bool eightChecks){
    State feedback{};
    blockFeedbackTile<0,Routed>(input,route,output,state,feedback,eightChecks);
    blockFeedbackTile<16,Routed>(input,route,output,state,feedback,eightChecks);
    blockFeedbackTile<32,Routed>(input,route,output,state,feedback,eightChecks);
    blockFeedbackTile<48,Routed>(input,route,output,state,feedback,eightChecks);
    mdsFinish<true,true>(state,feedback,mixer,t);
}

// Research dispatch for alternating feedback. The branch is once per
// 64-symbol packet; retaining the explicit kernels also supports regions
// with an odd packet count. It is included in the end-to-end timings.
template<unsigned Variant,bool Routed,class Index>
HDS_INLINE void stripedByStep(const Block* input,const Index* route,Block* output,State& state,
                              const Mixer& mixer,const DiagonalRefresh& t,std::size_t step){
    if constexpr(Variant>=41&&Variant<=45){
        constexpr bool globalEmission=Variant==42||Variant==44||Variant==45;
        constexpr bool prepared=Routed&&(Variant==43||Variant==44);
        const bool eight=Variant==45?(step&3)==3:bool(step&1);
        std::array<Index,prepared?64:0> reordered;
        if constexpr(prepared){
            for(unsigned j=0;j!=16;++j)for(unsigned lane=0;lane!=4;++lane)
                reordered[4*j+lane]=route[16*lane+j];
            route=reordered.data();
        }
        // Standalone packet tests include unaligned destinations. The full
        // hot loop dispatches once outside the loop, avoiding these checks.
        if(reinterpret_cast<std::uintptr_t>(output)&63U)
            mdsStriped<Routed,Index,true,globalEmission?0:8,globalEmission?16:4,0,!globalEmission,false,prepared>(input,route,output,state,mixer,t,eight);
        else{
            mdsStriped<Routed,Index,true,globalEmission?0:8,globalEmission?16:4,0,!globalEmission,true,prepared>(input,route,output,state,mixer,t,eight);
            _mm_sfence();
        }
    }else if constexpr(Variant>=38)mdsStriped<Routed,Index,true,8,4,0,true>(input,route,output,state,mixer,t);
    else if constexpr(Variant==37)mdsStriped<Routed,Index,true,0,4,0,true>(input,route,output,state,mixer,t,bool(step&1));
    else if constexpr(Variant==36)mdsStriped<Routed,Index,true,4,4,3>(input,route,output,state,mixer,t);
    else if constexpr(Variant==35)mdsStriped<Routed,Index,true,4,4,2>(input,route,output,state,mixer,t);
    else if constexpr(Variant==34)mdsStriped<Routed,Index,true,4,4,1>(input,route,output,state,mixer,t);
    else if constexpr(Variant==33)mdsStriped<Routed,Index,true,4,4>(input,route,output,state,mixer,t);
    else if constexpr(Variant==32)mdsStriped<Routed,Index,true,8,8>(input,route,output,state,mixer,t);
    else if constexpr(Variant==30)blockFeedbackPacket<Routed>(input,route,output,state,mixer,t,bool(step&1));
    else if constexpr(Variant==29||Variant==31)mdsStriped<Routed,Index,true,0>(input,route,output,state,mixer,t,bool(step&1));
    else if constexpr(Variant==28){
        if(step&1)mdsStriped<Routed,Index,true,8>(input,route,output,state,mixer,t);
        else mdsStriped<Routed,Index,true,4>(input,route,output,state,mixer,t);
    }else mdsStriped<Routed,Index,Variant>=24,Variant==26?8:Variant==27?4:16>(input,route,output,state,mixer,t);
}

// Scheduling control: a tile boundary bounds live SIMD temporaries without
// changing the map. Keep the monolithic kernel as the performance anchor.
template<unsigned Offset,bool Routed,class Index,bool Affine>
#if defined(_MSC_VER)
__declspec(noinline)
#else
__attribute__((noinline))
#endif
void mdsTileBounded(const Block* input,const Index* route,Block* output,
                    const State& q,State& feedback){
    mdsTile<Offset,Routed,Index,Affine>(input,route,output,q,feedback);
}

template<bool Routed,class Index,bool Affine=false,bool Bounded=false,bool Shared=false>
HDS_INLINE void mdsPacket(const Block* input,const Index* route,Block* output,State& state,
                         const Mixer& mixer,const DiagonalRefresh& t){
    State feedback{};
    if constexpr(Bounded){
        mdsTileBounded<0,Routed,Index,Affine>(input,route,output,state,feedback);
        mdsTileBounded<16,Routed,Index,Affine>(input,route,output,state,feedback);
        mdsTileBounded<32,Routed,Index,Affine>(input,route,output,state,feedback);
        mdsTileBounded<48,Routed,Index,Affine>(input,route,output,state,feedback);
    }else{
        mdsTile<0,Routed,Index,Affine>(input,route,output,state,feedback);
        mdsTile<16,Routed,Index,Affine>(input,route,output,state,feedback);
        mdsTile<32,Routed,Index,Affine>(input,route,output,state,feedback);
        mdsTile<48,Routed,Index,Affine>(input,route,output,state,feedback);
    }
    mdsFinish<Affine,Shared>(state,feedback,mixer,t);
}
#undef HDS_INLINE
} // namespace spin::research::char2::hdsumrank
