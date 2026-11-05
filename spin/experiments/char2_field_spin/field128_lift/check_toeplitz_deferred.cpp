// Correctness only; no benchmark. From this directory on an ISA-capable host:
// c++ -std=c++20 -O3 -mavx512f -mavx512bw -mgfni -mvpclmulqdq \
//   check_toeplitz_deferred.cpp -o /tmp/check_toeplitz_deferred
// /tmp/check_toeplitz_deferred
// Runtime dispatch cannot make a globally ISA-enabled binary safe on old CPUs.
#include "SparseByteField128ToeplitzDeferred.h"
#include <algorithm>
#include <array>
#include <cstring>
#include <iostream>
#include <stdexcept>

namespace f=spin::research::char2::field128_lift::sparse_byte;
namespace d=f::deferred;
using Matrix=std::array<f::Element,7>;
using Keys=std::array<f::Element,28>;
using Values=std::array<f::Element,16>;
using Raw=std::array<std::uint16_t,32>;
static void require(bool x,const char* s){if(!x)throw std::runtime_error(s);}
static std::uint64_t word(){static std::uint64_t s=0xd3fe44ed123987ULL;s^=s<<13;s^=s>>7;s^=s<<17;return s;}
static f::Element element(){return {word(),word()};}

// Independent bit-polynomial convolution, before any field reduction.
static Raw rawOracle(f::Element a,f::Element b){
    Raw raw{};
    for(unsigned i=0;i!=16;++i)for(unsigned j=0;j!=16;++j){
        const unsigned left=f::coefficient(a,i),right=f::coefficient(b,j);
        for(unsigned k=0;k!=8;++k)if((right>>k)&1U)raw[i+j]^=std::uint16_t(left<<k);
    }
    return raw;
}
static Keys keys(const std::array<Matrix,4>& t){
    Keys p{};
    for(unsigned lane=0;lane!=4;++lane){
        const auto& a=t[lane];
        p[lane]=a[3];p[4+lane]=a[3]^a[2];p[8+lane]=a[3]^a[4];
        p[12+lane]=a[5]^a[3];p[16+lane]=a[5]^a[3]^a[6]^a[4];
        p[20+lane]=a[1]^a[3];p[24+lane]=a[1]^a[3]^a[0]^a[2];
    }
    return p;
}
// Inverse of the compact linear key layout, used only when deliberately
// overlapping arbitrary key/input bytes. No invertibility assumption on T.
static std::array<Matrix,4> unpackKeys(const Keys& p){
    std::array<Matrix,4> t{};
    for(unsigned lane=0;lane!=4;++lane){
        auto& a=t[lane];
        a[3]=p[lane];a[2]=p[lane]^p[4+lane];a[4]=p[lane]^p[8+lane];
        a[5]=p[12+lane]^a[3];a[6]=p[16+lane]^p[12+lane]^a[4];
        a[1]=p[20+lane]^a[3];a[0]=p[24+lane]^p[20+lane]^a[2];
    }
    return t;
}
static Values oracle(const Values& x,const std::array<Matrix,4>& t){
    Values out{};
    for(unsigned i=0;i!=4;++i)for(unsigned j=0;j!=4;++j)for(unsigned lane=0;lane!=4;++lane)
        out[4*i+lane]^=f::multiplyOracle(t[lane][3+i-j],x[4*j+lane]);
    return out;
}

// Frozen nine-full-product arithmetic as an extra control, not the oracle.
static Values reducedControl(const Values& input,const Keys& p){
    alignas(64) __m512i x[4];std::memcpy(x,input.data(),sizeof(x));
    const auto a0=_mm512_loadu_si512(p.data()),a1=_mm512_loadu_si512(p.data()+4),a2=_mm512_loadu_si512(p.data()+8);
    const auto keyXor=_mm512_xor_si512(a1,a2),low=_mm512_xor_si512(x[0],x[2]),high=_mm512_xor_si512(x[1],x[3]);
    const auto sharedU=f::multiply4(_mm512_xor_si512(low,high),a0);
    const auto u0=_mm512_xor_si512(sharedU,f::multiply4(high,a1)),u1=_mm512_xor_si512(sharedU,f::multiply4(low,a2));
    const auto v0=_mm512_loadu_si512(p.data()+12),v2=_mm512_loadu_si512(p.data()+16);
    const auto sharedV=f::multiply4(_mm512_xor_si512(x[0],x[1]),v0);
    const auto outV0=_mm512_xor_si512(sharedV,f::multiply4(x[1],_mm512_xor_si512(v0,keyXor)));
    const auto outV1=_mm512_xor_si512(sharedV,f::multiply4(x[0],v2));
    const auto w0=_mm512_loadu_si512(p.data()+20),w1=_mm512_loadu_si512(p.data()+24);
    const auto sharedW=f::multiply4(_mm512_xor_si512(x[2],x[3]),w0);
    const auto outW0=_mm512_xor_si512(sharedW,f::multiply4(x[3],w1));
    const auto outW1=_mm512_xor_si512(sharedW,f::multiply4(x[2],_mm512_xor_si512(w0,keyXor)));
    x[0]=_mm512_xor_si512(u0,outW0);x[1]=_mm512_xor_si512(u1,outW1);
    x[2]=_mm512_xor_si512(u0,outV0);x[3]=_mm512_xor_si512(u1,outV1);
    Values out;std::memcpy(out.data(),x,sizeof(x));return out;
}
static void check(const Values& x,const std::array<Matrix,4>& t){
    const auto p=keys(t);
    const auto expected=oracle(x,t);
    require(unpackKeys(p)==t,"compact key layout round trip");
    Values out;d::toeplitz4Records(x.data(),p.data(),out.data());
    require(out==expected,"deferred versus schoolbook Toeplitz");
    require(out==reducedControl(x,p),"deferred versus nine reduced products");
}
static void rawCheck(){
    alignas(64) f::Element a[4],b[4],out[4];
    alignas(64) std::uint16_t slots[4][32];
    for(unsigned i=0;i!=128;++i)for(unsigned j=0;j!=128;j+=4){
        for(unsigned lane=0;lane!=4;++lane){a[lane]=f::binaryBasis(i);b[lane]=f::binaryBasis(j+lane);}
        const auto product=d::multiplyRaw4(_mm512_loadu_si512(a),_mm512_loadu_si512(b));
        _mm512_store_si512(slots[0],product.p0);_mm512_store_si512(slots[1],product.p1);
        _mm512_store_si512(slots[2],product.p2);_mm512_store_si512(slots[3],product.p3);
        _mm512_store_si512(out,d::reduce(product));
        for(unsigned lane=0;lane!=4;++lane){
            const auto expected=rawOracle(a[lane],b[lane]);
            for(unsigned k=0;k!=32;++k)require(slots[k/8][8*lane+k%8]==expected[k],"unreduced binary convolution");
            require(out[lane]==f::multiplyOracle(a[lane],b[lane]),"reduced binary product");
        }
    }
    for(unsigned trial=0;trial!=64;++trial){
        f::Element c[4],e[4],g[4],h[4];
        for(unsigned lane=0;lane!=4;++lane){a[lane]=element();b[lane]=element();c[lane]=element();e[lane]=element();g[lane]=element();h[lane]=element();}
        const auto p=d::multiplyRaw4(_mm512_loadu_si512(a),_mm512_loadu_si512(b));
        const auto q=d::multiplyRaw4(_mm512_loadu_si512(c),_mm512_loadu_si512(e));
        const auto r=d::multiplyRaw4(_mm512_loadu_si512(g),_mm512_loadu_si512(h));
        _mm512_store_si512(out,d::reduceSum3(p,q,r));
        for(unsigned lane=0;lane!=4;++lane)require(out[lane]==(f::multiplyOracle(a[lane],b[lane])^
            f::multiplyOracle(c[lane],e[lane])^f::multiplyOracle(g[lane],h[lane])),"linear deferred reduction");
    }
}
static void matrixCheck(){
    Values x{};std::array<Matrix,4> t{};check(x,t);
    for(auto& v:x)v=element();check(x,t);
    for(auto& m:t)m[3]={1,0};check(x,t);
    for(unsigned trial=0;trial!=128;++trial){
        for(auto& v:x)v=element();for(auto& m:t)for(auto& v:m)v=element();
        if(trial==0)x={};if(trial==1)for(auto& v:x)v={~std::uint64_t(0),~std::uint64_t(0)};
        check(x,t);
    }
    for(unsigned diagonal=0;diagonal!=7;++diagonal)for(unsigned row=0;row!=4;++row)
        for(unsigned bit=0;bit!=128;++bit){
            x={};t={};
            for(unsigned lane=0;lane!=4;++lane){
                t[lane][diagonal]=f::binaryBasis((37*bit+17*lane)%128);
                x[4*row+lane]=f::binaryBasis(bit);
            }
            check(x,t);
        }
    // Full-field coefficient/input superposition: no restriction to AES keys.
    for(unsigned bit=0;bit!=128;++bit){
        for(auto& m:t)for(auto& v:m)v=element();
        for(unsigned i=0;i!=16;++i)x[i]=f::binaryBasis((bit+13*i)%128);
        check(x,t);
    }
}
static void aliasCheck(){
    Values x;std::array<Matrix,4> t;
    for(auto& v:x)v=element();for(auto& m:t)for(auto& v:m)v=element();
    const auto p=keys(t);
    const auto expected=oracle(x,t);
    for(int displacement=-31;displacement<=31;++displacement){
        std::array<unsigned char,768> data;data.fill(0xa7);
        std::memcpy(data.data()+128,x.data(),sizeof(x));const auto before=data;
        const auto offset=128+displacement;
        d::toeplitz4Records(data.data()+128,p.data(),data.data()+offset);
        require(std::memcmp(data.data()+offset,expected.data(),sizeof(expected))==0,"partial input/output alias");
        for(unsigned i=0;i!=data.size();++i)if(i<unsigned(offset)||i>=unsigned(offset)+sizeof(x))
            require(data[i]==before[i],"alias output bounds");
    }
    for(unsigned offset=1;offset<=63;++offset){
        std::array<unsigned char,768> data;data.fill(0xb4);
        std::memcpy(data.data()+offset,p.data(),sizeof(p));const auto before=data;
        d::toeplitz4Records(x.data(),data.data()+offset,data.data()+offset+17);
        require(std::memcmp(data.data()+offset+17,expected.data(),sizeof(expected))==0,"unaligned key/output alias");
        for(unsigned i=0;i!=data.size();++i)if(i<offset+17||i>=offset+17+sizeof(x))
            require(data[i]==before[i],"key alias output bounds");
    }
    alignas(64) __m512i hot[4];std::memcpy(hot,x.data(),sizeof(hot));
    d::toeplitz4(hot,p.data());require(std::memcmp(hot,expected.data(),sizeof(hot))==0,"hot in-place interface");
    // All three ranges overlap. Expected values use the bytes present at call
    // entry, not the original disjoint x/p objects from earlier tests.
    for(unsigned displacement=0;displacement!=64;++displacement){
        std::array<unsigned char,768> data;
        for(auto& byte:data)byte=static_cast<unsigned char>(word());
        Values actualInput;Keys actualKeys;
        std::memcpy(actualInput.data(),data.data()+96,sizeof(actualInput));
        std::memcpy(actualKeys.data(),data.data()+32,sizeof(actualKeys));
        const auto combinedExpected=oracle(actualInput,unpackKeys(actualKeys));
        const auto before=data;
        d::toeplitz4Records(data.data()+96,data.data()+32,data.data()+64+displacement);
        require(std::memcmp(data.data()+64+displacement,combinedExpected.data(),sizeof(combinedExpected))==0,
                "simultaneous input/key/output overlap");
        for(unsigned i=0;i!=data.size();++i)if(i<64+displacement||i>=320+displacement)
            require(data[i]==before[i],"combined alias output bounds");
    }
    alignas(64) __m512i sharedHot[7];
    std::memcpy(sharedHot,p.data(),sizeof(p));
    Values sharedInput;std::memcpy(sharedInput.data(),sharedHot,sizeof(sharedInput));
    const auto sharedExpected=oracle(sharedInput,t);
    d::toeplitz4(sharedHot,sharedHot);
    require(std::memcmp(sharedHot,sharedExpected.data(),sizeof(sharedExpected))==0,"hot key/input alias");
    require(std::memcmp(reinterpret_cast<const unsigned char*>(sharedHot)+sizeof(sharedExpected),
                        reinterpret_cast<const unsigned char*>(p.data())+sizeof(sharedExpected),
                        sizeof(p)-sizeof(sharedExpected))==0,"hot key suffix preserved");
}
int main(){
    try{rawCheck();matrixCheck();aliasCheck();
        std::cout<<"PASS: raw basis/dense convolution, independent Toeplitz oracle, reduced control, zero/identity/full-field keys, aliases/guards\n";
        return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
