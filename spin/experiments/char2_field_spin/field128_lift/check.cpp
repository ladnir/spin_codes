// Standalone correctness check, not a benchmark. Portable build:
// c++ -std=c++20 -O2 check.cpp -o field128_lift_check
// Optional SIMD check: add -DFIELD128_LIFT_SIMD_CHECK=1 -mavx512f
// -mavx512bw -mvpclmulqdq, and run only on a host with that ISA.
#include "Field128.h"
#ifndef FIELD128_LIFT_SIMD_CHECK
#define FIELD128_LIFT_SIMD_CHECK 0
#endif
#if FIELD128_LIFT_SIMD_CHECK
#include "Field128Simd.h"
#endif
#include <algorithm>
#include <cstring>
#include <iostream>
#include <stdexcept>

namespace f=spin::research::char2::field128_lift;
static void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
static std::uint64_t randomWord(){
    static std::uint64_t state=0x128f13c0ffeeULL;
    state^=state<<13;state^=state>>7;state^=state<<17;return state;
}
static f::Scalar randomScalar(){return {randomWord(),randomWord()};}
static std::uint8_t byteMultiply(std::uint8_t a,std::uint8_t b){
    unsigned product=0;
    for(unsigned bit=0;bit!=8;++bit)if(b>>bit&1)product^=unsigned(a)<<bit;
    for(unsigned bit=15;bit-- >8;)if(product>>bit&1)product^=0x11bU<<(bit-8);
    return static_cast<std::uint8_t>(product);
}
static f::Scalar bytewise(f::Scalar x,std::uint8_t coefficient){
    f::Scalar result{};
    for(unsigned j=0;j!=8;++j){
        result.lo|=std::uint64_t(byteMultiply(static_cast<std::uint8_t>(x.lo>>(8*j)),coefficient))<<(8*j);
        result.hi|=std::uint64_t(byteMultiply(static_cast<std::uint8_t>(x.hi>>(8*j)),coefficient))<<(8*j);
    }
    return result;
}
static void checkScalar(){
    for(unsigned a=0;a!=128;++a)for(unsigned b=0;b!=128;++b)
        require(f::multiply(f::binaryBasis(a),f::binaryBasis(b))==
                f::multiplyOracle(f::binaryBasis(a),f::binaryBasis(b)),"native binary product basis");
    for(unsigned trial=0;trial!=128;++trial){
        const auto a=randomScalar(),b=randomScalar();
        require(f::multiply(a,b)==f::multiplyOracle(a,b),"native random product");
    }
    for(unsigned a=0;a!=256;++a)for(unsigned b=0;b!=256;++b){
        require((f::embeddedBytes[a]^f::embeddedBytes[b])==f::embeddedBytes[a^b],"embedding sum");
        require(f::multiplyOracle(f::embeddedBytes[a],f::embeddedBytes[b])==
                f::embeddedBytes[byteMultiply(a,b)],"embedding product");
    }
    const f::ByteBasis basis;
    for(unsigned bit=0;bit!=128;++bit){
        const auto value=f::binaryBasis(bit);
        require(basis.fromNative(basis.toNative(value))==value,"byte basis roundtrip");
        require(basis.toNative(basis.fromNative(value))==value,"native basis roundtrip");
        for(unsigned c=0;c!=8;++c){
            const auto coefficient=std::uint8_t(1U<<c);
            require(basis.toNative(bytewise(value,coefficient))==
                    f::multiplyOracle(basis.toNative(value),f::embeddedBytes[coefficient]),
                    "converted byte GFNI action");
        }
    }
    for(unsigned c=0;c!=256;++c){
        const auto value=randomScalar();
        require(basis.fromNative(f::multiplyOracle(basis.toNative(value),f::embeddedBytes[c]))==
                bytewise(value,static_cast<std::uint8_t>(c)),"dense converted byte action");
    }
    const f::Scalar x{1,0},a{0x80,0};
    require(bytewise(x,2)!=f::embeddedBytes[2],"negative control: raw byte is not embedded");
    require(bytewise(f::multiply(a,x),2)!=f::multiply(a,bytewise(x,2)),
            "negative control: raw byte GFNI is not GF128-linear");
    std::cout<<"scalar: 16384 binary products, all65536 embedded byte products, "
                 "rank128 basis/roundtrips, converted byte actions, negative controls pass\n";
}

#if FIELD128_LIFT_SIMD_CHECK
static void checkSimd(){
    alignas(64) f::Scalar input[4],coefficients[4],output[4],expected[4];
    for(unsigned left=0;left!=128;++left)for(unsigned right=0;right!=128;right+=4){
        for(unsigned lane=0;lane!=4;++lane){
            input[lane]=f::binaryBasis((left+17*lane)%128);
            coefficients[lane]=f::binaryBasis(right+lane);
            expected[lane]=f::multiplyOracle(input[lane],coefficients[lane]);
        }
        f::multiply4Records(input,coefficients,output);
        require(std::equal(output,output+4,expected),"SIMD binary product basis");
    }
    for(unsigned trial=0;trial!=128;++trial){
        for(unsigned lane=0;lane!=4;++lane){
            input[lane]=randomScalar();coefficients[lane]=randomScalar();
            if(trial<2)coefficients[lane]={trial,0};
            expected[lane]=f::multiplyOracle(input[lane],coefficients[lane]);
        }
        f::multiply4Records(input,coefficients,output);
        require(std::equal(output,output+4,expected),"SIMD independent random coefficients");
        f::multiply4Records(input,coefficients,input);
        require(std::equal(input,input+4,expected),"SIMD exact input alias");
    }
    // Vary every legal 8-byte Scalar alignment modulo64 and key alignment.
    // Byte storage plus memcpy permits aligned Scalar objects at each offset.
    alignas(64) std::array<std::uint8_t,192> source,keys,result;
    constexpr std::uint8_t guard=0xa7;
    for(unsigned offset=0;offset!=64;offset+=8)for(unsigned keyOffset=0;keyOffset!=64;keyOffset+=8){
        source.fill(guard);keys.fill(guard);result.fill(guard);
        for(unsigned lane=0;lane!=4;++lane){
            const auto value=randomScalar(),coefficient=randomScalar();
            std::memcpy(source.data()+64+offset+16*lane,&value,sizeof(value));
            std::memcpy(keys.data()+64+keyOffset+16*lane,&coefficient,sizeof(coefficient));
            expected[lane]=f::multiplyOracle(value,coefficient);
        }
        f::multiply4Records(reinterpret_cast<const f::Scalar*>(source.data()+64+offset),
                           reinterpret_cast<const f::Scalar*>(keys.data()+64+keyOffset),
                           reinterpret_cast<f::Scalar*>(result.data()+64+offset));
        std::memcpy(output,result.data()+64+offset,sizeof(output));
        require(std::equal(output,output+4,expected),"SIMD unaligned product");
        for(unsigned i=0;i!=result.size();++i)if(i<64+offset||i>=128+offset)
            require(result[i]==guard,"SIMD output guard");
    }
    std::cout<<"SIMD: all16384 binary products, independent lane coefficients, "
                 "zero/one, exact alias, alignments and output guards pass\n";
}
#endif

int main(){
    try{
#if FIELD128_LIFT_SIMD_CHECK && (defined(__GNUC__) || defined(__clang__))
        if(!__builtin_cpu_supports("avx512f")||!__builtin_cpu_supports("avx512bw")||
           !__builtin_cpu_supports("vpclmulqdq")){
            std::cout<<"SKIP: SIMD check requires AVX512F/BW and VPCLMULQDQ\n";return 77;
        }
#endif
        checkScalar();
#if FIELD128_LIFT_SIMD_CHECK
        checkSimd();
#endif
        return 0;
    }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
