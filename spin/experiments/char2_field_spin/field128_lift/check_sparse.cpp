// Standalone correctness check, no benchmark. Portable:
// c++ -std=c++20 -O2 check_sparse.cpp -o field128_sparse_check
// SIMD: add -DFIELD128_LIFT_SIMD_CHECK=1 -mavx512f -mavx512bw -mgfni
// -mvpclmulqdq. Run the SIMD build only on a host with that ISA.
#include "SparseByteField128.h"
#ifndef FIELD128_LIFT_SIMD_CHECK
#define FIELD128_LIFT_SIMD_CHECK 0
#endif
#if FIELD128_LIFT_SIMD_CHECK
#include "SparseByteField128Simd.h"
#endif
#include <algorithm>
#include <cstring>
#include <iostream>
#include <stdexcept>

namespace f=spin::research::char2::field128_lift::sparse_byte;
static void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
static std::uint64_t randomWord(){
    static std::uint64_t state=0xf1283216c0ffeeULL;
    state^=state<<13;state^=state>>7;state^=state<<17;return state;
}
static f::Element randomElement(){return {randomWord(),randomWord()};}
static f::Element bytewise(f::Element x,std::uint8_t c){
    std::array<std::uint8_t,16> result{};
    for(unsigned j=0;j!=16;++j)result[j]=f::byteMultiply(f::coefficient(x,j),c);
    return f::fromCoefficients(result);
}
static void checkScalar(){
    const f::Element relation{0x32ULL^(1ULL<<24)^(1ULL<<40),0};
    require(f::multiplyOracle(f::binaryBasis(120),f::binaryBasis(8))==relation,"defining relation");
    for(unsigned a=0;a!=128;++a)for(unsigned b=0;b!=128;++b)
        require(f::multiply(f::binaryBasis(a),f::binaryBasis(b))==
                f::multiplyOracle(f::binaryBasis(a),f::binaryBasis(b)),"scalar binary products");
    for(unsigned a=0;a!=256;++a)for(unsigned b=0;b!=256;++b)
        require(f::multiplyOracle(f::embedByte(a),f::embedByte(b))==f::embedByte(f::byteMultiply(a,b)),
                "all embedded byte products");
    for(unsigned bit=0;bit!=128;++bit)for(unsigned c=0;c!=8;++c){
        const auto value=f::binaryBasis(bit);const auto scalar=std::uint8_t(1U<<c);
        require(f::multiplyOracle(value,f::embedByte(scalar))==bytewise(value,scalar),"byte compatibility basis");
    }
    for(unsigned trial=0;trial!=128;++trial){
        const auto a=randomElement(),b=randomElement(),c=randomElement();
        require(f::multiply(a,b)==f::multiplyOracle(a,b),"scalar dense products");
        require(f::multiplyOracle(a,b^c)==(f::multiplyOracle(a,b)^f::multiplyOracle(a,c)),"field distributivity");
        require(f::multiplyOracle(f::multiplyOracle(a,b),c)==f::multiplyOracle(a,f::multiplyOracle(b,c)),
                "full-field multiplication linearity");
        require(bytewise(f::multiplyOracle(a,b),0x53)==f::multiplyOracle(bytewise(a,0x53),b),
                "fixed GF256 map is GF128-linear");
    }
    std::cout<<"sparse scalar: 16384 binary products, all65536 embedded byte products, "
                 "defining relation, fixed-byte compatibility, full-field linearity pass\n";
}

#if FIELD128_LIFT_SIMD_CHECK
static void checkSimd(){
    alignas(64) f::Element input[4],coefficients[4],output[4],expected[4],scalars[4],scaled[4],other[4];
    for(unsigned left=0;left!=128;++left)for(unsigned right=0;right!=128;right+=4){
        for(unsigned lane=0;lane!=4;++lane){
            input[lane]=f::binaryBasis((left+17*lane)%128);
            coefficients[lane]=f::binaryBasis(right+lane);
            expected[lane]=f::multiplyOracle(input[lane],coefficients[lane]);
        }
        f::multiply4Records(input,coefficients,output);
        require(std::equal(output,output+4,expected),"SIMD binary products");
    }
    // The literal GFNI instruction is compared with the complete multiplier.
    for(unsigned c=0;c!=256;++c){
        for(unsigned lane=0;lane!=4;++lane){input[lane]=randomElement();coefficients[lane]=f::embedByte(c);}
        f::multiply4Records(input,coefficients,output);
        _mm512_storeu_si512(expected,_mm512_gf2p8mul_epi8(_mm512_loadu_si512(input),_mm512_set1_epi8(char(c))));
        require(std::equal(output,output+4,expected),"literal byte GFNI compatibility");
    }
    for(unsigned trial=0;trial!=128;++trial){
        for(unsigned lane=0;lane!=4;++lane){
            input[lane]=randomElement();coefficients[lane]=randomElement();scalars[lane]=randomElement();
            if(trial<2)coefficients[lane]={trial,0};
            expected[lane]=f::multiplyOracle(input[lane],coefficients[lane]);
        }
        f::multiply4Records(input,coefficients,output);
        require(std::equal(output,output+4,expected),"SIMD dense products");
        f::multiply4Records(input,scalars,scaled);
        f::multiply4Records(scaled,coefficients,other);
        f::multiply4Records(output,scalars,scaled);
        require(std::equal(scaled,scaled+4,other),"SIMD full-field linearity");
        f::multiply4Records(input,coefficients,input);
        require(std::equal(input,input+4,expected),"SIMD exact input alias");
    }
    alignas(64) std::array<std::uint8_t,192> source,keys,result;
    constexpr std::uint8_t guard=0xa7;
    for(unsigned offset=0;offset!=64;offset+=8)for(unsigned keyOffset=0;keyOffset!=64;keyOffset+=8){
        source.fill(guard);keys.fill(guard);result.fill(guard);
        for(unsigned lane=0;lane!=4;++lane){
            const auto value=randomElement(),coefficient=randomElement();
            std::memcpy(source.data()+64+offset+16*lane,&value,sizeof(value));
            std::memcpy(keys.data()+64+keyOffset+16*lane,&coefficient,sizeof(coefficient));
            expected[lane]=f::multiplyOracle(value,coefficient);
        }
        f::multiply4Records(reinterpret_cast<const f::Element*>(source.data()+64+offset),
                           reinterpret_cast<const f::Element*>(keys.data()+64+keyOffset),
                           reinterpret_cast<f::Element*>(result.data()+64+offset));
        std::memcpy(output,result.data()+64+offset,sizeof(output));
        require(std::equal(output,output+4,expected),"SIMD alignments");
        for(unsigned i=0;i!=result.size();++i)if(i<64+offset||i>=128+offset)
            require(result[i]==guard,"SIMD output guards");
    }
    std::cout<<"sparse SIMD: all16384 binary products, literal GFNI compatibility, "
                 "dense products, GF128 linearity, zero/one, aliases, alignments/guards pass\n";
}
#endif

int main(){
    try{
#if FIELD128_LIFT_SIMD_CHECK && (defined(__GNUC__) || defined(__clang__))
        if(!__builtin_cpu_supports("avx512f")||!__builtin_cpu_supports("avx512bw")||
           !__builtin_cpu_supports("gfni")||!__builtin_cpu_supports("vpclmulqdq")){
            std::cout<<"SKIP: requires AVX512F/BW, GFNI, VPCLMULQDQ\n";return 77;
        }
#endif
        checkScalar();
#if FIELD128_LIFT_SIMD_CHECK
        checkSimd();
#endif
        return 0;
    }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
