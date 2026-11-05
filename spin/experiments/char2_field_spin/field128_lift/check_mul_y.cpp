// Standalone correctness only; no benchmark or production integration.
// Portable: c++ -std=c++20 -O2 check_mul_y.cpp -o field128_mul_y_check
// SIMD: add -DFIELD128_LIFT_MULY_SIMD_CHECK=1 -mavx512f -mavx512bw -mgfni
// No CLMUL flag is needed. Run the SIMD build only on a supporting host.
#include "SparseByteField128MulY.h"
#ifndef FIELD128_LIFT_MULY_SIMD_CHECK
#define FIELD128_LIFT_MULY_SIMD_CHECK 0
#endif
#if FIELD128_LIFT_MULY_SIMD_CHECK
#include "SparseByteField128MulYSimd.h"
#endif
#include <algorithm>
#include <cstring>
#include <iostream>
#include <stdexcept>

namespace f=spin::research::char2::field128_lift::sparse_byte;
static void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
static std::uint64_t randomWord(){
    static std::uint64_t state=0xf12848a15c0ffeeULL;
    state^=state<<13;state^=state>>7;state^=state<<17;return state;
}
static f::Element randomElement(){return {randomWord(),randomWord()};}
static f::Element reference(f::Element value){return f::multiplyOracle(value,f::towerY);}
static void checkScalar(){
    for(unsigned byte=0;byte!=256;++byte)
        require(f::multiplyByteBy32(static_cast<std::uint8_t>(byte))==f::byteMultiply(byte,0x32),"fixed AES scalar");
    require(f::multiplyY({})==f::Element{},"zero");
    require(f::multiplyY({1,0})==f::towerY,"one");
    require(f::multiplyY(f::binaryBasis(120))==f::Element{0x10001000032ULL,0},"defining relation");
    for(unsigned bit=0;bit!=128;++bit){
        const auto value=f::binaryBasis(bit);
        require(f::multiplyY(value)==reference(value),"scalar basis");
    }
    for(unsigned trial=0;trial!=1024;++trial){
        const auto a=randomElement(),b=randomElement();
        require(f::multiplyY(a)==reference(a),"scalar random oracle");
        require(f::multiplyY(a)==f::multiply(a,f::towerY),"scalar independent Horner");
        require(f::multiplyY(a^b)==(f::multiplyY(a)^f::multiplyY(b)),"scalar additivity");
        require(f::multiplyY(f::multiplyOracle(a,b))==f::multiplyOracle(f::multiplyY(a),b),"full-field linearity");
    }
    std::cout<<"mulY scalar: all256 byte constants,128 basis,1024 random/oracle/Horner/linearity cases pass\n";
}

#if FIELD128_LIFT_MULY_SIMD_CHECK
static void checkSimd(){
    alignas(64) f::Element input[4],output[4],expected[4];
    for(unsigned bit=0;bit!=128;++bit){
        for(unsigned lane=0;lane!=4;++lane){input[lane]=f::binaryBasis((bit+31*lane)%128);expected[lane]=reference(input[lane]);}
        f::multiplyY4Records(input,output);
        require(std::equal(output,output+4,expected),"SIMD basis/lane isolation");
    }
    for(unsigned top=0;top!=256;++top){
        for(unsigned lane=0;lane!=4;++lane){
            input[lane]=randomElement();
            input[lane].hi=(input[lane].hi&0x00ffffffffffffffULL)|(std::uint64_t((top+lane)%256)<<56);
            expected[lane]=reference(input[lane]);
        }
        f::multiplyY4Records(input,output);
        require(std::equal(output,output+4,expected),"SIMD all high-byte reductions");
        f::multiplyY4Records(input,input);
        require(std::equal(input,input+4,expected),"SIMD exact alias");
    }
    alignas(64) std::array<std::uint8_t,256> source,result;
    constexpr std::uint8_t guard=0xa7;
    for(unsigned offset=0;offset!=64;++offset){
        source.fill(guard);result.fill(guard);
        const unsigned inStart=64+offset,outStart=127-offset;
        for(unsigned lane=0;lane!=4;++lane){input[lane]=randomElement();expected[lane]=reference(input[lane]);}
        std::memcpy(source.data()+inStart,input,64);
        f::multiplyY4Records(source.data()+inStart,result.data()+outStart);
        std::memcpy(output,result.data()+outStart,64);
        require(std::equal(output,output+4,expected),"SIMD all byte alignments");
        for(unsigned i=0;i!=result.size();++i)if(i<outStart||i>=outStart+64)
            require(result[i]==guard,"SIMD unaligned output guards");
    }
    for(int delta=-63;delta<=63;++delta){
        source.fill(guard);
        for(unsigned lane=0;lane!=4;++lane){input[lane]=randomElement();expected[lane]=reference(input[lane]);}
        std::memcpy(source.data()+96,input,64);
        const auto before=source;
        const unsigned outStart=static_cast<unsigned>(96+delta);
        f::multiplyY4Records(source.data()+96,source.data()+outStart);
        std::memcpy(output,source.data()+outStart,64);
        require(std::equal(output,output+4,expected),"SIMD arbitrary partial overlap");
        for(unsigned i=0;i!=source.size();++i)if(i<outStart||i>=outStart+64)
            require(source[i]==before[i],"SIMD overlapping output guards");
    }
    std::cout<<"mulY SIMD:128 basis batches,256 high-byte/random batches,exact alias,all byte alignments,127 partial overlaps pass\n";
}
#endif

int main(){
    try{
#if FIELD128_LIFT_MULY_SIMD_CHECK && (defined(__GNUC__) || defined(__clang__))
        if(!__builtin_cpu_supports("avx512f")||!__builtin_cpu_supports("avx512bw")||!__builtin_cpu_supports("gfni")){
            std::cout<<"SKIP: requires AVX512F/BW and GFNI\n";return 77;
        }
#endif
        checkScalar();
#if FIELD128_LIFT_MULY_SIMD_CHECK
        checkSimd();
#endif
        return 0;
    }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
