#include "WideDiagonal16.h"
#include "WideDiagonal32.h"
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
#if defined(_MSC_VER)
#include <intrin.h>
#endif

namespace mf=spin::research::char2::medium32;
namespace w16=spin::research::char2::wide16;
namespace w32=spin::research::char2::wide32;
using Clock=std::chrono::steady_clock;
enum class Kernel { Gf8,Gf16,Tower32,Matrix32 };

static bool available() noexcept {
#if defined(_MSC_VER)
    int cpu[4]{};__cpuid(cpu,0);if(cpu[0]<7)return false;
    __cpuidex(cpu,1,0);if(!(cpu[2]&(1<<27))||(_xgetbv(0)&0xe6)!=0xe6)return false;
    __cpuidex(cpu,7,0);
    return (cpu[1]&(1<<16))&&(cpu[1]&(1<<30))&&(cpu[1]&(1U<<31))&&(cpu[2]&(1<<8));
#else
    __builtin_cpu_init();
    return __builtin_cpu_supports("avx512f")&&__builtin_cpu_supports("avx512bw")&&
        __builtin_cpu_supports("avx512vl")&&__builtin_cpu_supports("gfni");
#endif
}
static void require(bool condition,const char* message){if(!condition)throw std::runtime_error(message);}
static void compilerBarrier() noexcept {
#if defined(_MSC_VER)
    _ReadWriteBarrier();
#else
    asm volatile("" ::: "memory");
#endif
}

// Independent schoolbook quadratic oracle, using the existing polynomial
// long-division AES oracle rather than the kernel's shift/XOR multiplier.
static std::uint16_t oracle16(std::uint16_t x,std::uint16_t coefficient){
    const auto x0=std::uint8_t(x),x1=std::uint8_t(x>>8);
    const auto a=std::uint8_t(coefficient),b=std::uint8_t(coefficient>>8);
    const auto top=mf::byteOracle(x1,b);
    const auto low=mf::byteOracle(x0,a)^mf::byteOracle(mf::delta,top);
    const auto high=mf::byteOracle(x0,b)^mf::byteOracle(x1,a)^top;
    return std::uint16_t(low|(unsigned(high)<<8));
}
static std::uint16_t lane16(const mf::Record128& x,unsigned lane){
    return static_cast<std::uint16_t>(x.lane[lane/2]>>(16*(lane%2)));
}

template<Kernel K>
static void run(const mf::Record128* input,const std::uint8_t* c8,const std::uint32_t* c32,
    const w16::Coefficient* packed,const w32::Coefficient* matrix,mf::Record128* output,std::size_t count){
    if constexpr(K==Kernel::Gf16)w16::diagonal(input,packed,output,count);
    else if constexpr(K==Kernel::Matrix32)w32::diagonal(input,matrix,output,count);
    else if constexpr(K==Kernel::Gf8)mf::diagonal<mf::Kernel::Gf8>(input,c32,c8,output,count);
    else mf::diagonal<mf::Kernel::Tower32>(input,c32,c8,output,count);
}
template<Kernel K>
static mf::Record128 oracleRecord(const mf::Record128& input,std::uint8_t c8,
    std::uint16_t c16,std::uint32_t c32){
    mf::Record128 output{};
    for(unsigned j=0;j!=4;++j){
        const auto x=input.lane[j];
        if constexpr(K==Kernel::Gf16)output.lane[j]=oracle16(std::uint16_t(x),c16)|
            (std::uint32_t(oracle16(std::uint16_t(x>>16),c16))<<16);
        else if constexpr(K==Kernel::Tower32||K==Kernel::Matrix32)output.lane[j]=mf::towerMultiplyOracle(x,c32);
        else output.lane[j]=std::uint32_t(mf::byteOracle(std::uint8_t(x),c8))|
            (std::uint32_t(mf::byteOracle(std::uint8_t(x>>8),c8))<<8)|
            (std::uint32_t(mf::byteOracle(std::uint8_t(x>>16),c8))<<16)|
            (std::uint32_t(mf::byteOracle(std::uint8_t(x>>24),c8))<<24);
    }
    return output;
}

static void arithmeticChecks(){
    require(mf::byteTrace(mf::delta)==1,"GF216 quadratic trace criterion failed");
    require(oracle16(0x100,0x100)==(0x100|mf::delta),"GF216 u^2 relation failed");
    alignas(64) mf::Record128 input[4]{},output[4]{},byteOutput[4]{};
    w16::Coefficient packed[4]{};
    // Every scalar coefficient and every binary input basis vector. A
    // four-record SIMD batch tests four distinct scalars in parallel.
    for(unsigned first=0;first!=65536;first+=4){
        for(unsigned record=0;record!=4;++record)packed[record]=w16::prepareCoefficient(std::uint16_t(first+record));
        for(unsigned half=0;half!=2;++half){
            for(auto& record:input)for(unsigned j=0;j!=4;++j)
                record.lane[j]=(1U<<(8*half+2*j))|((1U<<(8*half+2*j+1))<<16);
            _mm512_store_si512(output,w16::multiply512(_mm512_load_si512(input),packed));
            for(unsigned record=0;record!=4;++record)for(unsigned j=0;j!=8;++j){
                const auto x=lane16(input[record],j),c=std::uint16_t(first+record);
                require(lane16(output[record],j)==oracle16(x,c),"GF216 SIMD coefficient/basis oracle mismatch");
                require(lane16(output[record],j)==mf::quadraticMultiply(x,c),"GF216 existing quadratic oracle mismatch");
            }
        }
    }
    std::mt19937_64 random(20261005);
    for(unsigned trial=0;trial!=1024;++trial){
        std::uint16_t coefficient[4]{};
        for(unsigned record=0;record!=4;++record){
            coefficient[record]=std::uint16_t(random());packed[record]=w16::prepareCoefficient(coefficient[record]);
            for(auto& word:input[record].lane)word=std::uint32_t(random());
        }
        _mm512_store_si512(output,w16::multiply512(_mm512_load_si512(input),packed));
        for(unsigned record=0;record!=4;++record)
            require(output[record]==oracleRecord<Kernel::Gf16>(input[record],0,coefficient[record],0),"GF216 SIMD random mismatch");
    }
    for(unsigned coefficient=0;coefficient!=256;++coefficient){
        std::uint8_t c8[4]{};
        for(unsigned record=0;record!=4;++record){
            packed[record]=w16::prepareCoefficient(std::uint16_t(coefficient));c8[record]=std::uint8_t(coefficient);
            for(auto& word:input[record].lane)word=std::uint32_t(random());
        }
        w16::diagonal(input,packed,output,4);
        mf::diagonal<mf::Kernel::Gf8>(input,nullptr,c8,byteOutput,4);
        require(std::equal(std::begin(output),std::end(output),std::begin(byteOutput)),"embedded GF256 action mismatch");
    }
}

static void matrix32Checks(){
    alignas(64) mf::Record128 input[4]{},output[4]{},tower[4]{},bytes[4]{};
    w32::Coefficient matrix[4]{};
    std::uint32_t scalar[4]{};
    const auto checkBatch=[&]{
        for(unsigned record=0;record!=4;++record){
            matrix[record]=w32::prepareCoefficient(scalar[record]);
            require(matrix[record].diagonal==w32::prepareCoefficientDirect(scalar[record]).diagonal,"GF232 setup table/oracle mismatch");
            require(w32::recoverScalar(matrix[record])==scalar[record],"GF232 matrix scalar recovery mismatch");
        }
        _mm512_store_si512(output,w32::multiply512(_mm512_load_si512(input),matrix));
        mf::diagonal<mf::Kernel::Tower32>(input,scalar,nullptr,tower,4);
        for(unsigned record=0;record!=4;++record){
            require(output[record]==oracleRecord<Kernel::Matrix32>(input[record],0,0,scalar[record]),"GF232 matrix independent oracle mismatch");
            require(output[record]==tower[record],"GF232 matrix/tower SIMD mismatch");
        }
    };
    for(unsigned coefficientBit=0;coefficientBit!=32;coefficientBit+=4)
        for(unsigned inputBit=0;inputBit!=32;inputBit+=4){
            for(unsigned record=0;record!=4;++record){
                scalar[record]=1U<<(coefficientBit+record);
                for(unsigned j=0;j!=4;++j)input[record].lane[j]=1U<<(inputBit+j);
            }
            checkBatch();
        }
    std::mt19937_64 random(320041);
    for(unsigned trial=0;trial!=1024;++trial){
        for(unsigned record=0;record!=4;++record){
            scalar[record]=std::uint32_t(random());
            for(auto& word:input[record].lane)word=std::uint32_t(random());
        }
        checkBatch();
    }
    for(unsigned coefficient=0;coefficient!=256;++coefficient){
        std::uint8_t c8[4]{};
        for(unsigned record=0;record!=4;++record){
            scalar[record]=coefficient;c8[record]=std::uint8_t(coefficient);
            for(auto& word:input[record].lane)word=std::uint32_t(random());
        }
        checkBatch();
        mf::diagonal<mf::Kernel::Gf8>(input,nullptr,c8,bytes,4);
        require(std::equal(std::begin(output),std::end(output),std::begin(bytes)),"GF232 matrix embedded GF256 mismatch");
    }
}

template<Kernel K>
static void diagonalChecks(){
    std::mt19937_64 random(518731);
    const mf::Record128 guard{{0x736baa91U,0xf152793dU,0xa0407b1eU,0x103abc47U}};
    for(std::size_t count:{0,1,2,3,4,5,7,8,9,15,16,17,31,63,64,65,257}){
        std::vector<mf::Record128> source(count),expected(count),input(count+8),output(count+8);
        std::vector<std::uint8_t> c8(count+4);
        std::vector<std::uint16_t> c16(count+4);
        std::vector<std::uint32_t> c32(count+4);
        std::vector<w16::Coefficient> packed(count+4);
        std::vector<w32::Coefficient> matrix(count+4);
        for(auto& record:source)for(auto& word:record.lane)word=std::uint32_t(random());
        for(std::size_t j=0;j!=count+4;++j){
            c8[j]=std::uint8_t(random());c16[j]=std::uint16_t(random());c32[j]=std::uint32_t(random());
        }
        w16::prepare(c16.data(),packed.data(),packed.size());
        w32::prepare(c32.data(),matrix.data(),matrix.size());
        const auto packedCopy=packed;
        const auto matrixCopy=matrix;
        for(unsigned coefficientOffset=0;coefficientOffset!=4;++coefficientOffset){
            for(std::size_t j=0;j!=count;++j)expected[j]=oracleRecord<K>(source[j],c8[j+coefficientOffset],c16[j+coefficientOffset],c32[j+coefficientOffset]);
            for(unsigned inputOffset=0;inputOffset!=4;++inputOffset){
                std::fill(input.begin(),input.end(),guard);
                std::copy(source.begin(),source.end(),input.begin()+inputOffset);
                const auto original=input;
                for(unsigned outputOffset=0;outputOffset!=4;++outputOffset){
                    std::fill(output.begin(),output.end(),guard);
                    run<K>(input.data()+inputOffset,c8.data()+coefficientOffset,c32.data()+coefficientOffset,
                        packed.data()+coefficientOffset,matrix.data()+coefficientOffset,output.data()+outputOffset,count);
                    require(std::equal(expected.begin(),expected.end(),output.begin()+outputOffset),"diagonal aligned/unaligned or tail mismatch");
                    require(input==original,"diagonal changed input");
                    for(std::size_t j=0;j!=output.size();++j)
                        if(j<outputOffset||j>=outputOffset+count)require(output[j]==guard,"diagonal output guard changed");
                }
                auto alias=input,aliasExpected=input;
                std::copy(expected.begin(),expected.end(),aliasExpected.begin()+inputOffset);
                run<K>(alias.data()+inputOffset,c8.data()+coefficientOffset,c32.data()+coefficientOffset,
                    packed.data()+coefficientOffset,matrix.data()+coefficientOffset,alias.data()+inputOffset,count);
                require(alias==aliasExpected,"diagonal in-place alias or guard mismatch");
            }
        }
        require(std::memcmp(packed.data(),packedCopy.data(),packed.size()*sizeof(w16::Coefficient))==0,"diagonal changed coefficient setup");
        require(std::memcmp(matrix.data(),matrixCopy.data(),matrix.size()*sizeof(w32::Coefficient))==0,"diagonal changed matrix setup");
        // Exact-size allocations let sanitizers catch full-batch/tail
        // overreads that a guard-only test could otherwise hide.
        std::vector<mf::Record128> exactOutput(count);
        std::vector<std::uint8_t> exact8(c8.begin(),c8.begin()+count);
        std::vector<std::uint32_t> exact32(c32.begin(),c32.begin()+count);
        std::vector<w16::Coefficient> exact16(packed.begin(),packed.begin()+count);
        std::vector<w32::Coefficient> exactMatrix(matrix.begin(),matrix.begin()+count);
        run<K>(source.data(),exact8.data(),exact32.data(),exact16.data(),exactMatrix.data(),exactOutput.data(),count);
        for(std::size_t j=0;j!=count;++j)
            require(exactOutput[j]==oracleRecord<K>(source[j],c8[j],c16[j],c32[j]),"exact-sized diagonal mismatch");
    }
}

static void check(){
    arithmeticChecks();matrix32Checks();diagonalChecks<Kernel::Gf8>();diagonalChecks<Kernel::Gf16>();
    diagonalChecks<Kernel::Tower32>();diagonalChecks<Kernel::Matrix32>();
    std::cout<<"GF216 delta=0x20 trace=1; all65536 coefficients x16 binary basis vectors; "
        "independent schoolbook and existing quadratic oracles; random SIMD; embedded GF256; "
        "GF232 matrix all1024 binary basis products, setup tables, scalar recovery, random/oracle/tower equality; "
        "GF8/GF216/tower32/matrix32 controls; exact-size tails, all record alignments modulo64, coefficient offsets, guards, aliases pass\n";
}

template<Kernel K>static constexpr const char* name(){
    if constexpr(K==Kernel::Gf8)return "gf8";
    else if constexpr(K==Kernel::Gf16)return "gf16-pairs";
    else if constexpr(K==Kernel::Matrix32)return "gf32-matrix4";
    else return "tower32";
}
template<Kernel K>
static std::uint64_t benchmark(std::size_t count,unsigned calls,unsigned passes,std::uint64_t seed){
    std::mt19937_64 random(seed);
    std::vector<mf::Record128> input(count),output(count);
    std::vector<std::uint8_t> c8(count);
    std::vector<std::uint32_t> c32(count);
    std::vector<w16::Coefficient> packed(count);
    std::vector<w32::Coefficient> matrix(K==Kernel::Matrix32?count:0);
    for(std::size_t i=0;i!=count;++i){
        for(auto& word:input[i].lane)word=std::uint32_t(random());
        std::uint16_t c16;
        do{c8[i]=std::uint8_t(random());}while(!c8[i]);
        do{c16=std::uint16_t(random());}while(!c16);
        do{c32[i]=std::uint32_t(random());}while(!c32[i]);
        packed[i]=w16::prepareCoefficient(c16);
    }
    if constexpr(K==Kernel::Matrix32)w32::prepare(c32.data(),matrix.data(),count);
    const auto apply=[&]{run<K>(input.data(),c8.data(),c32.data(),packed.data(),matrix.data(),output.data(),count);compilerBarrier();};
    for(unsigned warm=0;warm!=3;++warm)apply();
    std::vector<double> timings(calls);
    for(auto& timing:timings){
        const auto begin=Clock::now();
        for(unsigned pass=0;pass!=passes;++pass)apply();
        timing=std::chrono::duration<double,std::micro>(Clock::now()-begin).count()/passes;
    }
    std::sort(timings.begin(),timings.end());
    std::uint64_t checksum=0;
    for(const auto& record:output)for(const auto word:record.lane)checksum=(checksum^word)*0x100000001b3ULL;
    constexpr unsigned bits=K==Kernel::Gf8?8:K==Kernel::Gf16?16:32;
    constexpr unsigned coefficientBytes=K==Kernel::Gf8?1:K==Kernel::Matrix32?16:4;
    std::cout<<name<K>()<<','<<count<<",128,"<<128/bits<<','<<bits<<','<<coefficientBytes<<','
        <<calls<<','<<passes<<','<<seed<<','<<std::fixed<<std::setprecision(6)<<timings[calls/2]<<','
        <<timings[calls/10]<<','<<timings[9*calls/10]<<','<<1000*timings[calls/2]/count<<','
        <<coefficientBytes*count<<','<<32*count<<",1,"<<std::hex<<checksum<<std::dec<<'\n';
    return checksum;
}

int main(int argc,char** argv){try{
    if(!available()){std::cerr<<"wide_diagonal16 requires OS-enabled AVX512F/BW/VL and GFNI\n";return 77;}
    if(argc==2&&std::string(argv[1])=="--check"){check();return 0;}
    if(argc<2||argc>7||std::string(argv[1])!="--bench")
        throw std::invalid_argument("usage: wide_diagonal16 --check | --bench [all|records] [calls=101] [all|gf8|gf16|tower32|matrix32] [passes=auto] [seed=113]");
    const std::string size=argc>2?argv[2]:"all",kernel=argc>4?argv[4]:"all",passText=argc>5?argv[5]:"auto";
    const auto calls=argc>3?std::stoull(argv[3]):101,seed=argc>6?std::stoull(argv[6]):113;
    const auto explicitPasses=passText=="auto"?0:std::stoull(passText);
    if(!calls||calls>1000000||(passText!="auto"&&!explicitPasses)||explicitPasses>1000000||
        (kernel!="all"&&kernel!="gf8"&&kernel!="gf16"&&kernel!="tower32"&&kernel!="matrix32"))
        throw std::invalid_argument("invalid calls, passes, or kernel");
    const std::vector<std::size_t> sizes=size=="all"?std::vector<std::size_t>{4096,65536,262144,1048576}:
        std::vector<std::size_t>{std::stoull(size)};
    std::cout<<"kernel,records,payload_bits,field_lanes,coefficient_bits,stored_coefficient_bytes_per_record,calls,passes,seed,median_us,p10_us,p90_us,ns_per_record,coefficient_read_bytes,payload_io_bytes,setup_excluded,checksum\n";
    for(const auto count:sizes){
        if(!count||count>(1ULL<<28))throw std::invalid_argument("invalid record count");
        const auto passes=explicitPasses?unsigned(explicitPasses):count<=4096?31U:count<=65536?3U:1U;
        if(kernel=="all"||kernel=="gf8")benchmark<Kernel::Gf8>(count,unsigned(calls),passes,seed);
        if(kernel=="all"||kernel=="gf16")benchmark<Kernel::Gf16>(count,unsigned(calls),passes,seed);
        std::uint64_t towerChecksum=0,matrixChecksum=0;
        if(kernel=="all"||kernel=="tower32")towerChecksum=benchmark<Kernel::Tower32>(count,unsigned(calls),passes,seed);
        if(kernel=="all"||kernel=="matrix32")matrixChecksum=benchmark<Kernel::Matrix32>(count,unsigned(calls),passes,seed);
        if(kernel=="all")require(towerChecksum==matrixChecksum,"benchmark GF232 matrix/tower checksum mismatch");
    }
    return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
