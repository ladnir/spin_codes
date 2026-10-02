// Bounded screen only: polynomial construction, not a tuned algebraic encoder.
#define main bch_compare_reference_main
#include "bch_compare.cpp"
#undef main

static void algebraicBasis() {
    bc::Buffer input(1024),rawMask(512),expected(512),actual(512);
    std::memset(input.data,0,1024*16);
    std::size_t checked=0;
    for(unsigned coordinate=0;coordinate<1024;++coordinate) {
        // The dense raw oracle consists only of XORs of opaque elements.
        // Its all-ones image therefore gives this coordinate's binary column.
        input.data[coordinate]=k::block(~0ULL,~0ULL);
        k::bchTranspose4AlgebraicRawReference(input.data,rawMask.data);
        for(unsigned bit=0;bit<128;++bit) {
            const auto unit=bit<64?k::block(0,1ULL<<bit):k::block(1ULL<<(bit-64),0);
            input.data[coordinate]=unit;
            bc::bch<0>(input.data,expected.data);bc::bch<4>(input.data,actual.data);
            if(std::memcmp(actual.data,expected.data,512*16))throw std::runtime_error("exact algebraic basis mismatch");
            for(unsigned j=0;j<512;++j)expected.data[j]=k::block(_mm_and_si128(rawMask.data[j].mData,unit.mData));
            bc::bch<5>(input.data,actual.data);
            if(std::memcmp(actual.data,expected.data,512*16))throw std::runtime_error("raw algebraic basis mismatch");
            ++checked;
        }
        input.data[coordinate]=k::block{};
    }
    std::cout<<"Passed "<<checked<<" complete physical basis vectors: exact versus production; raw versus polynomial-basis dense oracle.\n";
    for(unsigned seed:{1U,17U,997U}) {
        bc::fillBlocks(input.data,1024,seed);
        k::bchTranspose4AlgebraicRawReference(input.data,expected.data);bc::bch<5>(input.data,actual.data);
        if(std::memcmp(actual.data,expected.data,512*16))throw std::runtime_error("raw algebraic dense mismatch");
    }
}
static void rawIsolated(bool hot,std::uint64_t seed,unsigned calls) {
    const std::size_t tiles=hot?1:2048,repeats=hot?2048:1;
    bc::Buffer input(tiles*tileStride),out(tiles*512),reference(512);
    bc::fillBlocks(input.data,tiles*tileStride,seed);
    bc::bulk<5>(input.data,out.data,tiles);
    for(std::size_t tile=0;tile<std::min<std::size_t>(tiles,2);++tile) {
        k::bchTranspose4AlgebraicRawReference(input.data+tile*tileStride,reference.data);
        if(std::memcmp(out.data+tile*512,reference.data,512*16))throw std::runtime_error("raw isolated dense mismatch");
    }
    if(!calls)return;
    bc::measure(hot?"raw-hot":"raw-bulk",hot?1:20,5,seed,calls,tiles*repeats,[&]{
        for(std::size_t r=0;r<repeats;++r)bc::bulk<5>(input.data,out.data,tiles);
    },out.data,tiles*512);
}
int main(int argc,char** argv) {try {
    if(!__builtin_cpu_supports("avx512f") || !__builtin_cpu_supports("avx512vl") ||
       !__builtin_cpu_supports("avx512bw") || !__builtin_cpu_supports("gfni"))throw std::runtime_error("AVX512/GFNI required");
    if(argc==2 && std::string(argv[1])=="check") {algebraicBasis();return 0;}
    if(argc!=5)throw std::invalid_argument("usage: bch_algebraic_perf hot|bulk|full exact|raw seed calls; or check");
    const std::string kind=argv[1],basis=argv[2];
    const auto seed=std::stoull(argv[3]);const unsigned calls=std::stoul(argv[4]);
    if(basis=="exact") {
        if(kind=="hot")bc::isolated<4>(true,1,seed,calls);
        else if(kind=="bulk")bc::isolated<4>(false,20,seed,calls);
        else if(kind=="full")bc::full<4>(20,seed,calls,false);
        else throw std::invalid_argument("bad kind");
    } else if(basis=="raw" && (kind=="hot" || kind=="bulk"))rawIsolated(kind=="hot",seed,calls);
    else throw std::invalid_argument("raw basis is isolated-only");
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
