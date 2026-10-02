// Research-only packed-domain mixers. Modes 6/7 relocate the exact old GF16 maps;
// modes 1--5 change the code, not just its implementation.
#define main packed_mixer_unused_compare_main
#include "bch_compare.cpp"
#undef main

namespace spin::detail::kernel {
void bchPackedDenseShared(const block*,block*,const std::uint64_t*);
void bchPackedFullShared(const block*,block*,const std::uint64_t*);
void bchPackedFullIndependent(const block*,block*,const std::uint64_t*);
void bchPackedFullGl32(const block*,block*,const std::uint64_t*);
void bchPackedFullLegacyAnd(const block*,block*,const std::uint64_t*);
}
namespace pm {
static unsigned rank(std::array<std::uint32_t,32> a,unsigned n) {
    unsigned r=0;
    for(unsigned c=0;c<n;++c) {
        unsigned p=r;while(p<n && !(a[p]&(1U<<c)))++p;
        if(p==n)continue;
        std::swap(a[r],a[p]);
        for(unsigned j=r+1;j<n;++j)if(a[j]&(1U<<c))a[j]^=a[r];
        ++r;
    }
    return r;
}
static std::array<std::uint32_t,32> draw(unsigned n,k::setup::Words& rng) {
    std::array<std::uint32_t,32> a{};
    do {for(unsigned j=0;j<n;++j)a[j]=std::uint32_t(rng())&(n==32?~0U:255U);}
    while(rank(a,n)!=n);
    return a;
}
template<unsigned Mode> struct Setup {
    static constexpr unsigned stride=Mode>=4?1024:(Mode==3?256:32);
    std::vector<std::uint64_t> coeff;
    Setup(std::size_t tiles,std::uint64_t seed):coeff(tiles*stride) {
        if constexpr(Mode>=6) {
            // Shared routing preserves the four lane identities. Reindex the
            // old packet transpose controls into canonical BCH columns.
            const Packets<4> packets(unsigned(4*tiles),seed,true,true);
            for(std::size_t p=0;p<packets.bases.size();++p) {
                const auto base=packets.bases[p],tile=base/tileStride,column=(base%tileStride)/4;
                const auto group=column/8,j=column%8;
                const auto control=std::uint32_t(packets.controls[p]);
                for(unsigned r=0;r<4;++r)for(unsigned d=0;d<4;++d) {
                    const unsigned q=(r+d)%4;
                    if((control>>(8*q+2*r))&1U) {
                        const auto bit=Mode==7?(0x0101010101010101ULL<<j):(std::uint64_t(1)<<(8*(7-j)+j));
                        coeff[tile*stride+32*group+8*d+2*r]|=bit;
                        coeff[tile*stride+32*group+8*d+2*r+1]|=bit;
                    }
                }
            }
            return;
        }
        k::setup::Words rng(seed^0x75a1dc09ULL);
        for(std::size_t t=0;t<tiles;++t)for(unsigned g=0;g<32;++g) {
            if constexpr(Mode>=4) {
                const auto a=draw(32,rng);
                for(unsigned d=0;d<4;++d)for(unsigned lane=0;lane<4;++lane) {
                    std::uint64_t m=0;
                    for(unsigned j=0;j<8;++j)m|=std::uint64_t((a[8*lane+j]>>(8*((lane+d)%4)))&255)<<(8*(7-j));
                    coeff[t*stride+g*32+d*8+2*lane]=m;
                    coeff[t*stride+g*32+d*8+2*lane+1]=m;
                }
            } else {
                for(unsigned lane=0;lane<(Mode==3?4U:1U);++lane) {
                    const auto a=draw(8,rng);std::uint64_t m=0;
                    for(unsigned j=0;j<8;++j)m|=std::uint64_t(a[j])<<(8*(7-j));
                    if constexpr(Mode==3) {
                        coeff[t*stride+g*8+2*lane]=m;coeff[t*stride+g*8+2*lane+1]=m;
                    } else coeff[t*stride+g]=m;
                }
            }
        }
    }
};
template<unsigned Mode> static std::uint32_t matrixRow(const std::uint64_t* coeff,unsigned group,unsigned row) {
    const unsigned lane=row/8,j=row%8;
    if constexpr(Mode>=4) {
        std::uint32_t result=0;
        for(unsigned d=0;d<4;++d) {
            const auto m=coeff[32*group+8*d+2*lane];
            const auto byte=Mode==7?((m&(1ULL<<j))?1U<<j:0U):std::uint32_t((m>>(8*(7-j)))&255);
            result|=byte<<(8*((lane+d)%4));
        }
        return result;
    } else {
        const auto m=coeff[Mode==3?8*group+2*lane:group];
        return std::uint32_t((m>>(8*(7-j)))&255)<<(8*lane);
    }
}
template<unsigned Mode> static void scalarMix(const k::block* in,k::block* out,const std::uint64_t* coeff) {
    for(unsigned group=0;group<32;++group)for(unsigned row=0;row<32;++row) {
        const unsigned at=4*(8*group+row%8)+row/8;
        if constexpr(Mode==1)if(group<16){out[at]=in[at];continue;}
        k::block sum{};auto mask=matrixRow<Mode>(coeff,group,row);
        while(mask){const auto c=std::countr_zero(mask);sum^=in[4*(8*group+c%8)+c/8];mask&=mask-1;}
        out[at]=sum;
    }
}
template<unsigned Mode> static SPIN_FORCEINLINE void kernel(const k::block* in,k::block* out,const std::uint64_t* coeff) {
    if constexpr(Mode==0)k::bchTranspose4GfniBlend(in,out);
    else if constexpr(Mode==1)k::bchPackedDenseShared(in,out,coeff);
    else if constexpr(Mode==2)k::bchPackedFullShared(in,out,coeff);
    else if constexpr(Mode==3)k::bchPackedFullIndependent(in,out,coeff);
    else if constexpr(Mode==7)k::bchPackedFullLegacyAnd(in,out,coeff);
    else k::bchPackedFullGl32(in,out,coeff);
}
template<unsigned Mode> static void reference(const k::block* in,k::block* out,const std::uint64_t* coeff) {
    if constexpr(Mode==0)k::bchTranspose4(in,out);
    else {alignas(64)k::block tmp[1024];scalarMix<Mode>(in,tmp,coeff);k::bchTranspose4(tmp,out);}
}
template<unsigned Mode,bool Reference=false> static SPIN_NOINLINE void bulk(const k::block* in,k::block* out,std::size_t tiles,const Setup<Mode>& setup) {
    for(std::size_t t=0;t<tiles;++t) {
        if constexpr(Reference)reference<Mode>(in+t*tileStride,out+t*512,setup.coeff.data()+t*Setup<Mode>::stride);
        else kernel<Mode>(in+t*tileStride,out+t*512,setup.coeff.data()+t*Setup<Mode>::stride);
    }
}
template<unsigned Mode,bool Reference=false> static SPIN_NOINLINE void fullEncode(k::block* input,k::block* scratch,std::size_t n,
    const unsigned* masks,const Packets<4>& packets,const Setup<Mode>& setup,double* phases) {
    const auto start=phases?Clock::now():Clock::time_point{};
    __m128i v1{},v2{},v3{};
    reverseInner<2>(input,n,masks,[&](std::size_t i,k::block v) {
        switch(i&3) {
            case 3:v3=v.mData;break;case 2:v2=v.mData;break;case 1:v1=v.mData;break;
            case 0:{
                auto x=_mm512_castsi128_si512(v.mData);
                x=_mm512_inserti32x4(x,v1,1);x=_mm512_inserti32x4(x,v2,2);x=_mm512_inserti32x4(x,v3,3);
                if constexpr(Mode!=4 && Mode<6)x=spin::research::gf16::transpose(x,std::uint32_t(packets.controls[i/4]));
                _mm512_stream_si512(reinterpret_cast<__m512i*>(scratch+packets.bases[i/4]),x);break;
            }
        }
    });
    _mm_sfence();const auto middle=phases?Clock::now():Clock::time_point{};
    bulk<Mode,Reference>(scratch,input,n/1024,setup);
    if(phases){phases[0]+=std::chrono::duration<double,std::milli>(middle-start).count();phases[1]+=std::chrono::duration<double,std::milli>(Clock::now()-middle).count();}
}
template<unsigned Mode> static void check(std::uint64_t seed) {
    Setup<Mode> setup(1,seed);bc::Buffer input(1024),actual(512),expected(512),images(1024*512);
    for(unsigned pattern=0;pattern<3;++pattern) {
        bc::fillBlocks(input.data,1024,seed+pattern);reference<Mode>(input.data,expected.data,setup.coeff.data());
        kernel<Mode>(input.data,actual.data,setup.coeff.data());
        if(std::memcmp(actual.data,expected.data,512*16))throw std::runtime_error("dense scalar/BCH mismatch");
    }
    std::memset(input.data,0,1024*16);
    for(unsigned c=0;c<1024;++c) {
        input.data[c]=k::block(~0ULL,~0ULL);reference<Mode>(input.data,images.data+c*512,setup.coeff.data());input.data[c]=k::block{};
    }
    for(unsigned c=0;c<1024;++c)for(unsigned bit=0;bit<128;++bit) {
        const auto unit=bit<64?k::block(0,1ULL<<bit):k::block(1ULL<<(bit-64),0);input.data[c]=unit;
        kernel<Mode>(input.data,actual.data,setup.coeff.data());
        for(unsigned j=0;j<512;++j)expected.data[j]=k::block(_mm_and_si128(images.data[c*512+j].mData,unit.mData));
        if(std::memcmp(actual.data,expected.data,512*16))throw std::runtime_error("complete physical basis mismatch");
        input.data[c]=k::block{};
    }
    k::setup::Words rng(seed);
    if constexpr(Mode>0)for(unsigned g=0;g<32;++g)for(unsigned trial=0;trial<128;++trial) {
        const auto x=std::uint32_t(rng()),y=std::uint32_t(rng());std::uint32_t mx=0,mty=0;
        for(unsigned row=0;row<32;++row) {
            const auto a=matrixRow<Mode>(setup.coeff.data(),g,row);
            mx|=(std::popcount(a&x)&1U)<<row;if((y>>row)&1)mty^=a;
        }
        if((std::popcount(mx&y)&1)!=(std::popcount(x&mty)&1))throw std::runtime_error("adjoint mismatch");
    }
    std::cout<<"checks: mode="<<Mode<<",seed="<<seed<<",131072 physical basis vectors, dense scalar reference, local matrix adjoint PASS\n";
}
template<unsigned Mode> static void run(const std::string& kind,unsigned exponent,std::uint64_t seed,unsigned calls) {
    const bool hot=kind=="hot",full=kind=="full"||kind=="profile";
    const std::size_t n=std::size_t(2)<<exponent,tiles=hot?1:n/1024;
    Setup<Mode> setup(tiles,seed);bc::Buffer input(full?n:tiles*tileStride),out(full?n:tiles*512),scratch(tiles*tileStride);
    if(!full) {
        bc::fillBlocks(input.data,tiles*tileStride,seed);bulk<Mode>(input.data,out.data,tiles,setup);
        bc::Buffer expected(tiles*512);bulk<Mode,true>(input.data,expected.data,tiles,setup);
        if(std::memcmp(out.data,expected.data,tiles*512*16))throw std::runtime_error("bulk reference mismatch");
        if(!calls)return;
        const unsigned repeat=hot?2048:1;
        bc::measure(kind,exponent,Mode,seed,calls,tiles*repeat,[&]{for(unsigned i=0;i<repeat;++i)bulk<Mode>(input.data,out.data,tiles,setup);},out.data,tiles*512);
        return;
    }
    const Packets<4> packets(unsigned(n/256),seed,true,true);std::vector<unsigned> masks(4*(n/128));
    k::setup::Words rng(seed^0x3f625a92ULL);
    for(std::size_t i=0;i<masks.size();i+=2){unsigned u;do{u=unsigned(rng())&((1U<<19)-1);}while(!u);unsigned v=unsigned(rng())&((1U<<19)-1);if(std::popcount(u&v)&1)v^=u&-u;masks[i]=u;masks[i+1]=v;}
    for(unsigned pattern=0;pattern<(calls?1U:3U);++pattern) {
        bc::fillBlocks(input.data,n,913+84*pattern);if(pattern==1){std::memset(input.data,0,n*16);input.data[n-3]=k::block(1,3);}
        std::memcpy(out.data,input.data,n*16);
        fullEncode<Mode,true>(out.data,scratch.data,n,masks.data(),packets,setup,nullptr);
        fullEncode<Mode>(input.data,scratch.data,n,masks.data(),packets,setup,nullptr);
        if(std::memcmp(input.data,out.data,n*16))throw std::runtime_error("full encoder/suffix reference mismatch");
        if constexpr(Mode>=6) {
            bc::fillBlocks(out.data,n,913+84*pattern);
            if(pattern==1){std::memset(out.data,0,n*16);out.data[n-3]=k::block(1,3);}
            ::encode<4,2,true,true,true>(out.data,scratch.data,n,masks.data(),packets,nullptr);
            if(std::memcmp(input.data,out.data,n*16))throw std::runtime_error("relocated GF16 differs from legacy encoder");
        }
    }
    std::cout<<"checks: full scalar-mixer/BCH reference and suffix PASS; coefficient_bytes="<<setup.coeff.size()*8<<"\n";
    if constexpr(Mode>=6)std::cout<<"checks: bit-for-bit original legacy shared-GF16 R2 encoder equality PASS\n";
    if(!calls)return;
    bc::fillBlocks(input.data,n,913);double phases[2]{};
    bc::measure(kind,exponent,Mode,seed,calls,tiles,[&]{fullEncode<Mode>(input.data,scratch.data,n,masks.data(),packets,setup,kind=="profile"?phases:nullptr);},input.data,n);
    if(kind=="profile")std::cout<<"phase_means_including_warmups,inner_route,"<<phases[0]/(calls+3)<<",packed_mixer_bch,"<<phases[1]/(calls+3)<<'\n';
}
template<unsigned Mode> static void dispatch(const std::string& kind,unsigned exponent,std::uint64_t seed,unsigned calls) {
    if(kind=="check")check<Mode>(seed);else run<Mode>(kind,exponent,seed,calls);
}
}
int main(int argc,char** argv){try {
    if(!__builtin_cpu_supports("avx512f")||!__builtin_cpu_supports("avx512vl")||!__builtin_cpu_supports("avx512bw")||!__builtin_cpu_supports("gfni"))throw std::runtime_error("AVX512/GFNI required");
    if(argc!=6)throw std::invalid_argument("usage: packed_mixer_perf check|hot|bulk|full|profile exponent mode seed calls");
    const std::string kind=argv[1];const unsigned exponent=std::stoul(argv[2]),mode=std::stoul(argv[3]),calls=std::stoul(argv[5]);const auto seed=std::stoull(argv[4]);
    if(exponent<14||exponent>20||mode>7||(kind!="check"&&kind!="hot"&&kind!="bulk"&&kind!="full"&&kind!="profile"))throw std::invalid_argument("unsupported case");
    switch(mode){case 0:pm::dispatch<0>(kind,exponent,seed,calls);break;case 1:pm::dispatch<1>(kind,exponent,seed,calls);break;case 2:pm::dispatch<2>(kind,exponent,seed,calls);break;case 3:pm::dispatch<3>(kind,exponent,seed,calls);break;case 4:pm::dispatch<4>(kind,exponent,seed,calls);break;case 5:pm::dispatch<5>(kind,exponent,seed,calls);break;case 6:pm::dispatch<6>(kind,exponent,seed,calls);break;case 7:pm::dispatch<7>(kind,exponent,seed,calls);break;}
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
