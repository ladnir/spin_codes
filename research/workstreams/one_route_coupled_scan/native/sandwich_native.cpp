// Research gate for GROUPED_45_CERTIFICATE.md. No production API changes.
// Reuse the frozen transpose-FFT implementation and its checks, not its main.
#define main retained_reset_main
#include "reset_native.cpp"
#undef main
#include "RsOuterOptimized.h"
#include "SparseByteField128Simd.h"
#include <memory>
#include <numeric>

namespace sandwich {
using E=spin::research::one_route::Element; // Storage only; NOT sparse-byte field multiplication.
namespace m=reset::m;
namespace rs=spin::research::one_route::rs_outer;
using Clock=std::chrono::steady_clock;
// J=GF(2^32) from MediumField32. Trace_J/F2(mu)=1 by the two quadratic
// trace identities. F=J[w,z], w^2+w+mu=0, z^2+z+mu*w=0.
// Storage order [1,w,z,wz], each coefficient in J's AES-byte tower basis.
constexpr std::uint32_t mu=0x20000000U;
static const reset::Packed muKey=reset::broadcast(mu);
static const reset::Packed mu2Key=reset::broadcast(m::towerMultiplyScalar(mu,mu));

template<unsigned A,unsigned B,unsigned C,unsigned D> FI __m512i select(__m512i x){
    constexpr unsigned imm=(A%4)|((B%4)<<2)|((C%4)<<4)|((D%4)<<6);
    constexpr unsigned mask=(A<4?1:0)|(B<4?2:0)|(C<4?4:0)|(D<4?8:0);
    auto y=_mm512_shuffle_epi32(x,static_cast<_MM_PERM_ENUM>(imm));
    if constexpr(mask!=15)y=_mm512_maskz_mov_epi32(mask*0x1111U,y);
    return y;
}
FI __m512i cross(__m512i x){
    const auto h=_mm512_xor_si512(x,select<2,3,0,1>(x));
    return _mm512_mask_blend_epi64(0xaa,h,_mm512_xor_si512(h,select<1,0,3,2>(h)));
}
template<std::uint32_t C,unsigned Rotation> FI __m512i constantTerm(__m512i x){
    constexpr auto diagonal=[] {
        std::uint32_t d=0;
        for(unsigned row=0;row<4;++row)
            d|=m::towerMultiplyScalar(C,1U<<(8*((row+Rotation)%4)))&(255U<<(8*row));
        return d;
    }();
    if constexpr(diagonal==0)return _mm512_setzero_si512();
    else {
        if constexpr(Rotation)x=_mm512_ror_epi32(x,8*Rotation);
        return _mm512_gf2p8mul_epi8(x,_mm512_set1_epi32(diagonal));
    }
}
template<std::uint32_t C> FI __m512i constantProduct(__m512i x){
    return _mm512_xor_si512(m::xor3(constantTerm<C,0>(x),constantTerm<C,1>(x),constantTerm<C,2>(x)),constantTerm<C,3>(x));
}
template<bool Specialized=true> FI __m512i mul128(__m512i a,__m512i b){
    const auto p=m::towerMultiply512(a,b);
    const auto q=m::towerMultiply512(_mm512_xor_si512(a,select<1,0,3,2>(a)),
                                    _mm512_xor_si512(b,select<1,0,3,2>(b)));
    const auto r=m::towerMultiply512(cross(a),cross(b));
    const auto base=_mm512_xor_si512(select<0,0,0,0>(p),
        m::xor3(select<4,0,4,0>(q),select<4,4,0,0>(r),select<4,4,4,2>(r)));
    const auto linear=m::xor3(select<1,4,1,4>(p),select<4,2,4,4>(q),select<4,4,1,4>(r));
    const auto quadratic=_mm512_xor_si512(select<2,3,4,4>(p),select<2,4,4,4>(q));
    if constexpr(Specialized)
        return m::xor3(base,constantProduct<mu>(linear),constantProduct<m::towerMultiplyScalar(mu,mu)>(quadratic));
    else return m::xor3(base,reset::product(linear,muKey),reset::product(quadratic,mu2Key));
}
using Pair=std::array<std::uint32_t,2>;
Pair qmul(Pair a,Pair b){
    const auto p=m::towerMultiplyOracle(a[0],b[0]);
    const auto q=m::towerMultiplyOracle(a[1],b[1]);
    const auto r=m::towerMultiplyOracle(a[0]^a[1],b[0]^b[1]);
    return {p^m::towerMultiplyOracle(mu,q),r^p};
}
E oracle(E aa,E bb){
    const auto a=std::bit_cast<std::array<std::uint32_t,4>>(aa);
    const auto b=std::bit_cast<std::array<std::uint32_t,4>>(bb);
    auto p=qmul({a[0],a[1]},{b[0],b[1]}),q=qmul({a[2],a[3]},{b[2],b[3]});
    auto r=qmul({a[0]^a[2],a[1]^a[3]},{b[0]^b[2],b[1]^b[3]});
    auto reduce=qmul({0,mu},q);
    return std::bit_cast<E>(std::array<std::uint32_t,4>{p[0]^reduce[0],p[1]^reduce[1],r[0]^p[0],r[1]^p[1]});
}
template<class T> struct Aligned {
    struct Free{void operator()(T* p)const{std::free(p);}};
    std::unique_ptr<T,Free> p;
    explicit Aligned(std::size_t n):p(static_cast<T*>(std::aligned_alloc(64,(n*sizeof(T)+63)&~std::size_t(63)))){
        if(!p)throw std::bad_alloc();std::memset(p.get(),0,n*sizeof(T));
    }
    T* data(){return p.get();}const T* data()const{return p.get();}
    T& operator[](std::size_t i){return p.get()[i];}
};
template<unsigned H> struct Map {
    using K=reset::Kernel<H>;
    static constexpr unsigned B=1U<<H;
    typename K::Fft zero=K::prepareFft(0),shift=K::prepareFft(reset::point(B));
    template<int I> static FI void inverseUpper(reset::E* x,const typename K::Fft& plan){
        if constexpr(I>=2){
            constexpr unsigned half=1U<<I;
            for(unsigned off=0;off<B;off+=2*half){
                const auto& key=plan[(B>>(I+1))-1+(off>>(I+1))];
                for(unsigned j=0;j<half;j+=4){
                    auto s=_mm512_loadu_si512(x+off+j),t=_mm512_loadu_si512(x+off+half+j);
                    auto b=_mm512_xor_si512(t,reset::product(s,key));
                    _mm512_storeu_si512(x+off+j,_mm512_xor_si512(s,b));
                    _mm512_storeu_si512(x+off+half+j,b);
                }
            }
            inverseUpper<I-1>(x,plan);
        }
    }
    template<unsigned I> static FI __m512i inverseSmall(__m512i x,const reset::Packed& p){
        auto s=_mm512_shuffle_i32x4(x,x,I==0?0xa0:0x44);
        auto t=_mm512_shuffle_i32x4(x,x,I==0?0xf5:0xee);
        auto b=_mm512_xor_si512(t,reset::product(s,p));
        return _mm512_mask_blend_epi64(I==0?0xcc:0xf0,_mm512_xor_si512(s,b),b);
    }
    static constexpr auto basis=[] {
        std::array<std::uint16_t,16> piv{},pre{},beta{};beta[0]=1;
        for(unsigned i=1;i<16;++i){
            auto x=std::uint16_t(1U<<i),y=std::uint16_t(m::quadraticMultiply(x,x)^x);
            while(y){auto j=std::bit_width(y)-1;if(!piv[j]){piv[j]=y;pre[j]=x;break;}y^=piv[j];x^=pre[j];}
        }
        for(unsigned i=1;i<16;++i){auto y=beta[i-1];std::uint16_t x=0;
            while(y){auto j=std::bit_width(y)-1;if(!piv[j])throw "Cantor root";y^=piv[j];x^=pre[j];}beta[i]=x;
        }
        return beta;
    }();
    static constexpr std::uint16_t pt(unsigned index){std::uint16_t c=0;for(unsigned j=0;index;index>>=1,++j)if(index&1)c^=basis[j];return c;}
    template<unsigned A,unsigned C,unsigned D,unsigned F> static FI __m512i scale(__m512i x){
        if constexpr((A|C|D|F)==0)return _mm512_setzero_si512();
        else if constexpr(A==1&&C==1&&D==1&&F==1)return x;
        else{
            constexpr auto lo=[](unsigned a){return std::uint32_t((a&255)|(((a&255)^(a>>8))<<8))*0x10001U;};
            const auto p=_mm512_setr_epi32(lo(A),lo(A),lo(A),lo(A),lo(C),lo(C),lo(C),lo(C),
                                          lo(D),lo(D),lo(D),lo(D),lo(F),lo(F),lo(F),lo(F));
            auto out=_mm512_gf2p8mul_epi8(x,p);
            if constexpr((A|C|D|F)>=256){
                constexpr auto hi=[](unsigned a){return std::uint32_t(m::byteMultiply(0x20,a>>8)|((a>>8)<<8))*0x10001U;};
                const auto q=_mm512_setr_epi32(hi(A),hi(A),hi(A),hi(A),hi(C),hi(C),hi(C),hi(C),
                                              hi(D),hi(D),hi(D),hi(D),hi(F),hi(F),hi(F),hi(F));
                const auto swap=_mm512_broadcast_i32x4(_mm_set_epi32(0x0e0f0c0d,0x0a0b0809,0x06070405,0x02030001));
                out=_mm512_xor_si512(out,_mm512_gf2p8mul_epi8(_mm512_shuffle_epi8(x,swap),q));
            }
            return out;
        }
    }
    template<unsigned I,unsigned A,unsigned C,bool Inverse> static FI __m512i small(__m512i x){
        auto a=_mm512_shuffle_i32x4(x,x,I==0?0xa0:0x44),b=_mm512_shuffle_i32x4(x,x,I==0?0xf5:0xee);
        if constexpr(Inverse){
            auto high=_mm512_xor_si512(b,scale<A,A,C,C>(a));
            return _mm512_mask_blend_epi64(I==0?0xcc:0xf0,_mm512_xor_si512(a,high),high);
        }else{
            auto sum=_mm512_xor_si512(a,b);
            return _mm512_mask_blend_epi64(I==0?0xcc:0xf0,sum,_mm512_xor_si512(b,scale<A,A,C,C>(sum)));
        }
    }
    template<unsigned L,unsigned Offset,bool Shift,bool Inverse> static FI void node(reset::E* x){
        if constexpr(L==2){
            constexpr auto a=pt(Offset)^(Shift?basis[H]:0),c=pt(Offset+2)^(Shift?basis[H]:0);
            constexpr auto b=pt(Offset>>1)^(Shift?basis[H-1]:0);
            auto z=_mm512_load_si512(x);
            if constexpr(Inverse){z=small<1,b,b,true>(z);z=small<0,a,c,true>(z);}
            else{z=small<0,a,c,false>(z);z=small<1,b,b,false>(z);}
            _mm512_store_si512(x,z);
        }else{
            constexpr unsigned half=1U<<(L-1),c=pt(Offset>>(L-1))^(Shift?basis[H-L+1]:0);
            if constexpr(!Inverse){node<L-1,Offset,Shift,Inverse>(x);node<L-1,Offset+half,Shift,Inverse>(x+half);}
            for(unsigned j=0;j<half;j+=4){auto a=_mm512_load_si512(x+j),b=_mm512_load_si512(x+half+j);
                if constexpr(Inverse){b=_mm512_xor_si512(b,scale<c,c,c,c>(a));a=_mm512_xor_si512(a,b);}
                else{a=_mm512_xor_si512(a,b);b=_mm512_xor_si512(b,scale<c,c,c,c>(a));}
                _mm512_store_si512(x+j,a);_mm512_store_si512(x+half+j,b);
            }
            if constexpr(Inverse){node<L-1,Offset,Shift,Inverse>(x);node<L-1,Offset+half,Shift,Inverse>(x+half);}
        }
    }
    template<bool Fast=true> FI void apply(reset::E* x)const{
        if constexpr(Fast){node<H,0,true,false>(x);node<H,0,false,true>(x);return;}
        // C0^-1 Cdelta = (Eval_delta Interp_0)^T, symmetric Cauchy here.
        K::fft(x,shift);inverseUpper<H-1>(x,zero);
        for(unsigned off=0;off<B;off+=4){
            auto v=_mm512_loadu_si512(x+off);
            v=inverseSmall<1>(v,zero[B/4-1+off/4]);
            v=inverseSmall<0>(v,zero[B/2-1+off/4]);
            _mm512_storeu_si512(x+off,v);
        }
    }
    NI void applyOutlined(reset::E* x)const{apply(x);}
};
E random(std::mt19937_64& rng){return {rng(),rng()};}
std::uint64_t digest(const E* p,std::size_t n){
    std::uint64_t v=1;for(std::size_t j=0;j<n;++j)v=(v^p[j].lo^p[j].hi)*1099511628211ULL;return v;
}
struct Plan {
    std::size_t K,N,blocks;
    Aligned<E> u,v;
    Aligned<std::uint32_t> route;
    explicit Plan(unsigned logK):K(std::size_t(1)<<logK),N(2*K),blocks(K/128),u(N),v(N),route(N){
        std::mt19937_64 rng(20261022);
        for(std::size_t j=0;j<N;++j){
            do{u[j]=random(rng);}while(!(u[j].lo|u[j].hi));
            do{v[j]=random(rng);}while(!(v[j].lo|v[j].hi));
        }
        for(unsigned j=0;j<256;++j){auto p=route.data()+j*blocks;std::iota(p,p+blocks,0U);std::shuffle(p,p+blocks,rng);}
    }
};
template<unsigned H,bool Routed,bool Fast=true,bool Outline=false,bool Stage=false> NI void inner(const Plan& p,const Map<H>& map,const E* x,E* y){
    constexpr unsigned B=1U<<H;
    alignas(64) std::array<reset::E,B> state{},work{},inputBatch;
    for(std::size_t off=0;off<p.N;off+=B){
        if constexpr(Routed&&Stage){
            const auto* source=x+(off/p.blocks)*p.blocks;
            for(unsigned j=0;j<B;j+=16){
                auto a=rs::gatherFour(source,p.route.data()+off+j);
                auto b=rs::gatherFour(source,p.route.data()+off+j+4);
                auto c=rs::gatherFour(source,p.route.data()+off+j+8);
                auto d=rs::gatherFour(source,p.route.data()+off+j+12);
                _mm512_store_si512(inputBatch.data()+j,a);_mm512_store_si512(inputBatch.data()+j+4,b);
                _mm512_store_si512(inputBatch.data()+j+8,c);_mm512_store_si512(inputBatch.data()+j+12,d);
            }
        }
        for(unsigned j=0;j<B;j+=4)_mm512_store_si512(work.data()+j,mul128(
            _mm512_load_si512(state.data()+j),_mm512_load_si512(p.v.data()+off+j)));
        if constexpr(Outline)map.applyOutlined(work.data());
        else map.template apply<Fast>(work.data());
        for(unsigned j=0;j<B;j+=4){
            auto feedback=mul128(_mm512_load_si512(work.data()+j),_mm512_load_si512(p.u.data()+off+j));
            __m512i input;
            if constexpr(Routed&&Stage)input=_mm512_load_si512(inputBatch.data()+j);
            else if constexpr(Routed){const auto base=(off/p.blocks)*p.blocks;
                input=rs::gatherFour(x+base,p.route.data()+off+j);
            }else input=_mm512_load_si512(x+off+j);
            auto output=_mm512_xor_si512(input,feedback);
            _mm512_store_si512(y+off+j,output);
            _mm512_store_si512(state.data()+j,_mm512_xor_si512(_mm512_load_si512(state.data()+j),output));
        }
    }
}
template<unsigned H,unsigned Mode,bool Fast=true,bool Specialized=true> NI void stream(const Plan& p,const Map<H>& map,const E* x,E* y){
    constexpr unsigned B=1U<<H;alignas(64) std::array<reset::E,B> work;
    for(std::size_t off=0;off<p.N;off+=B){
        for(unsigned j=0;j<B;j+=4){auto z=_mm512_load_si512(x+off+j);
            if constexpr(Mode!=1)z=mul128<Specialized>(z,_mm512_load_si512(p.v.data()+off+j));
            _mm512_store_si512(work.data()+j,z);
        }
        if constexpr(Mode!=0)map.template apply<Fast>(work.data());
        for(unsigned j=0;j<B;j+=4){auto z=_mm512_load_si512(work.data()+j);
            if constexpr(Mode==2)z=mul128(z,_mm512_load_si512(p.u.data()+off+j));
            _mm512_store_si512(y+off+j,z);
        }
    }
}
template<class Fn> double measure(unsigned reps,Fn&& fn){
    std::vector<double> times;
    for(unsigned i=0;i<reps+3;++i){auto start=Clock::now();fn();auto end=Clock::now();
        if(i>=3)times.push_back(std::chrono::duration<double,std::milli>(end-start).count());}
    std::sort(times.begin(),times.end());return times[times.size()/2];
}
// Representation diagnostic only. The sparse-byte basis cannot be substituted
// into the full tower encoder without compatible outer/MDS maps or conversions.
template<bool Sparse> NI void diagonalBasisControl(const Plan& p,const E* x,E* y){
    for(std::size_t j=0;j<p.N;j+=4){auto a=_mm512_load_si512(x+j),b=_mm512_load_si512(p.v.data()+j);
        if constexpr(Sparse)a=spin::research::char2::field128_lift::sparse_byte::multiply4(a,b);
        else a=mul128(a,b);
        _mm512_store_si512(y+j,a);
    }
}
template<unsigned H> void checkMap(){
    Map<H> map;constexpr unsigned B=1U<<H;alignas(64) typename Map<H>::K::Buffer x{};
    std::mt19937_64 rng(913+H);alignas(64) typename Map<H>::K::Buffer reference{};
    for(auto& row:x)for(auto& a:row.lane)a=rng();
    reference=x;map.apply(x.data());map.template apply<false>(reference.data());reset::require(x==reference);
    for(unsigned col:{0U,1U,B/2,B-1}){
        x={};for(auto& a:x[col].lane)a=1;map.apply(x.data());
        for(unsigned j=0;j<B;++j)for(auto a:x[j].lane)
            reset::require(a==reset::inverse(reset::point(B)^reset::point(j)^reset::point(col)));
        map.apply(x.data());for(unsigned j=0;j<B;++j)for(auto a:x[j].lane)reset::require(a==(j==col));
    }
    std::cout<<"Cauchy/involution B="<<B<<" PASS\n";
}
void checks(){
    std::uint32_t t=mu,trace=t;for(unsigned i=1;i<32;++i){t=m::towerMultiplyOracle(t,t);trace^=t;}reset::require(trace==1);
    std::mt19937_64 rng(713);alignas(64) E a[4],b[4],out[4];
    for(unsigned i=0;i<256;++i){
        for(unsigned j=0;j<4;++j){a[j]=random(rng);b[j]=random(rng);}
        if(i<128)a[0]=std::bit_cast<E>(std::array<std::uint32_t,4>{
            i/32==0?1U<<(i%32):0,i/32==1?1U<<(i%32):0,i/32==2?1U<<(i%32):0,i/32==3?1U<<(i%32):0});
        _mm512_store_si512(out,mul128(_mm512_load_si512(a),_mm512_load_si512(b)));
        for(unsigned j=0;j<4;++j)reset::require(out[j]==oracle(a[j],b[j]));
    }
    checkMap<8>();checkMap<9>();
    Plan p(16);Map<9> map;Aligned<E> x(p.N),y(p.N),z(p.N),planes(p.N);
    for(std::size_t i=0;i<p.N;++i)x[i]=random(rng);
    inner<9,false>(p,map,x.data(),y.data());
    alignas(64) std::array<reset::E,512> state{},work{};
    for(std::size_t off=0;off<p.N;off+=512){
        for(unsigned j=0;j<512;j+=4)_mm512_store_si512(work.data()+j,mul128(_mm512_load_si512(state.data()+j),_mm512_load_si512(p.v.data()+off+j)));
        map.apply(work.data());
        for(unsigned j=0;j<512;j+=4){auto out=_mm512_load_si512(y.data()+off+j);
            auto f=mul128(_mm512_load_si512(work.data()+j),_mm512_load_si512(p.u.data()+off+j));
            _mm512_store_si512(z.data()+off+j,_mm512_xor_si512(out,f));
            _mm512_store_si512(state.data()+j,_mm512_xor_si512(_mm512_load_si512(state.data()+j),out));}
    }
    reset::require(!std::memcmp(x.data(),z.data(),p.N*sizeof(E)));
    rs::routeRegional(p.blocks,p.route.data(),x.data(),planes.data());
    inner<9,false>(p,map,planes.data(),z.data());inner<9,true>(p,map,x.data(),y.data());
    reset::require(!std::memcmp(y.data(),z.data(),p.N*sizeof(E)));
    inner<9,true,true,true>(p,map,x.data(),y.data());reset::require(!std::memcmp(y.data(),z.data(),p.N*sizeof(E)));
    inner<9,true,true,false,true>(p,map,x.data(),y.data());reset::require(!std::memcmp(y.data(),z.data(),p.N*sizeof(E)));
    // Independent scalar stream for two batches, including actual nonzero state.
    std::array<E,512> scalar{},f{};
    inner<9,false>(p,map,x.data(),y.data());
    for(unsigned t=0;t<2;++t){
        for(unsigned j=0;j<512;++j){E sum{};
            for(unsigned k=0;k<512;++k){auto value=oracle(scalar[k],p.v[512*t+k]);
                auto coefficient=reset::inverse(reset::point(512)^reset::point(j)^reset::point(k));
                auto raw=std::bit_cast<std::array<std::uint32_t,4>>(value);
                for(auto& a:raw)a=m::towerMultiplyOracle(a,coefficient);
                sum^=std::bit_cast<E>(raw);}
            f[j]=sum;
        }
        for(unsigned j=0;j<512;++j){auto output=x[512*t+j]^oracle(p.u[512*t+j],f[j]);
            reset::require(output==y[512*t+j]);scalar[j]^=output;}
    }
    std::cout<<"1024 full-field products, trace/field, inverse, routed equivalence, scalar stream PASS\n";
}
void bench(unsigned logK,unsigned reps){
    Plan p(logK);Map<9> map;Aligned<E> x(p.N),planes(p.N),routed(p.N),y(p.N),ref(p.N);
    std::mt19937_64 rng(123);for(std::size_t i=0;i<p.N;++i)x[i]=random(rng);
    std::cout<<std::fixed<<std::setprecision(4)<<"K="<<p.K<<" N="<<p.N<<" outer=256 state=512 running_sum key_bytes="<<2*p.N*sizeof(E)<<"\n";
    auto show=[&](const char* name,auto&& fn){auto ms=measure(reps,fn);std::cout<<name<<"_ms="<<ms<<" checksum="<<digest(y.data(),p.N)<<"\n"<<std::flush;};
    show("one_diagonal",[&]{stream<9,0>(p,map,x.data(),y.data());});
    show("one_diagonal_control",[&]{stream<9,0,true,false>(p,map,x.data(),y.data());});
    show("tower_diagonal_flat",[&]{diagonalBasisControl<false>(p,x.data(),y.data());});
    show("sparse_basis_diagonal_diagnostic",[&]{diagonalBasisControl<true>(p,x.data(),y.data());});
    show("fixed_mds_control",[&]{stream<9,1,false>(p,map,x.data(),y.data());});
    show("fixed_mds",[&]{stream<9,1>(p,map,x.data(),y.data());});
    show("feedback",[&]{stream<9,2>(p,map,x.data(),y.data());});
    show("inner_contiguous",[&]{inner<9,false>(p,map,x.data(),y.data());});
    show("outer_only",[&]{rs::encodeRegional(p.blocks,x.data(),y.data());});
    show("outer_route",[&]{rs::encodeRegional(p.blocks,x.data(),planes.data());rs::routeRegional(p.blocks,p.route.data(),planes.data(),y.data());});
    rs::encodeRegional(p.blocks,x.data(),planes.data());rs::routeRegional(p.blocks,p.route.data(),planes.data(),routed.data());
    inner<9,false>(p,map,routed.data(),ref.data());inner<9,true>(p,map,planes.data(),y.data());
    reset::require(!std::memcmp(y.data(),ref.data(),p.N*sizeof(E)));
    show("inner_routed",[&]{inner<9,true>(p,map,planes.data(),y.data());});
    show("inner_routed_outlined",[&]{inner<9,true,true,true>(p,map,planes.data(),y.data());});
    show("inner_routed_staged",[&]{inner<9,true,true,false,true>(p,map,planes.data(),y.data());});
    reset::require(!std::memcmp(y.data(),ref.data(),p.N*sizeof(E)));
    show("full_fused",[&]{rs::encodeRegional(p.blocks,x.data(),planes.data());inner<9,true>(p,map,planes.data(),y.data());});
    show("full_materialized",[&]{rs::encodeRegional(p.blocks,x.data(),planes.data());rs::routeRegional(p.blocks,p.route.data(),planes.data(),routed.data());inner<9,false>(p,map,routed.data(),y.data());});
    show("full_staged",[&]{rs::encodeRegional(p.blocks,x.data(),planes.data());inner<9,true,true,false,true>(p,map,planes.data(),y.data());});
    reset::require(!std::memcmp(y.data(),ref.data(),p.N*sizeof(E)));
    std::cout<<"full_materialized_oracle=PASS\n";
}
}
int main(int argc,char** argv){
    try{
        if(argc==2&&std::string(argv[1])=="--check"){sandwich::checks();return 0;}
        if(argc!=4||std::string(argv[1])!="--bench")return 2;
        unsigned logK=std::stoul(argv[2]),reps=std::stoul(argv[3]);
        if(logK<16||logK>20||reps<3||reps>31||reps%2==0)return 2;
        sandwich::bench(logK,reps);
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
