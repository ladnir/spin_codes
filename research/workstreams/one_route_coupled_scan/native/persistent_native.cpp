// Research-only native gate for PERSISTENT_MDS.md. Four GF32 lanes per
// 128-bit record share one code. Setup is outside every timed operation.
#include "Persistent32.h"
#ifndef SPIN_OR_ADDITIVE
#define SPIN_OR_ADDITIVE 0
#endif
#ifndef SPIN_OR_FULL_CANTOR
#define SPIN_OR_FULL_CANTOR 0
#endif
#ifndef SPIN_OR_SINGLE_DIAGONAL
#define SPIN_OR_SINGLE_DIAGONAL 0
#endif
#if SPIN_OR_SINGLE_DIAGONAL != 0 && SPIN_OR_SINGLE_DIAGONAL != 1
#error "SPIN_OR_SINGLE_DIAGONAL must be 0 or 1"
#endif
#if SPIN_OR_SINGLE_DIAGONAL && !SPIN_OR_ADDITIVE
#error "experimental single-diagonal feedback requires additive mode"
#endif
#ifndef SPIN_OR_RAW_KEYS
#define SPIN_OR_RAW_KEYS 0
#endif
#if SPIN_OR_RAW_KEYS != 0 && SPIN_OR_RAW_KEYS != 1
#error "SPIN_OR_RAW_KEYS must be 0 or 1"
#endif
#if SPIN_OR_RAW_KEYS && !SPIN_OR_ADDITIVE
#error "compact random keys are instantiated only for additive profiles"
#endif
#ifndef SPIN_OR_PACKET
#define SPIN_OR_PACKET 1
#endif
#ifndef SPIN_OR_OUTER_TILE
#define SPIN_OR_OUTER_TILE 4
#endif
#ifndef SPIN_OR_OUTER_STREAM
#define SPIN_OR_OUTER_STREAM 0
#endif
#if SPIN_OR_OUTER_TILE != 4 && SPIN_OR_OUTER_TILE != 16 && SPIN_OR_OUTER_TILE != 32
#error "outer tile must contain 4, 16, or 32 rows"
#endif
#if SPIN_OR_OUTER_STREAM != 0 && SPIN_OR_OUTER_STREAM != 1
#error "outer streaming flag must be 0 or 1"
#endif
#if SPIN_OR_PACKET != 1 && SPIN_OR_PACKET != 2 && SPIN_OR_PACKET != 4
#error "SPIN_OR_PACKET must be 1, 2, or 4"
#endif
#if SPIN_OR_PACKET != 1 && !SPIN_OR_ADDITIVE
#error "packet routing is instantiated only for additive profiles"
#endif
#if SPIN_OR_ADDITIVE
#include "AdditivePersistent32.h"
#endif
#include <bit>
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace pn {
namespace tr=spin::research::one_route;
namespace gf=spin::research::char2::medium32;
namespace wide=spin::research::char2::wide32;
using E=tr::Element;
using Key=wide::Coefficient;
#if SPIN_OR_RAW_KEYS
using RandomKey=std::uint32_t;
#else
using RandomKey=Key;
#endif
using Clock=std::chrono::steady_clock;
inline constexpr unsigned packetWidth=SPIN_OR_PACKET;
inline constexpr unsigned randomDiagonals=SPIN_OR_SINGLE_DIAGONAL?1:2;
inline constexpr const char* profileLabel=SPIN_OR_SINGLE_DIAGONAL?"experimental uncertified":"certified";
#define NI __attribute__((noinline))

inline __attribute__((always_inline)) __m512i randomMultiply(__m512i value,const RandomKey* keys) noexcept {
#if SPIN_OR_RAW_KEYS
    return gf::towerMultiply512(value,gf::broadcastFour32(keys));
#else
    return wide::multiply512(value,keys);
#endif
}
inline __attribute__((always_inline)) std::uint32_t randomScalar(const RandomKey& key) noexcept {
#if SPIN_OR_RAW_KEYS
    return key;
#else
    return wide::recoverScalar(key);
#endif
}

void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
template<class T> struct Aligned {
    struct Release {void operator()(T* p)const noexcept{std::free(p);}};
    std::unique_ptr<T,Release> pointer;
    std::size_t count;
    explicit Aligned(std::size_t n):count(n){
        pointer.reset(static_cast<T*>(std::aligned_alloc(64,(n*sizeof(T)+63)&~std::size_t(63))));
        if(!pointer)throw std::bad_alloc();
        std::memset(pointer.get(),0,n*sizeof(T));
    }
    T* data(){return pointer.get();} const T* data()const{return pointer.get();}
    T& operator[](std::size_t j){return data()[j];}
    const T& operator[](std::size_t j)const{return data()[j];}
};
E random(std::mt19937_64& rng){return {rng(),rng()};}
std::uint32_t power32(std::uint32_t a,std::uint64_t n){
    std::uint32_t r=1;
    for(;n;n>>=1){if(n&1)r=gf::towerMultiplyOracle(r,a);a=gf::towerMultiplyOracle(a,a);}
    return r;
}
std::uint32_t primitive(){
    constexpr std::uint64_t order=0xffffffffULL;
    for(std::uint32_t a=0x10000;a<0x10100;++a){
        bool ok=power32(a,order)==1;
        for(unsigned prime:{3,5,17,257,65537})ok=ok&&power32(a,order/prime)!=1;
        if(ok)return a;
    }
    throw std::runtime_error("bounded primitive-element search failed");
}
E scale32(E x,std::uint32_t a){
    std::array<std::uint32_t,4> words;std::memcpy(words.data(),&x,sizeof(x));
    for(auto& word:words)word=gf::towerMultiplyOracle(word,a);
    std::memcpy(&x,words.data(),sizeof(x));return x;
}
E scale8(E x,std::uint8_t a){
    std::array<std::uint8_t,16> bytes;std::memcpy(bytes.data(),&x,sizeof(x));
    for(auto& byte:bytes)byte=gf::byteOracle(byte,a);
    std::memcpy(&x,bytes.data(),sizeof(x));return x;
}

// Bounded cross-representation check: all 32-by-32 binary basis products,
// then dense and zero fixtures. This does not change any Plan's RNG stream.
void checkRandomKeys(){
    std::mt19937_64 rng(929);
    alignas(64) E input[4],rawOutput[4],preparedOutput[4],selectedOutput[4],expected[4];
    alignas(64) std::uint32_t raw[4];
    alignas(64) Key prepared[4];
    alignas(64) RandomKey selected[4];
    for(unsigned fixture=0;fixture<290;++fixture){
        for(unsigned j=0;j<4;++j){
            if(fixture<256){
                const unsigned pair=4*fixture+j;
                const auto a=std::uint64_t{1}<<(pair/32);
                input[j]={a|(a<<32),a|(a<<32)};
                raw[j]=1U<<(pair%32);
            }else{
                input[j]=fixture==289?E{}:random(rng);
                raw[j]=fixture==288?0:std::uint32_t(rng());
            }
            prepared[j]=wide::prepareCoefficientDirect(raw[j]);
#if SPIN_OR_RAW_KEYS
            selected[j]=raw[j];
#else
            selected[j]=prepared[j];
#endif
            require(randomScalar(selected[j])==raw[j],"random coefficient recovery mismatch");
            expected[j]=scale32(input[j],raw[j]);
        }
        const auto value=_mm512_load_si512(input);
        _mm512_store_si512(rawOutput,gf::towerMultiply512(value,gf::broadcastFour32(raw)));
        _mm512_store_si512(preparedOutput,wide::multiply512(value,prepared));
        _mm512_store_si512(selectedOutput,randomMultiply(value,selected));
        require(std::memcmp(rawOutput,expected,sizeof(expected))==0,"raw random product oracle mismatch");
        require(std::memcmp(preparedOutput,expected,sizeof(expected))==0,"prepared random product oracle mismatch");
        require(std::memcmp(selectedOutput,expected,sizeof(expected))==0,"selected random product oracle mismatch");
    }
    std::cout<<"PASS raw/prepared random keys: 1024 basis products, 128 dense records, zero fixtures\n";
}

// Independent root-product basis and Lagrange interpolation, checking only.
template<unsigned P> struct Oracle {
    std::array<std::array<std::uint8_t,P>,P> basis{},g{},h{};
    Oracle(){
        std::array<std::array<std::uint8_t,P+1>,8> sub{};
        for(unsigned level=0;(1U<<level)<P;++level){
            const unsigned length=1U<<level;auto& s=sub[level];s[0]=1;
            for(unsigned root=0;root<length;++root){
                const auto a=tr::TransformCantor<P>::evaluationPoint(root);
                for(unsigned d=root+1;d>0;--d)s[d]=s[d-1]^gf::byteOracle(s[d],a);
                s[0]=gf::byteOracle(s[0],a);
            }
        }
        basis[0][0]=1;
        for(unsigned degree=1;degree<P;++degree){
            unsigned level=0;while((2U<<level)<=degree)++level;
            const auto half=1U<<level,low=degree-half;
            for(unsigned i=0;i<=low;++i)for(unsigned j=0;j<=half;++j)
                basis[degree][i+j]^=gf::byteOracle(basis[low][i],sub[level][j]);
        }
        for(unsigned point=0;point<P;++point){
            const auto a=tr::TransformCantor<P>::evaluationPoint(point);
            for(unsigned degree=0;degree<P;++degree){
                std::uint8_t value=0;
                for(unsigned d=degree+1;d-- >0;)value=gf::byteOracle(value,a)^basis[degree][d];
                g[degree][point]=value;
            }
            std::array<std::uint8_t,P> poly{};poly[0]=1;
            unsigned degree=0;std::uint8_t denominator=1;
            for(unsigned other=0;other<P;++other)if(other!=point){
                const auto z=tr::TransformCantor<P>::evaluationPoint(other);
                denominator=gf::byteOracle(denominator,a^z);
                for(unsigned d=degree+1;d>0;--d)poly[d]=poly[d-1]^gf::byteOracle(poly[d],z);
                poly[0]=gf::byteOracle(poly[0],z);++degree;
            }
            const auto inv=gf::bytePower(denominator,254);
            for(auto& c:poly)c=gf::byteOracle(c,inv);
            for(unsigned d=P;d-- >0;){
                const auto value=poly[d];require(basis[d][d]==1,"oracle basis not monic");
                h[d][point]=value;
                for(unsigned j=0;j<d;++j)poly[j]^=gf::byteOracle(value,basis[d][j]);
            }
        }
    }
};

template<unsigned B,unsigned b> struct Plan {
    std::size_t K,N,blocks,groups,batches;
    std::uint32_t alpha;
    Aligned<std::uint32_t> route;
    Aligned<RandomKey> keys;
    alignas(64) std::array<Key,b> fixed;
    std::array<std::uint32_t,b> powers;
#if SPIN_OR_ADDITIVE
    std::array<typename tr::AdditivePersistent32<b>::Shift,32> shifts;
#endif
    explicit Plan(std::size_t k):K(k),N(2*k),blocks(N/B),groups(blocks/packetWidth),batches(N/b),
        alpha(primitive()),route(N/packetWidth),keys(randomDiagonals*N){
        require(K>=256&&K%(2*B)==0&&N%b==0,"invalid research dimensions");
        require(blocks%packetWidth==0,"incomplete outer-row packet");
        require(batches<=16843009,"persistent MDS coset limit exceeded");
#if SPIN_OR_ADDITIVE
        for(unsigned k=std::countr_zero(b);k<=std::countr_zero(b)+std::bit_width(batches);++k)
            shifts[k]=tr::AdditivePersistent32<b>::prepare(k);
#endif
        std::uint32_t value=1;
        for(unsigned i=0;i<b;++i){powers[i]=value;fixed[i]=wide::prepareCoefficient(value);
            value=gf::towerMultiplyOracle(value,alpha);}
        std::mt19937_64 rng(113);
        // Reproducible sampling only; proofs concern independent ideal uniform coins.
        for(unsigned column=0;column<B;++column){
            // A label selects a fixed group of adjacent outer rows. Its lanes
            // retain their order; each region independently permutes groups.
            auto* p=route.data()+column*groups;std::iota(p,p+groups,0U);std::shuffle(p,p+groups,rng);
        }
#if SPIN_OR_RAW_KEYS
        for(std::size_t i=0;i<keys.count;++i)keys[i]=std::uint32_t(rng());
#else
        const auto& table=wide::byteTables();
        for(std::size_t i=0;i<keys.count;++i)keys[i]=wide::prepareFromTables(std::uint32_t(rng()),table);
#endif
    }
};

// Only used by the packet variants. The aligned outer-plane allocation and
// adjacent-row layout make every width-four source 64-byte aligned, and every
// width-two source 32-byte aligned. One label is stored per source packet.
inline __m512i packetFour(const E* source,const std::uint32_t* labels){
    if constexpr(packetWidth==4){
        return _mm512_load_si512(source+4*std::size_t(labels[0]));
    }else{
        const auto low=_mm256_load_si256(reinterpret_cast<const __m256i*>(source+2*std::size_t(labels[0])));
        const auto high=_mm256_load_si256(reinterpret_cast<const __m256i*>(source+2*std::size_t(labels[1])));
        return _mm512_inserti64x4(_mm512_castsi256_si512(low),high,1);
    }
}

template<unsigned B,unsigned Rows> inline __attribute__((always_inline))
void outerTile(const E* input,E* planes,std::size_t blocks,std::size_t row){
    alignas(64) E scratch[Rows*B];
    for(unsigned group=0;group<Rows;group+=4){
        auto* values=scratch+group*B;
        for(unsigned j=0;j<B/2;++j){
            values[4*j]=input[(row+group+0)*(B/2)+j];
            values[4*j+1]=input[(row+group+1)*(B/2)+j];
            values[4*j+2]=input[(row+group+2)*(B/2)+j];
            values[4*j+3]=input[(row+group+3)*(B/2)+j];
        }
        std::memset(values+2*B,0,2*B*sizeof(E));
        tr::TransformCantor<B>::evaluateNovelFour(values);
    }
    // Finish several adjacent cache lines in one plane before advancing
    // to the next power-of-two-strided plane. The scalar code is unchanged.
    for(unsigned j=0;j<B;++j)for(unsigned group=0;group<Rows;group+=4){
        const auto value=_mm512_load_si512(scratch+group*B+4*j);
#if SPIN_OR_OUTER_STREAM
        _mm512_stream_si512(reinterpret_cast<__m512i*>(planes+j*blocks+row+group),value);
#else
        _mm512_store_si512(planes+j*blocks+row+group,value);
#endif
    }
}
template<unsigned B,unsigned b> NI void outer(const Plan<B,b>& p,const E* input,E* planes){
#if SPIN_OR_OUTER_TILE == 4 && !SPIN_OR_OUTER_STREAM
    // Retain the original schedule exactly as the default A/B control.
    alignas(64) E scratch[4*B];
    for(std::size_t row=0;row<p.blocks;row+=4){
        for(unsigned j=0;j<B/2;++j){
            scratch[4*j]=input[(row+0)*(B/2)+j];scratch[4*j+1]=input[(row+1)*(B/2)+j];
            scratch[4*j+2]=input[(row+2)*(B/2)+j];scratch[4*j+3]=input[(row+3)*(B/2)+j];
        }
        std::memset(scratch+2*B,0,2*B*sizeof(E));
        tr::TransformCantor<B>::evaluateNovelFour(scratch);
        for(unsigned j=0;j<B;++j)_mm512_store_si512(planes+j*p.blocks+row,_mm512_load_si512(scratch+4*j));
    }
#else
    std::size_t row=0;
    for(;row+SPIN_OR_OUTER_TILE<=p.blocks;row+=SPIN_OR_OUTER_TILE)
        outerTile<B,SPIN_OR_OUTER_TILE>(input,planes,p.blocks,row);
    for(;row<p.blocks;row+=4)outerTile<B,4>(input,planes,p.blocks,row);
#if SPIN_OR_OUTER_STREAM
    // Order all streaming stores before the route reads the same allocation.
    _mm_mfence();
#endif
#endif
}
template<unsigned B,unsigned b,bool Routed> NI void inner(const Plan<B,b>& p,const E* input,E* output){
    alignas(64) E state[b]{},feedback[b],y[b];
#if SPIN_OR_ADDITIVE
    alignas(64) E shifted[b];
#if !SPIN_OR_SINGLE_DIAGONAL
    alignas(64) E weighted[b];
#endif
#endif
    std::size_t base=0,next=p.blocks;
    for(std::size_t t=0;t<p.batches;++t){
#if SPIN_OR_ADDITIVE
        const auto& shift=p.shifts[std::countr_zero(b)+std::countr_zero(t+1)];
        tr::AdditivePersistent32<b>::apply(shift,state,shifted);
#if SPIN_OR_SINGLE_DIAGONAL
        // Experimental recurrence: V=DZ, Y=X+r.*V, Z'=V+Y.
        // Its random coordinate set is support(V), not the whole live batch.
        for(unsigned j=0;j<b;j+=4)_mm512_store_si512(feedback+j,
            randomMultiply(_mm512_load_si512(shifted+j),p.keys.data()+b*t+j));
#else
        for(unsigned j=0;j<b;j+=4)_mm512_store_si512(weighted+j,
            randomMultiply(_mm512_load_si512(state+j),p.keys.data()+2*b*t+j));
        tr::AdditivePersistent32<b>::apply(shift,weighted,feedback);
        for(unsigned j=0;j<b;j+=4)_mm512_store_si512(feedback+j,_mm512_xor_si512(
            _mm512_load_si512(feedback+j),randomMultiply(_mm512_load_si512(shifted+j),p.keys.data()+2*b*t+b+j)));
#endif
#else
        tr::Persistent32<b>::feedback(state,p.keys.data()+2*b*t,feedback);
#endif
        const auto offset=b*t;
        if constexpr(Routed){
            if(p.blocks>=b){
                if(offset==next){base=next;next+=p.blocks;}
                if constexpr(packetWidth==1){
                    const auto* source=input+base;const auto* labels=p.route.data()+offset;
                    for(unsigned j=0;j<b;j+=4){
                        auto x=_mm512_castsi128_si512(_mm_loadu_si128(reinterpret_cast<const __m128i*>(source+labels[j])));
                        x=_mm512_inserti32x4(x,_mm_loadu_si128(reinterpret_cast<const __m128i*>(source+labels[j+1])),1);
                        x=_mm512_inserti32x4(x,_mm_loadu_si128(reinterpret_cast<const __m128i*>(source+labels[j+2])),2);
                        x=_mm512_inserti32x4(x,_mm_loadu_si128(reinterpret_cast<const __m128i*>(source+labels[j+3])),3);
                        const auto z=_mm512_xor_si512(x,_mm512_load_si512(feedback+j));
                        _mm512_store_si512(output+offset+j,z);_mm512_store_si512(y+j,z);
                    }
                }else{
                    const auto* source=input+base;const auto* labels=p.route.data()+offset/packetWidth;
                    for(unsigned j=0;j<b;j+=4){
                        const auto x=packetFour(source,labels+j/packetWidth);
                        const auto z=_mm512_xor_si512(x,_mm512_load_si512(feedback+j));
                        _mm512_store_si512(output+offset+j,z);_mm512_store_si512(y+j,z);
                    }
                }
            }else if constexpr(packetWidth==1){
                for(unsigned j=0;j<b;++j){const auto i=offset+j,r=(i/p.blocks)*p.blocks;
                    y[j]=output[i]=input[r+p.route[i]]^feedback[j];}
            }else{
                // Small check/benchmark sizes may put several whole regions
                // in one batch. Regions still contain a multiple of four rows.
                for(unsigned j=0;j<b;j+=4){const auto i=offset+j,r=(i/p.blocks)*p.blocks;
                    const auto x=packetFour(input+r,p.route.data()+i/packetWidth);
                    const auto z=_mm512_xor_si512(x,_mm512_load_si512(feedback+j));
                    _mm512_store_si512(output+i,z);_mm512_store_si512(y+j,z);
                }
            }
        }else for(unsigned j=0;j<b;j+=4){
            const auto z=_mm512_xor_si512(_mm512_load_si512(input+offset+j),_mm512_load_si512(feedback+j));
            _mm512_store_si512(output+offset+j,z);_mm512_store_si512(y+j,z);
        }
#if SPIN_OR_ADDITIVE
        for(unsigned j=0;j<b;j+=4)_mm512_store_si512(state+j,_mm512_xor_si512(
            _mm512_load_si512(shifted+j),_mm512_load_si512(y+j)));
#else
        tr::Persistent32<b>::update(state,y,p.fixed.data());
#endif
    }
}
template<unsigned B,unsigned b> NI void routeOnly(const Plan<B,b>& p,const E* input,E* output){
    if constexpr(packetWidth==1){
        for(unsigned r=0;r<B;++r){const auto base=r*p.blocks;
            for(std::size_t i=0;i<p.blocks;++i)output[base+i]=input[base+p.route[base+i]];}
    }else{
        for(unsigned r=0;r<B;++r){const auto base=r*p.blocks;const auto* labels=p.route.data()+r*p.groups;
            for(std::size_t i=0;i<p.blocks;i+=4)
                _mm512_store_si512(output+base+i,packetFour(input+base,labels+i/packetWidth));}
    }
}
template<unsigned B,unsigned b> NI void full(const Plan<B,b>& p,const E* input,E* workspace,E* output){
    outer(p,input,workspace);inner<B,b,true>(p,workspace,output);
}
// Independent scalar packet-index oracle: do not share packetFour's indexing.
template<unsigned B,unsigned b> void checkRoute(const Plan<B,b>& p,const E* source,const E* routed){
    for(unsigned r=0;r<B;++r){
        std::vector<bool> seen(p.groups);
        for(std::size_t group=0;group<p.groups;++group){
            const auto label=p.route[r*p.groups+group];
            require(label<p.groups&&!seen[label],"packet labels not bijective");seen[label]=true;
            for(unsigned lane=0;lane<packetWidth;++lane)
                require(routed[r*p.blocks+packetWidth*group+lane]==source[r*p.blocks+packetWidth*label+lane],
                    "routed packet source or lane order mismatch");
        }
    }
}
template<unsigned b> void feedbackOracle(const Oracle<2*b>& table,const E* state,const Key* keys,E* output){
    std::fill(output,output+b,E{});
    for(unsigned a=0;a<2*b;++a){
        E value{};for(unsigned j=0;j<b;++j)value^=scale8(state[j],table.g[j][a]);
        value=scale32(value,wide::recoverScalar(keys[a]));
        for(unsigned i=0;i<b;++i)output[i]^=scale8(value,table.h[b+i][a]);
    }
}
template<unsigned b> void momentOracle(const E* input,E* result){
    std::fill(result,result+b,E{});
    for(unsigned j=0;j<b;++j){const auto beta=tr::TransformCantor<b>::evaluationPoint(b+j);
        std::uint8_t coefficient=1;require(beta!=0,"zero moment point");
        for(unsigned i=0;i<b;++i){result[i]^=scale8(input[j],coefficient);coefficient=gf::byteOracle(coefficient,beta);}}
}
#if SPIN_OR_ADDITIVE
template<unsigned b> std::uint32_t shiftDelta(unsigned bit){
    return tr::AdditivePersistent32<b>::shiftDelta(bit);
}

template<unsigned b> void checkDirections(){
    std::array<std::uint32_t,32> rows{};
    for(unsigned bit=0;bit<32;++bit){
        const auto value=shiftDelta<b>(bit);
        if(bit<8)require(value==tr::TransformCantor<b>::evaluationPoint(1U<<bit),"changed base Cantor chain");
#if SPIN_OR_FULL_CANTOR
        if(bit)require((gf::towerMultiplyOracle(value,value)^value)==shiftDelta<b>(bit-1),"full Cantor oracle relation mismatch");
#endif
        rows[bit]=value;
    }
    unsigned rank=0;
    for(unsigned bit=0;bit<32;++bit){
        unsigned pivot=rank;while(pivot<32&&!(rows[pivot]&(1U<<bit)))++pivot;
        if(pivot==32)continue;
        std::swap(rows[pivot],rows[rank]);
        for(unsigned row=rank+1;row<32;++row)if(rows[row]&(1U<<bit))rows[row]^=rows[rank];
        ++rank;
    }
    require(rank==32,"additive directions lost binary independence");
}
template<unsigned b> std::array<std::uint32_t,b> shiftWeights(unsigned bit){
    const auto delta=shiftDelta<b>(bit);auto vanish=delta;
    for(unsigned i=0;i<std::countr_zero(b);++i)vanish=gf::towerMultiplyOracle(vanish,vanish)^vanish;
    require(vanish!=0,"translation lies in base subspace");
    std::array<std::uint32_t,b> weights;
    for(unsigned j=0;j<b;++j)weights[j]=gf::towerMultiplyOracle(vanish,
        power32(delta^tr::TransformCantor<b>::evaluationPoint(j),0xfffffffeULL));
    return weights;
}
template<unsigned b> void shiftOracle(const std::array<std::uint32_t,b>& weights,const E* input,E* output){
    for(unsigned i=0;i<b;++i){E value{};for(unsigned j=0;j<b;++j)value^=scale32(input[j],weights[i^j]);output[i]=value;}
}
template<unsigned b> void cosetMoment(std::uint32_t offset,const E* input,E* output){
    std::fill(output,output+b,E{});
    for(unsigned j=0;j<b;++j){const auto point=offset^tr::TransformCantor<b>::evaluationPoint(j);std::uint32_t a=1;
        for(unsigned i=0;i<b;++i){output[i]^=scale32(input[j],a);a=gf::towerMultiplyOracle(a,point);}}
}
template<unsigned B,unsigned b> void checkAdditiveStream(std::size_t k,const Oracle<B>& outerTable,
    std::mt19937_64& rng,bool direct){
    Plan<B,b> p(k);
    Aligned<E> message(p.K),planes(p.N),result(p.N),routed(p.N),reference(p.N);
    for(std::size_t i=0;i<p.K;++i)message[i]=random(rng);outer(p,message.data(),planes.data());
    if(direct)for(unsigned r=0;r<B;++r)
        for(std::size_t block=0;block<p.blocks;++block){
            E value{};for(unsigned j=0;j<B/2;++j)value^=scale8(message[block*(B/2)+j],outerTable.g[j][r]);
            require(value==planes[r*p.blocks+block],"RS outer oracle mismatch");}
    routeOnly(p,planes.data(),routed.data());checkRoute(p,planes.data(),routed.data());
    inner<B,b,true>(p,planes.data(),result.data());
    inner<B,b,false>(p,routed.data(),reference.data());
    require(std::memcmp(result.data(),reference.data(),p.N*sizeof(E))==0,"additive fused route mismatch");
    // Only the smallest fixture repeats the quadratic scalar field oracles.
    // Larger fixtures exercise nontrivial packet labels and multiple batches
    // per region without turning this correctness gate into a large job.
    if(!direct)return;
    alignas(64) E state[b]{},live[b],weighted[b],y[b],global[b]{},now[b],expected[b];std::uint32_t offset=0;
    for(std::size_t t=0;t<p.batches;++t){
        const unsigned bit=std::countr_zero(b)+std::countr_zero(t+1);
        const auto weights=shiftWeights<b>(bit);offset^=shiftDelta<b>(bit);
        shiftOracle<b>(weights,state,live);
#if SPIN_OR_SINGLE_DIAGONAL
        for(unsigned j=0;j<b;++j){
            expected[j]=scale32(live[j],randomScalar(p.keys[b*t+j]));
            y[j]=routed[b*t+j]^expected[j];
            require(y[j]==result[b*t+j],"experimental single-diagonal stream oracle mismatch");
            state[j]=live[j]^y[j];}
#else
        for(unsigned j=0;j<b;++j)weighted[j]=scale32(state[j],randomScalar(p.keys[2*b*t+j]));
        shiftOracle<b>(weights,weighted,expected);
        for(unsigned j=0;j<b;++j){
            y[j]=routed[b*t+j]^expected[j]^scale32(live[j],randomScalar(p.keys[2*b*t+b+j]));
            require(y[j]==result[b*t+j],"additive full stream oracle mismatch");state[j]=live[j]^y[j];}
#endif
        cosetMoment<b>(offset,y,now);for(unsigned j=0;j<b;++j)global[j]^=now[j];
        cosetMoment<b>(offset,state,now);
        require(std::memcmp(global,now,sizeof(now))==0,"additive global MDS unrolled oracle mismatch");
    }
}
template<unsigned B,unsigned b> NI void checkAdditive(){
    checkDirections<b>();
    std::mt19937_64 rng(227);Oracle<B> outerTable;
    alignas(64) E input[b],got[b],expected[b],back[b];
    for(unsigned bit=std::countr_zero(b);bit<=21;++bit){
        const auto shift=tr::AdditivePersistent32<b>::prepare(bit);const auto weights=shiftWeights<b>(bit);
        for(unsigned fixture=0;fixture<b+2;++fixture){
            for(unsigned j=0;j<b;++j)input[j]=(fixture<b?(fixture==j?random(rng):E{}):random(rng));
            tr::AdditivePersistent32<b>::apply(shift,input,got);
            if(fixture<b)for(unsigned j=0;j<b;++j)expected[j]=scale32(input[fixture],weights[j^fixture]);
            else shiftOracle<b>(weights,input,expected);
            require(std::memcmp(got,expected,sizeof(got))==0,"shift Cauchy oracle mismatch");
            tr::AdditivePersistent32<b>::apply(shift,got,back);
            require(std::memcmp(back,input,sizeof(back))==0,"shift involution mismatch");
        }
    }
    checkAdditiveStream<B,b>(256,outerTable,rng,true);
    checkAdditiveStream<B,b>(1024,outerTable,rng,false);
    checkAdditiveStream<B,b>(std::size_t(B)*b,outerTable,rng,false);
    std::cout<<"PASS additive B="<<B<<" b="<<b
        <<" packet="<<packetWidth<<" single_diagonal="<<SPIN_OR_SINGLE_DIAGONAL
        <<" profile="<<profileLabel<<" Cauchy/FFT shifts through bit21, involution, RS, packet bijections/order,"
        <<" fused/materialized routing at K256/K1024/K"<<B*b<<", full stream, global MDS identity\n";
}
#endif
template<unsigned B,unsigned b> NI void check(){
    Plan<B,b> p(256);Oracle<2*b> tables;Oracle<B> outerTable;std::mt19937_64 rng(227);
    alignas(64) E state[b],out[b],expected[b],y[b],saved[b];
    for(unsigned fixture=0;fixture<b+4;++fixture){
        for(unsigned j=0;j<b;++j)state[j]=(fixture<b?(fixture==j?E{rng(),rng()}:E{}):random(rng));
        tr::Persistent32<b>::moments(state,out);momentOracle<b>(state,expected);
        require(std::memcmp(out,expected,sizeof(out))==0,"moment map oracle mismatch");
        tr::Persistent32<b>::feedback(state,p.keys.data(),out);feedbackOracle<b>(tables,state,p.keys.data(),expected);
        require(std::memcmp(out,expected,sizeof(out))==0,"spectral feedback oracle mismatch");
    }
    for(auto& x:y)x=random(rng);std::memcpy(saved,state,sizeof(state));
    momentOracle<b>(y,expected);for(unsigned j=0;j<b;++j)expected[j]^=scale32(saved[j],p.powers[j]);
    tr::Persistent32<b>::update(state,y,p.fixed.data());
    require(std::memcmp(state,expected,sizeof(state))==0,"persistent update mismatch");
    // Independent direct check of the globally Vandermonde unrolled state.
    std::array<std::array<E,b>,4> history{};std::memset(state,0,sizeof(state));
    for(unsigned t=0;t<4;++t){
        for(auto& x:history[t])x=random(rng);
        tr::Persistent32<b>::update(state,history[t].data(),p.fixed.data());
        for(unsigned i=0;i<b;++i){
            E value{};
            for(unsigned h=0;h<=t;++h){
                const auto temporal=power32(p.powers[i],t-h);
                for(unsigned j=0;j<b;++j){const auto beta=tr::TransformCantor<b>::evaluationPoint(b+j);
                    value^=scale32(history[h][j],gf::towerMultiplyOracle(temporal,gf::bytePower(beta,i)));}
            }
            require(state[i]==value,"global MDS unrolled state mismatch");
        }
    }
    Aligned<E> message(p.K),planes(p.N),result(p.N),routed(p.N),contiguous(p.N);
    for(std::size_t i=0;i<p.K;++i)message[i]=random(rng);
    outer(p,message.data(),planes.data());
    for(unsigned r=0;r<B;++r){std::vector<bool> seen(p.blocks);
        for(std::size_t block=0;block<p.blocks;++block){
            E value{};for(unsigned j=0;j<B/2;++j)value^=scale8(message[block*(B/2)+j],outerTable.g[j][r]);
            require(value==planes[r*p.blocks+block],"RS outer oracle mismatch");
            const auto label=p.route[r*p.blocks+block];require(label<p.blocks&&!seen[label],"route not bijective");seen[label]=true;
        }}
    routeOnly(p,planes.data(),routed.data());inner<B,b,true>(p,planes.data(),result.data());
    inner<B,b,false>(p,routed.data(),contiguous.data());
    require(std::memcmp(result.data(),contiguous.data(),p.N*sizeof(E))==0,"fused and materialized routes disagree");
    std::memset(state,0,sizeof(state));
    for(std::size_t t=0;t<p.batches;++t){
        feedbackOracle<b>(tables,state,p.keys.data()+2*b*t,expected);
        for(unsigned j=0;j<b;++j){y[j]=routed[b*t+j]^expected[j];require(y[j]==result[b*t+j],"stream oracle mismatch");}
        momentOracle<b>(y,out);for(unsigned j=0;j<b;++j)state[j]=scale32(state[j],p.powers[j])^out[j];
    }
    std::cout<<"PASS B="<<B<<" b="<<b<<" alpha="<<p.alpha
        <<" primitive order, RS polynomial oracle, regional bijections, moments, uniform spectral map, update, full stream/free terminal\n";
}
std::uint64_t checksum(const E* data,std::size_t n){
    std::uint64_t sum=0x243f6a8885a308d3ULL;
    for(std::size_t i=0;i<n;++i)sum=(sum^(data[i].lo+data[i].hi))*0x9e3779b97f4a7c15ULL;return sum;
}
template<class F> double measure(unsigned repetitions,F&& fn){
    for(unsigned i=0;i<3;++i)fn();std::vector<double> times;times.reserve(repetitions);
    for(unsigned i=0;i<repetitions;++i){const auto start=Clock::now();fn();
        times.push_back(std::chrono::duration<double,std::milli>(Clock::now()-start).count());}
    std::sort(times.begin(),times.end());return times[times.size()/2];
}
template<unsigned B,unsigned b> void bench(unsigned logK,unsigned reps,const std::string& phase){
    require(logK>=12&&logK<=20&&reps>=3&&reps<=101,"benchmark cap logK12..20,repetitions3..101");
    require(phase=="all"||phase=="outer"||phase=="route"||phase=="inner"||phase=="full","invalid phase");
    const auto start=Clock::now();Plan<B,b> p(std::size_t(1)<<logK);
    std::cout<<std::fixed<<std::setprecision(6)<<"K="<<p.K<<" B="<<B<<" b="<<b<<" field_bits=32 alpha="<<p.alpha
        <<" setup_ms="<<std::chrono::duration<double,std::milli>(Clock::now()-start).count()
        <<" key_MiB="<<p.keys.count*sizeof(RandomKey)/1048576.0<<" repetitions="<<reps
        <<" additive="<<SPIN_OR_ADDITIVE<<" full_cantor="<<SPIN_OR_FULL_CANTOR<<" packet="<<packetWidth
        <<" raw_keys="<<SPIN_OR_RAW_KEYS
        <<" outer_tile="<<SPIN_OR_OUTER_TILE<<" outer_stream="<<SPIN_OR_OUTER_STREAM
        <<" single_diagonal="<<SPIN_OR_SINGLE_DIAGONAL<<" profile="<<profileLabel
        <<" route_MiB="<<p.route.count*sizeof(std::uint32_t)/1048576.0<<"\n"<<std::flush;
    Aligned<E> input(p.K),workspace(p.N),output(p.N);std::mt19937_64 rng(931);
    for(std::size_t i=0;i<p.K;++i)input[i]=random(rng);outer(p,input.data(),workspace.data());
    const auto show=[&](const char* name,double ms){std::cout<<name<<"_ms="<<ms<<" checksum="<<checksum(output.data(),p.N)<<"\n"<<std::flush;};
    if(phase=="all"||phase=="outer")show("outer",measure(reps,[&]{outer(p,input.data(),output.data());}));
    if(phase=="all"||phase=="route")show("route",measure(reps,[&]{routeOnly(p,workspace.data(),output.data());}));
    if(phase=="all"||phase=="inner"){
        show("inner_contiguous",measure(reps,[&]{inner<B,b,false>(p,workspace.data(),output.data());}));
        show("inner_routed",measure(reps,[&]{inner<B,b,true>(p,workspace.data(),output.data());}));}
    if(phase=="all"||phase=="full")show("full",measure(reps,[&]{full(p,input.data(),workspace.data(),output.data());}));
}
}
int main(int argc,char** argv){
    try{
        if(argc==2&&std::string(argv[1])=="--check"){
            pn::checkRandomKeys();
#if SPIN_OR_ADDITIVE
#if SPIN_OR_PACKET == 2
            pn::checkAdditive<128,64>();
#elif SPIN_OR_PACKET == 4
            pn::checkAdditive<128,128>();pn::checkAdditive<64,128>();
#else
            pn::checkAdditive<64,64>();pn::checkAdditive<128,32>();
#endif
#else
            pn::check<64,64>();pn::check<128,32>();
#endif
            return 0;}
        if(argc==7&&std::string(argv[1])=="--bench"){
            const auto B=std::stoul(argv[2]),b=std::stoul(argv[3]),logK=std::stoul(argv[4]),reps=std::stoul(argv[5]);
#if SPIN_OR_PACKET == 2
            if(B==128&&b==64)pn::bench<128,64>(logK,reps,argv[6]);
            else throw std::runtime_error(std::string("packet2: only ")+pn::profileLabel+" B128/b64 profile");
#elif SPIN_OR_PACKET == 4
            if(B==128&&b==128)pn::bench<128,128>(logK,reps,argv[6]);
            else if(B==64&&b==128)pn::bench<64,128>(logK,reps,argv[6]);
            else throw std::runtime_error(std::string("packet4: only ")+pn::profileLabel+" B128/b128 and B64/b128 profiles");
#else
            if(B==64&&b==64)pn::bench<64,64>(logK,reps,argv[6]);
            else if(B==128&&b==32)pn::bench<128,32>(logK,reps,argv[6]);
            else throw std::runtime_error(std::string("only ")+pn::profileLabel+" B64/b64 and B128/b32 profiles");
#endif
            return 0;
        }
        throw std::runtime_error("usage: --check | --bench B batch logK repetitions all|outer|route|inner|full");
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
