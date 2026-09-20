#pragma once
#include "Spin.h"
#include <array>
#include <bit>
#include <concepts>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include "GenericCircuits.h"

namespace bare_spin {
// Each value is an element of a characteristic-two additive group. XOR must
// operate componentwise with no hidden allocation. bool's vector specialization
// is deliberately excluded; use uint8_t for one byte (or a single bit).
template<class E> concept XorElement = std::is_trivially_copyable_v<E> &&
    std::default_initializable<E> && !std::same_as<E,bool> &&
    requires(const E& a,const E& b) {{a^b}->std::convertible_to<E>;};

class GenericTranspose {
public:
    template<XorElement E> struct Workspace {
        std::vector<E> values;
        explicit Workspace(const GenericTranspose& code):values(code.codeBlocks()) {}
        std::size_t bytes() const noexcept {return values.capacity()*sizeof(E);}
    };
    std::size_t messageBlocks() const noexcept {return mK;}
    std::size_t codeBlocks() const noexcept {return 2*mK;}
    std::size_t setupBytes() const noexcept {return sizeof(u32)*(mRoute.capacity()+mMasks.capacity());}
    // Separate buffers may not overlap. Use encodeInplace for an N-element
    // buffer whose first K elements are replaced and whose suffix is preserved.
    template<XorElement E> void encode(std::span<const E> in,std::span<E> out,Workspace<E>& w) const {
        validate(in.data(),in.size(),out.data(),out.size(),w);
        if(overlap(in.data(),in.size(),out.data(),out.size()))
            throw std::invalid_argument("generic input and output overlap; use encodeInplace");
        encodeUnchecked(in.data(),out.data(),w);
    }
    template<XorElement E> void encodeInplace(std::span<E> buffer,Workspace<E>& w) const {
        validate(buffer.data(),buffer.size(),buffer.data(),messageBlocks(),w);
        encodeUnchecked(buffer.data(),buffer.data(),w);
    }
    // No allocation or type erasure. Preconditions: correct geometry, workspace
    // disjoint from input/output, and disjoint buffers or exactly in-place.
    template<XorElement E> void encodeUnchecked(const E* in,E* out,Workspace<E>& w) const {
        switch(mConfig) {
            case Configuration::T128S19:run<generic_detail::Map128,1>(in,out,w);return;
            case Configuration::T64S12:run<generic_detail::Map64,1>(in,out,w);return;
            case Configuration::T64S12R2:run<generic_detail::Map64,2>(in,out,w);return;
            default:throw std::logic_error("unsupported generic map");
        }
    }
private:
    friend class Spin;
    Configuration mConfig;
    std::size_t mK;
    std::vector<u32> mRoute,mMasks;
    GenericTranspose(Configuration config,std::size_t k,std::vector<u32> route,std::vector<u32> masks)
        :mConfig(config),mK(k),mRoute(std::move(route)),mMasks(std::move(masks)) {}
    template<class E> static bool overlap(const E* a,std::size_t na,const E* b,std::size_t nb) {
        const auto x=reinterpret_cast<std::uintptr_t>(a),y=reinterpret_cast<std::uintptr_t>(b);
        return x<=y?y-x<na*sizeof(E):x-y<nb*sizeof(E);
    }
    template<XorElement E> void validate(const E* in,std::size_t ni,E* out,std::size_t no,const Workspace<E>& w) const {
        if(!in || !out || ni!=codeBlocks() || no!=messageBlocks() || w.values.size()!=codeBlocks())
            throw std::invalid_argument("generic input, output, or workspace geometry mismatch");
        if(overlap(in,ni,w.values.data(),w.values.size()) || overlap(out,no,w.values.data(),w.values.size()))
            throw std::invalid_argument("generic workspace overlaps input or output");
    }
    template<class E> static inline void imt(E* state,u32 u,u32 v) {
        E dot=generic_detail::vx(state[0],state[0]);
        while(u) {const auto j=std::countr_zero(u);u&=u-1;dot=generic_detail::vx(dot,state[j]);}
        while(v) {const auto j=std::countr_zero(v);v&=v-1;state[j]=generic_detail::vx(state[j],dot);}
    }
    template<unsigned T,unsigned D,unsigned Base=0,class E> static inline void zetaStage(E* x) {
        // Same compile-time pruning and paired loop shape as the SIMD kernel.
        if constexpr(std::popcount(Base/(2*D))<=2) {
            if constexpr(D==1) x[Base]=generic_detail::vx(x[Base],x[Base+1]);
            else for(unsigned j=0;j<D;j+=2) {
                x[Base+j]=generic_detail::vx(x[Base+j],x[Base+j+D]);
                x[Base+j+1]=generic_detail::vx(x[Base+j+1],x[Base+j+D+1]);
            }
        }
        if constexpr(Base+2*D<T) zetaStage<T,D,Base+2*D>(x);
    }
    template<unsigned T,unsigned D=T/2,class E> static inline void zeta(E* x) {
        zetaStage<T,D>(x);
        if constexpr(D>1) zeta<T,D/2>(x);
    }
    template<class Map,unsigned Rounds,XorElement E> void run(const E* in,E* out,Workspace<E>& w) const {
        std::array<E,Map::S> state,syndrome;
        std::array<E,Map::T> raw;
        const auto epochs=codeBlocks()/Map::T;
        auto emit=[&](std::size_t i,const E& v) {w.values[mRoute[i]]=v;};
        for(std::size_t epoch=epochs;epoch-->0;) {
            const auto base=epoch*Map::T;
            if(epoch+1==epochs) {
                for(unsigned j=Map::T;j-->0;) {raw[j]=in[base+j];emit(base+j,raw[j]);}
            } else Map::emitShared(in+base,raw.data(),state.data(),base,emit);
            if(epoch==0) break;
            zeta<Map::T>(raw.data());Map::finish(raw.data(),syndrome.data());
            if(epoch+1==epochs) state=syndrome;
            else {
                const auto* masks=mMasks.data()+2*Rounds*epoch;
                if constexpr(Rounds==2) imt(state.data(),masks[2],masks[3]);
                imt(state.data(),masks[0],masks[1]);
                for(unsigned j=0;j<Map::S;++j) state[j]=generic_detail::vx(state[j],syndrome[j]);
            }
        }
        // All input was consumed above: exactly in-place output is now safe.
        for(std::size_t row=0;row<mK/128;++row)
            generic_detail::bch(w.values.data()+256*row,out+128*row);
    }
};
}
