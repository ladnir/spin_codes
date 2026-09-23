#pragma once
#include "Code.h"
#include <array>
#include <bit>
#include <concepts>
#include <limits>
#include <utility>
#include <vector>
#include "detail/GenericCircuits.h"
#include "detail/BankState.h"

namespace spin {
// Componentwise XOR, with no hidden allocation. bool/vector<bool> is excluded.
template<class E> concept ValueElement = std::copyable<E> &&
    std::default_initializable<E> && !std::same_as<E,bool>;
template<class E> concept XorElement = ValueElement<E> &&
    requires(const E& a,const E& b) {{a^b}->std::convertible_to<E>;};
struct Xor {
    template<XorElement E> E operator()(const E& a,const E& b) const {return E(a^b);}
};

// Owned, portable XOR-element fallback for the same realized binary map.
// Code creates a fixed snapshot; PreparedEncoder creates a live view. Neither
// changes the optimized Code plan. Live-view refresh requires exclusive access.
class GenericTranspose {
    using u32=std::uint32_t;
    struct Data {
        Parameters config;
        std::size_t k;
        std::vector<u32> route,masks;
        std::shared_ptr<const detail::BankState> bank;
    };
public:
    template<ValueElement E> class Workspace {
    public:
        Workspace(Workspace&&) noexcept=default;
        Workspace& operator=(Workspace&&) noexcept=default;
        Workspace(const Workspace&)=delete;
        Workspace& operator=(const Workspace&)=delete;
        std::size_t bytes() const noexcept {return values_.capacity()*sizeof(E);}
    private:
        friend class GenericTranspose;
        std::shared_ptr<const Data> owner_;
        std::vector<E> values_;
        explicit Workspace(std::shared_ptr<const Data> d):owner_(std::move(d)),values_(2*owner_->k) {}
    };
    std::size_t message_size() const noexcept {return data_->k;}
    std::size_t code_size() const noexcept {return 2*data_->k;}
    std::size_t setup_bytes() const noexcept {return data_->bank?data_->bank->bytes():4*(data_->route.capacity()+data_->masks.capacity());}
    template<ValueElement E> Workspace<E> make_workspace() const {
        if(code_size()>std::numeric_limits<std::size_t>::max()/sizeof(E))
            throw std::length_error("SPIN generic workspace size overflow");
        return Workspace<E>(data_);
    }
    // Op implements characteristic-two addition. It is statically dispatched;
    // custom coefficient contexts need not provide operator^ on their values.
    template<ValueElement E,class Op=Xor> void transpose(std::span<const E> in,std::span<E> out,Workspace<E>& w,const Op& op={}) const {
        check(in.data(),in.size(),out.data(),out.size(),w);
        if(overlap(in.data(),in.size_bytes(),out.data(),out.size_bytes()))
            throw std::invalid_argument("SPIN generic buffers overlap");
        dispatch(in.data(),out.data(),w,op);
    }
    template<ValueElement E,class Op=Xor> void transpose_inplace(std::span<E> buf,Workspace<E>& w,const Op& op={}) const {
        check(buf.data(),buf.size(),buf.data(),message_size(),w);
        dispatch(buf.data(),buf.data(),w,op);
    }
private:
    friend class Code;
    friend class PreparedEncoder;
    std::shared_ptr<const Data> data_;
    GenericTranspose(Parameters c,std::size_t k,std::vector<u32> route,std::vector<u32> masks)
        :data_(std::make_shared<Data>(Data{c,k,std::move(route),std::move(masks),{}})) {}
    static bool overlap(const void* a,std::size_t na,const void* b,std::size_t nb) noexcept {
        const auto x=reinterpret_cast<std::uintptr_t>(a),y=reinterpret_cast<std::uintptr_t>(b);
        return x<=y?y-x<na:x-y<nb;
    }
    template<ValueElement E> void check(const E* in,std::size_t ni,E* out,std::size_t no,const Workspace<E>& w) const {
        if(!in || !out || ni!=code_size() || no!=message_size() || w.owner_!=data_ || w.values_.size()!=code_size())
            throw std::invalid_argument("SPIN generic geometry or workspace mismatch");
    }
    template<ValueElement E,class Op> void dispatch(const E* in,E* out,Workspace<E>& w,const Op& op) const {
        if(data_->bank) {
            const auto& b=*data_->bank;
            // Epochs visit complete regions in reverse order. Avoid division
            // per element, including at non-power-of-two natural lengths.
            struct BankRoute {
                const detail::BankState& bank;
                mutable std::size_t region=255,base;
                void begin_epoch(std::size_t i) const {
                    if(i<base){--region;base-=bank.rows;}
                }
                u32 operator()(std::size_t i) const {return bank.destination(region,i-base);}
            } route{b,255,255*b.rows};
            dispatchRoute(in,out,w,op,route,b.masks.data());
        } else {
            auto route=[&](std::size_t i){return data_->route[i];};
            dispatchRoute(in,out,w,op,route,data_->masks.data());
        }
    }
    template<ValueElement E,class Op,class Route> void dispatchRoute(const E* in,E* out,Workspace<E>& w,const Op& op,const Route& route,const u32* masks) const {
        switch(data_->config) {
        case Parameters::T128S19:run<detail::generic::Map128,1>(in,out,w,op,route,masks);return;
        case Parameters::T64S12:run<detail::generic::Map64,1>(in,out,w,op,route,masks);return;
        case Parameters::T64S12R2:run<detail::generic::Map64,2>(in,out,w,op,route,masks);return;
        }
        throw std::logic_error("SPIN invalid generic configuration");
    }
#include "detail/GenericMethods.h"
};
}
