#pragma once
#include "Gamma30Encoder.h"
#include "../BenchPayload.h"
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

// Shared correctness-only support, extracted from the frozen gamma30 dense
// oracle driver without changing its arithmetic. Not used by hot kernels or
// timing drivers. Future variants can reuse this without including a .cpp.
namespace spin::research::char2::field128_lift::gamma30_test {
namespace ref=spin::research::char2::field128_lift::gamma30;
namespace field=ref::field;
namespace c2=spin::research::char2;
using Element=ref::Element;
using Block=ref::Block;
using State30=std::array<Element,30>;
using Packet64=std::array<Element,64>;

static void require(bool condition,const char* label){
    if(!condition)throw std::runtime_error(label);
}
static Element element(Block x){return {x.lo,x.hi};}
static Block block(Element x){return {x.lo,x.hi};}
static Element randomElement(std::mt19937_64& random){return {random(),random()};}
static constexpr Element one{1,0},feedbackGamma{0x100,0};
static constexpr Block guard{0x2e71d6935a90bf24ULL,0xc031674f27e18ac5ULL};

static void same(const Block* actual,const Block* expected,std::size_t n,const char* label){
    for(std::size_t i=0;i!=n;++i)if(actual[i]!=expected[i])
        throw std::runtime_error(std::string(label)+": record "+std::to_string(i)+" differs");
}
template<std::size_t N>
static void same(const std::array<Element,N>& actual,const std::array<Element,N>& expected,const char* label){
    for(std::size_t i=0;i!=N;++i)if(actual[i]!=expected[i])
        throw std::runtime_error(std::string(label)+": coordinate "+std::to_string(i)+" differs");
}
static void produce(Block* output,std::size_t n,std::uint64_t seed){
    std::mt19937_64 random(seed);
    for(std::size_t i=0;i!=n;++i)output[i]=block(randomElement(random));
}

struct Oracle {
    std::array<std::array<std::uint8_t,256>,256> mul{};
    std::array<std::uint8_t,256> inverse{};
    std::array<std::array<std::uint8_t,32>,64> outer{};
    std::array<std::array<std::uint8_t,30>,64> a{};
    std::array<std::array<std::uint8_t,30>,30> h{};

    // Independent bit-serial AES multiplication, not field::byteMultiply.
    static std::uint8_t aes(unsigned left,unsigned right){
        unsigned result=0;
        for(unsigned bit=0;bit!=8;++bit){
            if(right&1U)result^=left;
            right>>=1;left=(left<<1)^((left&128U)?0x11bU:0U);
        }
        return static_cast<std::uint8_t>(result);
    }
    static std::uint8_t byte(Element x,unsigned i){
        return static_cast<std::uint8_t>((i<8?x.lo:x.hi)>>(8*(i%8)));
    }
    static Element pack(const std::array<std::uint8_t,16>& bytes){
        Element result{};
        for(unsigned i=0;i!=8;++i){
            result.lo|=std::uint64_t(bytes[i])<<(8*i);
            result.hi|=std::uint64_t(bytes[8+i])<<(8*i);
        }
        return result;
    }
    Oracle(){
        for(unsigned i=0;i!=256;++i)for(unsigned j=0;j!=256;++j){
            mul[i][j]=aes(i,j);
            if(mul[i][j]==1)inverse[i]=static_cast<std::uint8_t>(j);
        }
        for(unsigned i=0;i!=64;++i){
            for(unsigned j=0;j!=32;++j)outer[i][j]=novelBasis(j,i);
            for(unsigned j=0;j!=30;++j)a[i][j]=lagrange(128,30,j,i);
        }
        for(unsigned i=0;i!=30;++i)for(unsigned j=0;j!=30;++j)
            h[i][j]=lagrange(0,30,j,30+i);
    }
    std::uint8_t lagrange(unsigned first,unsigned count,unsigned column,unsigned point) const{
        std::uint8_t numerator=1,denominator=1;
        for(unsigned j=0;j!=count;++j)if(j!=column){
            numerator=mul[numerator][point^(first+j)];
            denominator=mul[denominator][(first+column)^(first+j)];
        }
        require(denominator!=0,"oracle Lagrange nodes must be distinct");
        return mul[numerator][inverse[denominator]];
    }
    std::uint8_t novelBasis(unsigned degree,unsigned point) const{
        std::uint8_t result=1;
        for(unsigned bit=0;bit!=5;++bit)if((degree>>bit)&1U){
            std::uint8_t value=1,normalizer=1;
            for(unsigned j=0;j!=(1U<<bit);++j){
                value=mul[value][point^j];normalizer=mul[normalizer][(1U<<bit)^j];
            }
            result=mul[result][mul[value][inverse[normalizer]]];
        }
        return result;
    }
    Element scale(Element x,std::uint8_t coefficient) const{
        std::array<std::uint8_t,16> bytes{};
        for(unsigned j=0;j!=16;++j)bytes[j]=mul[coefficient][byte(x,j)];
        return pack(bytes);
    }
    Element product(Element left,Element right) const{
        std::array<std::uint8_t,31> coefficients{};
        for(unsigned i=0;i!=16;++i)for(unsigned j=0;j!=16;++j)
            coefficients[i+j]^=mul[byte(left,i)][byte(right,j)];
        // Independent generic long division by the explicit tower modulus.
        // Do not call either the SIMD reduction or the multiply-by-Y helper.
        constexpr std::array<std::uint8_t,17> modulus{
            0x32,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,1};
        for(unsigned degree=31;degree-- >16;){
            const auto leading=coefficients[degree];
            for(unsigned j=0;j!=17;++j)coefficients[degree-16+j]^=mul[leading][modulus[j]];
        }
        std::array<std::uint8_t,16> bytes{};
        std::copy_n(coefficients.begin(),16,bytes.begin());return pack(bytes);
    }
    Packet64 emit(const State30& state) const{
        Packet64 result{};
        for(unsigned i=0;i!=64;++i)for(unsigned j=0;j!=30;++j)result[i]^=scale(state[j],a[i][j]);
        return result;
    }
    State30 feedback(const Packet64& input) const{
        State30 result{};
        for(unsigned j=0;j!=30;++j)for(unsigned i=0;i!=64;++i)result[j]^=scale(input[i],a[i][j]);
        return result;
    }
    State30 refresh(const State30& state) const{
        State30 result{};
        for(unsigned i=0;i!=30;++i)for(unsigned j=0;j!=30;++j)result[i]^=scale(state[j],h[i][j]);
        return result;
    }
    State30 finish(const State30& state,const State30& syndrome,const ref::Refresh& diagonal) const{
        State30 before{};
        for(unsigned j=0;j!=30;++j)
            before[j]=product(state[j]^product(feedbackGamma,syndrome[j]),diagonal.before[j]);
        auto result=refresh(before);
        for(unsigned j=0;j!=30;++j)result[j]=product(result[j],diagonal.after[j]);
        return result;
    }
    Packet64 step(const Packet64& input,State30& state,const ref::Refresh& diagonal) const{
        auto output=emit(state);
        for(unsigned i=0;i!=64;++i)output[i]^=input[i];
        // Consume C*x, deliberately independent of the hot path's C*y trick.
        state=finish(state,feedback(input),diagonal);
        return output;
    }
    Element determinant(const std::array<Element,7>& t) const{
        // Characteristic two: determinant equals the unsigned permutation sum.
        std::array<unsigned,4> permutation{0,1,2,3};Element result{};
        do{
            auto term=one;
            for(unsigned i=0;i!=4;++i)term=product(term,t[3+i-permutation[i]]);
            result^=term;
        }while(std::next_permutation(permutation.begin(),permutation.end()));
        return result;
    }
    void outerReference(const ref::Plan& plan,const Block* input,Block* planes) const{
        const auto& p=plan.outerPlan;
        for(std::size_t group=0;group!=p.rows;++group)for(unsigned c=0;c!=64;++c){
            std::array<Element,4> column{};
            for(unsigned r=0;r!=4;++r)for(unsigned d=0;d!=32;++d)
                column[r]^=scale(element(input[(group/4)*512+4*(r*32+d)+group%4]),outer[c][d]);
            const auto& t=p.toeplitz[c*p.rows+group];
            for(unsigned r=0;r!=4;++r){
                Element value{};
                for(unsigned j=0;j!=4;++j)value^=product(column[j],t[3+r-j]);
                planes[(4*c+r)*p.stride+group]=block(value);
            }
        }
        for(unsigned region=0;region!=256;++region){
            auto* plane=planes+region*p.stride;
            plane[p.rows]=plane[0];std::fill_n(plane+p.rows+1,3,Block{});
        }
    }
    void innerReference(const ref::Plan& plan,const Block* planes,Block* output) const{
        const auto& p=plan.outerPlan;
        State30 state{}; // One state for all packets and all 256 regions.
        for(std::size_t first=0;first!=p.n;first+=64){
            Packet64 input{};
            for(unsigned i=0;i!=64;++i)
                input[i]=element(planes[(first/p.rows)*p.stride+p.route[first+i]]);
            const auto y=step(input,state,plan.refresh[first/64]);
            for(unsigned i=0;i!=64;++i)output[first+4*(i%16)+i/16]=block(y[i]);
        }
    }
};

static ref::State packedState(const State30& state,bool poisonPadding=false){
    alignas(64) std::array<Element,32> storage{};
    std::copy(state.begin(),state.end(),storage.begin());
    if(poisonPadding){storage[30]={~0ULL,0x701239ca45b8ddedULL};storage[31]={0x528921a132456789ULL,~0ULL};}
    ref::State result;
    for(unsigned i=0;i!=8;++i)result.values[i]=_mm512_load_si512(storage.data()+4*i);
    return result;
}
static State30 naturalState(const ref::State& state,bool requirePadding=true){
    alignas(64) std::array<Element,32> storage{};
    for(unsigned i=0;i!=8;++i)_mm512_store_si512(storage.data()+4*i,state.values[i]);
    if(requirePadding)require(storage[30]==Element{}&&storage[31]==Element{},"state padding must be zero");
    State30 result{};std::copy_n(storage.begin(),30,result.begin());return result;
}
static ref::Packet packedPacket(const Packet64& input){
    ref::Packet result;
    for(unsigned i=0;i!=16;++i)result.values[i]=_mm512_loadu_si512(input.data()+4*i);
    return result;
}
static Packet64 naturalPacket(const ref::Packet& packet){
    Packet64 result{};
    for(unsigned i=0;i!=16;++i)_mm512_storeu_si512(result.data()+4*i,packet.values[i]);
    return result;
}


} // namespace spin::research::char2::field128_lift::gamma30_test
