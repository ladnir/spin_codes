#pragma once
#include <array>
#include <bit>
#include <cstdint>
#include <stdexcept>

namespace spin::research::char2::field128_lift {

// Native bit j is the coefficient of X^j modulo X^128+X^7+X^2+X+1.
// This is one field element, not sixteen independently scaled byte fields.
struct Scalar {
    std::uint64_t lo=0,hi=0;
    constexpr bool operator==(const Scalar&) const noexcept = default;
};
static_assert(sizeof(Scalar)==16);
constexpr Scalar operator^(Scalar a,Scalar b) noexcept {return {a.lo^b.lo,a.hi^b.hi};}
constexpr Scalar& operator^=(Scalar& a,Scalar b) noexcept {a=a^b;return a;}
constexpr bool nonzero(Scalar a) noexcept {return (a.lo|a.hi)!=0;}
constexpr Scalar binaryBasis(unsigned bit) noexcept {
    return bit<64?Scalar{std::uint64_t(1)<<bit,0}:Scalar{0,std::uint64_t(1)<<(bit-64)};
}
constexpr Scalar multiplyX(Scalar a) noexcept {
    return {(a.lo<<1)^(0x87ULL*(a.hi>>63)),(a.hi<<1)|(a.lo>>63)};
}

// Scalar reference; setup and checking only. The native SIMD kernel is separate.
constexpr Scalar multiply(Scalar a,Scalar b) noexcept {
    Scalar result{};
    for(unsigned bit=0;bit!=128;++bit){
        const auto mask=std::uint64_t(0)-(b.lo&1);
        result.lo^=a.lo&mask;result.hi^=a.hi&mask;
        b={b.lo>>1|(b.hi<<63),b.hi>>1};a=multiplyX(a);
    }
    return result;
}

// Independent full 256-bit carryless product and polynomial long division.
// Do not replace this with the shift/reduce recurrence used above.
inline Scalar multiplyOracle(Scalar a,Scalar b) noexcept {
    std::array<std::uint64_t,4> product{};
    const std::uint64_t left[2]={a.lo,a.hi},right[2]={b.lo,b.hi};
    for(unsigned limb=0;limb!=2;++limb)for(unsigned bit=0;bit!=64;++bit){
        if(!(right[limb]>>bit&1))continue;
        for(unsigned source=0;source!=2;++source){
            product[limb+source]^=left[source]<<bit;
            if(bit)product[limb+source+1]^=left[source]>>(64-bit);
        }
    }
    const auto toggle=[&](unsigned bit){product[bit/64]^=std::uint64_t(1)<<(bit%64);};
    for(unsigned bit=256;bit-- >128;){
        if(!(product[bit/64]>>(bit%64)&1))continue;
        const auto base=bit-128;
        toggle(bit);toggle(base+7);toggle(base+2);toggle(base+1);toggle(base);
    }
    return {product[0],product[1]};
}

// Deterministic pin from field128_oracle.py --constants. Let
// g=X^((2^128-1)/255). AES root alpha=g^95 is the first root of
// T^8+T^4+T^3+T+1 among g^0,...,g^254. Byte 0x02 embeds as alpha.
inline constexpr Scalar projectedGenerator{0xd9db18ed46787786ULL,0xe6114072b8ca57afULL};
inline constexpr unsigned aesRootPower=95;
inline constexpr Scalar aesRoot{0x7573da4a5f7710edULL,0x3d5bd35c94646a24ULL};
inline constexpr auto aesPowers=[] {
    std::array<Scalar,8> result{};result[0]={1,0};
    for(unsigned i=1;i!=8;++i)result[i]=multiply(result[i-1],aesRoot);
    return result;
}();
constexpr Scalar embedByte(std::uint8_t value) noexcept {
    Scalar result{};
    for(unsigned bit=0;bit!=8;++bit)if(value>>bit&1)result^=aesPowers[bit];
    return result;
}
inline constexpr auto embeddedBytes=[] {
    std::array<Scalar,256> result{};
    for(unsigned i=0;i!=256;++i)result[i]=embedByte(static_cast<std::uint8_t>(i));
    return result;
}();
static_assert(embedByte(0)==Scalar{}&&embedByte(1)==Scalar{1,0});
static_assert(embedByte(2)==aesRoot);
static_assert((multiply(aesPowers[4],aesPowers[4])^aesPowers[4]^aesPowers[3]^aesRoot^Scalar{1,0})==Scalar{});

// Reference coordinate conversion. Byte j is a coefficient in AES GF256,
// multiplying native X^j, j=0,...,15. These are polynomial extension
// coordinates; they are not the unrelated quadratic-tower coordinates.
// Keep these scalar conversions out of encoding kernels and include any
// required eventual optimized conversion in end-to-end encoding costs.
struct ByteBasis {
    std::array<Scalar,128> toNativeColumns{},fromNativeColumns{};

    static Scalar apply(Scalar value,const std::array<Scalar,128>& columns) noexcept {
        Scalar result{};
        for(auto bits=value.lo;bits;bits&=bits-1)result^=columns[std::countr_zero(bits)];
        for(auto bits=value.hi;bits;bits&=bits-1)result^=columns[64+std::countr_zero(bits)];
        return result;
    }

    ByteBasis(){
        for(unsigned bit=0;bit!=8;++bit){
            auto value=aesPowers[bit];
            for(unsigned byte=0;byte!=16;++byte){
                toNativeColumns[8*byte+bit]=value;value=multiplyX(value);
            }
        }
        std::array<Scalar,128> pivots{},preimages{};
        const auto highest=[](Scalar value){
            return value.hi?64U+unsigned(std::bit_width(value.hi))-1:
                              unsigned(std::bit_width(value.lo))-1;
        };
        for(unsigned bit=0;bit!=128;++bit){
            auto value=toNativeColumns[bit],preimage=binaryBasis(bit);
            while(nonzero(value)){
                const auto position=highest(value);
                if(!nonzero(pivots[position])){
                    pivots[position]=value;preimages[position]=preimage;break;
                }
                value^=pivots[position];preimage^=preimages[position];
            }
            if(!nonzero(value))throw std::logic_error("dependent GF128 byte basis");
        }
        for(unsigned bit=0;bit!=128;++bit){
            auto value=binaryBasis(bit);Scalar preimage{};
            while(nonzero(value)){
                const auto position=highest(value);
                value^=pivots[position];preimage^=preimages[position];
            }
            fromNativeColumns[bit]=preimage;
        }
    }
    Scalar toNative(Scalar byteCoordinates) const noexcept {return apply(byteCoordinates,toNativeColumns);}
    Scalar fromNative(Scalar native) const noexcept {return apply(native,fromNativeColumns);}
};

} // namespace spin::research::char2::field128_lift
