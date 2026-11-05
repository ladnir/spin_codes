#pragma once
#include <array>
#include <cstdint>

namespace spin::research::char2::field128_lift::sparse_byte {

// F=AES GF256[Y]/(Y^16+Y^5+Y^3+0x32). Byte j holds the coefficient
// of Y^j, with little-endian bytes in each word. The irreducibility check
// is in test_sparse_byte_field128.py. This type deliberately differs from
// the native binary-polynomial Scalar; raw bytes are not interchangeable.
struct Element {
    std::uint64_t lo=0,hi=0;
    constexpr bool operator==(const Element&) const noexcept = default;
};
static_assert(sizeof(Element)==16);
constexpr Element operator^(Element a,Element b) noexcept {return {a.lo^b.lo,a.hi^b.hi};}
constexpr Element& operator^=(Element& a,Element b) noexcept {a=a^b;return a;}
constexpr Element binaryBasis(unsigned bit) noexcept {
    return bit<64?Element{std::uint64_t(1)<<bit,0}:Element{0,std::uint64_t(1)<<(bit-64)};
}
constexpr std::uint8_t coefficient(Element a,unsigned index) noexcept {
    return static_cast<std::uint8_t>((index<8?a.lo:a.hi)>>(8*(index%8)));
}
constexpr Element fromCoefficients(const std::array<std::uint8_t,16>& a) noexcept {
    Element result{};
    for(unsigned i=0;i!=8;++i){result.lo|=std::uint64_t(a[i])<<(8*i);result.hi|=std::uint64_t(a[i+8])<<(8*i);}
    return result;
}
constexpr Element embedByte(std::uint8_t a) noexcept {return {a,0};}
inline constexpr std::array<std::uint8_t,17> modulus{
    0x32,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,1};

// Independent AES-byte schoolbook multiplication and polynomial division.
constexpr std::uint8_t byteMultiply(std::uint8_t a,std::uint8_t b) noexcept {
    unsigned product=0;
    for(unsigned bit=0;bit!=8;++bit)if(b>>bit&1)product^=unsigned(a)<<bit;
    for(unsigned bit=15;bit-- >8;)if(product>>bit&1)product^=0x11bU<<(bit-8);
    return static_cast<std::uint8_t>(product);
}

// Generic coefficient-polynomial oracle, independent of the SIMD spacing,
// CLMUL Karatsuba product, and two-fold sparse reduction. Setup/checking only.
inline Element multiplyOracle(Element a,Element b) noexcept {
    std::array<std::uint8_t,31> product{};
    for(unsigned i=0;i!=16;++i){
        const auto left=coefficient(a,i);if(!left)continue;
        for(unsigned j=0;j!=16;++j){
            const auto right=coefficient(b,j);
            if(right)product[i+j]^=byteMultiply(left,right);
        }
    }
    for(unsigned degree=31;degree-- >16;){
        const auto top=product[degree];if(!top)continue;
        for(unsigned j=0;j!=16;++j)if(modulus[j])product[degree-16+j]^=byteMultiply(top,modulus[j]);
        product[degree]=0;
    }
    std::array<std::uint8_t,16> result{};
    for(unsigned i=0;i!=16;++i)result[i]=product[i];
    return fromCoefficients(result);
}

// A second scalar reference uses Horner evaluation and multiplication by Y.
// In particular it does not construct or reduce the oracle's degree30 array.
inline Element multiply(Element a,Element b) noexcept {
    std::array<std::uint8_t,16> result{};
    for(unsigned i=16;i-- >0;){
        const auto top=result[15];
        for(unsigned j=15;j!=0;--j)result[j]=result[j-1];
        result[0]=byteMultiply(top,0x32);result[3]^=top;result[5]^=top;
        const auto c=coefficient(a,i);
        if(c)for(unsigned j=0;j!=16;++j)result[j]^=byteMultiply(c,coefficient(b,j));
    }
    return fromCoefficients(result);
}

} // namespace spin::research::char2::field128_lift::sparse_byte
