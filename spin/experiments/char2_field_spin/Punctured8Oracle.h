#pragma once
#include "Full.h"
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <random>
#include <stdexcept>
#include <vector>

namespace spin::research::char2 {

// Untimed, deliberately dense reference. No optimized packet, field helper,
// circulant circuit, or outer FFT is called. Plan is a template parameter so
// this header does not depend on Punctured8Code.h or Punctured8Inner.h.
inline std::uint8_t p8OracleMultiply(unsigned a,unsigned b)noexcept{
    unsigned result=0;
    for(unsigned bit=0;bit!=8;++bit){
        if(b&1U)result^=a;
        b>>=1;a<<=1;if(a&256U)a^=0x11bU;
    }
    return static_cast<std::uint8_t>(result);
}

struct Punctured8Oracle {
    static constexpr unsigned dimension=84,outer=128,width=8,state=8,outputs=7;
    std::array<std::array<std::uint8_t,256>,256> mul{};
    std::array<std::uint8_t,256> inverse{};
    std::array<std::array<std::uint8_t,84>,128> evaluation{};
    // a is the nominal A, so eta=2. effectiveA=2*A.
    std::uint8_t f[8][8]{},a[8][8]{},c[8][8]{},effectiveA[8][8]{};
    std::uint8_t extended[16][16]{};

    Punctured8Oracle(){
        for(unsigned x=0;x!=256;++x)for(unsigned y=0;y!=256;++y)
            mul[x][y]=p8OracleMultiply(x,y);
        for(unsigned x=1;x!=256;++x){
            unsigned y=1;while(y!=256&&mul[x][y]!=1)++y;
            if(y==256)throw std::runtime_error("punctured oracle field inverse");
            inverse[x]=static_cast<std::uint8_t>(y);
        }
        std::uint8_t q[16][16]{};
        for(unsigned i=0;i!=16;++i)for(unsigned j=0;j!=16;++j)q[i][j]=inverse[i^(128+j)];
        // Independently invert E by Gauss-Jordan; do not use the circulant
        // identity E^-1=E/(row sum)^2 used to design the fast circuit.
        std::uint8_t augmented[8][16]{};
        for(unsigned i=0;i!=8;++i)for(unsigned j=0;j!=8;++j){
            augmented[i][j]=q[8+i][8+j];augmented[i][8+j]=static_cast<std::uint8_t>(i==j);
        }
        for(unsigned column=0;column!=8;++column){
            unsigned pivot=column;
            while(pivot!=8&&!augmented[pivot][column])++pivot;
            if(pivot==8)throw std::runtime_error("punctured oracle singular E");
            for(unsigned j=0;j!=16;++j)std::swap(augmented[column][j],augmented[pivot][j]);
            const auto scale=inverse[augmented[column][column]];
            for(auto& value:augmented[column])value=mul[scale][value];
            for(unsigned i=0;i!=8;++i)if(i!=column){
                const auto coefficient=augmented[i][column];
                for(unsigned j=0;j!=16;++j)augmented[i][j]^=mul[coefficient][augmented[column][j]];
            }
        }
        std::uint8_t injection[8][8]{};
        for(unsigned i=0;i!=8;++i)for(unsigned j=0;j!=8;++j)
            for(unsigned h=0;h!=8;++h)injection[i][j]^=mul[augmented[i][8+h]][q[8+h][j]];
        std::uint8_t emissionSum=0,injectionSum=0;
        for(unsigned j=0;j!=8;++j){emissionSum^=q[0][8+j];injectionSum^=injection[0][j];}
        if(emissionSum!=161||injectionSum!=46)
            throw std::runtime_error("punctured oracle normalization geometry");
        const auto left=inverse[emissionSum],right=inverse[injectionSum];
        for(unsigned i=0;i!=8;++i)for(unsigned j=0;j!=8;++j){
            f[i][j]=mul[mul[left][right]][q[i][j]];
            effectiveA[i][j]=mul[left][q[i][8+j]];
            a[i][j]=mul[inverse[2]][effectiveA[i][j]];
            c[i][j]=mul[right][injection[i][j]];
            extended[i][j]=f[i][j];extended[i][8+j]=effectiveA[i][j];
            extended[8+i][j]=c[i][j];extended[8+i][8+j]=static_cast<std::uint8_t>(i==j);
        }
        // Direct evaluation of normalized additive polynomials. This does
        // not use interpolation, generated FFT constants, or FFT stages.
        for(unsigned point=0;point!=128;++point)for(unsigned degree=0;degree!=84;++degree){
            std::uint8_t value=1;
            for(unsigned bit=0;bit!=7;++bit)if(degree&(1U<<bit)){
                std::uint8_t numerator=1,denominator=1;
                for(unsigned label=0;label!=(1U<<bit);++label){
                    numerator=mul[numerator][point^label];
                    denominator=mul[denominator][(1U<<bit)^label];
                }
                value=mul[value][mul[numerator][inverse[denominator]]];
            }
            evaluation[point][degree]=value;
        }
    }
};

inline const Punctured8Oracle& p8Oracle(){static const Punctured8Oracle value;return value;}
inline std::uint8_t* p8OracleBytes(Block* p)noexcept{return reinterpret_cast<std::uint8_t*>(p);}
inline const std::uint8_t* p8OracleBytes(const Block* p)noexcept{return reinterpret_cast<const std::uint8_t*>(p);}
inline std::size_t p8OracleSourceIndex(std::size_t row,unsigned coordinate)noexcept{
    return (row/4)*336+4*coordinate+row%4;
}
inline std::size_t p8OracleOutputIndex(std::size_t step,unsigned coordinate)noexcept{
    // Forward packet order inside each coordinate-major four-packet tile.
    return 28*(step/4)+4*coordinate+step%4;
}

// Full eight-output packet. Both diagonals are included in this function.
// For a punctured transpose the caller supplies input[7]=0. The dense literal
// transpose includes F^T, C^T, (2A)^T, and both diagonal placements.
// Exact input/output alias and state/following alias are safe: all sixteen
// source coordinates for a byte are captured before any destination write.
inline void p8OraclePacket(const Block* input,const Block* state,
                           const std::uint8_t* inputDiagonal,const std::uint8_t* stateDiagonal,
                           Block* output,Block* following,bool transpose=false){
    const auto& o=p8Oracle();
    for(unsigned byte=0;byte!=16;++byte){
        std::uint8_t values[16]{},image[16]{};
        for(unsigned i=0;i!=8;++i){
            values[i]=p8OracleBytes(input+i)[byte];values[8+i]=p8OracleBytes(state+i)[byte];
            if(!transpose){
                values[i]=o.mul[inputDiagonal[i]][values[i]];
                values[8+i]=o.mul[stateDiagonal[i]][values[8+i]];
            }
        }
        for(unsigned i=0;i!=16;++i)for(unsigned j=0;j!=16;++j)
            image[i]^=o.mul[transpose?o.extended[j][i]:o.extended[i][j]][values[j]];
        for(unsigned i=0;i!=8;++i){
            p8OracleBytes(output+i)[byte]=transpose?o.mul[inputDiagonal[i]][image[i]]:image[i];
            p8OracleBytes(following+i)[byte]=transpose?o.mul[stateDiagonal[i]][image[8+i]]:image[8+i];
        }
    }
}

inline std::array<std::uint8_t,16> p8OracleDot(const Block* left,const Block* right,std::size_t count){
    const auto& o=p8Oracle();std::array<std::uint8_t,16> result{};
    for(std::size_t i=0;i!=count;++i)for(unsigned byte=0;byte!=16;++byte)
        result[byte]^=o.mul[p8OracleBytes(left+i)[byte]][p8OracleBytes(right+i)[byte]];
    return result;
}

template<class Plan>
inline void p8OracleValidatePlan(const Plan& p){
    if(!p.requested||!p.groups||p.groups%64||p.k!=84*p.groups||p.requested>p.k||
       p.inputN!=128*p.groups||p.outputN!=7*(p.inputN/8)||
       p.route.size()!=p.inputN||p.inverseRoute.size()!=p.inputN||
       p.inputDiagonal.size()!=p.inputN||p.stateDiagonal.size()!=p.inputN)
        throw std::invalid_argument("punctured oracle plan geometry");
    for(std::size_t i=0;i!=p.inputN;++i){
        const std::size_t target=p.route[i];
        if(target>=p.inputN||target/p.groups!=i/p.groups||p.inverseRoute[target]!=i||
           !p.inputDiagonal[i]||!p.stateDiagonal[i])
            throw std::invalid_argument("punctured oracle route or diagonal");
    }
}

// Plan-independent outer helpers for one scalar-code row of16-byte records.
inline void p8OracleOuterForward(const Block* input84,Block* output128){
    const auto& o=p8Oracle();
    for(unsigned point=0;point!=128;++point)for(unsigned byte=0;byte!=16;++byte){
        std::uint8_t value=0;
        for(unsigned degree=0;degree!=84;++degree)
            value^=o.mul[o.evaluation[point][degree]][p8OracleBytes(input84+degree)[byte]];
        p8OracleBytes(output128+point)[byte]=value;
    }
}
inline void p8OracleOuterTranspose(const Block* input128,Block* output84){
    const auto& o=p8Oracle();
    for(unsigned degree=0;degree!=84;++degree)for(unsigned byte=0;byte!=16;++byte){
        std::uint8_t value=0;
        for(unsigned point=0;point!=128;++point)
            value^=o.mul[o.evaluation[point][degree]][p8OracleBytes(input128+point)[byte]];
        p8OracleBytes(output84+degree)[byte]=value;
    }
}

template<class Plan>
inline std::vector<Block> p8OracleInnerForward(const Plan& p,const Block* planes){
    p8OracleValidatePlan(p);
    std::vector<Block> result(p.outputN);
    std::array<Block,8> state{},following{},packet{},emitted{};
    for(std::size_t step=0;step!=p.inputN/8;++step){
        for(unsigned j=0;j!=8;++j)packet[j]=planes[p.route[8*step+j]];
        p8OraclePacket(packet.data(),state.data(),p.inputDiagonal.data()+8*step,
                       p.stateDiagonal.data()+8*step,emitted.data(),following.data());
        for(unsigned j=0;j!=7;++j)result[p8OracleOutputIndex(step,j)]=emitted[j];
        state=following;
    }
    return result; // Forward terminal state is deliberately discarded.
}

template<class Plan>
inline std::vector<Block> p8OracleInnerTranspose(const Plan& p,const Block* input){
    p8OracleValidatePlan(p);
    std::vector<Block> planes(p.inputN);
    std::array<Block,8> state{},previous{},packet{},raw{};
    // Zero terminal dual is the transpose of discarding the terminal state.
    for(std::size_t step=p.inputN/8;step--;){
        for(unsigned j=0;j!=7;++j)packet[j]=input[p8OracleOutputIndex(step,j)];
        packet[7]=Block{};
        p8OraclePacket(packet.data(),state.data(),p.inputDiagonal.data()+8*step,
                       p.stateDiagonal.data()+8*step,raw.data(),previous.data(),true);
        for(unsigned j=0;j!=8;++j)planes[p.route[8*step+j]]=raw[j];
        state=previous;
    }
    return planes; // Initial state was fixed to zero; its dual is discarded.
}

template<class Plan>
inline std::vector<Block> p8OracleForward(const Plan& p,const Block* input){
    p8OracleValidatePlan(p);const auto& o=p8Oracle();
    std::vector<Block> planes(p.inputN);
    for(std::size_t row=0;row!=p.groups;++row)for(unsigned coordinate=0;coordinate!=84;++coordinate){
        const auto source=p8OracleSourceIndex(row,coordinate);if(source>=p.requested)continue;
        for(unsigned point=0;point!=128;++point)for(unsigned byte=0;byte!=16;++byte)
            p8OracleBytes(planes.data()+point*p.groups+row)[byte]^=
                o.mul[o.evaluation[point][coordinate]][p8OracleBytes(input+source)[byte]];
    }
    return p8OracleInnerForward(p,planes.data());
}

template<class Plan>
inline std::vector<Block> p8OracleTranspose(const Plan& p,const Block* input){
    const auto planes=p8OracleInnerTranspose(p,input);const auto& o=p8Oracle();
    std::vector<Block> result(p.requested);
    for(std::size_t source=0;source!=p.requested;++source){
        const auto row=(source/336)*4+source%4;
        const auto coordinate=static_cast<unsigned>((source%336)/4);
        for(unsigned point=0;point!=128;++point)for(unsigned byte=0;byte!=16;++byte)
            p8OracleBytes(result.data()+source)[byte]^=
                o.mul[o.evaluation[point][coordinate]][p8OracleBytes(planes.data()+point*p.groups+row)[byte]];
    }
    return result;
}

inline void p8OracleCheckLocal(){
    static_assert(sizeof(Block)==16,"punctured oracle requires16-byte records");
    const auto& o=p8Oracle();
    for(unsigned i=0;i!=8;++i)for(unsigned j=0;j!=8;++j){
        std::uint8_t product=0;
        for(unsigned h=0;h!=8;++h)product^=o.mul[o.effectiveA[i][h]][o.c[h][j]];
        if(o.f[i][j]!=o.mul[242][product])throw std::runtime_error("punctured oracle F=242AC");
    }
    std::mt19937_64 random(0x8a8407);
    for(unsigned trial=0;trial!=264;++trial){
        std::array<Block,8> input{},state{},dual{},terminal{},forward{},following{},backward{},previous{};
        if(trial<256){
            const unsigned coordinate=trial/16,base=16*(trial%16);
            auto* target=coordinate<8?input.data()+coordinate:state.data()+coordinate-8;
            for(unsigned byte=0;byte!=16;++byte)p8OracleBytes(target)[byte]=static_cast<std::uint8_t>(base+byte);
        }else{
            for(auto& value:input)value={random(),random()};
            for(auto& value:state)value={random(),random()};
        }
        for(auto& value:dual)value={random(),random()};
        for(auto& value:terminal)value={random(),random()};
        std::uint8_t dx[8],ds[8];
        for(unsigned j=0;j!=8;++j){dx[j]=static_cast<std::uint8_t>(1+random()%255);ds[j]=static_cast<std::uint8_t>(1+random()%255);}
        p8OraclePacket(input.data(),state.data(),dx,ds,forward.data(),following.data());
        p8OraclePacket(dual.data(),terminal.data(),dx,ds,backward.data(),previous.data(),true);
        auto left=p8OracleDot(forward.data(),dual.data(),8),tail=p8OracleDot(following.data(),terminal.data(),8);
        auto right=p8OracleDot(input.data(),backward.data(),8),head=p8OracleDot(state.data(),previous.data(),8);
        for(unsigned byte=0;byte!=16;++byte)if((left[byte]^tail[byte])!=(right[byte]^head[byte]))
            throw std::runtime_error("punctured oracle local literal adjoint");
        auto aliasInput=input,aliasState=state;
        p8OraclePacket(aliasInput.data(),aliasState.data(),dx,ds,aliasInput.data(),aliasState.data());
        if(aliasInput!=forward||aliasState!=following)throw std::runtime_error("punctured oracle local alias");
    }
}

// Deliberately expensive: use small natural plans in --check, not K20 timing.
template<class Plan>
inline void p8OracleCheckPlan(const Plan& p,std::uint64_t seed=0x7184){
    p8OracleValidatePlan(p);std::mt19937_64 random(seed);
    std::vector<Block> source(p.requested),dual(p.outputN);
    for(auto& value:source)value={random(),random()};
    for(auto& value:dual)value={random(),random()};
    const auto image=p8OracleForward(p,source.data());
    const auto back=p8OracleTranspose(p,dual.data());
    if(p8OracleDot(image.data(),dual.data(),p.outputN)!=p8OracleDot(source.data(),back.data(),p.requested))
        throw std::runtime_error("punctured oracle shortened full adjoint");
}

} // namespace spin::research::char2
