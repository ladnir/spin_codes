#pragma once
#include "Punctured8Oracle.h"
namespace spin::research::char2 {
// Independent dense4/4 map: construct Cauchy8x8 and Gauss-Jordan normalize E.
// No constants from the optimized4-point circuits are used by this oracle.
struct RateSpeed4Oracle {
    std::uint8_t matrix[8][8]{};
    RateSpeed4Oracle(){
        const auto& o=p8Oracle();std::uint8_t q[8][8]{},augmented[4][8]{};
        for(unsigned i=0;i!=8;++i)for(unsigned j=0;j!=8;++j)q[i][j]=o.inverse[i^(128+j)];
        for(unsigned i=0;i!=4;++i)for(unsigned j=0;j!=4;++j){augmented[i][j]=q[4+i][4+j];augmented[i][4+j]=static_cast<std::uint8_t>(i==j);}
        for(unsigned column=0;column!=4;++column){
            unsigned pivot=column;while(pivot!=4&&!augmented[pivot][column])++pivot;
            if(pivot==4)throw std::runtime_error("four-state oracle singular E");
            for(unsigned j=0;j!=8;++j)std::swap(augmented[pivot][j],augmented[column][j]);
            const auto scale=o.inverse[augmented[column][column]];
            for(auto& value:augmented[column])value=o.mul[value][scale];
            for(unsigned row=0;row!=4;++row)if(row!=column){const auto factor=augmented[row][column];for(unsigned j=0;j!=8;++j)augmented[row][j]^=o.mul[factor][augmented[column][j]];}
        }
        std::uint8_t injection[4][4]{},emissionSum=0,injectionSum=0;
        for(unsigned i=0;i!=4;++i)for(unsigned j=0;j!=4;++j)for(unsigned h=0;h!=4;++h)injection[i][j]^=o.mul[augmented[i][4+h]][q[4+h][j]];
        for(unsigned j=0;j!=4;++j){emissionSum^=q[0][4+j];injectionSum^=injection[0][j];}
        const auto left=o.inverse[emissionSum],right=o.inverse[injectionSum];
        for(unsigned i=0;i!=4;++i)for(unsigned j=0;j!=4;++j){
            matrix[i][j]=o.mul[o.mul[left][right]][q[i][j]];
            matrix[i][4+j]=o.mul[left][q[i][4+j]];matrix[4+i][j]=o.mul[right][injection[i][j]];
            matrix[4+i][4+j]=static_cast<std::uint8_t>(i==j);
        }
    }
};
inline const RateSpeed4Oracle& rs4Oracle(){static const RateSpeed4Oracle value;return value;}
inline void rs4OraclePacket(const Block* input,const Block* state,const std::uint8_t* dx,const std::uint8_t* ds,Block* output,Block* next,bool transpose){
    const auto& o=p8Oracle();const auto& map=rs4Oracle();
    for(unsigned byte=0;byte!=16;++byte){
        std::uint8_t values[8]{},image[8]{};
        for(unsigned i=0;i!=4;++i){values[i]=p8OracleBytes(input+i)[byte];values[4+i]=p8OracleBytes(state+i)[byte];if(!transpose){values[i]=o.mul[dx[i]][values[i]];values[4+i]=o.mul[ds[i]][values[4+i]];}}
        for(unsigned i=0;i!=8;++i)for(unsigned j=0;j!=8;++j)image[i]^=o.mul[transpose?map.matrix[j][i]:map.matrix[i][j]][values[j]];
        for(unsigned i=0;i!=4;++i){p8OracleBytes(output+i)[byte]=transpose?o.mul[dx[i]][image[i]]:image[i];p8OracleBytes(next+i)[byte]=transpose?o.mul[ds[i]][image[4+i]]:image[4+i];}
    }
}
template<class Plan>
inline std::vector<Block> rs4OracleTranspose(const Plan& p,const Block* input){
    const auto& o=p8Oracle();std::vector<Block> planes(p.inputN),result(p.requested);
    std::array<Block,4> state{},prior{},packet{},raw{};
    for(auto step=p.inputN/4;step--;){
        packet.fill(Block{});for(unsigned j=0;j!=3;++j)packet[j]=input[12*(step/4)+4*j+step%4];
        rs4OraclePacket(packet.data(),state.data(),p.inputDiagonal.data()+4*step,p.stateDiagonal.data()+4*step,raw.data(),prior.data(),true);
        for(unsigned j=0;j!=4;++j)planes[p.route[4*step+j]]=raw[j];state=prior;
    }
    for(std::size_t source=0;source!=p.requested;++source){
        const auto row=(source/(4*Plan::dimension))*4+source%4,coordinate=(source%(4*Plan::dimension))/4;
        for(unsigned point=0;point!=64;++point)for(unsigned byte=0;byte!=16;++byte)p8OracleBytes(result.data()+source)[byte]^=o.mul[o.evaluation[point][coordinate]][p8OracleBytes(planes.data()+point*p.groups+row)[byte]];
    }
    return result;
}
template<class Plan>
inline std::vector<Block> rs4OracleForward(const Plan& p,const Block* input){
    const auto& o=p8Oracle();std::vector<Block> planes(p.inputN),result(p.outputN);
    for(std::size_t source=0;source!=p.requested;++source){
        const auto row=(source/(4*Plan::dimension))*4+source%4,coordinate=(source%(4*Plan::dimension))/4;
        for(unsigned point=0;point!=64;++point)for(unsigned byte=0;byte!=16;++byte)p8OracleBytes(planes.data()+point*p.groups+row)[byte]^=o.mul[o.evaluation[point][coordinate]][p8OracleBytes(input+source)[byte]];
    }
    std::array<Block,4> state{},next{},packet{},raw{};
    for(std::size_t step=0;step!=p.inputN/4;++step){
        for(unsigned j=0;j!=4;++j)packet[j]=planes[p.route[4*step+j]];
        rs4OraclePacket(packet.data(),state.data(),p.inputDiagonal.data()+4*step,p.stateDiagonal.data()+4*step,raw.data(),next.data(),false);
        for(unsigned j=0;j!=3;++j)result[12*(step/4)+4*j+step%4]=raw[j];state=next;
    }
    return result;
}
} // namespace spin::research::char2
