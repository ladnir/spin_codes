#pragma once
#include "Punctured8Oracle.h"
#include <limits>

namespace spin::research::char2 {

// Untimed dense reference for the three rate-climb profiles. The field,
// directly evaluated outer columns, and dense 16-by-16 packet are shared
// with the frozen oracle. No native FFT, inner circuit, route kernel, or
// implementation-provided output index is used to compute a reference.

inline unsigned r8OracleRetained(unsigned dimension,std::size_t step,bool frozen75=false){
    if(frozen75){
        if(dimension!=84)throw std::invalid_argument("rate-climb oracle fixed-seven dimension");
        return 7;
    }
    switch(dimension){
    case 80:return 6;
    case 84:return step%4==3?7:6;
    case 74:return step%2==1?6:5;
    default:throw std::invalid_argument("rate-climb oracle outer dimension");
    }
}

template<class Plan>
inline constexpr bool r8OracleFrozen75(){
    static_assert((Plan::dimension==80&&Plan::period==1&&Plan::count==6)||
                  (Plan::dimension==84&&Plan::period==4&&Plan::count==25)||
                  (Plan::dimension==74&&Plan::period==2&&Plan::count==11)||
                  (Plan::dimension==84&&Plan::period==1&&Plan::count==7),
                  "rate-climb oracle requires a certified profile or frozen control");
    return Plan::dimension==84&&Plan::period==1&&Plan::count==7;
}

template<class Plan>
inline unsigned r8OracleRetained(const Plan&,std::size_t step){
    return r8OracleRetained(Plan::dimension,step,r8OracleFrozen75<Plan>());
}

inline std::size_t r8OracleSourceIndex(unsigned dimension,std::size_t row,
                                      unsigned coordinate){
    if(coordinate>=dimension)throw std::invalid_argument("rate-climb oracle source coordinate");
    return (row/4)*(4*dimension)+4*coordinate+row%4;
}

inline unsigned r8OracleTileRecords(unsigned dimension,bool frozen75=false){
    unsigned count=0;
    for(unsigned phase=0;phase!=4;++phase)count+=r8OracleRetained(dimension,phase,frozen75);
    return count;
}

inline std::size_t r8OracleOutputIndex(unsigned dimension,std::size_t step,
                                      unsigned coordinate,bool frozen75=false){
    if(coordinate>=r8OracleRetained(dimension,step,frozen75))
        throw std::invalid_argument("rate-climb oracle discarded output coordinate");
    // Independently enumerate a compact coordinate-major four-packet tile.
    // This intentionally does not reproduce the specialized native formulas.
    const auto tile=step/4;
    const auto phase=static_cast<unsigned>(step%4);
    std::size_t index=tile*r8OracleTileRecords(dimension,frozen75);
    for(unsigned j=0;j!=coordinate;++j)
        for(unsigned p=0;p!=4;++p)
            if(j<r8OracleRetained(dimension,4*tile+p,frozen75))++index;
    for(unsigned p=0;p!=phase;++p)
        if(coordinate<r8OracleRetained(dimension,4*tile+p,frozen75))++index;
    return index;
}

template<class Plan>
inline std::size_t r8OracleOutputIndex(const Plan&,std::size_t step,unsigned coordinate){
    return r8OracleOutputIndex(Plan::dimension,step,coordinate,r8OracleFrozen75<Plan>());
}

template<class Plan>
inline void r8OracleValidatePlan(const Plan& p){
    constexpr auto dimension=Plan::dimension;
    static_assert(dimension==74||dimension==80||dimension==84,
                  "rate-climb oracle supports the three certified profiles");
    constexpr auto maxRows=(std::numeric_limits<std::uint32_t>::max()/128/64)*64;
    constexpr std::size_t unit=dimension*64;
    if(!p.requested||!p.groups||p.groups>maxRows||p.groups%64||
       p.requested>std::numeric_limits<std::size_t>::max()-unit+1)
        throw std::invalid_argument("rate-climb oracle plan size");
    const auto expectedRows=((p.requested+unit-1)/unit)*64;
    const auto expectedInput=128*p.groups;
    const auto expectedOutput=(expectedInput/32)*r8OracleTileRecords(dimension,r8OracleFrozen75<Plan>());
    if(p.groups!=expectedRows||p.k!=dimension*p.groups||p.requested>p.k||
       p.inputN!=expectedInput||p.outputN!=expectedOutput||p.planeStride<p.groups||
       p.route.size()!=p.inputN||p.inverseRoute.size()!=p.inputN||
       p.inputDiagonal.size()!=p.inputN||p.stateDiagonal.size()!=p.inputN||
       (!p.inverseLocalRoute.empty()&&
        (p.groups>65536||p.inverseLocalRoute.size()!=p.inputN)))
        throw std::invalid_argument("rate-climb oracle plan geometry");
    for(std::size_t i=0;i!=p.inputN;++i){
        const std::size_t target=p.route[i];
        if(target>=p.inputN||target/p.groups!=i/p.groups||p.inverseRoute[target]!=i||
           !p.inputDiagonal[i]||!p.stateDiagonal[i])
            throw std::invalid_argument("rate-climb oracle route or diagonal");
        if(!p.inverseLocalRoute.empty()&&
           p.inverseLocalRoute[i]!=p.inverseRoute[i]-(i/p.groups)*p.groups)
            throw std::invalid_argument("rate-climb oracle compact inverse route");
    }
    if(p.packedInputDiagonal.size()!=p.inputN||p.foldedRoute.size()!=p.inputN||
       p.foldedInverseRoute.size()!=p.inputN||
       (!p.foldedInverseLocalRoute.empty()&&
        (p.groups>65536||p.foldedInverseLocalRoute.size()!=p.inputN)))
        throw std::invalid_argument("rate-climb oracle folded setup geometry");
    // Enumerate packet and coordinate labels directly rather than calling
    // foldedIndex. The state diagonal remains indexed by chronological step.
    for(std::size_t tile=0;tile!=p.inputN;tile+=32)
        for(unsigned coordinate=0;coordinate!=8;++coordinate)
            for(unsigned phase=0;phase!=4;++phase){
                const auto original=tile+8*phase+coordinate;
                const auto packed=tile+4*coordinate+phase;
                const auto target=p.route[original];
                if(packed/p.groups!=original/p.groups||
                   p.packedInputDiagonal[packed]!=p.inputDiagonal[original]||
                   p.foldedRoute[packed]!=target||p.foldedInverseRoute[target]!=packed)
                    throw std::invalid_argument("rate-climb oracle folded route or diagonal");
                if(!p.foldedInverseLocalRoute.empty()&&
                   p.foldedInverseLocalRoute[target]!=packed-(target/p.groups)*p.groups)
                    throw std::invalid_argument("rate-climb oracle folded compact inverse route");
            }
    std::size_t emitted=0;
    for(std::size_t step=0;step!=p.inputN/8;++step){
        const auto count=r8OracleRetained(p,step);
        if(p.retained(step)!=count)
            throw std::invalid_argument("rate-climb oracle retained-output pattern");
        for(unsigned j=0;j!=count;++j){
            const auto index=r8OracleOutputIndex(p,step,j);
            if(index>=p.outputN||p.outputIndex(step,j)!=index)
                throw std::invalid_argument("rate-climb oracle output layout");
            ++emitted;
        }
    }
    if(emitted!=p.outputN)throw std::invalid_argument("rate-climb oracle output count");
    for(std::size_t row=0;row!=p.groups;++row)
        for(unsigned j=0;j!=dimension;++j)
            if(p.sourceIndex(row,j)!=r8OracleSourceIndex(dimension,row,j))
                throw std::invalid_argument("rate-climb oracle source layout");
}

// Dense outer helpers act on one row in ordinary coordinate order, without
// the full encoder's four-row interleaving. All supported dimensions use a
// prefix of the independently evaluated additive-polynomial basis.
inline void r8OracleOuterForward(const Block* input,Block* output,unsigned dimension){
    (void)r8OracleRetained(dimension,0);
    const auto& o=p8Oracle();
    for(unsigned point=0;point!=128;++point)for(unsigned byte=0;byte!=16;++byte){
        std::uint8_t value=0;
        for(unsigned degree=0;degree!=dimension;++degree)
            value^=o.mul[o.evaluation[point][degree]][p8OracleBytes(input+degree)[byte]];
        p8OracleBytes(output+point)[byte]=value;
    }
}

inline void r8OracleOuterTranspose(const Block* input,Block* output,unsigned dimension){
    (void)r8OracleRetained(dimension,0);
    const auto& o=p8Oracle();
    for(unsigned degree=0;degree!=dimension;++degree)for(unsigned byte=0;byte!=16;++byte){
        std::uint8_t value=0;
        for(unsigned point=0;point!=128;++point)
            value^=o.mul[o.evaluation[point][degree]][p8OracleBytes(input+point)[byte]];
        p8OracleBytes(output+degree)[byte]=value;
    }
}

template<class Plan>
inline std::vector<Block> r8OracleInnerForward(const Plan& p,const Block* planes){
    r8OracleValidatePlan(p);
    std::vector<Block> result(p.outputN);
    std::array<Block,8> state{},following{},packet{},emitted{};
    for(std::size_t step=0;step!=p.inputN/8;++step){
        for(unsigned j=0;j!=8;++j)packet[j]=planes[p.route[8*step+j]];
        p8OraclePacket(packet.data(),state.data(),p.inputDiagonal.data()+8*step,
                       p.stateDiagonal.data()+8*step,emitted.data(),following.data());
        for(unsigned j=0;j!=r8OracleRetained(p,step);++j)
            result[r8OracleOutputIndex(p,step,j)]=emitted[j];
        state=following;
    }
    return result;
}

template<class Plan>
inline std::vector<Block> r8OracleInnerTranspose(const Plan& p,const Block* input){
    r8OracleValidatePlan(p);
    std::vector<Block> planes(p.inputN);
    std::array<Block,8> state{},previous{},packet{},raw{};
    // A zero terminal dual corresponds to discarding the forward state.
    for(std::size_t step=p.inputN/8;step--;){
        packet.fill(Block{});
        for(unsigned j=0;j!=r8OracleRetained(p,step);++j)
            packet[j]=input[r8OracleOutputIndex(p,step,j)];
        p8OraclePacket(packet.data(),state.data(),p.inputDiagonal.data()+8*step,
                       p.stateDiagonal.data()+8*step,raw.data(),previous.data(),true);
        for(unsigned j=0;j!=8;++j)planes[p.route[8*step+j]]=raw[j];
        state=previous;
    }
    return planes;
}

template<class Plan>
inline std::vector<Block> r8OracleForward(const Plan& p,const Block* input){
    r8OracleValidatePlan(p);const auto& o=p8Oracle();
    std::vector<Block> planes(p.inputN);
    for(std::size_t row=0;row!=p.groups;++row)
        for(unsigned coordinate=0;coordinate!=Plan::dimension;++coordinate){
            const auto source=r8OracleSourceIndex(Plan::dimension,row,coordinate);
            if(source>=p.requested)continue;
            for(unsigned point=0;point!=128;++point)for(unsigned byte=0;byte!=16;++byte)
                p8OracleBytes(planes.data()+point*p.groups+row)[byte]^=
                    o.mul[o.evaluation[point][coordinate]][p8OracleBytes(input+source)[byte]];
        }
    return r8OracleInnerForward(p,planes.data());
}

template<class Plan>
inline std::vector<Block> r8OracleTranspose(const Plan& p,const Block* input){
    const auto planes=r8OracleInnerTranspose(p,input);const auto& o=p8Oracle();
    std::vector<Block> result(p.requested);
    constexpr auto tile=4*Plan::dimension;
    for(std::size_t source=0;source!=p.requested;++source){
        const auto row=(source/tile)*4+source%4;
        const auto coordinate=static_cast<unsigned>((source%tile)/4);
        for(unsigned point=0;point!=128;++point)for(unsigned byte=0;byte!=16;++byte)
            p8OracleBytes(result.data()+source)[byte]^=
                o.mul[o.evaluation[point][coordinate]][p8OracleBytes(planes.data()+point*p.groups+row)[byte]];
    }
    return result;
}

inline void r8OracleCheckLocal(){
    p8OracleCheckLocal();
    // Check the complete independent tile enumeration, including the omitted
    // tail coordinates and consecutive tile boundaries for all profiles.
    for(unsigned dimension:{74U,80U,84U}){
        const auto tile=r8OracleTileRecords(dimension);
        std::vector<bool> seen(3*tile,false);
        for(std::size_t step=0;step!=12;++step){
            for(unsigned j=0;j!=r8OracleRetained(dimension,step);++j){
                const auto index=r8OracleOutputIndex(dimension,step,j);
                if(index>=seen.size()||seen[index])
                    throw std::runtime_error("rate-climb oracle tile is not bijective");
                seen[index]=true;
            }
        }
        if(std::find(seen.begin(),seen.end(),false)!=seen.end())
            throw std::runtime_error("rate-climb oracle tile has holes");
        for(std::size_t row=0;row!=12;++row)for(unsigned j=0;j!=dimension;++j){
            const auto index=r8OracleSourceIndex(dimension,row,j);
            if((index/(4*dimension))*4+index%4!=row||index%(4*dimension)/4!=j)
                throw std::runtime_error("rate-climb oracle source inverse");
        }
    }
    for(std::size_t step=0;step!=12;++step)for(unsigned j=0;j!=7;++j)
        if(r8OracleOutputIndex(84,step,j,true)!=p8OracleOutputIndex(step,j))
            throw std::runtime_error("rate-climb oracle frozen-control layout");
}

// This is an untimed dense check; call it on small natural plans, not K20.
template<class Plan>
inline void r8OracleCheckPlan(const Plan& p,std::uint64_t seed=0x8174){
    r8OracleValidatePlan(p);std::mt19937_64 random(seed);
    std::vector<Block> source(p.requested),dual(p.outputN);
    for(auto& value:source)value={random(),random()};
    for(auto& value:dual)value={random(),random()};
    const auto image=r8OracleForward(p,source.data());
    const auto back=r8OracleTranspose(p,dual.data());
    if(p8OracleDot(image.data(),dual.data(),p.outputN)!=p8OracleDot(source.data(),back.data(),p.requested))
        throw std::runtime_error("rate-climb oracle shortened full adjoint");
}

} // namespace spin::research::char2
