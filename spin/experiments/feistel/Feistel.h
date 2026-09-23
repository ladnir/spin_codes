#pragma once
// Experimental routing only. Not a cipher API or a certified SPIN profile.
#include "../../src/kernels/SetupRandom.h"
#include <array>
#include <cstdint>
#include <type_traits>
#include <vector>

namespace spin::experimental::feistel {
using kernelWords=detail::kernel::setup::Words;
template<unsigned Bits,unsigned Rounds> class Permutation {
    static_assert(Bits>=6 && Bits<=24);
    static_assert(Rounds==2 || Rounds==4 || Rounds==6 || Rounds==8);
    static constexpr unsigned Left=Bits/2,Right=Bits-Left;
    static constexpr unsigned NL=1U<<Left,NR=1U<<Right,Stride=NL+NR;
    using Entry=std::conditional_t<(Bits<=16),std::uint8_t,std::uint16_t>;
    std::array<Entry,Stride*(Rounds/2)> table_;
    template<unsigned Offset,unsigned Count,unsigned Mask>
    void fillRound(kernelWords& words) {
        constexpr unsigned Pack=8/sizeof(Entry);
        for(unsigned i=0;i<Count;i+=Pack) {
            const auto word=words();
            for(unsigned j=0;j<Pack;++j)table_[Offset+i+j]=Entry((word>>(8*sizeof(Entry)*j))&Mask);
        }
    }
    template<unsigned Pair=0> void fillPairs(kernelWords& words) {
        if constexpr(Pair<Rounds/2) {
            fillRound<Pair*Stride,NR,NL-1>(words);
            fillRound<Pair*Stride+NR,NL,NR-1>(words);
            fillPairs<Pair+1>(words);
        }
    }
    template<unsigned Pair=0> void forwardPairs(unsigned& l,unsigned& r) const {
        if constexpr(Pair<Rounds/2) {
            l^=table_[Pair*Stride+r];
            r^=table_[Pair*Stride+NR+l];
            forwardPairs<Pair+1>(l,r);
        }
    }
    template<unsigned Pair=Rounds/2> void inversePairs(unsigned& l,unsigned& r) const {
        if constexpr(Pair>0) {
            r^=table_[(Pair-1)*Stride+NR+l];
            l^=table_[(Pair-1)*Stride+r];
            inversePairs<Pair-1>(l,r);
        }
    }
public:
    // Deliberately leave storage uninitialized; fill() assigns every entry.
    Permutation() noexcept {}
    void fill(kernelWords& words) {fillPairs(words);}
    template<unsigned Pair> void inversePair(unsigned& l,unsigned& r) const {
        static_assert(Pair>0 && Pair<=Rounds/2);
        r^=table_[(Pair-1)*Stride+NR+l];
        l^=table_[(Pair-1)*Stride+r];
    }
    template<unsigned Pair=Rounds/2> void inverseBatchPairs(std::array<unsigned,16>& l,std::array<unsigned,16>& r) const {
        if constexpr(Pair>0) {
            for(unsigned i=0;i<16;++i)inversePair<Pair>(l[i],r[i]);
            inverseBatchPairs<Pair-1>(l,r);
        }
    }
    void inverseBatch(unsigned first,unsigned* out) const {
        std::array<unsigned,16> l,r;
        for(unsigned i=0;i<16;++i){l[i]=(first+i)&(NL-1);r[i]=(first+i)>>Left;}
        inverseBatchPairs(l,r);
        for(unsigned i=0;i<16;++i)out[i]=l[i]|(r[i]<<Left);
    }
    unsigned forward(unsigned x) const {
        unsigned l=x&(NL-1),r=x>>Left;forwardPairs(l,r);return l|(r<<Left);
    }
    unsigned inverse(unsigned x) const {
        unsigned l=x&(NL-1),r=x>>Left;inversePairs(l,r);return l|(r<<Left);
    }
};

// Each row and region receives its own round tables, without a shared bank.
// Two tagged streams separate the row and region schedules. This prototype
// retains the baseline's fast word generator; there is no secret-key claim.
template<unsigned RegionBits,unsigned Rounds> class Routing {
    static_assert(RegionBits>=6 && RegionBits<=23);
public:
    static constexpr unsigned Rows=1U<<RegionBits;
    static constexpr unsigned Size=256*Rows;
    using Row=Permutation<8,Rounds>;
    using Region=Permutation<RegionBits,Rounds>;
private:
    std::vector<Row> rows_;
    std::vector<Region> regions_;
    template<unsigned Pair=Rounds/2> void rowBatchPairs(const std::array<unsigned,16>& row,std::array<unsigned,16>& l,std::array<unsigned,16>& r) const {
        if constexpr(Pair>0) {
            for(unsigned i=0;i<16;++i)rows_[row[i]].template inversePair<Pair>(l[i],r[i]);
            rowBatchPairs<Pair-1>(row,l,r);
        }
    }
public:
    explicit Routing(std::uint64_t seed):rows_(Rows),regions_(256) {
        kernelWords rowWords(seed^0x726f772d66656973ULL),regionWords(seed^0x7265672d66656973ULL);
        for(auto& p:rows_)p.fill(rowWords);
        for(auto& p:regions_)p.fill(regionWords);
    }
    std::size_t bytes() const {return rows_.size()*sizeof(Row)+regions_.size()*sizeof(Region);}
    unsigned outer(unsigned inner) const {
        const auto region=inner>>RegionBits,position=inner&(Rows-1);
        const auto row=regions_[region].inverse(position);
        return 256*row+rows_[row].inverse(region);
    }
    unsigned inner(unsigned outer) const {
        const auto row=outer>>8,coordinate=outer&255;
        const auto region=rows_[row].forward(coordinate);
        return region*Rows+regions_[region].forward(row);
    }
    // first is 16-aligned, so this batch stays within one region.
    void outerBatch(unsigned first,unsigned* out) const {
        const unsigned region=first>>RegionBits;
        std::array<unsigned,16> row,l,r;
        regions_[region].inverseBatch(first&(Rows-1),row.data());
        l.fill(region&15);r.fill(region>>4);
        rowBatchPairs(row,l,r);
        for(unsigned i=0;i<16;++i)out[i]=(row[i]<<8)|l[i]|(r[i]<<4);
    }
    // Precompute only the tiny row permutations; keep the region round tables.
    // This separates setup speed from any change to the encoding kernels.
    std::vector<std::uint32_t> materialize() const {
        std::vector<std::uint8_t> coordinate(Size);
        for(unsigned row=0;row<Rows;++row)
            for(unsigned region=0;region<256;++region)
                coordinate[256*row+region]=std::uint8_t(rows_[row].inverse(region));
        std::vector<std::uint32_t> route(Size);
        for(unsigned region=0;region<256;++region)
            for(unsigned row=0;row<Rows;++row)
                route[region*Rows+regions_[region].forward(row)]=256*row+coordinate[256*row+region];
        return route;
    }
};
}
