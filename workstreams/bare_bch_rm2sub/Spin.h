#pragma once
#include <cryptoTools/Common/Defines.h>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace bare_spin {
using osuCrypto::block;
using u64=std::uint64_t;
using u32=std::uint32_t;
using u8=std::uint8_t;
enum class Configuration { T64S16, T64S20, T128S19, T256S14 };
enum class Layout { Packed24, Indices32 };

// Immutable after setup. Configuration dispatch happens once per encode.
class Spin {
public:
    struct Workspace {
        std::vector<block> buckets, tile;
        Workspace() = default;
        explicit Workspace(const Spin& code);
        std::size_t bytes() const noexcept;
    };
    Spin(Configuration configuration, unsigned messageExponent, u64 routeSeed=1,
         u64 coefficientSeed=2, unsigned tileRows=0);
    void encode(const block* input, std::size_t inputCount, block* output,
                std::size_t outputCount, Workspace& workspace, Layout layout=Layout::Packed24) const;
    void encodeUnchecked(const block* input, block* output, Workspace& workspace,
                         Layout layout=Layout::Packed24) const;
    // Independent dense oracle; intentionally not a performance path.
    void reference(const block* input, block* output) const;
    void validateSetup() const;
    // Discard the dense oracle and the unused routing representation after validation.
    // The chosen layout remains fixed; calling compact again with another layout fails.
    void compact(Layout layout=Layout::Packed24);
    unsigned step() const noexcept;
    unsigned state() const noexcept;
    std::size_t setupBytes() const noexcept;
    std::size_t messageBlocks() const noexcept { return mK; }
    std::size_t codeBlocks() const noexcept { return 2*mK; }
    std::size_t tileBlocks() const noexcept { return 256*mTileRows; }
    const char* name() const noexcept;
    u64 routeHash() const noexcept;
private:
    Configuration mConfig;
    std::size_t mK;
    unsigned mTileRows;
    bool mCompacted=false;
    Layout mRetainedLayout=Layout::Packed24;
    std::vector<u8> mSlots24, mOffsets24;
    std::vector<u32> mSlots32, mOffsets32, mRoute, mCoefficients, mFieldRows;
    template<class Map> void setupInner(u64 seed);
    template<class Map, bool Packed> void run(const block*,block*,Workspace&) const;
    template<class Map> void oracle(const block*,block*) const;
};
u64 splitmix(u64& state) noexcept;
}
