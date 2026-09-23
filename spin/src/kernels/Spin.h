#include "LengthGeometry.h"
#pragma once
#include "Block.h"
#include <cstddef>
#include <cstdint>
#include <vector>

namespace spin::detail::kernel {
using spin::detail::storage::block;
using u64=std::uint64_t;
using u32=std::uint32_t;
using u8=std::uint8_t;
enum class Configuration { T64S16, T64S20, T128S19, T256S14, T64S12, T64S12R2 };
enum class BchBackend { Auto, Avx2, Avx512 };
bool bchAvx512Available() noexcept;
struct MessageLength {std::size_t value;};
enum class Layout { Packed24, Indices32, Auto };

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
         u64 coefficientSeed=2, unsigned tileRows=0, BchBackend backend=BchBackend::Auto);
    Spin(Configuration configuration, MessageLength length, u64 routeSeed=1,
         u64 coefficientSeed=2, unsigned tileRows=0, BchBackend backend=BchBackend::Auto,
         bool compactSetup=false, std::vector<u32> preparedRoute={});
    void encode(const block* input, std::size_t inputCount, block* output,
                std::size_t outputCount, Workspace& workspace, Layout layout=Layout::Auto) const;
    void encodeUnchecked(const block* input, block* output, Workspace& workspace,
                         Layout layout=Layout::Auto) const;
    // Independent dense oracle; intentionally not a performance path.
    void reference(const block* input, block* output) const;
    // Expanding encoder E: K input blocks to 2K output blocks. Existing encode is E^T.
    void forward(const block* input, std::size_t inputCount, block* output,
                 std::size_t outputCount, Workspace& workspace, Layout layout=Layout::Auto) const;
    void forwardUnchecked(const block* input, block* output, Workspace& workspace,
                          Layout layout=Layout::Auto) const;
    void forwardReference(const block* input, block* output) const;
    // One binary row, bit-packed across coordinates. Scratch is codeBlocks()/64 words.
    void forwardBits(const std::uint64_t* input, std::size_t inputWords,
                     std::uint64_t* output, std::size_t outputWords,
                     std::uint64_t* scratch, std::size_t scratchWords) const;
    void encodeInplace(block*,std::size_t,Workspace&,Layout=Layout::Packed24) const;
    BchBackend bchBackend() const noexcept {return mBchBackend;}
    void validateSetup() const;
    // Discard the dense oracle and the unused routing representation after validation.
    // The chosen layout remains fixed; calling compact again with another layout fails.
    void compact(Layout layout=Layout::Auto);
    unsigned step() const noexcept;
    unsigned state() const noexcept;
    std::size_t setupBytes() const noexcept;
    static constexpr std::size_t maxMessageBlocks=length_geometry::maxMessageBlocks;
    bool packed24Available() const noexcept {return codeBlocks()<=(std::size_t{1}<<24);}
    Layout preferredLayout() const noexcept {
        return mCompacted?mRetainedLayout:(packed24Available()?Layout::Packed24:Layout::Indices32);
    }
    std::size_t messageBlocks() const noexcept { return mK; }
    std::size_t codeBlocks() const noexcept { return 2*mK; }
    std::size_t tileBlocks() const noexcept { return 256*mTileRows; }
    const char* name() const noexcept;
    u64 routeHash() const noexcept;
    struct WideView {
        std::size_t n, tile;
        const u8* slots;
        const u8* offsets;
        const u32* fieldRows; // IMT: (u,v) pairs, epoch then round.
        Configuration configuration=Configuration::T128S19;
        const u32* slots32=nullptr;
        const u32* offsets32=nullptr;
    };
    // Read-only schedules for the wide S19 implementations of the same map.
    WideView wideView() const;
    WideView wideForwardView() const;
private:
    friend struct Access;
    Configuration mConfig;
    std::size_t mK;
    unsigned mTileRows;
    bool mPartialTile=false;
    BchBackend mBchBackend=BchBackend::Avx2;
    std::vector<u8> mBchOffsets24;
    std::vector<u32> mBchOffsets32;
    bool mCompacted=false;
    Layout mRetainedLayout=Layout::Packed24;
    std::vector<u8> mSlots24, mOffsets24;
    std::vector<u32> mSlots32, mOffsets32, mRoute, mCoefficients, mFieldRows;
    std::vector<u32> mForwardFieldRows;
    template<bool Packed> void runK16R2Four(const block*,block*,Workspace&) const;
    template<bool Packed> void runK16R2FourTail(const block*,block*,Workspace&) const;
    template<bool Packed> void runK16Four(const block*,block*,Workspace&) const;
    template<bool Packed> void runK16FourTail(const block*,block*,Workspace&) const;
    template<bool Packed> void runFour(const block*,block*,Workspace&) const;
    template<bool Packed> void runFourTail(const block*,block*,Workspace&) const;
    std::vector<u32> mSmallRoute32;
    void runSmallK16R2(const block*,block*,Workspace&) const;
    void runSmallK16(const block*,block*,Workspace&) const;
    void runSmall(const block*,block*,Workspace&) const;
    template<class Map> void forwardBitsMap(const u64*,u64*,u64*) const;
    std::vector<u32> mK18Route32;
    void runK18(const block*,block*,Workspace&) const;
    template<class Map,bool Packed> void runTail(const block*,block*,Workspace&) const;
    template<class Map,bool Packed> void runForwardTail(const block*,block*,Workspace&) const;
    template<class Map> void setupInner(u64 seed);
    template<class Map, bool Packed> void run(const block*,block*,Workspace&) const;
    template<class Map> void oracle(const block*,block*) const;
    std::vector<u32> mRangeRoute32;
    template<class Map> void runRange(const block*,block*,Workspace&) const;
    std::vector<u32> mForwardDirect;
    template<class Map> void runForwardDirect(const block*,block*,Workspace&) const;
    template<class Map, bool Packed> void runForwardFour(const block*,block*,Workspace&) const;
    template<class Map, bool Packed> void runForwardFourTail(const block*,block*,Workspace&) const;
    template<class Map, bool Packed> void runForward(const block*,block*,Workspace&) const;
    template<class Map> void forwardOracle(const block*,block*) const;
};
u64 splitmix(u64& state) noexcept;
}
