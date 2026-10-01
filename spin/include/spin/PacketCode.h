#pragma once
#include "Code.h"
#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <stdexcept>
#include <type_traits>

namespace spin {
// Compatibility interface for four-row BCH packets with GL32 outer mixing and
// the t64/s16 inner with GL16 updates. New consumers use Code with PacketT64S16.
struct PacketSpec {
    std::size_t message_size;
    std::uint64_t seed=1;
    bool operator==(const PacketSpec&) const=default;
};
enum class PacketBackend { Automatic, Portable, Avx512Gfni };
using PacketMemory=MemoryPolicy;
using PacketBuffer=Buffer;
bool packet_fast_available() noexcept;
// Natural message-size unit: four BCH rows, each carrying 128 elements.
// Transposed buffers hold exactly 2K input elements and K output elements.
inline constexpr std::size_t packet_message_alignment() noexcept { return 512; }
// No allocation, setup, rounding, or certificate lookup. Includes representation
// bounds; a valid size can still exceed the caller's available memory.
bool valid_packet_message_size(std::size_t) noexcept;
namespace detail { struct PacketState; struct PacketScratch; }

class PacketCode {
public:
    class Workspace {
    public:
        ~Workspace();
        Workspace(Workspace&&) noexcept;
        Workspace& operator=(Workspace&&) noexcept;
        Workspace(const Workspace&)=delete;
        Workspace& operator=(const Workspace&)=delete;
        std::size_t bytes() const noexcept;
        Width width() const noexcept { return Width::Bits128; }
    private:
        friend class PacketCode;
        Workspace(std::shared_ptr<const detail::PacketState>, PacketMemory);
        std::unique_ptr<detail::PacketScratch> scratch_;
    };
    // Setup is immutable and shared by copied handles. K is a positive multiple
    // of packet_message_alignment(), not necessarily a power of two. No implicit
    // rounding. There is no distance claim for an arbitrary accepted K or seed.
    explicit PacketCode(PacketSpec, PacketBackend=PacketBackend::Automatic);
    std::size_t message_size() const noexcept;
    std::size_t code_size() const noexcept;
    PacketSpec specification() const noexcept;
    PacketBackend backend() const noexcept;
    std::size_t setup_bytes() const noexcept;
    bool supports_forward(Width=Width::Bits128) const noexcept { return false; }
    bool supports_transpose(Width width=Width::Bits128) const noexcept {
        return message_size()!=0 && width==Width::Bits128;
    }
    bool supports_generic_transpose() const noexcept { return false; }
    // Little-endian SPKP descriptor: version1, family1, reserved0, K, seed.
    // Backend and memory policy never change the represented binary map.
    std::array<std::byte,32> descriptor() const noexcept;
    // Scratch is 64-byte aligned with best-effort advice on owned pages.
    // PreferHugePages additionally requests 2 MiB-aligned backing storage.
    Workspace make_workspace(PacketMemory=PacketMemory::Normal) const;
    // Zero-initialized storage for exactly 2K 128-bit elements.
    PacketBuffer make_buffer(PacketMemory=PacketMemory::Normal) const;
    // No allocation during encoding. One workspace per concurrent call.
    // Exact sizes: 2K inputs, K outputs, 16 bytes per element. No overlap.
    void transpose_bytes(std::span<const std::byte>,std::span<std::byte>,Workspace&) const;
    // Overwrite the first K elements; preserve the last K elements.
    void transpose_inplace_bytes(std::span<std::byte>,Workspace&) const;
    template<class E> void transpose(std::span<const E> in,std::span<E> out,Workspace& w) const {
        static_assert(std::is_trivially_copyable_v<E> && sizeof(E)==16);
        if(in.size()!=code_size() || out.size()!=message_size())
            throw std::invalid_argument("SPIN packet transpose element count mismatch");
        transpose_bytes(std::as_bytes(in),std::as_writable_bytes(out),w);
    }
    template<class E> void transpose_inplace(std::span<E> in,Workspace& w) const {
        static_assert(std::is_trivially_copyable_v<E> && sizeof(E)==16);
        if(in.size()!=code_size()) throw std::invalid_argument("SPIN packet inplace element count mismatch");
        transpose_inplace_bytes(std::as_writable_bytes(in),w);
    }
private:
    std::shared_ptr<const detail::PacketState> state_;
    const detail::PacketState& checked_state() const;
    detail::PacketScratch& checked_workspace(Workspace&) const;
};
}
