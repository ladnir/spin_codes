#pragma once
#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <stdexcept>
#include <type_traits>

namespace spin {
// Stable identifiers for supplied maps, not arbitrary (t,s) synthesis.
// Value 4 identified the retired BCH packet construction; it is not reused.
enum class Parameters : std::uint32_t { T128S19=1, T64S12=2, T64S12R2=3, PacketRsT64S20=5 };
enum class Backend { Automatic, Avx2, Avx512, Portable };
enum class Width : unsigned { Bits128=16, Bits256=32, Bits512=64 };
enum class MemoryPolicy { Normal, PreferHugePages, Automatic };
struct CodeSpec {
    std::size_t message_size;
    Parameters parameters=Parameters::T128S19;
    // For PacketRsT64S20, route_seed also generates the GF(2^32) outer maps;
    // inner_seed generates GL20 updates. Equal seeds reproduce PacketCode(seed).
    std::uint64_t route_seed=1;
    std::uint64_t inner_seed=2;
};
struct ExecutionOptions {
    Backend backend=Backend::Automatic;
    // 0 selects the family default: 256 for IMT, 8 for PacketRsT64S20.
    // PacketRsT64S20 accepts only 0 or 8; its eight-row geometry is fixed.
    unsigned tile_rows=0;
};
struct Capabilities { bool avx2, avx512_bch, forward256, forward512; };
Capabilities capabilities() noexcept;
std::size_t message_alignment(Parameters);
bool valid_message_size(Parameters, std::size_t) noexcept;

namespace detail { struct Plan; struct Scratch; struct BufferAllocation; }
class GenericTranspose;
class Code;

// Optional zero-initialized owned storage, independent of the selected family.
// Automatic (default): PreferHugePages for allocations >=16 MiB on supported
// platforms (currently Linux), Normal otherwise. Explicit overrides stay fixed.
// Normal: 64-byte alignment. PreferHugePages: 2 MiB alignment/rounding and
// best-effort Linux page advice during allocation, never during encoding.
class Buffer {
public:
    explicit Buffer(std::size_t bytes,MemoryPolicy=MemoryPolicy::Automatic);
    ~Buffer();
    Buffer(Buffer&&) noexcept;
    Buffer& operator=(Buffer&&) noexcept;
    Buffer(const Buffer&)=delete;
    Buffer& operator=(const Buffer&)=delete;
    std::span<std::byte> bytes() noexcept;
    std::span<const std::byte> bytes() const noexcept;
    std::size_t allocation_bytes() const noexcept;
private:
    std::unique_ptr<detail::BufferAllocation> allocation_;
};

// Move-only scratch. Retains its plan; moving/destroying Code handles is safe.
// One workspace per concurrent call. No allocation occurs during encoding.
class Workspace {
public:
    ~Workspace();
    Workspace(Workspace&&) noexcept;
    Workspace& operator=(Workspace&&) noexcept;
    Workspace(const Workspace&)=delete;
    Workspace& operator=(const Workspace&)=delete;
    std::size_t bytes() const noexcept;
    Width width() const noexcept;
private:
    friend class Code;
    Workspace(std::shared_ptr<const detail::Plan>, Width, MemoryPolicy);
    std::unique_ptr<detail::Scratch> scratch_;
};

// Immutable shared setup. Copying a handle shares the plan, not its scratch.
// Supplied half-rate codes: forward K -> 2K; transpose 2K -> K.
class Code {
public:
    using Workspace=spin::Workspace;
    explicit Code(CodeSpec, ExecutionOptions={});
    std::size_t message_size() const noexcept;
    std::size_t code_size() const noexcept;
    CodeSpec specification() const noexcept;
    Backend backend() const noexcept;
    std::size_t setup_bytes() const noexcept;
    // Canonical little-endian descriptor, version 1. Includes map and both seeds;
    // excludes backend, tile size, and route packing. Consumers choose the hash.
    std::array<std::byte,40> descriptor() const noexcept;
    bool supports_forward(Width=Width::Bits128) const noexcept;
    bool supports_transpose(Width=Width::Bits128) const noexcept;
    bool supports_generic_transpose() const noexcept;
    Workspace make_workspace(Width=Width::Bits128,MemoryPolicy=MemoryPolicy::Automatic) const;
    // Optional storage for exactly code_size() records at the selected width.
    Buffer make_buffer(Width=Width::Bits128,MemoryPolicy=MemoryPolicy::Automatic) const;
    // Rebind compatible 128-bit scratch without clearing or allocating; otherwise
    // recreate it at the same width (128 bits if moved-from). No concurrent use.
    void prepare_workspace(Workspace&) const;
    GenericTranspose generic_transpose() const;

    // Byte views are the zero-copy foreign-buffer interface. Records are 16/32/64
    // bytes according to the workspace; buffers require 16-byte alignment.
    // Exact sizes only, no overlap. Wider records apply the same map lane-wise.
    // PacketRsT64S20 supports 128-bit forward/transpose. A 64-byte-aligned
    // forward output (as supplied by Buffer) enables its fastest large-size path.
    void forward_bytes(std::span<const std::byte>, std::span<std::byte>, Workspace&) const;
    void transpose_bytes(std::span<const std::byte>, std::span<std::byte>, Workspace&) const;
    // Exactly 2K 128-bit records; overwrite the first K and preserve the suffix.
    void transpose_inplace_bytes(std::span<std::byte>, Workspace&) const;

    template<class E> void forward(std::span<const E> in, std::span<E> out, Workspace& w) const {
        check_record<E>(w);
        if(in.size()!=message_size() || out.size()!=code_size())
            throw std::invalid_argument("SPIN forward element count mismatch");
        forward_bytes(std::as_bytes(in), std::as_writable_bytes(out), w);
    }
    template<class E> void transpose(std::span<const E> in, std::span<E> out, Workspace& w) const {
        check_record<E>(w);
        if(in.size()!=code_size() || out.size()!=message_size())
            throw std::invalid_argument("SPIN transpose element count mismatch");
        transpose_bytes(std::as_bytes(in), std::as_writable_bytes(out), w);
    }
    template<class E> void transpose_inplace(std::span<E> buffer, Workspace& w) const {
        check_record<E>(w);
        if(buffer.size()!=code_size())throw std::invalid_argument("SPIN inplace element count mismatch");
        transpose_inplace_bytes(std::as_writable_bytes(buffer), w);
    }
    // One bit-packed binary message: word bit j is coordinate 64*i+j.
    // Caller scratch is exactly code_size()/64 words. All three ranges disjoint.
    void forward_bits(std::span<const std::uint64_t>, std::span<std::uint64_t>,
                      std::span<std::uint64_t> scratch) const;
private:
    std::shared_ptr<const detail::Plan> plan_;
    template<class E> static void check_record(const Workspace& w) {
        static_assert(std::is_trivially_copyable_v<E>);
        static_assert(sizeof(E)==16 || sizeof(E)==32 || sizeof(E)==64);
        if(sizeof(E)!=static_cast<unsigned>(w.width()))
            throw std::invalid_argument("SPIN record width does not match workspace");
    }
    detail::Scratch& check_workspace(Workspace&) const;
};
}
