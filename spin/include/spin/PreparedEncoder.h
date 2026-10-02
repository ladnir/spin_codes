#pragma once
#include "Generic.h"

namespace spin {
enum class SetupMode : std::uint32_t { Full=0, BankedHeuristic=1 };
struct CodeSeed {
    std::uint64_t route=1, inner=2;
    bool operator==(const CodeSeed&) const=default;
};
struct SetupOptions {
    SetupMode mode=SetupMode::Full;
    std::uint64_t bank_seed=0x5350494e42414e4bULL;
};
namespace detail { struct PreparedState; struct PreparedScratch; }

// Mutable, move-only setup controller for transposed encoding. Code remains the
// immutable forward/transpose API. Refresh requires exclusive access, including
// all generic views; encoding is concurrent with distinct workspaces only.
class PreparedEncoder {
public:
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
        friend class PreparedEncoder;
        Workspace(std::shared_ptr<detail::PreparedState>,Width,MemoryPolicy);
        std::unique_ptr<detail::PreparedScratch> scratch_;
    };
    explicit PreparedEncoder(CodeSpec,SetupOptions={},ExecutionOptions={});
    ~PreparedEncoder();
    PreparedEncoder(PreparedEncoder&&) noexcept;
    PreparedEncoder& operator=(PreparedEncoder&&) noexcept;
    PreparedEncoder(const PreparedEncoder&)=delete;
    PreparedEncoder& operator=(const PreparedEncoder&)=delete;
    // Same seed is a no-op. Full replaces sampled setup; BankedHeuristic updates
    // wrappers/masks in place without allocation or rebuilding the bank.
    void setCodeSeed(CodeSeed);
    CodeSpec specification() const noexcept;
    SetupOptions setup_options() const noexcept;
    std::size_t message_size() const noexcept;
    std::size_t code_size() const noexcept;
    std::size_t setup_bytes() const noexcept;
    Backend backend() const noexcept;
    bool supports_transpose(Width=Width::Bits128) const noexcept;
    bool supports_generic_transpose() const noexcept;
    // Version 2 includes setup mode, bank seed, and bank-family revision.
    // Execution backend and workspace layout do not identify the binary map.
    std::array<std::byte,56> descriptor() const noexcept;
    Workspace make_workspace(Width=Width::Bits128,MemoryPolicy=MemoryPolicy::Automatic) const;
    Buffer make_buffer(Width=Width::Bits128,MemoryPolicy=MemoryPolicy::Automatic) const;
    // After setCodeSeed, rebind Full scratch and release its previous setup
    // outside a critical path. Compatible scratch does not allocate. Encoding
    // also prepares automatically for existing callers. Banked scratch is unchanged.
    void prepare_workspace(Workspace&) const;
    // Live view where supported: follows seed updates; scratch remains valid.
    GenericTranspose generic_transpose() const;
    bool owns(const GenericTranspose&) const noexcept;
    void transpose_bytes(std::span<const std::byte>,std::span<std::byte>,Workspace&) const;
    void transpose_inplace_bytes(std::span<std::byte>,Workspace&) const;
    template<class E> void transpose(std::span<const E> in,std::span<E> out,Workspace& w) const {
        static_assert(std::is_trivially_copyable_v<E> && sizeof(E)==16);
        if(in.size()!=code_size() || out.size()!=message_size())
            throw std::invalid_argument("SPIN prepared transpose element count mismatch");
        transpose_bytes(std::as_bytes(in),std::as_writable_bytes(out),w);
    }
    template<class E> void transpose_inplace(std::span<E> x,Workspace& w) const {
        static_assert(std::is_trivially_copyable_v<E> && sizeof(E)==16);
        if(x.size()!=code_size())
            throw std::invalid_argument("SPIN prepared inplace element count mismatch");
        transpose_inplace_bytes(std::as_writable_bytes(x),w);
    }
private:
    std::shared_ptr<detail::PreparedState> state_;
};
}
