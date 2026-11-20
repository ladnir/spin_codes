#include <spin/Code.h>
#include <spin/PreparedEncoder.h>
#include <array>
#include <vector>
// Model a consumer-owned type with the name formerly leaked by the import.
namespace osuCrypto { struct alignas(16) block { std::array<unsigned long long,2> words{}; }; }
static_assert(static_cast<std::uint32_t>(spin::Parameters::PacketRs8T64)==8);
int main() {
    constexpr auto family=spin::Parameters::PacketRs8T64;
    const auto k=3*spin::message_alignment(family);
    if(k!=384 || !spin::valid_message_size(family,k))return 1;
    spin::Code code({k,family,17,17});
    if(!code.supports_transpose() || !code.supports_forward(spin::Width::Bits512))return 1;
    spin::Code::Workspace work=code.make_workspace();
    std::vector<osuCrypto::block> input(code.code_size()),output(code.message_size());
    code.transpose<osuCrypto::block>(input,output,work);
    code.transpose_inplace<osuCrypto::block>(input,work);
    auto owned=code.make_buffer();code.transpose_inplace_bytes(owned.bytes(),work);
    code.forward<osuCrypto::block>(output,input,work);
    spin::PreparedEncoder prepared({k,family,17,17},{spin::SetupMode::BankedHeuristic});
    auto preparedWork=prepared.make_workspace();
    prepared.transpose_inplace<osuCrypto::block>(input,preparedWork);
    prepared.setCodeSeed({18,19});
    prepared.transpose_inplace<osuCrypto::block>(input,preparedWork);
}
