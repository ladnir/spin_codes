#include <spin/Code.h>
#include <spin/Generic.h>
#include <array>
#include <vector>
// Model a consumer-owned type with the name formerly leaked by the import.
namespace osuCrypto { struct alignas(16) block { std::array<unsigned long long,2> words{}; }; }
int main() {
    if(!spin::capabilities().avx2)return 0;
    spin::Code code({49152,spin::Parameters::T128S19});
    auto work=code.make_workspace();
    std::vector<osuCrypto::block> input(code.message_size()),output(code.code_size());
    code.forward<osuCrypto::block>(input,output,work);
    code.transpose_inplace<osuCrypto::block>(output,work);
    auto generic=code.generic_transpose();auto bytes=generic.make_workspace<unsigned char>();
    std::vector<unsigned char> bits(code.code_size());
    generic.transpose_inplace<unsigned char>(bits,bytes);
}
