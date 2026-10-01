#include <spin/Code.h>
#include <spin/Generic.h>
#include <array>
#include <vector>
// Model a consumer-owned type with the name formerly leaked by the import.
namespace osuCrypto { struct alignas(16) block { std::array<unsigned long long,2> words{}; }; }
int main() {
    for(auto family:{spin::Parameters::T128S19,spin::Parameters::T64S12,
                    spin::Parameters::T64S12R2,spin::Parameters::PacketT64S16}) {
        if(family!=spin::Parameters::PacketT64S16 && !spin::capabilities().avx2)continue;
        const auto k=3*spin::message_alignment(family);
        if(!spin::valid_message_size(family,k))return 1;
        spin::Code code({k,family,17,17});
        if(!code.supports_transpose())return 1;
        spin::Code::Workspace work=code.make_workspace();
        std::vector<osuCrypto::block> input(code.code_size()),output(code.message_size());
        code.transpose<osuCrypto::block>(input,output,work);
        code.transpose_inplace<osuCrypto::block>(input,work);
        auto owned=code.make_buffer();code.transpose_inplace_bytes(owned.bytes(),work);
        if(code.supports_forward())code.forward<osuCrypto::block>(output,input,work);
        if(code.supports_generic_transpose()) {
            auto generic=code.generic_transpose();auto bytes=generic.make_workspace<unsigned char>();
            std::vector<unsigned char> bits(code.code_size());
            generic.transpose_inplace<unsigned char>(bits,bytes);
        }
    }
}
