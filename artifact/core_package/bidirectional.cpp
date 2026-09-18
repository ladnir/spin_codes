#include "Spin.h"
#include <algorithm>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace bare_spin;
static void require(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}
static bool equal(const std::vector<block>& a, const std::vector<block>& b) {
    return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size()*sizeof(block)) == 0;
}
static __m128i product(const std::vector<block>& a, const std::vector<block>& b) {
    auto result = _mm_setzero_si128();
    for (std::size_t i=0; i<a.size(); ++i)
        result = _mm_xor_si128(result, _mm_and_si128(a[i].mData, b[i].mData));
    return result;
}
int main() {
    try {
        Spin code(Configuration::T128S19, 16, 1, 2);
        code.validateSetup();
        Spin::Workspace work(code);
        std::vector<block> x(code.messageBlocks()), q(code.codeBlocks());
        std::vector<block> y(q.size()), expectedY(q.size()), z(x.size()), expectedZ(x.size());
        u64 seed=123;
        for (auto* values : {&x, &q}) for (auto& v : *values) {
            const auto lo=splitmix(seed), hi=splitmix(seed);
            v=block(hi,lo);
        }
        code.forwardReference(x.data(), expectedY.data());
        code.reference(q.data(), expectedZ.data());
        for (auto layout : {Layout::Packed24, Layout::Indices32}) {
            code.forward(x.data(), x.size(), y.data(), y.size(), work, layout);
            code.encode(q.data(), q.size(), z.data(), z.size(), work, layout);
            require(equal(y, expectedY), "forward dense-oracle mismatch");
            require(equal(z, expectedZ), "transpose dense-oracle mismatch");
            require(_mm_movemask_epi8(_mm_cmpeq_epi8(product(y,q), product(x,z)))==0xffff,
                    "adjoint identity mismatch");
        }
        code.compact();
        code.forward(x.data(), x.size(), y.data(), y.size(), work);
        code.encode(q.data(), q.size(), z.data(), z.size(), work);
        require(equal(y, expectedY) && equal(z, expectedZ), "compaction changed the map");
        std::cout << "forward, transpose, layouts, adjoint, and compaction passed\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
