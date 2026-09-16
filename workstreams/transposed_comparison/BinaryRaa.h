#pragma once
#include "libOTe/Tools/Ppcg/Permutation.h"
#include <array>
#include <vector>

namespace comparison {
using namespace osuCrypto;

// Binary transpose R Pi_0 Acc Pi_1 Acc. Each contiguous group of R
// coordinates is one repetition block. No field multiplications or diagonals.
// Eight-way scan/scatter follows the batching of libOTe's chosen-block BAA.
template<unsigned R> class BinaryRaa {
    static_assert(R == 3 || R == 4);
    u64 k_, n_;
    std::array<Feistel2KPerm, 2> pi_;
    std::vector<block> mid_, state_;

    void stage(const block* in, block* out, unsigned s) const {
        block sum = ZeroBlock;
        u64 i = 0;
        for (; i + 8 <= n_; i += 8) {
            const block v0 = sum ^ in[i];
            const block v1 = v0 ^ in[i+1];
            const block v2 = v1 ^ in[i+2];
            const block v3 = v2 ^ in[i+3];
            const block v4 = v3 ^ in[i+4];
            const block v5 = v4 ^ in[i+5];
            const block v6 = v5 ^ in[i+6];
            const block v7 = v6 ^ in[i+7];
            sum = v7;
            const block p = pi_[s].feistelBijection(block(_mm_set_epi32(i+3,i+2,i+1,i)));
            const block q = pi_[s].feistelBijection(block(_mm_set_epi32(i+7,i+6,i+5,i+4)));
            out[p.get<u32>(0)] = v0; out[p.get<u32>(1)] = v1;
            out[p.get<u32>(2)] = v2; out[p.get<u32>(3)] = v3;
            out[q.get<u32>(0)] = v4; out[q.get<u32>(1)] = v5;
            out[q.get<u32>(2)] = v6; out[q.get<u32>(3)] = v7;
        }
        for (; i < n_; ++i) { sum ^= in[i]; out[pi_[s].feistelBijection(u32(i))] = sum; }
    }
public:
    explicit BinaryRaa(u64 k) : k_(k), n_(k*R), mid_(n_), state_(n_) {
        if (!k || k > (1ull<<32)/R) throw std::invalid_argument("RAA size out of range");
        PRNG master(block(34123421, 2134123));
        pi_[0].init(n_, master.get<block>(), 4);
        pi_[1].init(n_, master.get<block>(), 4);
    }
    void encode(const block* in, block* out) {
        stage(in, mid_.data(), 1);
        stage(mid_.data(), state_.data(), 0);
        for (u64 i=0; i<k_; ++i) {
            auto x = state_[R*i] ^ state_[R*i+1] ^ state_[R*i+2];
            if constexpr (R == 4) x ^= state_[R*i+3];
            out[i] = x;
        }
    }
    // Forward map, for an independent bilinear transpose test. Acc^T is
    // a suffix scan and a transposed permutation gathers instead of scatters.
    void forward(const block* in, block* out) const {
        std::vector<block> current(n_), next(n_);
        for (u64 i=0; i<k_; ++i) for (unsigned j=0; j<R; ++j) current[R*i+j]=in[i];
        for (unsigned s=0; s<2; ++s) {
            block sum = ZeroBlock;
            for (u64 i=n_; i--;) {
                sum ^= current[pi_[s].feistelBijection(u32(i))];
                next[i] = sum;
            }
            current.swap(next);
        }
        std::copy(current.begin(), current.end(), out);
    }
};
}
