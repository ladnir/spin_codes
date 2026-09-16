#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <stdexcept>

struct Mask { std::uint64_t lo=0, hi=0; };

// Each completed group is a disjoint basis. Keep original columns for a
// separate low-pivot rank check; discovery uses high pivots.
unsigned bases(const std::array<std::uint32_t,128>& columns, unsigned t,
               unsigned s, Mask selected, bool reverse) {
    std::uint32_t pivots[128][19]{};
    std::uint32_t originals[128][19]{};
    unsigned ranks[128]{};
    unsigned complete=0;
    const unsigned groups=t/s;
    for (unsigned j=0;j<t;++j) {
        const unsigned p=reverse?t-1-j:j;
        if (!((p<64?selected.lo:selected.hi)>>(p%64)&1)) continue;
        for (unsigned g=0;g<groups;++g) {
            if (ranks[g]==s) continue;
            auto v=columns[p];
            while (v) {
                const auto bit=31u-std::countl_zero(v);
                if (pivots[g][bit]) v^=pivots[g][bit];
                else {
                    pivots[g][bit]=v;
                    originals[g][ranks[g]++]=columns[p];
                    if (ranks[g]==s) {
                        std::uint32_t low[19]{};
                        for (unsigned k=0;k<s;++k) {
                            auto x=originals[g][k];
                            while (x) {
                                const auto bit2=std::countr_zero(x);
                                if (low[bit2]) x^=low[bit2];
                                else { low[bit2]=x; break; }
                            }
                            if (!x) throw std::runtime_error("dependent saved basis");
                        }
                        ++complete;
                    }
                    break;
                }
            }
            if (v) break; // This column was used once; groups stay disjoint.
        }
    }
    return complete;
}

int main() {
    unsigned t,s;
    if (!(std::cin>>t>>s) || !s || s>19 || t<s || t>128) return 2;
    std::array<std::uint32_t,128> a{},b{};
    for (unsigned i=0;i<t;++i) if (!(std::cin>>a[i]) || a[i]>=(1u<<s)) return 2;
    for (unsigned i=0;i<t;++i) if (!(std::cin>>b[i]) || b[i]>=(1u<<s)) return 2;
    std::array<Mask,19> generators{};
    for (unsigned j=0;j<s;++j) for (unsigned i=0;i<t;++i) if (b[i]>>j&1) {
        if (i<64) generators[j].lo|=std::uint64_t(1)<<i;
        else generators[j].hi|=std::uint64_t(1)<<(i-64);
    }
    const Mask full{t>=64?~std::uint64_t(0):(std::uint64_t(1)<<t)-1,
                    t==128?~std::uint64_t(0):t>64?(std::uint64_t(1)<<(t-64))-1:0};
    // At most floor(128/s) bases. Fixed storage avoids per-state allocations.
    static std::uint32_t hist[129][129][129]{};
    Mask word;
    for (std::uint32_t i=1;i<(1u<<s);++i) {
        const auto& changed=generators[std::countr_zero(i)];
        word.lo^=changed.lo; word.hi^=changed.hi;
        const Mask other{word.lo^full.lo,word.hi^full.hi};
        const auto d=std::max(bases(a,t,s,word,false),bases(a,t,s,word,true));
        const auto e=std::max(bases(a,t,s,other,false),bases(a,t,s,other,true));
        ++hist[std::popcount(word.lo)+std::popcount(word.hi)][d][e];
    }
    std::cout<<"{\"states\":"<<((1u<<s)-1)<<",\"groups\":[";
    bool first=true;
    for (unsigned w=0;w<=t;++w) for (unsigned d=0;d<=t/s;++d) for (unsigned e=0;e<=t/s;++e)
        if (hist[w][d][e]) {
            if (!first) std::cout<<',';
            first=false;
            std::cout<<'['<<w<<','<<d<<','<<e<<','<<hist[w][d][e]<<']';
        }
    std::cout<<"]}\n";
}
