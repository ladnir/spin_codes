"""Instrument a separate control binary; never use its timings as final benchmarks."""
from pathlib import Path
PARENT=Path(__file__).resolve().parent.parent
source=PARENT/'generated/asymmetric_greedy3_2_sparse'
code=(source/'CandidateSpin.cpp').read_text()
code=code.replace('#include <algorithm>','#include <algorithm>\n#include <chrono>\n#include <cstdio>')
code=code.replace('namespace bare_spin {','''namespace bare_spin {
namespace {
using Clock=std::chrono::steady_clock;
struct Profile {
    double inner=0,scatter=0,outer=0;unsigned calls=0;
    ~Profile(){if(calls) std::fprintf(stderr,"PROFILE calls=%u inner_route_ms=%.6f scatter_ms=%.6f outer_ms=%.6f\\n",calls,inner/calls,scatter/calls,outer/calls);}
} profile;
double millis(Clock::time_point start) {return std::chrono::duration<double,std::milli>(Clock::now()-start).count();}
}
''',1)
code=code.replace('    block* values=w.buckets.data();','    const auto profileStart=Clock::now();\n    block* values=w.buckets.data();',1)
code=code.replace('    block* tile=w.tile.data();','    profile.inner+=millis(profileStart);++profile.calls;\n    block* tile=w.tile.data();',1)
code=code.replace('        for(std::size_t j=0;j<tileSize;++j) {','        const auto scatterStart=Clock::now();\n        for(std::size_t j=0;j<tileSize;++j) {',1)
code=code.replace('        if constexpr(Quarter) {\n            for(std::size_t j=0;j<tileSize;j+=256)','        profile.scatter+=millis(scatterStart);\n        const auto outerStart=Clock::now();\n        if constexpr(Quarter) {\n            for(std::size_t j=0;j<tileSize;j+=256)',1)
code=code.replace('    }\n}\nvoid Spin::encodeUnchecked','        profile.outer+=millis(outerStart);\n    }\n}\nvoid Spin::encodeUnchecked',1)
dest=PARENT/'generated/routeopt_profile'
dest.mkdir(exist_ok=True)
(dest/'CandidateSpin.cpp').write_text(code,encoding='utf-8',newline='\n')
(dest/'AsymmetricMap.h').write_text((source/'AsymmetricMap.h').read_text(),encoding='utf-8',newline='\n')
