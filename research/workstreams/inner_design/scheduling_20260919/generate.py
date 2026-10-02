"""Generate isolated, exact-map scheduling candidates; leave pinned sources alone."""
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BASE = ROOT / 'workstreams/inner_design/asymmetric/bch256/weight5/implementation'


def replace_once(text, old, new):
    assert text.count(old) == 1, (old[:100], text.count(old))
    return text.replace(old, new)


def main():
    output = HERE / 'generated'
    output.mkdir(exist_ok=True)
    cpp = (BASE / 'Weight5Spin.cpp').read_text()
    inner = (BASE / 'Weight5Inner.h').read_text()
    inner = replace_once(inner, 'const auto base=epoch*Map::T;',
                         'const auto base=epoch*Map::T;\n'
                         '            if constexpr(requires {emit.prepare(base);}) emit.prepare(base);')
    (output / 'ScheduledInner.h').write_text(inner, newline='\n')
    cpp = replace_once(cpp, '#include "Weight5Inner.h"',
                       '#include "ScheduledInner.h"\n#include "Schedule.h"')
    cpp = replace_once(cpp, 'namespace bare_spin {', '''#if W5_PROFILE
#include <chrono>
#include <cstdio>
namespace {
using Clock=std::chrono::steady_clock;
double ms(Clock::time_point start) {return std::chrono::duration<double,std::milli>(Clock::now()-start).count();}
struct StageProfile {
    double inner=0,scatter=0,bch=0;
    unsigned calls=0;
    ~StageProfile() { if(calls) std::fprintf(stderr,"PROFILE calls=%u inner_ms=%.6f scatter_ms=%.6f bch_ms=%.6f\\n",calls,inner/calls,scatter/calls,bch/calls); }
} profile;
}
#endif
namespace bare_spin {''')
    cpp = replace_once(cpp, '    block* values=w.buckets.data();',
        '    block* values=w.buckets.data();\n#if W5_PROFILE\n    auto stageStart=Clock::now();\n#endif')
    cpp = replace_once(cpp, '    block* tile=w.tile.data();',
        '#if W5_PROFILE\n    profile.inner+=ms(stageStart); ++profile.calls;\n#endif\n'
        '    block* tile=w.tile.data();')
    start = cpp.index('        for(std::size_t j=0;j<tileSize;++j) {', cpp.index('void Spin::run('))
    end = cpp.index('        if constexpr(Quarter)', start)
    old = cpp[start:end]
    cpp = replace_once(cpp, old, '''#if W5_PROFILE
        stageStart=Clock::now();
#endif
#if SCHEDULE_GROUP
        schedule::scatter<Packed,SCHEDULE_GROUP,W5_PREFETCH,W5_WRITE_PREFETCH>(
            values,tile,mOffsets24.data(),mOffsets32.data(),base,tileSize);
#else
''' + old + '''#endif
#if W5_PROFILE
        profile.scatter+=ms(stageStart);stageStart=Clock::now();
#endif
''')
    cpp = replace_once(cpp, '''        }
    }
}
void Spin::encodeUnchecked''', '''        }
#if W5_PROFILE
        profile.bch+=ms(stageStart);
#endif
    }
}
void Spin::encodeUnchecked''')
    old = '''    asymmetricReverse<Map,true,W5_MASKED>(in,codeBlocks(),mFieldRows.data(),[&](std::size_t i,block v) {
        if constexpr(Packed) values[unpack(mSlots24.data()+3*i)]=v;
        else values[mSlots32[i]]=v;
    });'''
    cpp = replace_once(cpp, old, '''#if SCHEDULE_DECODE || SCHEDULE_EMIT_LEAD
    schedule::Emitter<Packed,SCHEDULE_DECODE,SCHEDULE_EMIT_LEAD> emitter{
        values,mSlots24.data(),mSlots32.data()};
    asymmetricReverse<Map,true,W5_MASKED>(in,codeBlocks(),mFieldRows.data(),emitter);
#else
''' + old + '\n#endif')
    (output / 'ScheduledSpin.cpp').write_text(cpp, newline='\n')
    bench = (ROOT / 'workstreams/bare_bch_rm2sub/benchmark.cpp').read_text()
    bench = replace_once(bench, 'const int lock=open("/tmp/bare-spin-benchmark.lock",O_CREAT|O_RDWR,0600);',
        'const int common=open("/tmp/prindal-addition-encoder-benchmark.lock",O_CREAT|O_RDWR,0600);\n'
        '        if(common<0 || flock(common,LOCK_EX|LOCK_NB)) throw std::runtime_error("cross-project benchmark lock unavailable");\n'
        '        const int lock=open("/tmp/bare-spin-benchmark.lock",O_CREAT|O_RDWR,0600);')
    bench = replace_once(bench, 'const unsigned trials=argc>1?',
        'const u64 routeSeed=argc>8?std::stoull(argv[8]):1;\n'
        '        const unsigned trials=argc>1?')
    bench = replace_once(bench, 'm,1,2,tile,quarter?', 'm,routeSeed,2,tile,quarter?')
    bench = replace_once(bench, 'std::sort(samples.begin(),samples.end());',
        'const auto rawSamples=samples;\n                std::sort(samples.begin(),samples.end());')
    bench = replace_once(bench, '<<std::dec<<"\\\"}\\n"<<std::flush;',
        '<<std::dec<<"\\\",\\\"route_seed\\\":"<<routeSeed<<",\\\"samples_ms\\\":[";\n'
        '                for(std::size_t j=0;j<rawSamples.size();++j) std::cout<<(j?",":"")<<rawSamples[j];\n'
        '                std::cout<<"]}\\n"<<std::flush;')
    (HERE / 'benchmark.cpp').write_text(bench, newline='\n')
    inputs = [BASE / 'Weight5Spin.cpp', BASE / 'Weight5Inner.h', BASE / 'AsymmetricMap.h']
    (output / 'source_pins.json').write_text(json.dumps({
        str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs
    }, indent=2) + '\n', newline='\n')
    print('Generated scheduling copies; certified source files unchanged.')


if __name__ == '__main__':
    main()
