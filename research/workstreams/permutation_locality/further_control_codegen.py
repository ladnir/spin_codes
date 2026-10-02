"""Controlled exact-map experiments on the retained creative winner.

The alignment variant changes only owned benchmark storage. The outer variants
change a fixed loop schedule or the order of disjoint output tiles.
"""
import argparse
from pathlib import Path

import outer_transpose_codegen as transposed
from packed_coeff_codegen import definition


def replace_once(source, old, new):
    if source.count(old) != 1:
        raise ArithmeticError(f'expected one source occurrence: {old[:100]}')
    return source.replace(old, new)


def aligned_driver(path):
    source = path.read_text()
    begin = source.index('template<unsigned Mode> static void experiment(')
    end = source.index('static std::uint64_t unsignedArgument', begin)
    body = source[begin:end]
    body = replace_once(body,
        'std::vector<k::block> input(n),expected(n),initial(n),guarded((n/1024)*tileStride+12);',
        '''static_assert(SPIN_INPUT_OFFSET==0 || SPIN_INPUT_OFFSET==16 ||
                      SPIN_INPUT_OFFSET==32 || SPIN_INPUT_OFFSET==48);
    std::vector<k::block> ownedInput(n+12),expected(n),initial(n),guarded((n/1024)*tileStride+12);
    auto* inputBase=reinterpret_cast<k::block*>(((reinterpret_cast<std::uintptr_t>(ownedInput.data()+4)+63)&~std::uintptr_t(63))+SPIN_INPUT_OFFSET);
    std::span<k::block> input(inputBase,n);''')
    body = replace_once(body, 'expected=input;initial=input;',
        'std::copy(input.begin(),input.end(),expected.begin());\n'
        '        std::copy(input.begin(),input.end(),initial.begin());')
    body = replace_once(body,
        'const k::block canary(0x781aca29312357bfULL,0x928fafcac352697bULL);',
        '''const k::block canary(0x781aca29312357bfULL,0x928fafcac352697bULL);
    inputBase[-1]=canary;inputBase[n]=canary;''')
    body = replace_once(body, 'checkScratchCanaries<(Mode >= 7)>(scratch,n/1024,canary);',
        '''checkScratchCanaries<(Mode >= 7)>(scratch,n/1024,canary);
        equal(inputBase-1,&canary,16,"leading aligned input canary");
        equal(inputBase+n,&canary,16,"trailing aligned input canary");''')
    body = replace_once(body, '<<" prepared matrices; 1/2/3/4-epoch basis; whole output/suffix/canaries/padding"',
        '<<" input_mod64="<<(reinterpret_cast<std::uintptr_t>(inputBase)%64)'
        '<<" prepared matrices; 1/2/3/4-epoch basis; whole output/suffix/canaries/padding"')
    helpers = '''static void fill(std::span<k::block> x,std::uint64_t seed) {
    k::setup::Words rng(seed);for(auto& v:x){auto a=rng(),b=rng();v=k::block(a,b);}
}
static std::uint64_t hash(std::span<k::block> x) {
    std::uint64_t h=0;for(const auto& v:x){std::uint64_t halves[2];std::memcpy(halves,&v,16);
        for(auto a:halves)h=(h^a)*0x100000001b3ULL;}return h;
}
'''
    return source[:begin] + helpers + body + source[end:]


def outer(header, kind):
    source = transposed.generate(header)
    if kind == 'reverse':
        old = '\n'.join(f'packedCoeffTile<{i}>(src,out);' for i in range(0,16,2))
        new = '\n'.join(f'packedCoeffTile<{i}>(src,out);' for i in range(14,-1,-2))
        if source.count(old) != 3:
            raise ArithmeticError('expected three outer entry points')
        source = source.replace(old,new)
    elif kind == 'unroll':
        old = definition(source, 'for(unsigned input=0;input<16;input+=2){')
        body = old[old.index('{')+1:old.rfind('}')]
        new = '\n'.join('{constexpr unsigned input='+str(i)+';'+body+'}'
                        for i in range(0,16,2))
        source = replace_once(source,old,new)
    return source


def wide_driver(path):
    source=path.read_text()
    source=replace_once(source,'#include "InnerComposed.h"',
        '#include "InnerComposed.h"\n#include "further_mapWide.h"\n'
        '#ifndef SPIN_FURTHER_PRUNED\n#define SPIN_FURTHER_PRUNED 0\n#endif\n'
        '#ifndef SPIN_FURTHER_INCREMENTAL\n#define SPIN_FURTHER_INCREMENTAL 0\n#endif')
    source=replace_once(source,'    else\n        ip::transposeInnerComposedWide<',
        '    else if constexpr(Mode == 4)\n'
        '        ip::transposeInnerFurtherWide<bool(SPIN_FURTHER_PRUNED),bool(SPIN_FURTHER_INCREMENTAL)>(input,n,composed.rows(),emit);\n'
        '    else\n        ip::transposeInnerComposedWide<')
    source=replace_once(source,'        sizeprobe::checkInnerBoundaries();',
        '        sizeprobe::checkInnerBoundaries();\n'
        '        if(!ip::furtherMapWideSelfCheck())throw std::runtime_error("wide emission basis");')
    return source


def raw_driver(path):
    source=path.read_text()
    source=replace_once(source,'#include <limits>', '#include <limits>\n#include "InnerComposed.h"')
    old='const ip::PreparedPacketUpdates16& updates'
    if source.count(old)!=2: raise ArithmeticError('expected two raw update parameters')
    source=source.replace(old,'const ip::PreparedComposedUpdates16& updates')
    source=replace_once(source,
        'ip::transposeInner<ip::InnerRecipe::StreamVbmi>(input, n, updates.conjugatedRows(), emit);',
        'ip::transposeInnerComposedWide<ip::ComposedRecipe::Stream>(input,n,updates.rows(),emit);')
    source=replace_once(source,'const ip::PreparedPacketUpdates16 transformed(setup.reverse);',
        'const ip::PreparedPacketUpdates16 transformed(setup.reverse);\n'
        '    const ip::PreparedComposedUpdates16 composed(transformed.conjugatedRows());')
    source=replace_once(source,'route, coeff, transformed, phases);','route, coeff, composed, phases);')
    source=replace_once(source,'// StreamVbmi/cached/no-prefetch recipe and allocates nothing.',
        '// ComposedWide Stream/cached/no-prefetch recipe and allocates nothing.')
    source=replace_once(source,'int main(int argc, char** argv) { try {',
        'int main(int argc, char** argv) { try {\n'
        '    if(!__builtin_cpu_supports("vpclmulqdq"))throw std::runtime_error("VPCLMULQDQ required");')
    return source


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path)
    parser.add_argument('--kind',choices=('align','reverse','unroll','wide-driver','raw-driver'),required=True)
    args=parser.parse_args()
    if args.kind=='align': source=aligned_driver(args.input)
    elif args.kind=='wide-driver': source=wide_driver(args.input)
    elif args.kind=='raw-driver': source=raw_driver(args.input)
    else: source=outer(args.input,args.kind)
    print(source)


if __name__=='__main__':
    main()
