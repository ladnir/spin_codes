"""Research-only packed GL8/BCH kernels, keeping the existing generator intact."""
from pathlib import Path
import subprocess
import sys


def main():
    generator = Path(__file__).with_name("gfni_bch.py")
    print(subprocess.check_output([sys.executable, str(generator), sys.argv[1]], text=True))
    print('namespace spin::detail::kernel {')
    # The matrix rows use Intel GFNI's reversed byte order; affine constant is zero.
    print('template<bool Independent> static SPIN_FORCEINLINE __m512i mixMatrix(const std::uint64_t* m,unsigned group) {')
    print('if constexpr(Independent)return _mm512_loadu_si512(m+8*group);')
    print('else return _mm512_set1_epi64(m[group]);\n}')
    print('template<unsigned output> static SPIN_NOINLINE void packedFullTile(const __m512i* __restrict src,block* __restrict out) {')
    for p in range(2):
        for j in range(8):
            systematic = f'_mm512_load_si512(src+8*(output+{p})+{j})'
            print(f'auto y{8*p+j}={systematic};')
            print(f'if constexpr(output+{p}==15)y{8*p+j}=_mm512_and_si512(y{8*p+j},_mm512_set1_epi8(0x7f));')
        print(f'constexpr auto parity{p}=[] {{std::uint64_t m=0;for(unsigned j=0;j<8;++j)if(BchRows[8*(output+{p})+j][1]>>63)m|=std::uint64_t(0x80)<<(8*(7-j));return m;}}();')
        for j in range(8):
            print(f'y{8*p+j}=_mm512_xor_si512(y{8*p+j},_mm512_gf2p8affine_epi64_epi8(_mm512_load_si512(src+120+{j}),_mm512_set1_epi64(parity{p}),0));')
    print('for(unsigned input=0;input<16;input+=2) {')
    for p in range(2):
        print(f'const auto m{p}=_mm512_set1_epi64(matrices[output+{p}][input]);')
        print(f'const auto n{p}=_mm512_set1_epi64(matrices[output+{p}][input+1]);')
    for j in range(8):
        print(f'const auto x{j}=_mm512_load_si512(src+128+8*input+{j});')
        print(f'const auto z{j}=_mm512_load_si512(src+128+8*input+8+{j});')
        for p in range(2):
            print(f'y{8*p+j}=_mm512_ternarylogic_epi64(y{8*p+j},_mm512_gf2p8affine_epi64_epi8(x{j},m{p},0),_mm512_gf2p8affine_epi64_epi8(z{j},n{p},0),0x96);')
    print('}')
    for p in range(2):
        print('orthoBlend('+','.join(f'y{8*p+j}' for j in range(8))+');')
        for j in range(8):
            for lane in range(4):
                print(f'out[{128*lane}+8*(output+{p})+{j}]=block(_mm512_extracti32x4_epi32(y{8*p+j},{lane}));')
    print('}')
    for mode, name in ((1,'DenseShared'),(2,'FullShared'),(3,'FullIndependent'),(4,'FullGl32'),(7,'FullLegacyAnd')):
        full=mode>1
        print(f'SPIN_NOINLINE void bchPacked{name}(const block* __restrict a,block* __restrict out,const std::uint64_t* __restrict coeff) {{')
        print(f'alignas(64) __m512i src[{256 if full else 128}];')
        print(f'for(unsigned group={0 if full else 16};group<32;++group) {{')
        for j in range(8):
            print(f'auto v{j}=_mm512_loadu_si512(a+4*(8*group+{j}));')
        print('orthoBlend('+','.join(f'v{j}' for j in range(8))+');')
        if mode>=4:
            for d in range(4):
                print(f'const auto m{d}=_mm512_loadu_si512(coeff+32*group+8*{d});')
        else:
            print(f'const auto m=mixMatrix<{"true" if mode>2 else "false"}>(coeff,group);')
        for j in range(8):
            if mode>=4:
                term=f'_mm512_and_si512(z{j},m0)' if mode==7 else f'_mm512_gf2p8affine_epi64_epi8(z{j},m0,0)'
                print(f'const auto z{j}=v{j};v{j}={term};')
                for d,imm in ((1,'0x39'),(2,'0x4e'),(3,'0x93')):
                    rotated=f'_mm512_shuffle_i32x4(z{j},z{j},{imm})'
                    term=f'_mm512_and_si512({rotated},m{d})' if mode==7 else f'_mm512_gf2p8affine_epi64_epi8({rotated},m{d},0)'
                    print(f'v{j}=_mm512_xor_si512(v{j},{term});')
            else:
                print(f'v{j}=_mm512_gf2p8affine_epi64_epi8(v{j},m,0);')
            print(f'_mm512_store_si512(src+8*(group-{0 if full else 16})+{j},v{j});')
        print('}')
        for output in range(0,16,2):
            print(f'packedFullTile<{output}>(src,out);' if full else f'tileBlend<{output}>(src,a,out);')
        print('}')
    print('}')


if __name__ == '__main__':
    main()
