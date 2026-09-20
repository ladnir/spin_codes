"""Isolate K20 tiled streaming stores; preserve the mathematical encoder."""
from pathlib import Path
import sys

p=Path(sys.argv[1])
mode=int(sys.argv[2]) if len(sys.argv)>2 else 5
assert mode in (5,6,7)
def edit(text,old,new):
    assert text.count(old)==1,old
    return text.replace(old,new)
h=(p/'Spin.h').read_text();s=(p/'Spin.cpp').read_text();f=(p/'Fast.cpp').read_text()
bidir='forwardBits(' in h
signature='template<bool Packed>' if bidir else 'template<class Map,bool Packed>'
start=f.index(signature+' void Spin::runFour(')
brace=f.index('{',start)
depth=1;end=brace+1
while depth:
    depth+=(f[end]=='{')-(f[end]=='}');end+=1
kernel=f[start:end].replace('Spin::runFour(', 'Spin::runK20Tiled(')
if mode==5:
    kernel=edit(kernel,'values[unpack(mSlots24.data()+3*i)]=v;',
        '_mm_stream_si128(reinterpret_cast<__m128i*>(values+unpack(mSlots24.data()+3*i)),v.mData);')
    kernel=edit(kernel,'values[mSlots32[i]]=v;',
        '_mm_stream_si128(reinterpret_cast<__m128i*>(values+mSlots32[i]),v.mData);')
else:
    # The setup schedule visits each bucket in strictly descending slot order.
    # With four fixed buckets of 2^19 blocks, each group arrives at lanes 3,2,1,0.
    # All pending lanes are initialized before lane 0 triggers the four stores.
    kernel=edit(kernel,'    block* values=w.buckets.data();',
        '    alignas(64) block pending[4][4];\n    block* values=w.buckets.data();')
    begin=kernel.index('        if constexpr(Packed) values[')
    finish=kernel.index('\n    });',begin)
    writes=''.join(f'            _mm_stream_si128(reinterpret_cast<__m128i*>(values+slot+{j}),pending[bucket][{j}].mData);\n' if mode==6 else
                   f'            values[slot+{j}]=pending[bucket][{j}];\n' for j in range(4))
    kernel=kernel[:begin]+'''        const u32 slot=Packed?unpack(mSlots24.data()+3*i):mSlots32[i];
        const u32 bucket=slot>>19;
        pending[bucket][slot&3]=v;
        if(!(slot&3)) {
'''+writes+'        }'+kernel[finish:]
# Keep the same barrier in the cached-store control to isolate the store policy.
kernel=edit(kernel,'    block* tile=w.tile.data();',
    '    // Complete bucket writes before consuming them.\n'
    '    _mm_sfence();\n    block* tile=w.tile.data();')
h=edit(h,'    template<class Map> void setupInner',
    '    '+signature+' void runK20Tiled(const block*,block*,Workspace&) const;\n'
    '    template<class Map> void setupInner')
anchor='    if(mBchBackend==BchBackend::Avx512) {\n'
at=s.index(anchor,s.index('void Spin::encodeUnchecked'))+len(anchor)
args='' if bidir else 'AsymmetricMap,'
condition='mK==(1U<<20)'+(' && tileBlocks()==(1U<<19)' if mode!=5 else '')
s=s[:at]+f'''        if({condition}) {{
            if(layout==Layout::Packed24) runK20Tiled<{args}true>(in,out,w);
            else runK20Tiled<{args}false>(in,out,w);
            return;
        }}
'''+s[at:]
f=f[:-2]+kernel+'\n'+''.join(
    f'template void Spin::runK20Tiled<{args}{packed}>(const block*,block*,Workspace&) const;\n'
    for packed in ('true','false'))+'}\n'
for name,text in [('Spin.h',h),('Spin.cpp',s),('Fast.cpp',f)]:
    (p/name).write_text(text,newline='\n')
print('Generated K20 tiled store mode',mode)
