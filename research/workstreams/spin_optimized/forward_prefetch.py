"""Optional compile-time gather lookahead; changes no routing or arithmetic."""
from pathlib import Path
import sys

root = Path(sys.argv[1])
distance = int(sys.argv[2])
assert distance in (16, 32, 64, 128)
for name in ('Spin.cpp', 'ForwardFast.cpp'):
    path = root / name
    if not path.exists():
        continue
    text = path.read_text()
    anchor = '''    innerForward<Map>(n,mForwardFieldRows.data(),[&](std::size_t i) {
        if constexpr(Packed) return values[unpack(mSlots24.data()+3*i)];'''
    assert text.count(anchor) == 1, name
    replacement = f'''    innerForward<Map>(n,mForwardFieldRows.data(),[&](std::size_t i) {{
        if(i+{distance}<n) {{
            const auto next=Packed?unpack(mSlots24.data()+3*(i+{distance})):mSlots32[i+{distance}];
            _mm_prefetch(reinterpret_cast<const char*>(values+next),_MM_HINT_T0);
        }}
        if constexpr(Packed) return values[unpack(mSlots24.data()+3*i)];'''
    path.write_text(text.replace(anchor, replacement), newline='\n')
path = root / 'wide/WideKernel.h'
if path.exists():
    text = path.read_text()
    anchor = '        return values[unpack(view.slots+3*i)];'
    assert text.count(anchor) == 1
    replacement = f'''        if(i+{distance}<view.n)
            _mm_prefetch(reinterpret_cast<const char*>(values+unpack(view.slots+3*(i+{distance}))),_MM_HINT_T0);
{anchor}'''
    path.write_text(text.replace(anchor, replacement), newline='\n')
