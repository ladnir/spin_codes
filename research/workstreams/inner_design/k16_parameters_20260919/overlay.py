"""Modify a fresh generated experimental build, never the certified sources."""
from pathlib import Path
import sys
root,header,t,s,shared = Path(sys.argv[1]),Path(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]),int(sys.argv[5])

def replace(text,a,b):
    assert text.count(a)==1,a
    return text.replace(a,b)

p=root/'generated'
(p/'AsymmetricMap.h').write_bytes(header.read_bytes())
text=(p/'Spin.cpp').read_text()
text=replace(text,'if constexpr(Map::T==128 && Map::S==19)','if constexpr(true)')
text=replace(text,'case Configuration::T128S19:return 128;',f'case Configuration::T128S19:return {t};')
text=replace(text,'case Configuration::T128S19:return 19;',f'case Configuration::T128S19:return {s};')
text=replace(text,'(u>>19) || (v>>19)',f'(u>>{s}) || (v>>{s})')
text=text.replace('t128_s19_weight5_seed0_r1',f't{t}_s{s}_adaptive_r1')
(p/'Spin.cpp').write_text(text)
text=(p/'Weight5Inner.h').read_text()
text=replace(text,'static_assert(Map::T==128 && Map::S==19);',f'static_assert(Map::T=={t} && Map::S=={s});')
text=replace(text,'#if W5_SHARED_EMISSION==2',f'#if 0')
text=replace(text,'#elif W5_SHARED_EMISSION==1',f'#elif {shared}')
(p/'Weight5Inner.h').write_text(text)
# Restrict the retained dense-oracle test to the intended K16 regime.
source=Path(sys.argv[6])
(root/'correctness.cpp').write_text(source.read_text().replace('{16U,18U,20U}','{16U}'))
