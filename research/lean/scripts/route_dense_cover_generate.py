import pathlib,re
out=['import SpinCodes.Structured.ConcreteOuterTailDenseSound','import SpinCodes.Cover.DenseTail','','namespace Spin.Structured.ConcreteOuter','open Spin.Numeric Spin.Numeric.Fix Spin.Cover','']
for n in range(44):
 s=pathlib.Path(f'SpinCodes/Cover/Part{n}.lean').read_text(encoding='utf8')
 body=s[s.index(f'theorem part{n}_claim'):s.index('end Spin.Cover')]
 body=body.replace(f'part{n}_claim',f'closed_part{n}_claim').replace('DenseClaim','ClosedDenseClaim').replace('denseTail_of_checkTree','denseTail_closed_of_checkTree')
 out.append(body)
s=pathlib.Path('SpinCodes/Cover/DenseTail.lean').read_text(encoding='utf8')
body=s[s.index('theorem denseTail '):s.index('end Spin.Cover')]
body=body.replace('theorem denseTail ','theorem denseTail_closed ').replace('DenseClaim','ClosedDenseClaim')
body=re.sub(r'part(\d+)_claim',r'closed_part\1_claim',body)
out.append(body)
out.append('end Spin.Structured.ConcreteOuter\n')
pathlib.Path('SpinCodes/Structured/ConcreteOuterTailDenseCover.lean').write_text('\n'.join(out),encoding='utf8')
