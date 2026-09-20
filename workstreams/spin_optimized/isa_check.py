"""Check GCC/Clang object-code isolation: EVEX must stay in the AVX-512 object."""
from pathlib import Path
import re
import subprocess
import sys
build = Path(sys.argv[1])
paths=[
    ('CMakeFiles/spin_half_transpose.dir/generated/Spin.cpp.o',False),
    ('CMakeFiles/spin_half_transpose.dir/generated/generated/BchCircuit.cpp.o',False),
    ('CMakeFiles/spin_bch512.dir/generated/Fast.cpp.o',True),
    ('CMakeFiles/spin_bch512.dir/generated/generated/BchAvx512.cpp.o',True)]
for path in ('CMakeFiles/spin_half_transpose.dir/generated/GenericSpin.cpp.o',
             'CMakeFiles/spin_generic_test.dir/generic_test.cpp.o'):
    if (build/path).exists(): paths.append((path,False))
if (build/'bidirectional').exists():
    paths += [
        ('CMakeFiles/spin_half_bidirectional.dir/bidirectional/Spin.cpp.o',False),
        ('CMakeFiles/spin_half_bidirectional.dir/bidirectional/generated/BchCircuit.cpp.o',False),
        ('CMakeFiles/spin_half_bidirectional.dir/bidirectional/generated/BchForward.cpp.o',False),
        ('CMakeFiles/spin_bidir512.dir/bidirectional/Fast.cpp.o',True),
        ('CMakeFiles/spin_bidir512.dir/bidirectional/generated/BchAvx512.cpp.o',True)]
if (build/'CMakeFiles/spin_wide.dir').exists():
    paths += [
        ('CMakeFiles/spin_wide.dir/Wide.cpp.o',False),
        ('CMakeFiles/spin_wide256.dir/bidirectional/wide/Wide256.cpp.o',False),
        ('CMakeFiles/spin_wide512.dir/bidirectional/wide/Wide512.cpp.o',True)]
    if (build/'CMakeFiles/spin_column512.dir').exists():
        paths += [('CMakeFiles/spin_column512.dir/ColumnAssembly.cpp.o',True),
                  ('CMakeFiles/spin_column512.dir/bidirectional/Column128.cpp.o',True),
                  ('CMakeFiles/spin_wide_benchmark.dir/wide_benchmark.cpp.o',False)]
for path, expected in paths:
    assembly = subprocess.check_output(['objdump','-d',str(build / path)],text=True)
    evex = bool(re.search(r'^\s*[0-9a-f]+:\s+62\s',assembly,re.M))
    if evex != expected:
        raise RuntimeError(f'unexpected EVEX presence ({evex}) in {path}')
    print(path, 'EVEX confined' if expected else 'no EVEX', 'PASS')
