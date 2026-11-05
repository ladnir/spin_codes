"""Read-only code-generation comparison for the frozen Compact7 experiment."""
import re
import subprocess

for name in ('control', 'compact7'):
    raw = subprocess.check_output(['objdump', '-Cd', '-Mintel', name], text=True)
    chunks = re.split(r'\n([0-9a-f]+) <([^\n]+)>:\n', raw)
    for j in range(1, len(chunks), 3):
        address, label, body = chunks[j:j+3]
        if 'hdgroup::coset<4u' not in label:
            continue
        instructions = []
        for line in body.splitlines():
            columns = line.split('\t')
            if len(columns) >= 3:
                instructions.append(columns[-1].strip())
        print(name, address, label)
        print('instructions', len(instructions),
              'gfni', sum('vgf2p8mulb' in s for s in instructions),
              'key_shuffles', sum('vpshufb' in s for s in instructions),
              'stack_vector_refs', sum(bool(re.search(r'\b[zy]mm', s)) and
                  re.search(r'\[(?:r[bs]p)[^]]*\]', s) is not None for s in instructions),
              'calls', sum(s.startswith('call') for s in instructions))
        print('\n'.join(s for s in instructions if s.startswith('call')))
