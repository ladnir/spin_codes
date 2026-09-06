"""Wrap LP tokens to avoid QSopt_ex's line-length limit; no algebraic edits."""
import json
import textwrap
from pathlib import Path
from run_higher_endpoint_preflight import sha,write_new
ROOT=Path(__file__).resolve().parents[1]
FOLDER=ROOT/'generated/split38_probe'
source=FOLDER/'h38.lp'
output=FOLDER/'h38_wrapped.lp'
text=source.read_text()
wrapped='\n'.join('\n'.join(textwrap.wrap(line,width=2000,break_long_words=False,break_on_hyphens=False,
                                        subsequent_indent=' ')) for line in text.splitlines())+'\n'
assert text.split()==wrapped.split()
with output.open('x',encoding='ascii') as stream:
    stream.write(wrapped)
write_new(FOLDER/'wrapping.json',dict(tokens_identical=True,
          source_sha256={str(p.relative_to(ROOT)):sha(p) for p in (Path(__file__),source,output)}))
print('LP wrapped with exactly identical token stream')
