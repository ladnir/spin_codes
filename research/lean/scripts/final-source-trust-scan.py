"""Read-only comment-aware project source scan, independent of Lean compilation."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, importlib.util, json, re

ROOT = Path(__file__).resolve().parents[1]
PARSER = ROOT / 'scripts/replay-native-closure.py'
OUT = ROOT / 'scripts/map_data/final_source_trust_scan.json'
SPEC = importlib.util.spec_from_file_location('closure_comment_parser', PARSER)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)
TOKEN = re.compile(r'\b(?:sorry|sorryAx|admit|axiom|native_decide|ofReduceBool|trustCompiler)\b')
STRING = re.compile(r'"(?:\\.|[^"\\])*"')

def sha_bytes(content):
    return hashlib.sha256(content).hexdigest()

def files():
    return sorted([ROOT / 'SpinCodes.lean', *(ROOT / 'SpinCodes').rglob('*.lean')],
                  key=lambda p: p.relative_to(ROOT).as_posix())

def main():
    script_sha = sha_bytes(Path(__file__).read_bytes())
    parser_sha = sha_bytes(PARSER.read_bytes())
    inventory, code_hits, quoted_hits = {}, [], []
    source_files = files()
    for path in source_files:
        relative = path.relative_to(ROOT).as_posix()
        raw = path.read_bytes()
        inventory[relative] = sha_bytes(raw)
        cleaned = MODULE.uncomment(raw.decode('utf-8-sig'))
        # Retain and separately report quoted tokens. Interpolated strings can
        # contain executable Lean expressions, so such hits are not dismissed.
        string_matches = list(STRING.finditer(cleaned))
        string_index = 0
        for match in TOKEN.finditer(cleaned):
            while string_index < len(string_matches) and string_matches[string_index].end() <= match.start():
                string_index += 1
            quote = (string_matches[string_index] if string_index < len(string_matches)
                     and string_matches[string_index].start() <= match.start() else None)
            record = {'source': relative, 'line': cleaned.count('\n', 0, match.start()) + 1,
                      'token': match.group(),
                      'context': cleaned[max(0, match.start()-65):min(len(cleaned), match.end()+65)]}
            if quote:
                prefix = cleaned[max(0, quote.start()-12):quote.start()]
                record['interpolated_string_review_required'] = bool(re.search(r'[A-Za-z_][A-Za-z_0-9]*!\s*$', prefix))
                quoted_hits.append(record)
                if record['interpolated_string_review_required']:
                    code_hits.append(dict(record, classification='interpolated string requires review'))
            else:
                code_hits.append(record)

    # Bind the result to a deterministic inventory and reject a moving snapshot.
    assert [p.relative_to(ROOT).as_posix() for p in files()] == list(inventory), 'Source set changed'
    for path in source_files:
        relative = path.relative_to(ROOT).as_posix()
        assert sha_bytes(path.read_bytes()) == inventory[relative], f'Source changed: {relative}'
    assert sha_bytes(Path(__file__).read_bytes()) == script_sha
    assert sha_bytes(PARSER.read_bytes()) == parser_sha
    inventory_text = ''.join(relative + '\0' + digest + '\n' for relative, digest in inventory.items())
    report = {'status': 'PASS' if not code_hits else 'REVIEW_REQUIRED',
              'checked_at_utc': datetime.now(timezone.utc).isoformat(),
              'scope': 'All SpinCodes/**/*.lean plus SpinCodes.lean; nested comments stripped using the existing replay parser; quoted text separately classified. No compiler run and no proof-source changes. This lexical scan complements recursive axiom audits.',
              'source_count': len(inventory),
              'inventory_sha256': sha_bytes(inventory_text.encode('utf-8')),
              'inventory_encoding': 'UTF-8 concatenation of sorted relative path, NUL, source SHA256, newline',
              'script': Path(__file__).relative_to(ROOT).as_posix(), 'script_sha256': script_sha,
              'comment_parser': PARSER.relative_to(ROOT).as_posix(), 'comment_parser_sha256': parser_sha,
              'tokens': ['sorry', 'sorryAx', 'admit', 'axiom', 'native_decide', 'ofReduceBool', 'trustCompiler'],
              'code_hits': code_hits, 'quoted_text_hits': quoted_hits,
              'note': 'Token matching catches private/protected axiom declarations and qualified uses such as Lean.ofReduceBool and Lean.trustCompiler.'}
    tmp = OUT.with_suffix('.tmp')
    tmp.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    tmp.replace(OUT)
    print(json.dumps(report, indent=2))

if __name__ == '__main__':
    main()
