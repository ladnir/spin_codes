"""Verify that a completed v5 sparse run resumes without changing its evidence."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

import study_bch_dominance_v1 as study

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cover', type=Path, required=True)
    args = parser.parse_args()
    cover = args.cover.resolve(); data = json.loads(cover.read_text())
    if data['status'] != 'BINARY64_COMPLETE_SPARSE_INTERVAL':
        raise ValueError('a completed interval is required')
    producer = HERE/'close_bch_sparse_tail_v5.py'
    key = producer.relative_to(study.ROOT).as_posix()
    if data['source_sha256'].get(key) != study.sha(producer):
        raise ValueError('expected the current v5 producer')
    files = [cover, *(cover.parent/f'occupation_{q}.json'
                     for q in range(data['occupation_min'], data['occupation_max']+1))]
    before = {str(path): study.sha(path) for path in files}
    command = [sys.executable, '-u', str(producer)]
    for name, value in data['arguments'].items():
        command += ['--'+name.replace('_', '-'), str(value)]
    with (HERE/'bch_sparse_resume_v1.log').open('w') as log:
        subprocess.run(command, cwd=HERE, stdout=log, stderr=subprocess.STDOUT, check=True)
    after = {str(path): study.sha(path) for path in files}
    if before != after:
        raise ArithmeticError('resume changed a completed occupation or full interval')
    result = dict(status='VERIFIED_BYTE_IDENTICAL_SPARSE_RESUME', occupations=len(files)-1,
                  source_sha256={str(Path(__file__)): study.sha(Path(__file__)),
                                 str(producer): study.sha(producer), **after})
    (HERE/'bch_sparse_resume_verification_v1.json').write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], result['occupations'], flush=True)


if __name__ == '__main__':
    main()
