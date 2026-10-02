"""Match the measured header and generated map manifest to the full bound."""
import argparse
import hashlib
import json
from pathlib import Path
import re


def run(a):
    proof=json.loads(a.proof.read_text())
    manifest=json.loads(a.manifest.read_text())
    inner=proof['instance']['inner']
    assert proof['full_distance_proved'] and proof['margin_bits']>40
    assert (inner['t'],inner['s'])==(64,12)
    for key in ('t','s','expansion_columns','feedback_columns','transvection_rounds'):
        assert manifest['inner'][key]==inner[key]
    text=a.header.read_text()
    # Transpose code calls the feedback map's columns "columns". Its
    # "feedbackColumns" feed the transpose of the forward expansion map.
    for name,key in [('columns','feedback_columns'),('feedbackColumns','expansion_columns')]:
        match=re.search(r'\b'+name+r'\{([^}]+)\};',text)
        assert match
        assert list(map(int,match[1].split(',')))==inner[key]
    digest=hashlib.sha256(a.header.read_bytes()).hexdigest()
    matches=[line for line in a.sources.read_text().splitlines()
             if line.endswith('/AsymmetricMap.h')]
    assert len(matches)==2
    assert all(line.split()[0]==digest for line in matches)
    print('Full-bound expansion/feedback columns match both measured circuit variants.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('proof','manifest','header','sources'):p.add_argument('--'+name,type=Path,required=True)
    run(p.parse_args())
