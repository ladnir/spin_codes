"""Apply the reviewed setup-only overlay after historical kernel generation.

Exact context matching deliberately fails if the input generator changes.
No external patch executable or Python package is required.
"""
from pathlib import Path


def apply_setup_overlay(name, text):
    patch = Path(__file__).with_suffix('.patch').read_text().splitlines(True)
    selected = False
    old, new = [], []

    def flush():
        nonlocal text, old, new
        if old:
            before, after = ''.join(old), ''.join(new)
            if text.count(before) != 1:
                raise ValueError(f'setup overlay context changed: {name}')
            text = text.replace(before, after, 1)
        old, new = [], []

    for line in patch:
        if line.startswith('diff --git '):
            flush()
            selected = line.rstrip().endswith('/' + name)
        elif selected and line.startswith('@@'):
            flush()
        elif selected and not line.startswith(('index ', '--- ', '+++ ')):
            if line.startswith((' ', '-')):
                old.append(line[1:])
            if line.startswith((' ', '+')):
                new.append(line[1:])
    flush()
    return text
