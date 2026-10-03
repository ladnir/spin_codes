"""Preserve the package's per-call IMT output-store policy on kernel import.

The historical research generator owns the XOR schedules. This overlay owns
only the explicit cached/streaming output writers and their dispatch plumbing.
Exact hunk matching rejects source drift instead of dropping the policy.
"""
from pathlib import Path


FILES = frozenset((
    'Spin.h', 'Spin.cpp', 'ForwardFast.cpp', 'Inner.h',
    'wide/WideKernel.h', 'wide/generated/WideCircuit.h',
))


def apply_output_stores_overlay(name, text):
    if name not in FILES:
        return text
    patch = Path(__file__).with_suffix('.patch').read_text().splitlines(True)
    selected = False
    in_hunk = False
    applied = 0
    old, new = [], []

    def flush():
        nonlocal text, old, new, applied
        if old:
            before, after = ''.join(old), ''.join(new)
            if text.count(before) != 1:
                raise ValueError(f'output-store overlay context changed: {name}, hunk {applied + 1}')
            text = text.replace(before, after, 1)
            applied += 1
        old, new = [], []

    target = f'diff --git a/spin/src/kernels/{name} b/spin/src/kernels/{name}'
    for line in patch:
        if line.startswith('diff --git '):
            flush()
            selected = line.rstrip() == target
            in_hunk = False
        elif selected and line.startswith('@@'):
            flush()
            in_hunk = True
        elif selected and in_hunk:
            if line.startswith((' ', '-')):
                old.append(line[1:])
            if line.startswith((' ', '+')):
                new.append(line[1:])
            if not line.startswith((' ', '-', '+')):
                raise ValueError(f'unsupported output-store patch line: {name}: {line!r}')
    flush()
    if applied == 0:
        raise ValueError(f'output-store overlay has no hunks for required file: {name}')
    return text
