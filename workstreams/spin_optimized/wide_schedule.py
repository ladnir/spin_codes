"""Verify and reorder the wide BCH XOR circuit; never alter the linear map."""
import argparse
from functools import reduce
import hashlib
import json
from operator import xor
from pathlib import Path
import re


def require(condition, message):
    # These checks must also execute under python -O.
    if not condition:
        raise ValueError(message)


def xor_inputs(expression):
    """Accept only signal references or binary vx expressions."""
    tokens = re.findall(r"[ogz]\d+|vx|[(),]", expression)
    require("".join(tokens) == expression, f"Unsupported XOR expression: {expression}")
    position = 0
    signals = []

    def consume(expected=None):
        nonlocal position
        require(position < len(tokens), "Truncated XOR expression")
        token = tokens[position]
        position += 1
        require(expected is None or token == expected, "Malformed XOR expression")
        return token

    def parse():
        token = consume()
        if re.fullmatch(r"[ogz]\d+", token):
            signals.append(token)
        else:
            require(token == "vx", "Expected a signal or binary vx")
            consume("(")
            parse()
            consume(",")
            parse()
            consume(")")

    parse()
    require(position == len(tokens), "Trailing tokens in XOR expression")
    return signals


def reorder(source, matrix):
    """Return (rewritten header, XOR count) after checking all 256 outputs."""
    anchor = "static void bchForward(const Vec* a,Vec* x) {"
    require(source.count(anchor) == 1, "Wide BCH function anchor changed")
    start = source.index(anchor) + len(anchor)
    end = source.index("\n}", start)
    definitions, stores = {}, []
    for line in source[start:end].splitlines():
        line = line.strip()
        if not line:
            continue
        declaration = re.fullmatch(r"const auto ([ogz]\d+)=(.*);", line)
        store = re.fullmatch(r"Ops::store\(x\+(\d+),z(\d+)\);", line)
        if declaration:
            name, expression = declaration.groups()
            require(name not in definitions, f"Duplicate signal: {name}")
            definitions[name] = expression
        elif store:
            index, signal = map(int, store.groups())
            require(index == signal, "BCH store index does not match output")
            stores.append(index)
        else:
            raise ValueError(f"Unrecognized BCH statement: {line}")
    inputs = {f"o{i}" for i in range(128)}
    outputs = {f"z{i}" for i in range(256)}
    require({v for v in definitions if v.startswith("o")} == inputs, "Unexpected BCH inputs")
    require({v for v in definitions if v.startswith("z")} == outputs, "Unexpected BCH outputs")
    require(sorted(stores) == list(range(256)), "Missing or duplicate BCH stores")
    for i in range(128):
        require(definitions[f"o{i}"] == f"Ops::load(a+{i})", "Unexpected BCH input load")
    dependencies = {name: [] if name in inputs else xor_inputs(expr)
                    for name, expr in definitions.items()}
    forms = {f"o{i}": 1 << i for i in range(128)}
    pending = set()

    def form(name):
        require(name in definitions, f"Unknown signal: {name}")
        require(name not in pending, f"Cycle in BCH circuit at {name}")
        if name not in forms:
            pending.add(name)
            forms[name] = reduce(xor, (form(dep) for dep in dependencies[name]), 0)
            pending.remove(name)
        return forms[name]

    words = [int(x, 16) for x in re.findall(r"0x([0-9a-fA-F]+)ULL", matrix)]
    require(len(words) == 512, "Expected 128 by 256 BCH matrix")
    targets = [sum(((words[4*i+j//64] >> (j % 64)) & 1) << i for i in range(128))
               for j in range(256)]
    require([form(f"z{j}") for j in range(256)] == targets, "BCH output differs from matrix")
    emitted, code = set(), []

    def visit(name):
        if name not in emitted:
            for dep in dependencies[name]:
                visit(dep)
            code.append(f"const auto {name}={definitions[name]};")
            emitted.add(name)

    for j in range(256):
        visit(f"z{j}")
        code.append(f"Ops::store(x+{j},z{j});")
    return (source[:start] + "\n" + "\n".join(code) + source[end:],
            sum(definitions[name].count("vx(") for name in emitted))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build_source", type=Path, help="Generated bidirectional source directory")
    root = parser.parse_args().build_source
    path = root / "wide/generated/WideCircuit.h"
    result, count = reorder(path.read_text(), (root / "generated/BchCircuit.h").read_text())
    manifest = root / "wide/generated/WIDE_MANIFEST.json"
    record = json.loads(manifest.read_text())
    record.update(bch_schedule="dfs", bch_xors=count, bch_outputs_verified=256,
                  output_sha256=hashlib.sha256(result.encode()).hexdigest())
    path.write_text(result, newline="\n")
    manifest.write_text(json.dumps(record, indent=2) + "\n", newline="\n")
    print(f"Wide BCH DFS: {count} XORs; all 256 outputs verified")


if __name__ == "__main__":
    main()
