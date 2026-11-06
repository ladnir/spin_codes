"""Check the exact s22 feedback interface used in the pair-return proof."""

from pathlib import Path
from collections import Counter
import argparse
import re


SOURCE = Path(__file__).resolve().parents[3] / "spin/src/certified22/productionMap.h"
PAPER = Path(__file__).resolve().parents[3] / "research/paper/structured_uniform_appendix.tex"


def literal_array(source: str, name: str, expected: int) -> list[int]:
    marker = name + "{{"
    text = source.split(marker, 1)[1].split("}}", 1)[0]
    values = [int(word, 16) for word in re.findall(r"0x([0-9a-fA-F]+)U(?:LL)?", text)]
    assert len(values) == expected, (name, len(values))
    return values


def binary_rank(columns: list[int]) -> int:
    pivots: dict[int, int] = {}
    for column in columns:
        value = column
        while value:
            bit = value.bit_length() - 1
            if bit not in pivots:
                pivots[bit] = value
                break
            value ^= pivots[bit]
    return len(pivots)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", action="store_true",
                        help="enumerate all nonzero s22 state output weights")
    args = parser.parse_args()
    source = SOURCE.read_text(encoding="utf-8")
    rows = literal_array(source, "oldRows", 22)
    columns = literal_array(source, "oldColumns", 64)
    paper = PAPER.read_text(encoding="utf-8")
    table = paper.split(r"\begin{tabular}{@{}llll@{}}", 1)[1].split(
        r"\end{tabular}", 1)[0]
    paper_rows = [int(word, 16) for word in re.findall(r"\b[0-9a-f]{16}\b", table)]
    assert paper_rows == rows, "the manuscript's s22 map differs from production"
    reconstructed = [
        sum(((row >> bit) & 1) << index for index, row in enumerate(rows))
        for bit in range(64)
    ]
    assert reconstructed == columns, "feedback columns differ from map rows"
    ranks = [binary_rank(columns[4 * slot:4 * slot + 4]) for slot in range(16)]
    assert ranks == [4] * 16, ranks
    print("s22 map: manuscript matches production; 64 feedback columns match rows; all 16 packet slots have rank four")
    if args.spectrum:
        output = 0
        weights: Counter[int] = Counter()
        for gray_index in range(1, 1 << 22):
            output ^= rows[(gray_index & -gray_index).bit_length() - 1]
            weights[output.bit_count()] += 1
        assert sum(weights.values()) == (1 << 22) - 1
        assert 0 not in weights
        assert sum(weight * count for weight, count in weights.items()) == 64 * (1 << 21)
        print("nonzero-state output spectrum:", sorted(weights.items()))


if __name__ == "__main__":
    main()
