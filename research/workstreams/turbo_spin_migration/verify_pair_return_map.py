"""Check the exact s22 feedback interface used in the pair-return proof."""

from pathlib import Path
import re


SOURCE = Path(__file__).resolve().parents[3] / "spin/src/certified22/productionMap.h"


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
    source = SOURCE.read_text(encoding="utf-8")
    rows = literal_array(source, "oldRows", 22)
    columns = literal_array(source, "oldColumns", 64)
    reconstructed = [
        sum(((row >> bit) & 1) << index for index, row in enumerate(rows))
        for bit in range(64)
    ]
    assert reconstructed == columns, "feedback columns differ from map rows"
    ranks = [binary_rank(columns[4 * slot:4 * slot + 4]) for slot in range(16)]
    assert ranks == [4] * 16, ranks
    print("s22 feedback: 64 columns match rows; all 16 packet slots have rank four")


if __name__ == "__main__":
    main()
