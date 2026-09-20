#!/usr/bin/env python3
"""Query or export the finite-SPIN landscape database."""

from __future__ import annotations

import argparse
import csv
import pathlib
import sqlite3
import sys


HERE = pathlib.Path(__file__).resolve().parent
DEFAULT_DB = HERE / "spin_landscape.sqlite3"


def print_rows(columns: list[str], rows: list[tuple[object, ...]]) -> None:
    widths = [len(column) for column in columns]
    rendered = [["" if value is None else str(value) for value in row] for row in rows]
    for row in rendered:
        widths = [max(width, len(value)) for width, value in zip(widths, row)]
    print("  ".join(column.ljust(width) for column, width in zip(columns, widths)))
    print("  ".join("-" * width for width in widths))
    for row in rendered:
        print("  ".join(value.ljust(width) for value, width in zip(row, widths)))


def summary(db: sqlite3.Connection) -> None:
    query = """
        SELECT result_class, transfer_review_status, outer_family, coverage_kind, COUNT(*) AS rows
        FROM landscape
        GROUP BY result_class, transfer_review_status, outer_family, coverage_kind
        ORDER BY result_class, outer_family, coverage_kind
    """
    cursor = db.execute(query)
    print_rows([item[0] for item in cursor.description], cursor.fetchall())


def curves(db: sqlite3.Connection, family: str | None) -> None:
    conditions = ["coverage_kind='single_occupation'", "occupation_min=1"]
    parameters: list[object] = []
    if family:
        conditions.append("outer_family=?")
        parameters.append(family)
    query = f"""
        SELECT study,outer_label,message_exponent,step_bits,state_bits,
               ROUND(persistence_offset,3) AS persistence_offset,
               ROUND(margin_bits,6) AS margin_bits,result_class,transfer_review_status,comparison_eligible
        FROM q1_curves
        WHERE {' AND '.join(conditions)}
        ORDER BY outer_family,block_bits,outer_label,study,message_bits,state_bits
    """
    cursor = db.execute(query, parameters)
    print_rows([item[0] for item in cursor.description], cursor.fetchall())


def export(db: sqlite3.Connection, output: pathlib.Path) -> None:
    cursor = db.execute("SELECT * FROM landscape ORDER BY study,message_bits,occupation_min")
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([item[0] for item in cursor.description])
        writer.writerows(cursor)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--database", type=pathlib.Path, default=DEFAULT_DB)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("summary")
    curve_parser = subparsers.add_parser("curves")
    curve_parser.add_argument("--family", choices=["bch", "rm", "random"])
    export_parser = subparsers.add_parser("export")
    export_parser.add_argument("output", type=pathlib.Path)
    sql_parser = subparsers.add_parser("sql")
    sql_parser.add_argument("query")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    db = sqlite3.connect(args.database)
    try:
        if args.command == "summary":
            summary(db)
        elif args.command == "curves":
            curves(db, args.family)
        elif args.command == "export":
            export(db, args.output)
        elif args.command == "sql":
            cursor = db.execute(args.query)
            print_rows([item[0] for item in cursor.description], cursor.fetchall())
        else:
            raise AssertionError(args.command)
    finally:
        db.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
