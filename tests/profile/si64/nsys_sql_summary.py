#!/usr/bin/env python3
"""Summarize Nsight Systems SQLite exports for CP-PAW traces."""

from __future__ import annotations

import argparse
import os
import sqlite3
from collections.abc import Iterable


def has_table(conn: sqlite3.Connection, name: str) -> bool:
    row = conn.execute(
        "select 1 from sqlite_master where type='table' and name=?", (name,)
    ).fetchone()
    return row is not None


def enum_labels(conn: sqlite3.Connection, table: str) -> dict[int, str]:
    if not has_table(conn, table):
        return {}
    labels: dict[int, str] = {}
    for row in conn.execute(f"select id, coalesce(label, name) from {table}"):
        labels[int(row[0])] = str(row[1])
    return labels


def print_rows(title: str, rows: Iterable[tuple[object, ...]]) -> None:
    print(title)
    for row in rows:
        print("  " + "  ".join(str(value) for value in row))
    print()


def top_stringid_rows(
    conn: sqlite3.Connection, table: str, id_column: str, limit: int
) -> list[tuple[str, int, str]]:
    if not has_table(conn, table) or not has_table(conn, "StringIds"):
        return []
    return conn.execute(
        f"""
        select s.value, count(*) as calls,
               printf('%.6f s', sum(t.end - t.start) / 1e9) as seconds
          from {table} t
          join StringIds s on s.id = t.{id_column}
         group by s.value
         order by sum(t.end - t.start) desc
         limit ?
        """,
        (limit,),
    ).fetchall()


def summarize_memcpy(conn: sqlite3.Connection) -> list[tuple[str, int, str, str]]:
    if not has_table(conn, "CUPTI_ACTIVITY_KIND_MEMCPY"):
        return []
    labels = enum_labels(conn, "ENUM_CUDA_MEMCPY_OPER")
    rows = conn.execute(
        """
        select copyKind, count(*) as calls, sum(bytes) as bytes,
               sum(end - start) / 1e9 as seconds
          from CUPTI_ACTIVITY_KIND_MEMCPY
         group by copyKind
         order by sum(end - start) desc
        """
    ).fetchall()
    return [
        (
            labels.get(int(kind), str(kind)),
            int(calls),
            f"{int(byte_count) / 1e9:.6f} GB",
            f"{float(seconds):.6f} s",
        )
        for kind, calls, byte_count, seconds in rows
    ]


def summarize_sync(conn: sqlite3.Connection) -> list[tuple[str, int, str]]:
    if not has_table(conn, "CUPTI_ACTIVITY_KIND_SYNCHRONIZATION"):
        return []
    labels = enum_labels(conn, "ENUM_CUPTI_SYNC_TYPE")
    rows = conn.execute(
        """
        select syncType, count(*) as calls, sum(end - start) / 1e9 as seconds
          from CUPTI_ACTIVITY_KIND_SYNCHRONIZATION
         group by syncType
         order by sum(end - start) desc
        """
    ).fetchall()
    return [
        (labels.get(int(kind), str(kind)), int(calls), f"{float(seconds):.6f} s")
        for kind, calls, seconds in rows
    ]


def summarize(path: str, limit: int) -> None:
    print(f"## {path}")
    with sqlite3.connect(path) as conn:
        print_rows(
            "Top CUDA kernels",
            top_stringid_rows(conn, "CUPTI_ACTIVITY_KIND_KERNEL", "shortName", limit),
        )
        print_rows(
            "Top CUDA runtime calls",
            top_stringid_rows(conn, "CUPTI_ACTIVITY_KIND_RUNTIME", "nameId", limit),
        )
        print_rows("CUDA memcpy summary", summarize_memcpy(conn))
        print_rows("CUDA synchronization summary", summarize_sync(conn))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Summarize CUDA activity from Nsight Systems SQLite files."
    )
    parser.add_argument("sqlite", nargs="+", help="Nsight Systems .sqlite file(s)")
    parser.add_argument("--limit", type=int, default=12, help="Top rows to print")
    args = parser.parse_args()

    for index, path in enumerate(args.sqlite):
        if index:
            print()
        summarize(os.path.abspath(path), args.limit)


if __name__ == "__main__":
    main()
