#!/usr/bin/env python3

import argparse
import math
import struct
from pathlib import Path


def read_records(data):
    records = []
    offset = 0
    while offset < len(data):
        if offset + 8 > len(data):
            raise ValueError("truncated Fortran record")
        size = struct.unpack_from("<I", data, offset)[0]
        end = offset + 4 + size
        if end + 4 > len(data):
            raise ValueError("record payload extends beyond the file")
        if struct.unpack_from("<I", data, end)[0] != size:
            raise ValueError("Fortran record markers do not match")
        records.append(bytearray(data[offset + 4 : end]))
        offset = end + 4
    return records


def write_records(records):
    output = bytearray()
    for record in records:
        marker = struct.pack("<I", len(record))
        output.extend(marker)
        output.extend(record)
        output.extend(marker)
    return output


def record_name(record):
    if len(record) < 36:
        return ""
    return bytes(record[4:36]).decode("ascii", errors="ignore").strip()


def displace(data, atom, axis, delta):
    """Shift both coordinate time levels without changing electronic records."""
    if axis not in (1, 2, 3) or not math.isfinite(delta):
        raise ValueError("Need a Cartesian axis and finite displacement")
    records = read_records(data)
    atom_header = next(
        (index for index, record in enumerate(records) if record_name(record) == "ATOMS"),
        None,
    )
    if atom_header is None:
        raise ValueError("ATOMS section not found")
    natom = struct.unpack("<i", records[atom_header + 1])[0]
    if not 1 <= atom <= natom:
        raise ValueError(f"atom index must be between 1 and {natom}")

    coordinate_records = (atom_header + 3, atom_header + 4)
    expected_size = 3 * natom * 8
    if any(len(records[index]) != expected_size for index in coordinate_records):
        raise ValueError("unexpected ATOMS coordinate record size")
    component = 3 * (atom - 1) + axis - 1
    for index in coordinate_records:
        offset = 8 * component
        value = struct.unpack_from("<d", records[index], offset)[0]
        if not math.isfinite(value):
            raise ValueError("Nonfinite starting coordinate")
        struct.pack_into("<d", records[index], offset, value + delta)
    return write_records(records)


def main():
    parser = argparse.ArgumentParser(description="Displace one atom in a CP-PAW restart.")
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("atom", type=int, help="one-based atom index")
    parser.add_argument("axis", type=int, choices=(1, 2, 3))
    parser.add_argument("delta", type=float, help="Cartesian displacement in bohr")
    args = parser.parse_args()

    args.output.write_bytes(displace(args.input.read_bytes(), args.atom, args.axis, args.delta))


if __name__ == "__main__":
    main()
