#!/usr/bin/env python3

import argparse
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
        records.append(bytearray(data[offset + 4:end]))
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


def matrix_product(left, right):
    return [
        [sum(left[i][k] * right[k][j] for k in range(3)) for j in range(3)]
        for i in range(3)
    ]


def transform_vector(matrix, vector):
    return [sum(matrix[i][j] * vector[j] for j in range(3)) for i in range(3)]


def unpack_fortran_matrix(record, offset):
    values = struct.unpack_from("<9d", record, offset)
    return [[values[i + 3 * j] for j in range(3)] for i in range(3)]


def pack_fortran_matrix(record, offset, matrix):
    values = [matrix[i][j] for j in range(3) for i in range(3)]
    struct.pack_into("<9d", record, offset, *values)


def main():
    parser = argparse.ArgumentParser(
        description="Apply an affine strain to a CP-PAW restart cell and atoms."
    )
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "strain",
        nargs=9,
        type=float,
        metavar=(
            "E11", "E12", "E13", "E21", "E22", "E23", "E31", "E32", "E33"
        ),
        help="row-major displacement-gradient components; F = I + strain",
    )
    args = parser.parse_args()

    strain = [args.strain[3 * i : 3 * i + 3] for i in range(3)]
    deformation = [
        [strain[i][j] + (1.0 if i == j else 0.0) for j in range(3)]
        for i in range(3)
    ]
    records = read_records(args.input.read_bytes())

    cell_header = next(
        (index for index, record in enumerate(records) if record_name(record) == "CELL"),
        None,
    )
    if cell_header is None or len(records[cell_header + 1]) != 27 * 8:
        raise ValueError("unexpected or missing CELL section")
    cell_record = records[cell_header + 1]
    for offset in (0, 9 * 8, 18 * 8):
        cell = unpack_fortran_matrix(cell_record, offset)
        pack_fortran_matrix(cell_record, offset, matrix_product(deformation, cell))

    atom_header = next(
        (index for index, record in enumerate(records) if record_name(record) == "ATOMS"),
        None,
    )
    if atom_header is None:
        raise ValueError("ATOMS section not found")
    natom = struct.unpack("<i", records[atom_header + 1])[0]
    expected_size = 3 * natom * 8
    for index in (atom_header + 3, atom_header + 4):
        if len(records[index]) != expected_size:
            raise ValueError("unexpected ATOMS coordinate record size")
        coordinates = list(struct.unpack(f"<{3 * natom}d", records[index]))
        for atom in range(natom):
            begin = 3 * atom
            coordinates[begin : begin + 3] = transform_vector(
                deformation, coordinates[begin : begin + 3]
            )
        struct.pack_into(f"<{3 * natom}d", records[index], 0, *coordinates)

    args.output.write_bytes(write_records(records))


if __name__ == "__main__":
    main()
