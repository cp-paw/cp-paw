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


def deform(data, components):
    """Strain all cell/coordinate time levels, preserving electronic records."""
    if len(components) != 9 or not all(math.isfinite(x) for x in components):
        raise ValueError("Need nine finite row-major strain components")
    strain = [components[3 * i : 3 * i + 3] for i in range(3)]
    deformation = [
        [strain[i][j] + (1.0 if i == j else 0.0) for j in range(3)]
        for i in range(3)
    ]
    a, b, c = deformation
    determinant = (a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0])
                   + a[2]*(b[0]*c[1]-b[1]*c[0]))
    if not math.isfinite(determinant) or determinant <= 0:
        raise ValueError("Deformation must preserve positive volume")
    records = read_records(data)
    for name in ("CELL", "ATOMS"):
        if sum(record_name(row) == name for row in records) != 1:
            raise ValueError(f"Missing or duplicate {name} section")

    cell_header = next(
        (index for index, record in enumerate(records) if record_name(record) == "CELL"),
        None,
    )
    if cell_header + 1 >= len(records) or len(records[cell_header + 1]) != 27 * 8:
        raise ValueError("unexpected or missing CELL section")
    cell_record = records[cell_header + 1]
    for offset in (0, 9 * 8, 18 * 8):
        cell = unpack_fortran_matrix(cell_record, offset)
        if not all(math.isfinite(x) for row in cell for x in row):
            raise ValueError("Nonfinite starting cell")
        transformed = matrix_product(deformation, cell)
        if not all(math.isfinite(x) for row in transformed for x in row):
            raise ValueError("Nonfinite deformed cell")
        pack_fortran_matrix(cell_record, offset, transformed)

    atom_header = next(
        (index for index, record in enumerate(records) if record_name(record) == "ATOMS"),
        None,
    )
    if atom_header + 4 >= len(records) or len(records[atom_header + 1]) != 4:
        raise ValueError("Truncated ATOMS section")
    natom = struct.unpack("<i", records[atom_header + 1])[0]
    if natom < 1:
        raise ValueError("Invalid atom count")
    expected_size = 3 * natom * 8
    for index in (atom_header + 3, atom_header + 4):
        if len(records[index]) != expected_size:
            raise ValueError("unexpected ATOMS coordinate record size")
        coordinates = list(struct.unpack(f"<{3 * natom}d", records[index]))
        if not all(math.isfinite(x) for x in coordinates):
            raise ValueError("Nonfinite starting coordinate")
        for atom in range(natom):
            begin = 3 * atom
            coordinates[begin : begin + 3] = transform_vector(
                deformation, coordinates[begin : begin + 3]
            )
        if not all(math.isfinite(x) for x in coordinates):
            raise ValueError("Nonfinite deformed coordinate")
        struct.pack_into(f"<{3 * natom}d", records[index], 0, *coordinates)

    return write_records(records)


def main():
    parser = argparse.ArgumentParser(
        description="Apply an affine strain to a CP-PAW restart cell and atoms."
    )
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("strain", nargs=9, type=float,
                        help="row-major displacement-gradient components; F = I + strain")
    args = parser.parse_args()
    args.output.write_bytes(deform(args.input.read_bytes(), args.strain))


if __name__ == "__main__":
    main()
