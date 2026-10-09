"""Shared conventional CP-PAW restart and trajectory helpers."""

import hashlib
import json
import math
import struct


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


def wave_records(data):
    records = read_records(data)
    indices = [i for i, row in enumerate(records) if record_name(row) == "WAVES"]
    if len(indices) != 1:
        raise ValueError("Missing or duplicate WAVES restart section")
    start = indices[0]
    count, = struct.unpack_from("<i", records[start])
    if count < 1 or start + count >= len(records):
        raise ValueError("Invalid WAVES record count")
    return records[start + 1:start + count + 1]


def wave_difference(a, b, *, cell_tolerance=0.):
    """Compare numeric wave/Lambda payloads, requiring identical binary metadata."""
    if not math.isfinite(cell_tolerance) or cell_tolerance < 0:
        raise ValueError("Cell tolerance must be finite and nonnegative")
    if (len(a) != len(b) or not a or len(a[0]) != 84 or len(b[0]) != 84
            or a[0][:8] != b[0][:8] or a[0][80:] != b[0][80:]):
        raise ValueError("Changed electronic restart metadata")
    cells = [struct.unpack_from("<9d", header, 8) for header in (a[0], b[0])]
    if not all(math.isfinite(x) for cell in cells for x in cell):
        raise ValueError("Nonfinite electronic restart cell")
    if max(abs(x-y) for x, y in zip(*cells)) > cell_tolerance:
        raise ValueError("Changed electronic restart cell")
    nk, ns = struct.unpack_from("<2i", a[0])
    nw, = struct.unpack_from("<i", a[0], 80)
    if min(nk, ns, nw) < 1:
        raise ValueError("Invalid electronic restart dimensions")
    index, maximum = 1, 0.

    def metadata(size):
        nonlocal index
        if index >= len(a) or len(a[index]) != size or a[index] != b[index]:
            raise ValueError("Changed electronic restart metadata")
        record = a[index]
        index += 1
        return record

    def numeric(count):
        nonlocal index, maximum
        if index >= len(a) or len(a[index]) != 8 * count or len(b[index]) != 8 * count:
            raise ValueError("Invalid electronic restart payload")
        x, y = (struct.unpack(f"<{count}d", records[index]) for records in (a, b))
        if not all(math.isfinite(v) for v in (*x, *y)):
            raise ValueError("Nonfinite electronic restart")
        maximum = max(maximum, max(abs(u - v) for u, v in zip(x, y)))
        index += 1

    for _ in range(nk):
        key, ng, ndim, nb, real = struct.unpack("<8s4i", metadata(24))
        if key.strip() != b"PSI" or min(ng, ndim, nb) < 1 or real not in (0, 1):
            raise ValueError("Invalid PSI restart header")
        metadata(24 + 12 * ng)
        nbh = (nb + 1) // 2 if real else nb
        for _ in range(nw * ns * nbh):
            numeric(2 * ng * ndim)
        key, n, nd, nspin, nl = struct.unpack("<8s4i", metadata(24))
        if key.strip() != b"LAMBDA" or (n, nd, nspin) != (nb, ndim, ns) or nl < 1:
            raise ValueError("Invalid LAMBDA restart header")
        for _ in range(ns * nl):
            numeric(2 * nb * nb)
    if index != len(a):
        raise ValueError("Unparsed electronic restart records")
    return maximum


def trajectory_difference(first, second):
    a, b = read_records(first), read_records(second)
    if len(a) != len(b) or not a:
        raise ValueError("Changed or empty trajectory")
    difference = 0.
    for x, y in zip(a, b):
        if len(x) < 16 or x[:16] != y[:16]:
            raise ValueError("Changed trajectory step, time or size")
        _, time, count = struct.unpack("<idi", x[:16])
        if not math.isfinite(time) or count < 1 or len(x) != 16 + 8 * count or len(y) != len(x):
            raise ValueError("Invalid trajectory payload")
        u, v = (struct.unpack(f"<{count}d", row[16:]) for row in (x, y))
        if not all(math.isfinite(value) for value in (*u, *v)):
            raise ValueError("Nonfinite trajectory")
        difference = max(difference, max(abs(p - q) for p, q in zip(u, v)))
    return difference


def geometry(data):
    records = read_records(data)
    sections = {}
    for name in ("ATOMS", "CELL"):
        indices = [i for i, row in enumerate(records) if record_name(row) == name]
        if len(indices) != 1:
            raise ValueError(f"Missing or duplicate {name} section")
        sections[name] = indices[0]
    atom, cell = sections["ATOMS"], sections["CELL"]
    natom, = struct.unpack("<i", records[atom + 1])
    if natom < 1 or len(records[cell + 1]) != 27 * 8:
        raise ValueError("Invalid restart geometry layout")
    positions = [list(struct.unpack(f"<{3 * natom}d", records[atom + offset]))
                 for offset in (3, 4)]
    cells = list(struct.unpack("<27d", records[cell + 1]))
    if not all(math.isfinite(x) for x in [*cells, *positions[0], *positions[1]]):
        raise ValueError("Nonfinite restart geometry")
    return {"natom": natom, "positions": positions, "cells": cells}


def displace_geometry(data, atom, axis, delta, rigid_translation=False):
    """Translate all nuclei, or one nucleus, at both coordinate time levels."""
    if not rigid_translation:
        return displace(data, atom, axis, delta)
    for index in range(1, geometry(data)["natom"] + 1):
        data = displace(data, index, axis, delta)
    return data


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def digest(path):
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()
