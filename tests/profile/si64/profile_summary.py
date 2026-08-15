#!/usr/bin/env python3
import collections
import csv
import glob
import sys


def transfer_kind(op):
    if "OFFDEN" in op:
        return "offden"
    if "DENMAT" in op:
        return "denmat"
    if any(token in op for token in ("PROPSI", "PRO_CACHE", "THIS_PROJ")):
        return "proj"
    if any(
        token in op
        for token in (
            "_PSI0",
            "_PSI1",
            "_PSI2",
            "_PSI_IN",
            "_PSI_OUT",
            "_PSIM",
            "_OPSI",
            "_HPSI",
        )
    ):
        return "wave"
    if any(token in op for token in ("PROJ", "ADDPRO")):
        return "proj"
    return "other"


def transfer_bucket(op, prefix):
    return f"{prefix} {transfer_kind(op)}"


def copy_bucket(op):
    return transfer_bucket(op, "ACC copy")


def update_bucket(op):
    return transfer_bucket(op, "ACC update")


def is_copy_detail_row(op):
    return op in (
        "ACC_COPY_SERIAL3D_ACC_INPUT",
        "ACC_COPY_SERIAL3D_ACC_OUTPUT",
        "ACC_COPY_SERIAL3D_ACC_MAP_META",
    )


def is_fft_kernel_detail(op):
    return op.startswith(("FFT1D_", "FFT3D_", "CUFFT1D_", "CUFFT3D_"))


def category(op):
    if op.startswith("SKALA_"):
        return "Skala detail"
    if op.startswith("PHASE_"):
        return "Phase trace"
    if op.startswith("ACC_COPY"):
        return copy_bucket(op)
    if op.startswith("ACC_UPDATE"):
        return update_bucket(op)
    if op.startswith("ACC_PRESENT"):
        return "ACC residency"
    if op.startswith("ACC_SETUP"):
        return "ACC setup"
    if op.startswith("PAW_"):
        return "PAW envelope"
    if op.startswith("PW_") and not op.startswith("PW_FFT"):
        if "MPE_TRANSPOSE" in op:
            return "PW MPI envelope"
        return "PW local trace"
    if op.startswith("MPI_ALLTOALL"):
        return "MPI alltoall"
    if is_fft_kernel_detail(op):
        return "FFT kernel detail"
    if op.startswith("FFT") or op.startswith("PW_FFT") or op.startswith("CUFFT"):
        return "FFT envelope"
    if op.startswith("LAPACK") or op.startswith("CUSOLVER"):
        return "LAPACK"
    if "GEMM" in op or "HERK" in op or "SYRK" in op:
        return "BLAS-like"
    return "other"


def main(argv):
    files = []
    for arg in argv[1:]:
        files.extend(glob.glob(arg))
    files = sorted(set(files))
    if not files:
        print("No profile CSV files found.", file=sys.stderr)
        return 1

    per_op = collections.defaultdict(
        lambda: {"calls": 0, "seconds": 0.0, "gflop": 0.0, "gbyte": 0.0}
    )
    per_shape = collections.defaultdict(
        lambda: {"calls": 0, "seconds": 0.0, "max": 0.0, "gflop": 0.0}
    )

    for path in files:
        with open(path, newline="") as handle:
            for row in csv.DictReader(handle):
                op = row["op"]
                calls = int(row["calls"])
                seconds = float(row["total_seconds"])
                gflop = float(row["gflop"])
                gbyte = float(row["gbyte"])
                key = (op, row["n1"], row["n2"], row["n3"], row["n4"])

                op_data = per_op[op]
                op_data["calls"] += calls
                op_data["seconds"] += seconds
                op_data["gflop"] += gflop
                op_data["gbyte"] += gbyte

                shape_data = per_shape[key]
                shape_data["calls"] += calls
                shape_data["seconds"] += seconds
                shape_data["max"] = max(shape_data["max"], float(row["max_seconds"]))
                shape_data["gflop"] += gflop

    by_category = collections.defaultdict(float)
    for op, data in per_op.items():
        by_category[category(op)] += data["seconds"]
    total = sum(by_category.values())
    primary_total = sum(
        data["seconds"] for op, data in per_op.items()
        if not op.startswith("ACC_")
        and not op.startswith("PAW_")
        and not op.startswith("SKALA_")
        and not is_fft_kernel_detail(op)
        and not (op.startswith("PW_") and not op.startswith("PW_FFT"))
        and not op.startswith("PHASE_")
    )
    phase_total = sum(
        data["seconds"] for op, data in per_op.items()
        if op.startswith("PHASE_")
    )
    paw_total = sum(
        data["seconds"] for op, data in per_op.items()
        if op.startswith("PAW_")
    )
    skala_total = sum(
        data["seconds"] for op, data in per_op.items()
        if op.startswith("SKALA_")
    )
    fft_kernel_total = sum(
        data["seconds"] for op, data in per_op.items()
        if is_fft_kernel_detail(op)
    )
    trace_total = sum(
        data["seconds"] for op, data in per_op.items()
        if op.startswith("PW_") and not op.startswith("PW_FFT")
    )
    setup_total = sum(
        data["seconds"] for op, data in per_op.items()
        if op.startswith("ACC_SETUP")
    )
    copy_gbyte = sum(
        data["gbyte"] for op, data in per_op.items()
        if op.startswith("ACC_COPY") and not is_copy_detail_row(op)
    )
    update_gbyte = sum(
        data["gbyte"] for op, data in per_op.items()
        if op.startswith("ACC_UPDATE")
    )
    copy_bucket_gbyte = collections.defaultdict(float)
    update_bucket_gbyte = collections.defaultdict(float)
    for op, data in per_op.items():
        if op.startswith("ACC_COPY") and not is_copy_detail_row(op):
            copy_bucket_gbyte[copy_bucket(op)] += data["gbyte"]
        elif op.startswith("ACC_UPDATE"):
            update_bucket_gbyte[update_bucket(op)] += data["gbyte"]

    print("Profile files: {}".format(len(files)))
    print("Instrumented rank-seconds: {:.6f}".format(primary_total))
    if paw_total:
        print("PAW envelope rank-seconds (nested): {:.6f}".format(paw_total))
    if skala_total:
        print("Skala detail rank-seconds (nested): {:.6f}".format(skala_total))
    if fft_kernel_total:
        print("FFT kernel rank-seconds (nested): {:.6f}".format(fft_kernel_total))
    if phase_total:
        print("Diagnostic phase rank-seconds: {:.6f}".format(phase_total))
    if trace_total:
        print("Diagnostic PW trace rank-seconds: {:.6f}".format(trace_total))
    if setup_total or copy_gbyte or update_gbyte:
        print(
            "Accelerator setup seconds: {:.6f}  transfer estimate: {:.6f} GB  copy estimate: {:.6f} GB  update estimate: {:.6f} GB".format(
                setup_total, copy_gbyte + update_gbyte, copy_gbyte, update_gbyte
            )
        )
    if copy_bucket_gbyte:
        print("Copy bucket estimates")
        for name, gbyte in sorted(
            copy_bucket_gbyte.items(), key=lambda item: -item[1]
        ):
            print("  {:<16s} {:10.4f} GB".format(name, gbyte))
    if update_bucket_gbyte:
        print("Update bucket estimates")
        for name, gbyte in sorted(
            update_bucket_gbyte.items(), key=lambda item: -item[1]
        ):
            print("  {:<16s} {:10.4f} GB".format(name, gbyte))
    print("")
    print("Category summary")
    for name, seconds in sorted(by_category.items(), key=lambda item: -item[1]):
        percent = 100.0 * seconds / total if total else 0.0
        print("  {:<14s} {:10.4f}s {:6.2f}%".format(name, seconds, percent))

    print("")
    print("Top operations")
    for op, data in sorted(per_op.items(), key=lambda item: -item[1]["seconds"])[:15]:
        seconds = data["seconds"]
        rate = data["gflop"] / seconds if seconds else 0.0
        bw = data["gbyte"] / seconds if seconds else 0.0
        print(
            "  {:<24s} calls={:8d} seconds={:9.4f} GF/s={:8.2f} GB/s={:8.2f}".format(
                op, data["calls"], seconds, rate, bw
            )
        )

    copy_ops = [
        (op, data) for op, data in per_op.items() if op.startswith("ACC_COPY")
    ]
    if copy_ops:
        print("")
        print("Top copy estimates")
        for op, data in sorted(copy_ops, key=lambda item: -item[1]["gbyte"])[:10]:
            print(
                "  {:<24s} calls={:8d} GB={:10.4f}".format(
                    op, data["calls"], data["gbyte"]
                )
            )

    present_ops = [
        (op, data) for op, data in per_op.items() if op.startswith("ACC_PRESENT")
    ]
    if present_ops:
        print("")
        print("OpenACC present observations")
        for op, data in sorted(present_ops, key=lambda item: item[0])[:20]:
            print("  {:<24s} calls={:8d}".format(op, data["calls"]))

    phase_ops = [
        (op, data) for op, data in per_op.items() if op.startswith("PHASE_")
    ]
    if phase_ops:
        print("")
        print("Top high-level phases")
        for op, data in sorted(phase_ops, key=lambda item: -item[1]["seconds"])[:15]:
            print(
                "  {:<24s} calls={:8d} seconds={:9.4f}".format(
                    op, data["calls"], data["seconds"]
                )
            )

    print("")
    print("Top shapes")
    for key, data in sorted(per_shape.items(), key=lambda item: -item[1]["seconds"])[:15]:
        op, n1, n2, n3, n4 = key
        seconds = data["seconds"]
        rate = data["gflop"] / seconds if seconds else 0.0
        print(
            "  {:<24s} n=({},{},{},{}) calls={:7d} seconds={:9.4f} GF/s={:8.2f}".format(
                op, n1, n2, n3, n4, data["calls"], seconds, rate
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
