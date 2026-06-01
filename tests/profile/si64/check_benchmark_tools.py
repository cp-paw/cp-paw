#!/usr/bin/env python3
import csv
import os
import subprocess
import sys
import tempfile


HERE = os.path.dirname(os.path.abspath(__file__))


PROFILE_HEADER = (
    "op,n1,n2,n3,n4,calls,total_seconds,max_seconds,avg_seconds,"
    "gflop,gbyte,measured_gflop_per_s,measured_gbyte_per_s"
)


def run_tool(*args):
    return subprocess.run(
        [sys.executable, os.path.join(HERE, args[0]), *args[1:]],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    ).stdout


def write(path, text):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def assert_equal(actual, expected, label):
    if actual != expected:
        raise AssertionError(f"{label}: expected {expected!r}, got {actual!r}")


def assert_contains(text, needle, label):
    if needle not in text:
        raise AssertionError(f"{label}: missing {needle!r}")


def check_summary_and_markdown(tmpdir):
    run_dir = os.path.join(tmpdir, "case", "rep1")
    os.makedirs(run_dir)
    write(
        os.path.join(run_dir, "test_profile.csv"),
        "\n".join(
            [
                PROFILE_HEADER,
                "ACC_COPY_HPSI_PSI0_IN,1,1,1,0,1,0,0,0,0,1.5,0,0",
                "ACC_COPY_PROJ_SETUP0_PSI_IN,1,1,1,0,1,0,0,0,0,0.25,0,0",
                "ACC_COPY_PROJ_PROPSI_OUT,1,1,1,0,1,0,0,0,0,2.5,0,0",
                "ACC_COPY_ADDPRO_HPSI_PROPSI_IN,1,1,1,0,1,0,0,0,0,0.5,0,0",
                "ACC_COPY_OFFDEN_DPACK_PROJ_IN,1,1,1,0,1,0,0,0,0,3.5,0,0",
                "ACC_COPY_DENMAT_LAGR_IN,1,1,1,0,1,0,0,0,0,4.5,0,0",
                "ACC_COPY_MISC_IN,1,1,1,0,1,0,0,0,0,0.75,0,0",
                "PW_GTOR_TOTAL,1,1,1,0,1,1,1,1,0,0,0,0",
                "PW_RTOG_TOTAL,1,1,1,0,1,2,2,2,0,0,0,0",
                "PW_FFT_RTOG_TOTAL,1,1,1,0,1,4,4,4,0,0,0,0",
            ]
        )
        + "\n",
    )
    write(os.path.join(run_dir, "time.txt"), "real 10\n")
    write(
        os.path.join(run_dir, "out.log"),
        "CONSTANT ENERGY 302.280854\nNORMAL STOP\n",
    )
    write(
        os.path.join(run_dir, "run.env"),
        "ranks=1\nnsteps=1\nexpected_energy=302.280854\n",
    )

    tsv_path = os.path.join(tmpdir, "benchmark.tsv")
    write(tsv_path, run_tool("benchmark_summary.py", tmpdir))
    rows = read_tsv(tsv_path)
    assert_equal(len(rows), 1, "benchmark row count")
    row = rows[0]
    assert_equal(row["pw_gtor_s"], "1", "pw_gtor_s")
    assert_equal(row["pw_rtog_s"], "2", "pw_rtog_s")
    assert_equal(row["fft_s"], "4", "fft_s")
    assert_equal(row["copy_gb"], "13.5", "copy_gb")
    assert_equal(row["copy_wave_gb"], "1.75", "copy_wave_gb")
    assert_equal(row["copy_proj_gb"], "3", "copy_proj_gb")
    assert_equal(row["copy_offden_gb"], "3.5", "copy_offden_gb")
    assert_equal(row["copy_denmat_gb"], "4.5", "copy_denmat_gb")

    markdown = run_tool("benchmark_markdown.py", tsv_path)
    assert_contains(markdown, "copy_wave_gb", "benchmark markdown header")
    assert_contains(markdown, "|  | case | 1 | yes | 10.00 |", "benchmark markdown row")

    profile_summary = run_tool(
        "profile_summary.py", os.path.join(run_dir, "test_profile.csv")
    )
    assert_contains(profile_summary, "ACC copy wave", "profile copy wave bucket")
    assert_contains(profile_summary, "ACC copy proj", "profile copy proj bucket")
    assert_contains(profile_summary, "ACC copy offden", "profile copy offden bucket")
    assert_contains(profile_summary, "ACC copy denmat", "profile copy denmat bucket")


def check_compare(tmpdir):
    path = os.path.join(tmpdir, "combined.tsv")
    header = (
        "suite\tcase\trepeat\tnsteps\tranks\tok\twall_s\twall_rank_s\t"
        "rank_s\tgap_s\tcoverage_pct\tpaw_s\tblas_s\tlapack_s\tfft_s\tmpi_s\t"
        "pw_trace_s\tpw_gtor_s\tpw_rtog_s\tphase_s\tphase_gap_s\tsetup_s\t"
        "copy_gb\tcopy_wave_gb\tcopy_proj_gb\tcopy_offden_gb\t"
        "copy_denmat_gb\tenergy\tenergy_delta\tenergy_ok\tenv"
    )
    rows = [
        "gpu_acc_1steps_1rank\tgpu_resident_stack\trep1\t1\t1\tyes\t25\t25\t20\t5\t80\t0\t0\t0\t0\t0\t0\t0\t0\t20\t5\t0\t10\t4\t3\t2\t1\t302.280854\t0\tyes\t",
        "gpu_acc_1steps_1rank\tnvhpc_cpu\trep1\t1\t1\tyes\t100\t100\t80\t20\t80\t0\t0\t0\t0\t0\t0\t0\t0\t80\t20\t0\t0\t0\t0\t0\t0\t302.280854\t0\tyes\t",
        "cpu_ref_1steps_8ranks\tnvhpc_cpu\trep1\t1\t8\tyes\t40\t320\t260\t60\t81\t0\t0\t0\t0\t0\t0\t0\t0\t260\t60\t0\t0\t0\t0\t0\t0\t302.280854\t0\tyes\t",
        "si64_bands_empty128_1steps_1ranks_gpu\tgpu_resident_stack\trep1\t1\t1\tyes\t30\t30\t24\t6\t80\t0\t0\t0\t0\t0\t0\t0\t0\t24\t6\t0\t11\t5\t3\t2\t1\t302.280854\t0\tyes\t",
        "si64_bands_empty128_1steps_1rank_cpu\tnvhpc_cpu\trep1\t1\t1\tyes\t90\t90\t70\t20\t78\t0\t0\t0\t0\t0\t0\t0\t0\t70\t20\t0\t0\t0\t0\t0\t0\t302.280854\t0\tyes\t",
        "si64_bands_empty128_1steps_8ranks_cpu_ref\tnvhpc_cpu\trep1\t1\t8\tyes\t45\t360\t280\t80\t78\t0\t0\t0\t0\t0\t0\t0\t0\t280\t80\t0\t0\t0\t0\t0\t0\t302.280854\t0\tyes\t",
        "empty128_1rank_cusolver\tcusolver_generalized\trep1\t1\t1\tyes\t25\t25\t20\t5\t80\t0\t0\t0\t0\t0\t0\t0\t0\t20\t5\t0\t0\t0\t0\t0\t0\t302.280854\t0\tyes\t",
        "empty128_1rank_cpu\tnvhpc_cpu\trep1\t1\t1\tyes\t100\t100\t80\t20\t80\t0\t0\t0\t0\t0\t0\t0\t0\t80\t20\t0\t0\t0\t0\t0\t0\t302.280854\t0\tyes\t",
        "empty128_8rank_cpu_ref\tnvhpc_cpu\trep1\t1\t8\tyes\t40\t320\t260\t60\t81\t0\t0\t0\t0\t0\t0\t0\t0\t260\t60\t0\t0\t0\t0\t0\t0\t302.280854\t0\tyes\t",
        "empty512-nstep2-1r\tgpu_resident\trep1\t2\t1\tyes\t20\t20\t16\t4\t80\t0\t0\t0\t0\t0\t0\t0\t0\t16\t4\t0\t8\t3\t2\t2\t1\t302.280854\t0\tyes\t",
        "empty512-nstep2-1r\tgpu_resident_stack\trep1\t2\t1\tyes\t10\t10\t8\t2\t80\t0\t0\t0\t0\t0\t0\t0\t0\t8\t2\t0\t4\t1.5\t1\t1\t0.5\t302.280854\t0\tyes\t",
        "gpu_1rank\tgpu_resident_stack\trep1\t1\t1\tyes\t20\t20\t16\t4\t80\t0\t0\t0\t0\t0\t0\t0\t0\t16\t4\t0\t8\t3\t2\t2\t1\t302.280854\t0\tyes\t",
        "cpu_1rank\tnvhpc_cpu\trep1\t1\t1\tyes\t80\t80\t60\t20\t75\t0\t0\t0\t0\t0\t0\t0\t0\t60\t20\t0\t0\t0\t0\t0\t0\t302.280854\t0\tyes\t",
        "cpu_8rank_ref\tnvhpc_cpu\trep1\t1\t8\tyes\t32\t256\t210\t46\t82\t0\t0\t0\t0\t0\t0\t0\t0\t210\t46\t0\t0\t0\t0\t0\t0\t302.280854\t0\tyes\t",
    ]
    write(path, "\n".join([header, *rows]) + "\n")

    compare = run_tool("benchmark_compare.py", path)
    assert_contains(compare, "| gpu_acc_1steps_1rank | gpu_resident_stack | 1 | yes | 25.00 | gpu_resident_stack | 1.00 | 4.00 | 1.60 |", "gpu_acc speedups")
    assert_contains(compare, "| si64_bands_empty128_1steps_1ranks_gpu | gpu_resident_stack | 1 | yes | 30.00 | gpu_resident_stack | 1.00 | 3.00 | 1.50 |", "band speedups")
    assert_contains(compare, "| empty128_1rank_cusolver | cusolver_generalized | 1 | yes | 25.00 | cusolver_generalized | 1.00 | 4.00 | 1.60 |", "cusolver speedups")
    assert_contains(compare, "| empty=512, nsteps=2 | gpu_resident_stack | 1 | yes | 10.00 | gpu_resident | 2.00 |  |  |", "gpu-only base speedup")
    assert_contains(compare, "| gpu_1rank | gpu_resident_stack | 1 | yes | 20.00 | gpu_resident_stack | 1.00 | 4.00 | 1.60 |", "standard speedups")


def main():
    with tempfile.TemporaryDirectory() as tmpdir:
        check_summary_and_markdown(tmpdir)
        check_compare(tmpdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
