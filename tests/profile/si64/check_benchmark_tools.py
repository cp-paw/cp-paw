#!/usr/bin/env python3
import csv
import os
import re
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


def markdown_table(text):
    rows = [line for line in text.splitlines() if line.startswith("| ")]
    if len(rows) < 3:
        raise AssertionError("markdown table: no data rows")
    header = [cell.strip() for cell in rows[0].strip("|").split("|")]
    data = [
        [cell.strip() for cell in row.strip("|").split("|")]
        for row in rows[2:]
    ]
    return header, data


def case_env_cases(script_name):
    path = os.path.join(HERE, script_name)
    with open(path, encoding="utf-8") as handle:
        text = handle.read()

    in_case_env = False
    cases = set()
    for line in text.splitlines():
        if line.startswith("case_env()"):
            in_case_env = True
            continue
        if in_case_env and line.strip() == "esac":
            break
        if not in_case_env:
            continue

        match = re.match(r"\s*([A-Za-z0-9_|*]+)\)\s+echo", line)
        if match:
            cases.update(
                case for case in match.group(1).split("|") if "*" not in case
            )
    return cases


def check_nsys_case_coverage():
    benchmark_cases = case_env_cases("run_benchmark.sh")
    nsys_cases = case_env_cases("run_nsys.sh")
    missing_nsys = sorted(benchmark_cases - nsys_cases)
    if missing_nsys:
        raise AssertionError(
            "benchmark case_env cases missing from Nsight case_env: "
            + ", ".join(missing_nsys)
        )


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
                "ACC_COPY_SKALA_GRID_BACK_RESIDENCY_IN,1,1,1,0,1,0,0,0,0,0.25,0,0",
                "ACC_COPY_SKALA_GRID_BACK_INPUT,1,1,1,0,1,0,0,0,0,0.5,0,0",
                "ACC_COPY_SKALA_GRID_BACK_RESIDENCY_OUT,1,1,1,0,1,0,0,0,0,0.75,0,0",
                "ACC_COPY_MISC_IN,1,1,1,0,1,0,0,0,0,0.75,0,0",
                "ACC_COPY_SERIAL3D_ACC_MAP,1,1,1,0,1,0,0,0,0,4.0,0,0",
                "ACC_COPY_SERIAL3D_ACC_INPUT,1,1,1,0,1,0,0,0,0,1.0,0,0",
                "ACC_COPY_SERIAL3D_ACC_OUTPUT,1,1,1,0,1,0,0,0,0,2.0,0,0",
                "ACC_COPY_SERIAL3D_ACC_MAP_META,1,1,1,0,1,0,0,0,0,0.5,0,0",
                "ACC_UPDATE_VPSI_HPSI_IN,1,1,1,0,1,0,0,0,0,1.25,0,0",
                "ACC_PRESENT_ADDOPSI_PSIM,1,1,1,0,2,0,0,0,0,0,0,0",
                "PW_GTOR_TOTAL,1,1,1,0,1,1,1,1,0,0,0,0",
                "PW_RTOG_TOTAL,1,1,1,0,1,2,2,2,0,0,0,0",
                "PW_FFT_GTOR_TOTAL,1,1,1,0,1,3,3,3,0,0,0,0",
                "PW_FFT_RTOG_TOTAL,1,1,1,0,1,4,4,4,0,0,0,0",
                "FFT1D_C8,1,1,1,0,3,6.5,3,2.1666667,0,0,0,0",
                "PAW_VPSI_FFT_GTOR,1,1,1,0,1,5,5,5,0,0,0,0",
                "PAW_VPSI_FFT_RTOG,1,1,1,0,1,6,6,6,0,0,0,0",
                "PAW_VPSI_TOTAL,1,1,1,0,1,12,12,12,0,0,0,0",
                "SKALA_PARTITION,1,1,1,0,1,0.5,0.5,0.5,0,0,0,0",
                "SKALA_MODEL,1,1,1,0,2,1.5,1,0.75,0,0,0,0",
                "CUBLAS_OZAKI_DGEMM_USED,64,1,2048,0,3,0,0,0,0,0,0,0",
                "CUBLAS_OZAKI_DGEMM_FALLBACK,-1,1,2048,0,2,0,0,0,0,0,0,0",
                "CUBLAS_OZAKI_ZGEMM_USED,72,1,2048,0,4,0,0,0,0,0,0,0",
                "CUBLAS_OZAKI_ZHERK_FALLBACK,-1,1,2048,0,5,0,0,0,0,0,0,0",
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
        os.path.join(run_dir, "case.prot"),
        "HYBRID-GRID ROWS 1348944\n"
        "PARTITION TRANSLATION CLASSES 2\n"
        "MODEL XC ENERGY -1.21801589257683E+03\n"
        "TOTAL ENERGY : -8.483135182E+02 H\n",
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
    assert_equal(row["fft_s"], "7", "fft_s")
    assert_equal(row["fft_kernel_s"], "6.5", "fft_kernel_s")
    assert_equal(row["pw_fft_gtor_s"], "3", "pw_fft_gtor_s")
    assert_equal(row["pw_fft_rtog_s"], "4", "pw_fft_rtog_s")
    assert_equal(row["vpsi_s"], "12", "vpsi_s")
    assert_equal(row["vpsi_gtor_s"], "5", "vpsi_gtor_s")
    assert_equal(row["vpsi_rtog_s"], "6", "vpsi_rtog_s")
    assert_equal(row["rank_s"], "7", "nested Skala exclusion")
    assert_equal(row["skala_s"], "2", "skala_s")
    assert_equal(row["skala_partition_s"], "0.5", "skala_partition_s")
    assert_equal(row["skala_model_s"], "1.5", "skala_model_s")
    assert_equal(row["hybrid_grid_rows"], "1348944", "hybrid_grid_rows")
    assert_equal(row["partition_classes"], "2", "partition_classes")
    assert_equal(row["model_xc_energy"], "-1218.01589", "model_xc_energy")
    assert_equal(row["static_total_energy"], "-848.313518", "static_total_energy")
    assert_equal(row["transfer_gb"], "20.25", "transfer_gb")
    assert_equal(row["copy_gb"], "19", "copy_gb")
    assert_equal(row["copy_wave_gb"], "1.75", "copy_wave_gb")
    assert_equal(row["copy_proj_gb"], "3", "copy_proj_gb")
    assert_equal(row["copy_offden_gb"], "3.5", "copy_offden_gb")
    assert_equal(row["copy_denmat_gb"], "4.5", "copy_denmat_gb")
    assert_equal(row["copy_skala_gb"], "1.5", "copy_skala_gb")
    assert_equal(row["update_gb"], "1.25", "update_gb")
    assert_equal(row["update_wave_gb"], "1.25", "update_wave_gb")
    assert_equal(row["update_proj_gb"], "0", "update_proj_gb")
    assert_equal(row["ozaki_dgemm_used"], "3", "ozaki_dgemm_used")
    assert_equal(row["ozaki_dgemm_fallback"], "2", "ozaki_dgemm_fallback")
    assert_equal(row["ozaki_dgemm_max_bits"], "64", "ozaki_dgemm_max_bits")
    assert_equal(row["ozaki_zgemm_used"], "4", "ozaki_zgemm_used")
    assert_equal(row["ozaki_zgemm_max_bits"], "72", "ozaki_zgemm_max_bits")
    assert_equal(row["ozaki_zherk_fallback"], "5", "ozaki_zherk_fallback")

    markdown = run_tool("benchmark_markdown.py", tsv_path)
    assert_contains(markdown, "transfer_gb", "benchmark markdown transfer header")
    assert_contains(markdown, "copy_wave_gb", "benchmark markdown header")
    assert_contains(markdown, "copy_skala_gb", "benchmark markdown Skala copy header")
    assert_contains(markdown, "update_wave_gb", "benchmark markdown update header")
    assert_contains(markdown, "skala_atom_grid_s", "benchmark markdown Skala header")
    assert_contains(markdown, "model_xc_energy", "benchmark markdown energy header")
    assert_contains(markdown, "ozaki_zgemm_used", "benchmark markdown Ozaki header")
    assert_contains(markdown, "|  | case | 1 | yes | 10.00 |", "benchmark markdown row")

    profile_summary = run_tool(
        "profile_summary.py", os.path.join(run_dir, "test_profile.csv")
    )
    assert_contains(profile_summary, "ACC copy wave", "profile copy wave bucket")
    assert_contains(profile_summary, "ACC copy proj", "profile copy proj bucket")
    assert_contains(profile_summary, "ACC copy offden", "profile copy offden bucket")
    assert_contains(profile_summary, "ACC copy denmat", "profile copy denmat bucket")
    assert_contains(profile_summary, "ACC copy skala", "profile copy Skala bucket")
    assert_contains(profile_summary, "ACC update wave", "profile update wave bucket")
    assert_contains(profile_summary, "transfer estimate: 20.250000 GB", "profile transfer total")
    assert_contains(
        profile_summary,
        "Skala detail rank-seconds (nested): 2.000000",
        "profile nested Skala detail",
    )
    assert_contains(
        profile_summary,
        "FFT kernel rank-seconds (nested): 6.500000",
        "profile nested FFT kernel detail",
    )
    assert_contains(
        profile_summary,
        "cuBLAS FP64 Ozaki dispatch",
        "profile Ozaki dispatch section",
    )

    copy_rows = run_tool(
        "profile_copy_rows.py",
        "--op-prefix",
        "ACC_COPY",
        "--op-prefix",
        "ACC_UPDATE",
        "--op-prefix",
        "ACC_PRESENT",
        "--include-zero",
        "--top",
        "20",
        os.path.join(run_dir, "test_profile.csv"),
    )
    assert_contains(
        copy_rows,
        "suite\tcase\trepeat\top\tcalls\tgbyte\tseconds\tfiles",
        "profile row TSV header",
    )
    assert_contains(
        copy_rows,
        "ACC_UPDATE_VPSI_HPSI_IN\t1\t1.25",
        "profile row update bytes",
    )
    assert_contains(
        copy_rows,
        "ACC_COPY_SERIAL3D_ACC_INPUT\t1\t1",
        "profile row detail bytes",
    )
    assert_contains(
        copy_rows,
        "ACC_PRESENT_ADDOPSI_PSIM\t2\t0",
        "profile row present zero bytes",
    )

    copy_rows_md = run_tool(
        "profile_copy_rows.py",
        "--op-prefix",
        "ACC_UPDATE",
        "--markdown",
        os.path.join(run_dir, "test_profile.csv"),
    )
    assert_contains(copy_rows_md, "| suite | case | repeat | op | calls | gbyte | seconds | files |", "profile row markdown header")


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
    assert_contains(compare, "| gpu_acc_1steps_1rank | gpu_resident_stack | 1 | yes | 25.00 | gpu_resident_stack | 1.00 | 4.00 | 1.60 | 20.00 | 10.00 | 10.00 |", "gpu_acc speedups and legacy transfer")
    assert_contains(compare, "| si64_bands_empty128_1steps_1ranks_gpu | gpu_resident_stack | 1 | yes | 30.00 | gpu_resident_stack | 1.00 | 3.00 | 1.50 |", "band speedups")
    assert_contains(compare, "| empty128_1rank_cusolver | cusolver_generalized | 1 | yes | 25.00 | cusolver_generalized | 1.00 | 4.00 | 1.60 |", "cusolver speedups")
    assert_contains(compare, "| empty=512, nsteps=2 | gpu_resident_stack | 1 | yes | 10.00 | gpu_resident | 2.00 |  |  |", "gpu-only base speedup")
    assert_contains(compare, "| gpu_1rank | gpu_resident_stack | 1 | yes | 20.00 | gpu_resident_stack | 1.00 | 4.00 | 1.60 |", "standard speedups")

    markdown = run_tool("benchmark_markdown.py", path)
    header_cells, data_rows = markdown_table(markdown)
    transfer_index = header_cells.index("transfer_gb")
    copy_index = header_cells.index("copy_gb")
    assert_equal(data_rows[0][transfer_index], "10.00", "legacy markdown transfer_gb")
    assert_equal(data_rows[0][copy_index], "10.00", "legacy markdown copy_gb")


def check_ozaki_validate(tmpdir):
    reference = os.path.join(tmpdir, "ozaki-reference")
    candidate = os.path.join(tmpdir, "ozaki-candidate")
    os.makedirs(reference)
    os.makedirs(candidate)
    common = (
        "CPPAW FP64 ORTHO RESIDUAL 5.0000000000000000E-09\n"
        "CPPAW FP64 FORCE ACTIVE T\n"
        "CPPAW FP64 FORCE      1 1.0000000000000000E-03 2.0000000000000000E-03 3.0000000000000000E-03\n"
        "CPPAW FP64 STRESS ACTIVE T\n"
        "CPPAW FP64 STRESS  1 1.0000000000000000E-02 2.0000000000000000E-02 3.0000000000000000E-02\n"
        "CPPAW FP64 STRESS  2 4.0000000000000000E-02 5.0000000000000000E-02 6.0000000000000000E-02\n"
        "CPPAW FP64 STRESS  3 7.0000000000000000E-02 8.0000000000000000E-02 9.0000000000000000E-02\n"
    )
    write(
        os.path.join(reference, "out.log"),
        "CONSTANT ENERGY 1.0000000000000000E+00 0.0\n" + common,
    )
    write(
        os.path.join(candidate, "out.log"),
        "CONSTANT ENERGY 1.0000000010000000E+00 0.0\n" + common,
    )
    protocol = "NUMBER OF K-POINTS.....................................: 2\n"
    write(os.path.join(reference, "test.prot"), protocol)
    write(os.path.join(candidate, "test.prot"), protocol)
    write(
        os.path.join(candidate, "test_profile.csv"),
        PROFILE_HEADER
        + "\nCUBLAS_OZAKI_ZGEMM_USED,64,2,8192,0,4,0,0,0,0,0,0,0\n",
    )
    result = run_tool(
        "ozaki_validate.py",
        reference,
        candidate,
        "--expected-kpoints",
        "2",
        "--require-force",
        "--require-stress",
        "--require-ozaki",
        "ZGEMM",
    )
    assert_contains(result, "ozaki_validation=pass", "Ozaki validation")
    assert_contains(result, "kernel=ZGEMM used=4", "Ozaki dispatch count")


def main():
    with tempfile.TemporaryDirectory() as tmpdir:
        check_summary_and_markdown(tmpdir)
        check_compare(tmpdir)
        check_ozaki_validate(tmpdir)
    check_nsys_case_coverage()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
