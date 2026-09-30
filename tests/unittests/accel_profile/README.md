# Profiling report regression

Use a build with `CPPVAR_ACCEL_PROFILE`, such as `profile` or
`nvhpc_gpu_profile`. No GPU or Skala model is required for this test.

```
make -C bin/Build_profile -f Makefile \
  -f ../../tests/unittests/accel_profile/driver.mk accel-profile-report-test
```

The test checks the default CSV name for an unset, empty or blank
`CPPAW_ACCEL_PROFILE_FILE`, explicit names, disabled output and exact
accumulated counters. Each case runs in its own temporary directory.
For an MPI profiling build, compile the same driver and pass its executable
to `check_report.py --mpi-ranks 2` to check rank suffixes as well.
