#!/usr/bin/env python3
"""Check profile output names and counters using a compiled profiling build."""
import argparse
import csv
import math
import os
from pathlib import Path
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('--mpi-ranks', type=int, default=1)
    args = parser.parse_args()
    if args.mpi_ranks < 1:
        parser.error('--mpi-ranks must be positive')
    command = [str(args.executable.resolve())]
    if args.mpi_ranks > 1:
        command = ['mpirun', '-np', str(args.mpi_ranks), *command]
    cases = [('unset', None, 'cppaw_accel_profile'),
             ('empty', '', 'cppaw_accel_profile'),
             ('blank', '   ', 'cppaw_accel_profile'),
             ('explicit', 'custom_profile', 'custom_profile'),
             ('padded', '  custom profile  ', 'custom profile'),
             ('disabled', 'disabled_profile', None)]
    environment = {key: value for key, value in os.environ.items()
                   if not key.startswith('CPPAW_')}
    environment.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    for name, value, expected in cases:
        selected = dict(environment, CPPAW_ACCEL_PROFILE='0' if name == 'disabled' else '1')
        if value is not None:
            selected['CPPAW_ACCEL_PROFILE_FILE'] = value
        with tempfile.TemporaryDirectory(prefix='cppaw-profile-') as directory:
            work = Path(directory)
            run = subprocess.run(command, cwd=work, env=selected, check=True,
                                 capture_output=True, text=True, timeout=60)
            expected_files = set()
            if expected is not None:
                expected_files = ({expected + '.csv'} if args.mpi_ranks == 1 else
                                  {f'{expected}.rank{rank:05d}.csv'
                                   for rank in range(1, args.mpi_ranks + 1)})
            actual_files = {path.name for path in work.glob('*.csv')}
            if actual_files != expected_files:
                raise AssertionError((name, actual_files, expected_files, run.stdout, run.stderr))
            for filename in sorted(expected_files):
                with (work / filename).open(newline='') as stream:
                    rows = list(csv.DictReader(stream))
                if len(rows) != 1 or rows[0]['op'] != 'PROFILE_REPORT_TEST':
                    raise AssertionError((name, rows))
                for field, number in dict(n1=2, n2=3, n3=4, n4=5, calls=2,
                                          total_seconds=0.75, max_seconds=0.5,
                                          avg_seconds=0.375, gflop=48e-9, gbyte=96e-9).items():
                    if not math.isclose(float(rows[0][field]), number, rel_tol=1e-8, abs_tol=0):
                        raise AssertionError((name, field, rows[0][field], number))
            print(f'PROFILE REPORT {name}: passed ({len(expected_files)} files)')


if __name__ == '__main__':
    main()
