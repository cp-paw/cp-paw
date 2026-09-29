# Periodic Si2 checks

`skala_si2.strc` is the two-atom primitive diamond cell, not an isolated Si2
molecule. It uses eight k points. `run_isolated.sh` performs a cold one-step
integration test in a fresh directory; this is not an SCF convergence test.

`RADIALPOINTS` selects the model atom-grid radial rule; `LEBEDEVEXACTNESS`
selects the minimum angular polynomial order. The maximum is 65, with 1454
directions; a request for 64 selects 65. The default remains 53. Orientations
average rotated rules without increasing their polynomial exactness.

## Electronic stationarity

Enable `!DFT!SKALA CHECK=T` and inspect the completed protocol:

```sh
python3 tests/fulltests/skala_si2/stationarity.py /path/to/si2.prot
python3 tests/fulltests/skala_si2/stationarity.py /path/to/si2.prot \
  --residual-tolerance 1e-6 --commutator-tolerance 1e-6 --last 3
python3 -m unittest discover -s tests/fulltests/skala_si2 -p test_stationarity.py
```

Without the two explicit convergence tolerances, the script checks diagnostic
consistency only. The convergence mode limits the **maximum occupied residual**
and occupation commutator, in Hartree, over the requested final window; the RMS
is also reported. Overlap and Hermiticity limits apply to every recorded step.
The example tolerances are not universal force/stress accuracy guarantees.
Neither mode certifies the global electronic minimum or grid convergence.
`APPLY=F` measures the conventional-XC Hamiltonian; `APPLY=T` measures Skala.

A PBE restart must have the same cell, positions, k mesh, band layout and
Fourier basis as the Skala probe. A small time step preserves that starting
state but cannot make it stationary for a different functional. A one-k-point
restart cannot be used for this eight-k-point test.

## PAW setup grid

The partial-wave setup grid is independent of the Skala atom-grid quadrature.
The built-in `SI_.75_6.0` uses `DMIN=1.E-6 DMAX=0.1 RMAX=20.`. To refine it,
provide a complete `!SPECIES!AUGMENT` block, including its `!GRID` block, or use
the current external `AUGPARMS` format (`!ACNTL!AUGMENT`). Check the resolved
`*.strc_out` file to verify which settings were actually used.

The historical `PARMS_STP` file attachment in older test controls does **not**
override the built-in setup. Changing that `stp.cntl` alone therefore does not
constitute a radial convergence test. The explicit Si audit used identical HBS
parameters with `DMIN/DMAX` halved and quartered in the structure's `AUGMENT`
block. See [the validation checkpoint](../skala_crystals/VALIDATION.md) for
results and unresolved physical acceptance gates.
