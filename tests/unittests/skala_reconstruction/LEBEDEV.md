# Higher-order Lebedev rules

`LEBEDEVEXACTNESS` is a minimum algebraic exactness, not a point count.
CP-PAW selects the smallest supported rule satisfying it:

| Requested minimum | Actual exactness | Points per radial shell |
| --- | ---: | ---: |
| 48-53 | 53 | 974 |
| 54-59 | 59 | 1202 |
| 60-65, including 64 | 65 | 1454 |

The existing lower-order rules and default 53 are unchanged. Values greater
than 65 fail explicitly. The Skala report prints the requested minimum, the
selected rule's actual exactness and its point count. Multiple orientations
remain averages of rotated copies of the selected rule, not higher-degree
rules themselves.

## Coefficient provenance

The additional numerical tables are the published Lebedev-Laikov 1202- and
1454-point rules, as distributed in PySCF's `pyscf/dft/LebedevGrid.py` at commit
`e47d1279f72ffd50cbad43fbfed50511d8415b7e`:

https://github.com/pyscf/pyscf/blob/e47d1279f72ffd50cbad43fbfed50511d8415b7e/pyscf/dft/LebedevGrid.py

That distribution credits Dmitri Laikov's original coefficients and routines,
Christoph van Wuellen's Fortran translation, and Gerald Knizia's subsequent
implementation. The numerical constants are retained at their published
16-digit precision. They are mapped to CP-PAW's existing CP2K-derived
octahedral-orbit generator; PySCF is not a build or runtime dependency. The
PySCF distribution's Apache-2.0 license accompanies these tables in
`LEBEDEV-LICENSE.txt`. The pre-existing CP2K code retains its original notice.

Original reference: V. I. Lebedev and D. N. Laikov, "A quadrature formula for
the sphere of the 131st algebraic order of accuracy", Doklady Mathematics
59(3), 477-481 (1999). The order-59 formula is also described by V. I. Lebedev,
Russian Academy of Sciences Doklady Mathematics 50, 283-286 (1995).

## Tests

`lebedev_exactness.f90` verifies point counts, positive weights, unit radii,
normalization, all even Cartesian monomials through each rule's exactness,
and odd-parity moments. Analytic moments use double-factorial ratios, without
reference to the source tables. Tests cover the minimum-order rounding
(`54 -> 59`, `60..65 -> 65`) and reject an unsupported request of 66.

GNU CPU-only on Terok and NVHPC on Spark pass these tests. Maximum relative
even-moment errors are `1.61e-14`, `1.25e-14`, and `1.59e-14` for orders 53,
59, and 65, respectively (tolerance `5e-13`). Normalized odd-moment errors
are below `2.7e-17`. The existing source, primitive, and periodic-partition
derivative tests also pass with the extended library.

## Periodic Si2 volume probe

For the compact periodized Becke partition, shell 1, 200 radial points,
one orientation, and exact primitive-cell volume `270.011394 bohr^3`:

| Actual order | Integrated volume / bohr^3 | Relative volume error | First cosine integral / bohr^3 |
| ---: | ---: | ---: | ---: |
| 53 | 270.0265316420 | 5.60630e-5 | -0.00681261 |
| 59 | 270.0176379699 | 2.31248e-5 | -0.00281980 |
| 65 (requested 64) | 270.0140090452 | 9.68494e-6 | -0.00118577 |

These model-free integrals were evaluated on Terok. The exact first cosine
integral is zero. The higher angular rules improve both measures in this
case, without rescaling weights or density. Polynomial exactness and volume
accuracy alone do not guarantee converged neural-functional energies or
forces; radial resolution, periodic layouts and native-grid interpolation
remain independent convergence parameters.

The two-shell control at 200 radial points and one orientation also improves:
order 53 gives relative volume error `1.809252458e-4` and cosine integral
`-0.021658174220 bohr^3`; order 65 (requested 64) gives `4.166476066e-5`
and `-0.005028693319 bohr^3`, respectively. The latter passes the explicitly
requested `1e-4` relative tolerance on both integrals. This does not establish
image-shell convergence of the Skala energy or its descriptor window.

The higher-rule Terok volume probes used executable SHA-256
`e14a4fd24055361dd483591d2797163b244b3ebada5a289b562fd9249e94b342`.
The applied Si2 Skala result is recorded in
[`VALIDATION.md`](../../fulltests/skala_crystals/VALIDATION.md#higher-lebedev-rules).
