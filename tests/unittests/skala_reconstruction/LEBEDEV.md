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

Polynomial exactness alone does not guarantee converged functional energies
or forces. The separate `partition_measure.x` probe tests the periodic
constant-field volume and first reciprocal cosine integral; its invocation
is documented in [README.md](README.md#model-free-integration-probe).
