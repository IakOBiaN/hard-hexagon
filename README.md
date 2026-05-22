# Baxter's exact solution for hard hexagon model

The Python implementation of an exact solution for the hard hexagon model, proposed by Baxter in 1980. It is a 2D lattice model of a gas, where particles are allowed to be on the vertices of a triangular lattice but no two particles may be adjacent.

## Description

A schematic representation of the model is shown in the figure. Alternatively, the model can be formulated as a model of hard disks of size <img src="images/disk_formula.png" height="19" /> on a triangular lattice, where <i>d</i> is a hard disk diameter and <i>a</i> is a lattice spacing.

<p align="center"><img src="images/hardhex.png" width="600" /></p>

The hard hexagon model occurs within the framework of the grand canonical ensemble, where the total number of particles (the "hexagons") is allowed to vary naturally, and is fixed by a chemical potential. In the hard hexagon model, all valid states have zero energy, and so the only important thermodynamic control variable is the ratio of chemical potential to temperature <i>mu/(kT)</i>. The exponential of this ratio, <i>z = exp(mu/(kT))</i> is called the activity and larger values correspond roughly to denser configurations.

For a triangular lattice with N sites, the grand partition function is

<img src="images/partition_formula.png" />

where g(n, N) is the number of ways of placing n particles on distinct lattice sites such that no 2 are adjacent. The function kappa is defined by

<img src="images/k_formula.png" />

so that log(kappa) is the free energy per unit site. Solving the hard hexagon model means (roughly) finding an exact expression for kappa as a function of z.

## Solution

The solution is different for the areas before and after the critical point. For <i>z</i> less than critical:

<img src="images/z_formula_less.png" />
<img src="images/k_formula_less.png" />
<img src="images/rho_formula_less.png" />

For <i>z</i> more than critical the solution is given by:

<img src="images/z_formula_more.png" />
<img src="images/k_formula_more.png" />
<img src="images/rho1_formula_more.png" />
<img src="images/rho23_formula_more.png" />

where

<img src="images/G_formula.png" />
<img src="images/H_formula.png" />
<img src="images/P_formula.png" />
<img src="images/Q_formula.png" />

## Python usage

The formulas are parameterized by Baxter's internal parameter `x`, while the physically convenient input is the chemical potential `mu`. The project exposes `solve_by_mu(mu)`, which finds the matching `x` numerically and returns the thermodynamic quantities:

```python
import functions as hardhex

solution = hardhex.solve_by_mu(2.437)
print(solution.activity)
print(solution.free_energy)
print(solution.density)
print(solution.parameter)  # Baxter's internal x parameter
print(solution.phase)
```

The critical activity and chemical potential are available as constants:

```python
hardhex.CRITICAL_ACTIVITY
hardhex.CRITICAL_CHEMICAL_POTENTIAL
```
