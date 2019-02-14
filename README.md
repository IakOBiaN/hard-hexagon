# Baxter's exact solution for hard hexagon model
The Python implementation of an exact solution for the hard hexagon model, proposed by Baxter in 1980. It is a 2D lattice model of a gas, where particles are allowed to be on the vertices of a triangular lattice but no two particles may be adjacent.

## Description

A schematic representation of the model is shown in the figure. Alternatively, the model can be formulated as a model of hard disks of size <img src="images/disk_formula.png" height="19" /> on a triangular lattice, where <i>d</i> is a hard disk diameter and <i>a</i> is a lattice spacing.

<p align="center"><img src="images/hardhex.png" width="600" /></p>

The hard hexagon model occurs within the framework of the grand canonical ensemble, where the total number of particles (the "hexagons") is allowed to vary naturally, and is fixed by a chemical potential. In the hard hexagon model, all valid states have zero energy, and so the only important thermodynamic control variable is the ratio of chemical potential to temperature <i>µ/(kT)</i>. The exponential of this ratio, <i>z = exp(µ/(kT))</i> is called the activity and larger values correspond roughly to denser configurations.

For a triangular lattice with N sites, the grand partition function is 

<img src="images/partition_formula.svg" />

where g(n, N) is the number of ways of placing n particles on distinct lattice sites such that no 2 are adjacent. The function κ is defined by 

<img src="images/k_formula.svg" />

so that log(κ) is the free energy per unit site. Solving the hard hexagon model means (roughly) finding an exact expression for κ as a function of z.

## Solution

The solution is different for the areas before and after the critical point. For <i>z</i> less than critical:

<img src="images/z_formula_less.svg" />
<img src="images/k_formula_less.svg" />
<img src="images/rho_formula_less.svg" />

For <i>z</i> more than critical the solution is given by:

<img src="images/z_formula_more.svg" />
<img src="images/k_formula_more.svg" />
<img src="images/rho1_formula_more.svg" />
<img src="images/rho23_formula_more.svg" />

,where

<img src="images/G_formula.svg" />
<img src="images/H_formula.svg" />
<img src="images/P_formula.svg" />
<img src="images/Q_formula.svg" />
