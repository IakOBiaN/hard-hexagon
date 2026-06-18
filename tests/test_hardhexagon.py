"""Correctness tests for Baxter's exact hard hexagon solution.

The hard hexagon model has a number of analytically known quantities that the
implementation in :mod:`hardhexagon` must reproduce.  Each test below pins one of
them, so that any change to the numerics is checked against the physics rather
than against a previously recorded output.

Known exact results used here:

* critical activity        z_c = (11 + 5*sqrt(5)) / 2          ~= 11.0902
* critical density         rho_c = (5 - sqrt(5)) / 10          ~= 0.276393
* close-packing density    rho_max = 1/3   (one full sublattice)
* low-density expansion    rho = z - 7 z^2 + O(z^3)
"""
from fractions import Fraction
from math import exp, log, sqrt

import pytest

import hardhexagon as fn

CRITICAL_DENSITY = (5.0 - sqrt(5.0)) / 10.0
CLOSE_PACKING_DENSITY = 1.0 / 3.0
SECOND_VIRIAL_COEFFICIENT = -7  # coefficient of z^2 in rho(z)


@pytest.fixture(scope="module")
def default_solutions():
    """Solve the default chemical-potential grid once and share it.

    ``get_solution`` inverts mu -> x by bisection at every grid point, so it is
    the most expensive call in the suite; computing it once keeps the run fast.
    """
    return fn.get_solution()


# --------------------------------------------------------------------------- #
# Critical point
# --------------------------------------------------------------------------- #
def test_critical_activity_closed_form():
    """z_c equals the golden-ratio expression (11 + 5*sqrt(5)) / 2."""
    assert fn.CRITICAL_ACTIVITY == pytest.approx((11.0 + 5.0 * sqrt(5.0)) / 2.0)


def test_critical_chemical_potential_is_log_activity():
    assert fn.CRITICAL_CHEMICAL_POTENTIAL == pytest.approx(log(fn.CRITICAL_ACTIVITY))


@pytest.mark.parametrize("offset", [-1e-7, +1e-7])
def test_critical_density_approached_from_both_phases(offset):
    """Both phases meet at rho_c = (5 - sqrt(5)) / 10 (continuous transition)."""
    solution = fn.solve_by_mu(fn.CRITICAL_CHEMICAL_POTENTIAL + offset)
    assert solution.density == pytest.approx(CRITICAL_DENSITY, abs=1e-4)


def test_density_continuous_across_transition():
    eps = 1e-7
    low = fn.solve_by_mu(fn.CRITICAL_CHEMICAL_POTENTIAL - eps)
    high = fn.solve_by_mu(fn.CRITICAL_CHEMICAL_POTENTIAL + eps)
    assert low.density == pytest.approx(high.density, abs=1e-4)


def test_free_energy_continuous_across_transition():
    eps = 1e-7
    low = fn.solve_by_mu(fn.CRITICAL_CHEMICAL_POTENTIAL - eps)
    high = fn.solve_by_mu(fn.CRITICAL_CHEMICAL_POTENTIAL + eps)
    assert low.free_energy == pytest.approx(high.free_energy, abs=1e-4)


def test_exact_critical_call_labels_phase_and_density():
    solution = fn.solve_by_mu(fn.CRITICAL_CHEMICAL_POTENTIAL)
    assert solution.phase == "critical"
    assert solution.density == pytest.approx(CRITICAL_DENSITY, abs=1e-5)


# --------------------------------------------------------------------------- #
# Limiting regimes
# --------------------------------------------------------------------------- #
def test_low_density_leading_order():
    """As z -> 0 the density tends to the activity: rho ~ z, with rho < z."""
    solution = fn.solve_by_mu(-12.0)
    z, rho = solution.activity, solution.density
    assert rho == pytest.approx(z, rel=1e-3)
    assert rho < z  # the first correction (-7 z^2) is negative


def test_second_virial_coefficient_from_exact_solution():
    """rho = z - 7 z^2 + ... ; recover the -7 from the exact solution."""
    solution = fn.solve_by_mu(-13.0)
    z, rho = solution.activity, solution.density
    b2 = (rho - z) / z ** 2
    assert b2 == pytest.approx(SECOND_VIRIAL_COEFFICIENT, abs=1e-2)


def test_close_packing_limit():
    """As mu -> +infinity the density saturates one sublattice: rho -> 1/3."""
    solution = fn.solve_by_mu(14.0)
    assert solution.density == pytest.approx(CLOSE_PACKING_DENSITY, abs=1e-5)
    assert solution.density < CLOSE_PACKING_DENSITY


@pytest.mark.parametrize("mu", [-6.0, -2.0, 0.0, 2.0, 2.5, 4.0, 8.0])
def test_density_within_physical_bounds(mu):
    solution = fn.solve_by_mu(mu)
    assert 0.0 < solution.density < CLOSE_PACKING_DENSITY


def test_density_strictly_increasing_in_chemical_potential(default_solutions):
    densities = [s.density for s in default_solutions]
    assert all(a < b for a, b in zip(densities, densities[1:]))


# --------------------------------------------------------------------------- #
# Inversion mu <-> x and the public API
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("mu", [-5.0, -1.0, 0.0, 1.0, 2.0, 2.4, 2.5, 3.0, 5.0])
def test_inversion_round_trip(mu):
    """solve_by_mu must return the chemical potential it was asked for."""
    solution = fn.solve_by_mu(mu)
    assert solution.chemical_potential == pytest.approx(mu, abs=1e-9)
    assert solution.activity == pytest.approx(exp(mu), rel=1e-9)


@pytest.mark.parametrize(
    "mu, expected_phase",
    [
        (-3.0, "low_activity"),
        (2.0, "low_activity"),
        (3.0, "high_activity"),
        (6.0, "high_activity"),
    ],
)
def test_phase_labels(mu, expected_phase):
    assert fn.solve_by_mu(mu).phase == expected_phase


def test_solution_iterates_as_mu_free_energy_density():
    solution = fn.solve_by_mu(1.0)
    assert tuple(solution) == (
        solution.chemical_potential,
        solution.free_energy,
        solution.density,
    )


def test_get_solution_covers_default_grid(default_solutions):
    assert len(default_solutions) == len(fn.DEFAULT_CHEMICAL_POTENTIALS)
    assert [s.chemical_potential for s in default_solutions] == list(
        fn.DEFAULT_CHEMICAL_POTENTIALS
    )


# --------------------------------------------------------------------------- #
# Independent cross-check: brute-force enumeration on a small torus
# --------------------------------------------------------------------------- #
def _triangular_torus_edges(length):
    """Undirected edges of a ``length`` x ``length`` triangular lattice on a torus.

    On the triangular lattice every site (i, j) has six neighbours: the four
    square-lattice neighbours plus the two along the (1, -1) diagonal.  Listing
    only the three "forward" directions yields each undirected edge exactly once.
    """

    def index(i, j):
        return (i % length) * length + (j % length)

    edges = set()
    for i in range(length):
        for j in range(length):
            a = index(i, j)
            for di, dj in ((1, 0), (0, 1), (1, -1)):
                edges.add(frozenset((a, index(i + di, j + dj))))
    return [tuple(edge) for edge in edges]


def _independent_set_counts(length):
    """g[n] = number of ways to place n mutually non-adjacent particles."""
    sites = length * length
    edges = _triangular_torus_edges(length)
    counts = [0] * (sites + 1)
    for state in range(1 << sites):
        if all(not ((state >> a) & 1 and (state >> b) & 1) for a, b in edges):
            counts[bin(state).count("1")] += 1
    return counts


def test_torus_enumeration_matches_lattice_combinatorics():
    """A 4x4 torus has 3N edges, N single-particle and C(N,2)-3N two-particle states."""
    length = 4
    sites = length * length
    edges = _triangular_torus_edges(length)
    counts = _independent_set_counts(length)

    assert len(edges) == 3 * sites
    assert counts[0] == 1
    assert counts[1] == sites
    assert counts[2] == sites * (sites - 1) // 2 - 3 * sites


def test_second_virial_coefficient_from_enumeration():
    """The brute-force torus reproduces the -7 coefficient exactly (rational arithmetic).

    Expanding rho = (1/N) * z d/dz log(sum_n g[n] z^n) gives
    rho = z + (2 g2 - g1^2)/N * z^2 + ...  with the linear term normalised to z.
    """
    length = 4
    sites = length * length
    counts = _independent_set_counts(length)
    g1, g2 = counts[1], counts[2]

    linear = Fraction(g1, sites)
    quadratic = Fraction(2 * g2 - g1 * g1, sites)

    assert linear == 1
    assert quadratic == SECOND_VIRIAL_COEFFICIENT
