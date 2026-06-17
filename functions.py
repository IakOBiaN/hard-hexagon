from dataclasses import dataclass
from math import exp, isfinite, log, sqrt

DEFAULT_ITERATIONS = 1000
DEFAULT_CHEMICAL_POTENTIALS = tuple(value / 10.0 for value in range(-80, 81))
CRITICAL_ACTIVITY = (11.0 + 5.0 * sqrt(5.0)) / 2.0
CRITICAL_CHEMICAL_POTENTIAL = log(CRITICAL_ACTIVITY)

def G(x, iterations=DEFAULT_ITERATIONS):
    prod = 1.
    for i in range(1, iterations):
        prod = prod*1./((1.-x**(5.*i-4.))*(1.-x**(5.*i-1.)))
    return prod

def H(x, iterations=DEFAULT_ITERATIONS):
    prod = 1.
    for i in range(1, iterations):
        prod = prod*1./((1.-x**(5.*i-3.))*(1.-x**(5.*i-2.)))
    return prod

def P(x, iterations=DEFAULT_ITERATIONS):
    prod = 1.
    for i in range(1, iterations):
        prod = prod*(1.-x**(2.*i-1.))
    return prod

def Q(x, iterations=DEFAULT_ITERATIONS):
    prod = 1.
    for i in range(1, iterations, 1):
        prod = prod*(1.-x**i)
    return prod

@dataclass(frozen=True)
class Solution:
    chemical_potential: float
    activity: float
    free_energy: float
    density: float
    parameter: float
    phase: str

    def __iter__(self):
        yield self.chemical_potential
        yield self.free_energy
        yield self.density

def calculate(x, iterations=DEFAULT_ITERATIONS):
    if x < 0:
        z = (-x*H(x, iterations)**5)/G(x, iterations)**5
        rho = -x*G(x, iterations)*H(x**6, iterations)*P(x**3, iterations)/P(x, iterations)
        prod = 1.0
        for i in range(1, iterations, 1):
            prod = prod*((1.0-x**(6.0*i-4.0))*(1.0-x**(6.0*i-3.0))**2*(1.0-x**(6.0*i-2.0)))/((1.0-x**(6.0*i-5.0))*(1.0-x**(6.0*i-1.0))*(1.0-x**(6.0*i))**2)
        btp = H(x, iterations)**3*Q(x**5, iterations)**2/G(x, iterations)**2*prod
    else:
        z = G(x, iterations)**5/((x)*H(x, iterations)**5)
        rho23 = x**2*H(x, iterations)*Q(x, iterations)*H(x**9, iterations)*Q(x**9, iterations)/Q(x**3, iterations)**2
        rho1 = H(x, iterations)*Q(x, iterations)*(G(x, iterations)*Q(x, iterations)+x**2*H(x**9, iterations)*Q(x**9, iterations))/Q(x**3, iterations)**2
        rho = (rho1 + rho23*2.0)/3.0
        prod = 1.0
        for i in range(1, iterations, 1):
            prod = prod*((1.0-x**(3.0*i-2.0))*(1.0-x**(3.0*i-1.0)))/(1.0-x**(3.0*i))**2
        btp = x**(-1.0/3.0)*G(x, iterations)**3*Q(x**5, iterations)**2/H(x, iterations)**2*prod
    return log(z), -log(btp),rho

def solve_by_parameter(x):
    chemical_potential, free_energy, density = calculate(x)
    phase = "low_activity" if x < 0 else "high_activity"
    return Solution(
        chemical_potential=chemical_potential,
        activity=exp(chemical_potential),
        free_energy=free_energy,
        density=density,
        parameter=x,
        phase=phase,
    )

def _chemical_potential_error(x, target_mu):
    return calculate(x)[0] - target_mu

def _find_parameter_by_bisection(target_mu, left, right, tolerance=1e-12, max_steps=80):
    left_error = _chemical_potential_error(left, target_mu)
    right_error = _chemical_potential_error(right, target_mu)
    if left_error == 0.0:
        return left
    if right_error == 0.0:
        return right
    if left_error * right_error > 0.0:
        raise ValueError("Could not bracket the Baxter parameter for this chemical potential.")

    for _ in range(max_steps):
        middle = (left + right) / 2.0
        middle_error = _chemical_potential_error(middle, target_mu)
        if abs(middle_error) <= tolerance or abs(right - left) <= 1e-15:
            return middle
        if left_error * middle_error <= 0.0:
            right = middle
            right_error = middle_error
        else:
            left = middle
            left_error = middle_error
    return (left + right) / 2.0

def _low_activity_parameter_for_mu(target_mu):
    left = -0.99
    right = -1e-12
    while _chemical_potential_error(right, target_mu) > 0.0:
        right *= 1e-2
        if right == 0.0:
            raise ValueError("Chemical potential is too small for double-precision inversion.")
    return _find_parameter_by_bisection(target_mu, left, right)

def _high_activity_parameter_for_mu(target_mu):
    left = 1e-12
    right = 0.99
    while _chemical_potential_error(left, target_mu) < 0.0:
        left *= 1e-2
        if left == 0.0:
            raise ValueError("Chemical potential is too large for double-precision inversion.")
    return _find_parameter_by_bisection(target_mu, left, right)

def parameter_from_chemical_potential(chemical_potential):
    if not isfinite(chemical_potential):
        raise ValueError("Chemical potential must be a finite number.")
    if abs(chemical_potential - CRITICAL_CHEMICAL_POTENTIAL) <= 1e-12:
        # The product formulas are singular at the endpoint; this is the critical limit.
        return -0.99
    if chemical_potential < CRITICAL_CHEMICAL_POTENTIAL:
        return _low_activity_parameter_for_mu(chemical_potential)
    return _high_activity_parameter_for_mu(chemical_potential)

def solve_by_mu(chemical_potential):
    parameter = parameter_from_chemical_potential(chemical_potential)
    solution = solve_by_parameter(parameter)
    phase = "critical" if abs(chemical_potential - CRITICAL_CHEMICAL_POTENTIAL) <= 1e-12 else solution.phase
    return Solution(
        chemical_potential=chemical_potential,
        activity=exp(chemical_potential),
        free_energy=solution.free_energy,
        density=solution.density,
        parameter=solution.parameter,
        phase=phase,
    )

def get_solution(chemical_potentials=DEFAULT_CHEMICAL_POTENTIALS):
    return [solve_by_mu(mu) for mu in chemical_potentials]
