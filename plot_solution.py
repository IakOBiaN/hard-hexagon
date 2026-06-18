import matplotlib.pyplot as plt

import hardhexagon as hh

solutions = hh.get_solution()
chemical_potential = [solution.chemical_potential for solution in solutions]
density = [solution.density for solution in solutions]

critical_mu = hh.CRITICAL_CHEMICAL_POTENTIAL
critical_density = hh.solve_by_mu(critical_mu).density

fig, ax = plt.subplots()
ax.plot(chemical_potential, density, color="#3b82c4", label="density")
ax.axvline(critical_mu, color="#e3a008", linestyle="--", linewidth=1)
ax.plot(critical_mu, critical_density, "o", color="#e3a008")
ax.annotate(
    "critical point",
    xy=(critical_mu, critical_density),
    xytext=(critical_mu + 1.0, critical_density - 0.06),
    color="#e3a008",
    arrowprops=dict(arrowstyle="->", color="#e3a008"),
)
ax.set_xlabel(r"chemical potential  $\mu / kT$")
ax.set_ylabel(r"density  $\rho$")
ax.set_title("Hard hexagon density vs. chemical potential")
ax.legend()
plt.show()
