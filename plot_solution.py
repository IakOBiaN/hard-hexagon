import matplotlib.pyplot as plt
import hardhexagon as fn

chemical_potential = []
density = []
for result in fn.get_solution():
	chemical_potential.append(result.chemical_potential)
	density.append(result.density)
	print(
		"z =", result.activity,
		" mu =", result.chemical_potential,
		" free_energy =", result.free_energy,
		" density =", result.density,
		" x =", result.parameter,
		" phase =", result.phase,
	)
plt.plot(chemical_potential, density)
plt.xlabel('Chemical potential')
plt.ylabel('Density')
plt.show()
