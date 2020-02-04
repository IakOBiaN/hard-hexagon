from math import exp
import matplotlib.pyplot as plt
import functions as fn

chemical_potential = []
density = []
for pair in reversed(fn.compressibility_to_chemical_potential):
	result = fn.calculate(pair[1])
	chemical_potential.append(result[0])
	density.append(result[2])
	print("z =",exp(result[0])," mu =",result[0], " free_energy =",result[1], " density =",result[2])
plt.plot(chemical_potential, density)
plt.xlabel('Chemical potential')
plt.ylabel('Density')
plt.show()
