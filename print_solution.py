import functions as fn

solution = fn.get_solution()
for result in solution:
	print(
		"z =", result.activity,
		" mu =", result.chemical_potential,
		" free_energy =", result.free_energy,
		" density =", result.density,
		" x =", result.parameter,
		" phase =", result.phase,
	)
