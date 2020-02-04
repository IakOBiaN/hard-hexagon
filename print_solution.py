from math import exp
import functions as fn

solution = fn.get_solution()
for result in solution:
	print("z =",exp(result[0])," mu =",result[0], " free_energy =",result[1], " density =",result[2])
