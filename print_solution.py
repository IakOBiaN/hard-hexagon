from math import exp
import functions

solution = functions.get_solution()
for result in solution:
	print("z =",exp(result[0])," mu =",result[0], " free_energy =",result[1], " density =",result[2])
