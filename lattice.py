import numpy as np
from ase import Atoms
from ase.io import read, write

even_numbers = np.arange(0, 6, 2)

valid_combinations = []

for x in even_numbers:
        print(x)
        for y in even_numbers:
                for z in even_numbers:
                        if (x + y + z) % 4 == 0:   
                                valid_combinations.append([x, y, z])
odds = []
for combo in valid_combinations:
    odds.append([combo[i] + 1 for i in range(3)])

pos = valid_combinations + odds
pos = np.array(pos)
print(pos)
atoms = Atoms(f'C{len(pos)}', pos) 
write('test.xyz', atoms)
