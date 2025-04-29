import random
import numpy as np
from ase import Atoms
from ase.io import write
import os

def add(vec1, vec2):
    return tuple([vec1[i] + vec2[i] for i in range(3)])

def minus(vec1, vec2):
    return tuple([vec1[i] - vec2[i] for i in range(3)])

moves_even = [[1, 1, 1], [-1, 1, -1], [-1, -1, 1], [1, -1, -1]] 
moves_odd = [[-1, -1, -1], [1, -1, 1], [1, 1, -1], [-1, 1, 1]] 

kB = 8.617333262e-5

class Defect:
    def __init__(self, pos, lattice=None):
        self.pos = tuple(pos) # has to be a tuple as you can't hash lists
        self.lattice = lattice
        if lattice != None:
            self.lattice.add_defect(self)
        self.rate = 0 # by default any defect has a diffusion rate of 0

    def __str__(self):
        return f"{self.__class__.__name__} {str(self.pos)}"

    def __repr__(self):
        return f"{self.__class__.__name__} {str(self.pos)}"

    def move_to(self, vec):
        old_pos = self.pos
        del self.lattice.defects[self.pos] # remove defect from lattice
        self.pos = tuple(vec)
        if self.lattice.add_defect(self) != 0:
            self.pos = old_pos # if we don't move successfully then reinstate the old position
            self.lattice.add_defect(self)
        
    def move_by(self, vec):
        del self.lattice.defects[self.pos] # remove defect from lattice
        self.pos = add(self.pos, vec)
        if self.lattice.add_defect(self) != 0: # if we didn't move successfully then just stay put
            self.pos = minus(self.pos, vec)
            self.lattice.add_defect(self)

    def get_neighbours(self, radius=3):
        x, y, z = self.pos
        neighbours = []
        for dx in range(-radius, radius+1):
            for dy in range(-radius, radius+1):
                for dz in range(-radius, radius+1):
                    if dx == dy == dz == 0:
                        continue
                    neighbour_pos = (x+dx, y+dy, z+dz)
                    if neighbour_pos in self.lattice.defects:
                        neighbours.append(self.lattice.defects[neighbour_pos])
        return neighbours
    
    def wrap(self):
        # this works nicely
        pos = [0, 0, 0]
        for i in range(3):
            pos[i] = self.pos[i] % self.lattice.box[i] 
        self.pos = tuple(pos)

    def available_moves(self):
        global moves_even
        global moves_odd
        if self.pos[0] % 2 == 0:
            return moves_even
        else:
            return moves_odd


class Vacancy(Defect):
    def __init__(self, pos, lattice=None):
        self.pos = tuple(pos) # has to be a tuple as you can't hash lists
        self.lattice = lattice
        if lattice != None:
            self.lattice.add_defect(self)
        self.rate = 40e12 * np.exp(-2.3/(kB*T))
        self.species = 'C'

    def vacancy_merge(self):
        neighbours = self.get_neighbours(radius=3)
        pos = self.pos
        if len(neighbours) > 1:
            print("WE HAVE MANY NEIGHBOURS")
        if len(neighbours) > 0:
            for neighbour in neighbours:
                self.lattice.remove_defect(neighbour.pos)
            self.lattice.remove_defect(pos)
            if type(neighbours[0]) == Vacancy or type(neighbours[0]) == VacancyCluster:
                self.lattice.add_defect(VacancyCluster(pos)) 
            elif type(neighbours[0]) == Nitrogen or type(neighbours[0]) == NitrogenVacancy:
                self.lattice.add_defect(NitrogenVacancy(pos)) 



class VacancyCluster(Defect):
    def __init__(self, pos, lattice=None): 
        self.pos = tuple(pos)
        self.lattice = lattice
        if lattice != None:
            self.lattice.add_defect(self)
        self.rate = 0
        self.species = 'O'

class Divacancy(Vacancy):
    # for the divacancy it will need two positions, and a vector for which 100 plane it's pointing in -> vacancy chain will have the same but many positions
    def __init__(self, pos1, pos2, vec, lattice=None): 
        self.pos = tuple(pos)
        self.lattice = lattice
        if lattice != None:
            self.lattice.add_defect(self)
        self.rate = 0

class Nitrogen(Defect):
    def __init__(self, pos, lattice=None): 
        self.pos = tuple(pos)
        self.lattice = lattice
        if lattice != None:
            self.lattice.add_defect(self)
        self.rate = 0
        self.species = 'N'

class NitrogenVacancy(Defect):
    def __init__(self, pos, lattice=None): 
        self.pos = tuple(pos)
        self.lattice = lattice
        if lattice != None:
            self.lattice.add_defect(self)
        self.rate = 0
        self.species = 'B'
        

class Lattice():
    """
    Contains all the defects, moves the simulation forward
    """
    def __init__(self, box):
        assert sum(box[i] % 4 for i in range(3)) == 0, "Every lattice length must be divisible by 4"
        self.box = box # in order to make sure that any input passed is a valid shape,
        # it is probably a good idea to make box in terms of 8-atom blocks
        self.defects = {}
        self.time = 0
        try:
            os.remove("output.extxyz") 
        except:
            0

    def __str__(self):
        return f"Lattice Defects {self.defects}"

    def add_defect(self, defect):
        # currently this both initialises defects *and* is responsible for adding them after moving, should maybe turn into two separate functions
        if defect != type(Divacancy):
            try:
                print(f"A defect already exists on site {self.defects[defect.pos]}")
            except:
                # this allows us to add defects via lattice.add_defect directly, could be removed?
                if defect.lattice == None:
                    self.defects[defect.pos] = defect
                    self.defects[defect.pos].lattice = self

                # actually may be worth making one giant array with every valid position in it and then comparing against that -> may be faster

                # I know that this looks stupid but it's much faster than if any(a) < 0:
                #if lattice.box[0]-1 < defect.pos[0] < 0 or lattice.box[1]-1 < defect.pos[1] < 0 or lattice.box[2]-1 < defect.pos[2] < 0:
                    # do something special if we care about hitting the boundary.
                
                # wrap first then check if we're cute and valid    
                defect.wrap()
                # checks if the position is *valid* i.e. could exist in an infinite diamond lattice
                if defect.pos[0] % 2 == 0 and defect.pos[1] % 2 == 0 and defect.pos[2] % 2 == 0 and sum(defect.pos) % 4 == 0:
                    self.defects[defect.pos] = defect
                    self.defects[defect.pos].lattice = self
                    return 0 
                elif defect.pos[0] % 2 != 0 and defect.pos[1] % 2 != 0 and defect.pos[2] % 2 != 0 and (sum(defect.pos) + 1) % 4 == 0:
                    self.defects[defect.pos] = defect
                    self.defects[defect.pos].lattice = self
                    return 0 
                else: 
                    del self.defects[defect.pos].lattice
                    print(f"Not a valid position: {defect.pos}")
                    return 1

    def remove_defect(self, pos):
        try:
            del self.defects[tuple(pos)]
        except:
            print(f"No defect exists at site {pos}")

    def write_atoms(self):
        spec = ''
        pos = []
        for defect in list(self.defects.values()):
            spec += defect.species
            pos.append(defect.pos)
        atoms = Atoms(spec, pos, cell=self.box)

        #atoms.cell = self.box  
        write('output.extxyz', atoms, append=True)

    def choose_defect(self):
        # maybe make it so that we only change rates in the array induvidually instead of always
        # recreating this array? it's currently the bottleneck imo
        # will have to maybe rethink this when it comes to dissociating
        rates = []
        for defect in list(self.defects.values()):
            rates.append(defect.rate)
        tot_rate = sum(rates)
        rates = np.array(rates) / tot_rate
        self.time += 1/tot_rate * np.log(1/np.random.uniform())
        return np.random.choice(list(self.defects.values()), 1, p=rates)[0]

    def random_pos(self):
        even_numbers = np.arange(0, self.box[0], 2) 
        valid_combinations = []

        for x in even_numbers:
                for y in even_numbers:
                        for z in even_numbers:
                                if (x + y + z) % 4 == 0:   
                                        valid_combinations.append([int(x), int(y), int(z)])
        odds = []
        for combo in valid_combinations:
            odds.append([combo[i] + 1 for i in range(3)])
        print(tuple(random.choice(valid_combinations)))
        return tuple(random.choice(valid_combinations))
    
    def get_grid(self):
        """Returns a basically unordered array of every valid grid point.
            This gets prohibitively expensive for anything above 512 basically."""
        valid_combinations = []

        for x in np.arange(0, self.box[0], 2) :
                for y in np.arange(0, self.box[1], 2) :
                        for z in np.arange(0, self.box[2], 2) :
                                if (x + y + z) % 4 == 0:   
                                        valid_combinations.append([int(x), int(y), int(z)])
        odds = []
        for combo in valid_combinations:
            odds.append([combo[i] + 1 for i in range(3)])
        print(tuple(random.choice(valid_combinations)))
        return valid_combinations + odds

    def random_gaussian_pos(self, sigma=64):
        # this is biased towards being placed on even cells currently
        arr = [4 * round(np.random.normal(self.box[i]/2, sigma)/4) for i in range(3)]
        return tuple(arr)

    def random_uniform_pos(self):
        # this is biased towards being placed on even cells currently
        return tuple([4 * round(np.random.uniform(0, self.box[i])/4) for i in range(3)])

    def get_num_type(self, typee):
        num = 0
        for defect in list(self.defects.values()):
            if type(defect) == typee:
                num += 1
        return num

        
###### SIMULATION PARAMETERS
T = 1100
num_steps = int(1e5)
N_ppm = 40 
V_ppm = 10


#64 is just over 5nm -> 7.1074 / 8 * 64 * 10e-10
lattice = Lattice([1024, 1024, 1024])
#Vacancy((20, 20, 20), lattice)
#lattice.add_defect(Vacancy([0, 0, 0]))
print("INITIALISING SYSTEM")

V_num = round(V_ppm * np.prod(lattice.box) / 80e6)
print(V_num)
for i in range(V_num):
    lattice.add_defect(Vacancy(lattice.random_gaussian_pos()))

N_num = round(N_ppm * np.prod(lattice.box) / 80e6)
for _ in range(N_num):
    lattice.add_defect(Nitrogen(lattice.random_uniform_pos()))

for i in range(num_steps):
    defect = lattice.choose_defect()
    #defect = list(lattice.defects.values())[0]
    defect.move_by(random.choice(defect.available_moves()))
    defect.vacancy_merge()
    if i % 100 == 0:
        print("Iteration:", f"{i}/{num_steps}", "Time:", lattice.time, "seconds", "V:", lattice.get_num_type(Vacancy), "NV:", lattice.get_num_type(NitrogenVacancy), "Vn:", lattice.get_num_type(VacancyCluster))
        lattice.write_atoms()

print("Total simulation time:", lattice.time, "seconds")