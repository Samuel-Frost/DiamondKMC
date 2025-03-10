##### SOME TRUTHS
# Vacancy clusters cannot move, only the single vacancy
# Intertitials of any species can move
# NVx doesn't move, but V can split off NV2 and move elsewhere
# Ns doesn't move and only has one size
from itertools import permutations
import numpy as np
from ase import Atoms
from ase.io import read, write
from scipy.spatial import cKDTree
import random
from ase.spacegroup import crystal 
s = 400 # not too sure what this number actually is 
a = 3.567 
nu = 40e12
kB = 8.6173303e-5
T = 800
kBT = kB * T

# really the rate should be inherently linked to the species?
# also need a strain field component
# EITHER IT DISSOCIATES OR IT MIGRATES OR IT DOES NOTHING (Ns)
# these are the base rates, will need to be modified for strain fields
diffusion_dict = {'V': nu * np.exp(-2.3/kBT), 'I': nu * np.exp(-1.8/kBT)}
moves_even = [[1, 1, 1], [-1, 1, -1], [-1, -1, 1], [1, -1, -1]] # if on an odd site then we can only move the neg version of these
moves_odd = [[-1, -1, -1], [1, -1, 1], [1, 1, -1], [-1, 1, 1]] 

class defect:
    def __init__(self, species, pos, pairs=[]):
        assert species in ['V', 'I', 'Ns', 'NVx'], "Defect must be of valid species" 
        assert type(pos) == list or type(pos) == type(np.array(0)), "Defect must be a list or numpy array"
        self.species = species
        self.pos = pos
        self.migration = 0
        self.dissociation = 0
        self.pairs = pairs
        if len(pairs) == 0:
            self.size = 1
        else:
            self.size = len(pairs)
        # self.maxsize = 4 # to stop things growing forever
        if self.size == 1:
            try:
                self.migration = diffusion_dict[self.species]
            except:
                self.migration = 0 

        # can also just go in a dictionary?
        if self.size > 2 and self.species == 'V':    
            # random small number
            self.migration = 0
            self.dissociation = 0.0001
        if self.size > 2 and self.species == 'NVx':
            self.dissociation = 3e-6
        
    def update(self):
        self.size = len(self.pairs)
        print("size", self.size)
        if self.size == 0:
            print("WHAT THE FUCK IS GOING ON")
            quit()
        if self.size == 1:
            try:
                self.migration = diffusion_dict[self.species]
                self.dissociation = 0
            except:
                self.migration = 0
                self.dissociation = 0

        # can also just go in a dictionary?
        if  self.size >= 2 and self.species == 'V':    
            print("does this ever happen?")
            self.migration = 0
            self.dissociation = 0.0001
        if self.size > 2 and self.species == 'NVx':
            self.migration = 0
            self.dissociation = 3e-6
        

    def __str__(self):
        return f"species: {self.species}, pos: {self.pos}, migration: {self.migration}, dissociation: {self.dissociation}"
    
    def move(self, vec):
        """
        Just moves the atom in the prescribed direction, 
        I think choosing where to move it should be its own function
        """
        self.pos = self.wrap(self.pos + vec)
        return self.pos.copy()

    def wrap(self, pos):
        limit = s / 4 * a
        for i in range(len(pos)):
            if pos[i] > limit:
                if pos[i] % 2 == 0:
                    pos[i] = 0
                else:
                    pos[i] = 1

            if pos[i] < 0:
                if pos[i] % 2 == 0:
                    pos[i] = limit - 1
                else:
                    pos[i] = limit
        return pos 
 

    def get_pos(self):
        return self.pos.copy()

    def random_walk(self):
        if self.pos[0] % 2 != 0:
            self.move(random.choice(moves_odd))
        else:
            self.move(random.choice(moves_even))
        return self.pos.copy()
    
    def get_nearest_neighbour(self, kdtree, r=100): 
        # for now k = 2 by default because we're not strain fields, so just find nearest neighbour
        query = kdtree.query([self.get_pos()], k=2, distance_upper_bound=r)
        return [query[0][0][1], query[1][0][1]]

def gauss():
        r = np.random.normal(s/2, 6, 1)[0]
        return 2 * round(r / 2)

def gen_pos(func): 
    pos = np.array([func() for _ in range(3)])
    if sum(pos) % 4 != 0:
    # might have a slight bias for the pos direction
        pos[random.randint(0, 2)] += 2
    pos += random.choice([np.array([0, 0, 0]), np.array([1, 1, 1])])
    return pos

def round_to_nearest_site(pos):
    for i in range(len(pos)):
       pos[i] = round(pos[i]) 
    if sum(pos) % 4 != 0:
    # might have a slight bias for the pos direction
        pos[random.randint(0, 2)] += 2
    pos += random.choice([np.array([0, 0, 0]), np.array([1, 1, 1])])
    return pos  

class defect_list:
    """
    List which contains every defect, make sure to do some species checking that everything is a defect
    """
    def __init__(self, defects):
        self.defects = defects 
        for i in range(len(self.defects)):
            self.defects[i].pairs = [i]
            

    def __getitem__(self, indices):
        return self.defects[indices] 

    def __len__(self):
        return len(self.defects)

    def get_positions(self):
        pos = []
        for defect in self.defects:
            pos.append(defect.get_pos())
        return pos

    def get_migrations(self):
        arr = []
        for defect in self.defects:
            arr.append(defect.migration)
        return arr
                

    def get_dissociations(self):
        arr = []
        for defect in self.defects:
            arr.append(defect.dissociation)
        return arr

    def append(self, appendix):
        self.defects.append(appendix)
    
    def atoms(self):
        # need to be able to export as certain kinds
        #pos = self.get_positions()
        #atoms = Atoms(f'C{len(self)}', pos) 
        atoms = Atoms()
        atoms.cell[0][0] = s/4 * a 
        atoms.cell[1][1] = s/4 * a 
        atoms.cell[2][2] = s/4 * a 
        for defect in self.defects:
            #print(defect.size)
            if defect.species == 'V':
                if defect.size == 1:
                    # could get this as a pointer instead and then update things quicker?
                    d = Atoms(f'C', [defect.get_pos()])
                if defect.size >= 2:
                    d = Atoms(f'N', [defect.get_pos()])

            atoms += d

        return atoms 
    
    def create_kdtree(self):
        # no boxsize if we don't want PBC
        # has to be rebuilt every time we move something, otherwise positions aren't updated
        return cKDTree(self.get_positions())#, boxsize = [s/4 * a]*3)
       
    def check_neighbours(self, index):
        pos = self.get_positions()
        kdtree = cKDTree(pos, boxsize = [s/4 * a]*3)
        query = kdtree.query([pos[index]], k=2)
        # nearest neighbour
        d = query[0][0][1]
        i = query[1][0][1]
        # if it's closer than 3 then it has basically been captured (if applicable)
        if d < 3:
            # index is the mobile one, so it should do all the moving
            self.combine(index, i)
        return [d, i]

    def combine(self, i, j):
        # i is moving into j, so j should keep its position
        print("AAAAAAAAAah, I'm combining")
        pair = set((self.defects[j].species, self.defects[i].species))
        pos = self.defects[j].get_pos() 
        if pos[0] % 2 != 0:
            possible_pos = np.array(pos - moves_even) 
        else:
            possible_pos = np.array(pos - moves_odd) 
        # finds all the possible neighbours of j, then takes the one that is closest to i to move i to
        self.defects[i].pos = possible_pos[np.argmin(np.sum(np.abs(possible_pos - self.defects[i].get_pos()), axis=1))]
        

        if pair == set(('V', 'V')):
            print("COMBINING")
            print("i", self.defects[i].pairs)
            print("j", self.defects[j].pairs)
            self.defects[i].pairs = self.defects[j].pairs # they share the same place in memory
            self.defects[j].pairs.append(i)
            self.defects[i].update()
            self.defects[j].update()
            print("i updated", self.defects[i].pairs)
            print("j updated", self.defects[j].pairs)


    def dissociate(self, i):
        # has to be a defects func and not a defect func as it needs to know the list that it is in
        print("AAAAAaah I'm dissociating!")
        pairs = self.defects[i].pairs.copy()
        print(i, "is dissociating")
        print("all the pairs are ", pairs)
        print("j pairs", self.defects[j].pairs)
        self.defects[j].pairs.remove(i)
        print("j pairs", self.defects[j].pairs)
        self.defects[i].pairs = [i]
        print("j pairs", self.defects[j].pairs)
        # make this random or something
        self.defects[i].pos = self.defects[i].pos + random.choice([[-3, -3, -3], [3, -3, 3], [3, 3, -3], [-3, 3, 3]]) 
        for i in pairs:
            self.defects[i].update()
        # needs to move away from its nearest neighbour in plane

# initialise system
# move defect by one 
# build kdtree and check if it's near anything to recombine (one function)
#  

defects = defect_list([defect('V', gen_pos(gauss))])
for i in range(20):
    defects.append(defect('V', gen_pos(gauss), pairs=[i+1]))
arr = [] 

#while len(defects) > 5:
i = 0
while i < 2000:
    if i % 1 == 0:
        #print(i)
        arr.append(defects.atoms())
    if i % 1000 == 0:
        #print(i)
        0
    total = sum(defects.get_migrations()) + sum(defects.get_dissociations())
    print(sum(defects.get_migrations()), sum(defects.get_dissociations()))
    if random.random() * total < sum(defects.get_migrations()):
        j = np.random.choice(range(len(defects)), p=defects.get_migrations()/sum(defects.get_migrations()))
        defects[j].random_walk()
        defects.check_neighbours(j)
    else:
        j = np.random.choice(range(len(defects)), p=np.array(defects.get_dissociations())/sum(defects.get_dissociations()))
        defects.dissociate(j)

    i += 1
    
arr.append(defects.atoms())
write('test.extxyz', arr)
