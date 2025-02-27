### TODO: REPLACE KDTREE WITH RTREE OR SOMETHING ELSE QUICKER 

import numpy as np
import random
from scipy.spatial import cKDTree
# only constants needed to define our cell
a = 3.567
# one million atoms is 50x50x50
s = 50
limit = 3*s + (s-1) + 0.01
cell = [limit, limit, limit]
# really don't know where to put this
moves = [[1, 1, 1], [-1, 1, -1], [-1, -1, 1], [1, -1, -1]]
rate_dict = {'vacancy': 2.3, 'interstitial': 1.8} 

nu = 40e12
kB = 8.6173303e-5
T = 1100
kBT = kB * T 

N = 100000
max_time = 3600
num_vacancies = 5
num_interstitials = 5

class defect:
    def __init__(self, type, coords=None, mobile=True):
        assert type in ['vacancy', 'interstitial', 'divacancy',
                        'di-interstitial', 'recombined'], "Defect must be of valid type"
        self.type = type
        if self.type != 'vacancy' or self.type != 'interstitial':
            mobile = False
        self.mobile = mobile
        if coords == None:
            self.coords = defect.gen_coords(s)
            self.wrap()
        if mobile:
            self.rate = rate_dict[self.type]            

    def __str__(self):
        return f"type: {self.type}, coords: {self.coords}, mobile: {self.mobile}"

    def gen_coords(s):
        """
        Generates random coordinates for a defect according to size of cell, s
        """
        coords = np.array([(random.randrange(0, 4 * s, 2)) for _ in range(3)])
        if sum(coords) % 4 != 0:
            # currently this can put defects 'outside' the box, so maybe rethink how this is done 
            coords[random.randint(0, 2)] += random.choice([1, -1]) * 2
            coords = abs(coords) 
        coords += random.choice([np.array([0, 0, 0]), np.array([1, 1, 1])])
        return coords 

    def wrap(self):
        # maybe worth having a site_type variable to avoid doing % 2 every time
        for i in range(len(self.coords)):
            if self.coords[i] > limit:
                if self.coords[i] % 2 == 0:
                   self.coords[i] = 0
                else:
                    self.coords[i] = 1 

            if self.coords[i] < 0:
                if self.coords[i] % 2 == 0:
                    self.coords[i] = limit - 1
                else:
                    self.coords[i] = limit

    def hop(self):
        """
        Hops one lattice site, returns the new coordinates
        TODO: find a proper place to put the moves array,
        different movement for interstitials (need to probe first)
        """
        move = random.choice(moves)
        # if we're on an 'odd' site then minus the moves array
        if self.coords[0] % 2 != 0:
            self.coords -= move
            self.wrap()
        else:
            self.coords += move
            self.wrap()
        return self.coords.copy()

    def get_coords(self):
        # since coords are a pointer needs own function to get them
        return self.coords.copy()

    def neighbour_index(self, positions):
        kdtree = cKDTree(positions, boxsize=cell)
        query = kdtree.query([self.coords], k=2)
        if query[0][0][1] < 2:
            return query[1][0][1]
        else:
            return None

def combine(array, i, j):
    n = max(i, j)
    m = min(i, j)
    del(array[n])
    del(array[m])

defects = []
for _ in range(num_vacancies):
    defects.append(defect('vacancy'))
for _ in range(num_interstitials):
    defects.append(defect('interstitial'))

# array of pointers to coordinates

coords = [defect.coords for defect in defects]
i = 0
while len(defects) > 0 and i < 1000000:
    i += 1
    index = random.randint(0, len(defects)-1)
    defect = defects[index]
    defect.hop()
    neigh = defect.neighbour_index(coords)  
    if neigh != None:
        combine(defects, index, neigh) 
        coords = [defect.coords for defect in defects]
        print(i, neigh, len(coords))
