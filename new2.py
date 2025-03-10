"""
Where do I store the rates information? I suppose the base rate should be stored in the defect class,
also not every defect should have a rate surely as some are immobile so should be kept separate?
"""

from scipy.spatial import KDTree
import numpy as np

class Defect:
    def __init__(self, pos, lattice):
        self.pos = pos
        self.lattice = lattice
    
    def move(self, vec):
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

class Vacancy(Defect):
    def __init__(self, pos):
        super().__init__(self.pos)

    def check_for_merging(self):
        ####
        # 1. search kdtree for things nearby under like 3
        # 2. need to check if other vacancy is nearby then put it in 110 plane
        # 2.5 implement checking for different angles 
        # 3. just call VacancyChain and delete the old ones from the lattice 

class VacancyChain(Defect):
    """
    1D chain of vacancies, only available in [110] for now, but looks like they can form in other directions.
    It's not storing the Vacancy class because they can't migrate unless they break off from the chain.
    """
    def __init__(self, vacancies, lattice):
        # sorted by their global position in 110 direction for now
        # vacancies here are Defects, but not Vacancies
        self.vacancies = sorted(vacancies, key=lambda v: np.sqrt(v.pos[0]**2 + v.pos[1]**2))
        # this may give weird memory erros, do get_pos() instead?
        self.positions = [v.pos for v in vacancies]
        # initialises the parent Defect class with the lower vacancy in the chain
        super().__init__(self.positions[0])
        self.lattice = lattice

        # basically stores loads of pointers to the vacancy chain for each induvidual position in the 
        # vacancy chain, making it possible for the kdtree to see each of them
        for pos in self.positions:
            self.lattice.defects[tuple(pos)] = self

    def extend_chain(self, vacancy, lattice):
        # takes in the vacancy which is being added (so it can be deleted)
        # and the lattice because we need to update what's going on there

        # can't have chains of longer than 4 to match up with Kai's simulation
        if len(self.vacancies) <= 4:
            # check if the vacancy is close to one of the side things
            if np.array_equal(self.pos[-1] + [2, 2, 0], vacancy.pos) or np.array_equal(self.pos[0] - [2, 2, 0], vacancy.pos):

                # Remove the vacancy from the lattice
                del lattice.defects[tuple(vacancy.pos)]

                # Add it to the chain
                # IMPLEMENT TACKING THE VACANCY ONTO THE END OF THE CHAIN PROPERLY
                self.pos.append(vacancy.pos)
                self.pos.sort(key=lambda p: p[0]) 


class Lattice():
    """
    Contains all the defects, moves the simulation forward
    """
    def __init__(self, box):
        self.box = box # should probably assert that this is a 1x3 array of ints or something
        self.defects = {}
        self.time = 0
        self.kdtree = None
    
    def add_defect(self, defect):
        # store in a dictionary for faster indexing and removing (I think)
        self.defects[tuple(defect.pos)] = defect

    def remove_defect(self, defect):
        del self.defects[tuple(defect.pos)]

    def update_kdtree(self):
        # as the vacancy chain only has the initial vacancy stored as the position,
        # only that will be interacted with here, MUST FIX
        positions = list(self.defects.keys())
        return cKDTree(positions, boxsize = box)

        
lattice = Lattice([100, 100, 100])
print(lattice)
# when it's initialised you can place it on the lattice 
vac = VacancyChain([Defect([10, 10, 10]), Defect([11, 11, 11])], lattice)
print(vac.positions)

print(lattice.defects[(10, 10, 10)])
print(lattice.defects[(11, 11, 11)])
print(lattice.defects[(11, 11, 11)] == lattice.defects[(10, 10, 11)])

lattice.add_defect(vac)
#print(lattice.defects[(10, 10, 10)].pos)
