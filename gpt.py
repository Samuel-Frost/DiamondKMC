import numpy as np
import random

class Defect:
    """Base class for all defects"""
    def __init__(self, position):
        self.position = np.array(position)

    def get_possible_moves(self, lattice):
        """Each defect should define its own movement rules."""
        raise NotImplementedError

class Vacancy(Defect):
    """Vacancy that can merge into chains in the [100] direction."""
    def __init__(self, position):
        super().__init__(position)

    def check_for_merging(self, lattice):
        """Vacancies can only merge in the [100] direction."""
        x, y, z = self.position
        allowed_neighbors = [(x + 1, y, z), (x - 1, y, z)]  # Only ±x direction

        for n in allowed_neighbors:
            if tuple(n) in lattice.defects and isinstance(lattice.defects[tuple(n)], Vacancy):
                neighbor = lattice.defects[tuple(n)]
                
                # Remove both vacancies from the lattice
                del lattice.defects[tuple(self.position)]
                del lattice.defects[tuple(neighbor.position)]
                
                # Create a new vacancy chain
                new_chain = VacancyChain([self, neighbor])
                
                # Add it to the lattice
                lattice.defects[tuple(self.position)] = new_chain
                return  # Only merge one neighbor at a time

    def get_possible_moves(self, lattice):
        """Isolated vacancies move freely; chains stay in ±x direction."""
        moves = [
            self.position + np.array([1, 0, 0]),
            self.position + np.array([-1, 0, 0])
        ]
        return [tuple(m) for m in moves if tuple(m) not in lattice.defects]

class VacancyChain(Defect):
    """Represents a 1D chain of vacancies along [100]."""
    def __init__(self, vacancies):
        self.vacancies = sorted(vacancies, key=lambda v: v.position[0])  # Ordered by x

    def get_possible_moves(self, lattice):
        """Chains move together in the ±x direction."""
        head_x = self.vacancies[-1].position[0]
        tail_x = self.vacancies[0].position[0]

        moves = [
            (head_x + 1, self.vacancies[-1].position[1], self.vacancies[-1].position[2]),
            (tail_x - 1, self.vacancies[0].position[1], self.vacancies[0].position[2])
        ]

        return [m for m in moves if m not in lattice.defects]

    def split_chain(self, lattice):
        """Allow breaking chains into two if energy permits."""
        if len(self.vacancies) > 2 and random.random() < 0.1:  # 10% chance of splitting
            split_index = random.randint(1, len(self.vacancies) - 1)
            left_chain = VacancyChain(self.vacancies[:split_index])
            right_chain = VacancyChain(self.vacancies[split_index:])

            # Remove old chain from lattice
            del lattice.defects[tuple(self.vacancies[0].position)]
            
            # Add the two new chains
            lattice.defects[tuple(left_chain.vacancies[0].position)] = left_chain
            lattice.defects[tuple(right_chain.vacancies[0].position)] = right_chain

class NV_Center(Defect):
    """Nitrogen-Vacancy center that interacts with vacancies."""
    def __init__(self, position):
        super().__init__(position)

    def check_for_interaction(self, lattice):
        """If a vacancy is nearby, form an NV-Vacancy complex."""
        x, y, z = self.position
        possible_vacancy_sites = [
            (x + 1, y, z), (x - 1, y, z),
            (x, y + 1, z), (x, y - 1, z),
            (x, y, z + 1), (x, y, z - 1)
        ]
        
        for n in possible_vacancy_sites:
            if tuple(n) in lattice.defects and isinstance(lattice.defects[tuple(n)], Vacancy):
                vacancy = lattice.defects[tuple(n)]
                
                # Remove both defects
                del lattice.defects[tuple(self.position)]
                del lattice.defects[tuple(vacancy.position)]
                
                # Create a new NV+Vacancy complex
                nv_complex = NVComplex(self, vacancy)
                
                # Add it to the lattice
                lattice.defects[tuple(self.position)] = nv_complex
                return

class NVComplex(Defect):
    """NV center bound to a vacancy"""
    def __init__(self, nv, vacancy):
        self.nv = nv
        self.vacancy = vacancy
        super().__init__(nv.position)

    def get_possible_moves(self, lattice):
        """The NV-complex may be stationary or migrate."""
        return []

class Lattice:
    """Represents the simulation grid and contains all defects."""
    def __init__(self, size):
        self.size = size
        self.defects = {}

    def add_defect(self, defect):
        self.defects[tuple(defect.position)] = defect

    def step(self):
        """Perform a single KMC step: move defects or merge them."""
        defect_list = list(self.defects.values())
        random.shuffle(defect_list)

        for defect in defect_list:
            if isinstance(defect, Vacancy):
                defect.check_for_merging(self)
            
            if isinstance(defect, VacancyChain):
                defect.split_chain(self)

            possible_moves = defect.get_possible_moves(self)
            if possible_moves:
                new_position = random.choice(possible_moves)
                del self.defects[tuple(defect.position)]
                defect.position = np.array(new_position)
                self.defects[tuple(new_position)] = defect

    def print_lattice(self):
        """Prints the defect positions for debugging."""
        print(f"Defects: {list(self.defects.keys())}")

# === RUNNING THE SIMULATION ===
lattice = Lattice(size=10)

# Add initial defects
lattice.add_defect(Vacancy((2, 2, 2)))
lattice.add_defect(Vacancy((3, 2, 2)))
lattice.add_defect(Vacancy((5, 2, 2)))
lattice.add_defect(NV_Center((7, 2, 2)))

# Run the simulation for 10 steps
for _ in range(10):
    lattice.step()
    lattice.print_lattice()

