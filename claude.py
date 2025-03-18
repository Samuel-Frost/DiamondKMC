import numpy as np
from collections import defaultdict
from scipy.spatial import cKDTree

class DefectKMC:
    def __init__(self, box_dimensions, cutoff_radius, temperature):
        """
        Initialize the KMC simulation.
        
        Parameters:
        -----------
        box_dimensions : tuple or list
            Dimensions of the simulation box (Lx, Ly, Lz)
        cutoff_radius : float
            Interaction radius for defect-defect interactions
        temperature : float
            Simulation temperature in Kelvin
        """
        self.box = np.array(box_dimensions)
        self.cutoff = cutoff_radius
        self.temperature = temperature
        self.kb = 8.617333262e-5  # Boltzmann constant in eV/K
        
        # Store defects as a dictionary: {defect_id: {'position': array, 'type': str, 'properties': dict}}
        self.defects = {}
        self.defect_counter = 0
        
        # Neighbor list for efficient interaction calculations
        self.neighbor_list = None
        
    def add_defect(self, position, defect_type, properties=None):
        """Add a defect to the simulation."""
        if properties is None:
            properties = {}
            
        # Apply periodic boundary conditions
        position = np.array(position) % self.box
        
        # Add the defect to our dictionary
        self.defects[self.defect_counter] = {
            'position': position,
            'type': defect_type,
            'properties': properties
        }
        
        # Increment counter and return the defect ID
        defect_id = self.defect_counter
        self.defect_counter += 1
        
        # Neighbor list is now outdated
        self.neighbor_list = None
        
        return defect_id
    
    def update_neighbor_list(self):
        """Update the neighbor list for efficient interaction calculations."""
        if not self.defects:
            self.neighbor_list = None
            return
            
        # Get positions and IDs
        positions = np.array([self.defects[defect_id]['position'] for defect_id in self.defects])
        ids = list(self.defects.keys())
        
        # Build KD-tree for efficient nearest neighbor search
        tree = cKDTree(positions, boxsize=self.box)
        
        # Find all pairs within cutoff radius
        pairs = tree.query_pairs(self.cutoff)
        
        # Build neighbor dictionary
        self.neighbor_list = defaultdict(list)
        for i, j in pairs:
            self.neighbor_list[ids[i]].append(ids[j])
            self.neighbor_list[ids[j]].append(ids[i])
    
    def get_migration_energy(self, defect_id, destination):
        """
        Calculate migration energy for a defect to move to a new position,
        accounting for interactions with nearby defects.
        """
        if self.neighbor_list is None:
            self.update_neighbor_list()
            
        defect = self.defects[defect_id]
        defect_type = defect['type']
        current_pos = defect['position']
        
        # Base migration energy (depends on defect type)
        base_energy = self._get_base_migration_energy(defect_type)
        
        # Calculate energy modification due to nearby defects
        energy_modifier = 0.0
        
        for neighbor_id in self.neighbor_list.get(defect_id, []):
            neighbor = self.defects[neighbor_id]
            neighbor_type = neighbor['type']
            neighbor_pos = neighbor['position']
            
            # Calculate distance with periodic boundary conditions
            vec = neighbor_pos - current_pos
            vec = vec - self.box * np.round(vec / self.box)
            distance = np.linalg.norm(vec)
            
            # Apply interaction model (this is where your physics model comes in)
            interaction = self._get_interaction_energy(defect_type, neighbor_type, distance)
            energy_modifier += interaction
            
        return base_energy + energy_modifier
    
    def _get_base_migration_energy(self, defect_type):
        """
        Return the base migration energy for a specific defect type.
        Override this method with your physical model.
        """
        # Example values - replace with your actual migration energies
        energy_dict = {
            'vacancy': 0.5,  # eV
            'interstitial': 0.2,  # eV
            'impurity': 1.0   # eV
        }
        return energy_dict.get(defect_type, 0.5)
    
    def _get_interaction_energy(self, defect_type1, defect_type2, distance):
        """
        Calculate interaction energy between two defects.
        Override this method with your physical model.
        """
        # Simple example - inverse square law with type-dependent prefactor
        if distance < 1e-10:
            return 0.0  # Avoid division by zero
            
        # Example interaction parameters - replace with your physical model
        interaction_matrix = {
            ('vacancy', 'vacancy'): 0.1,
            ('interstitial', 'interstitial'): 0.2,
            ('vacancy', 'interstitial'): -0.3,
            ('impurity', 'vacancy'): 0.4,
            ('impurity', 'interstitial'): 0.2
        }
        
        # Get parameter or use default
        key = (defect_type1, defect_type2)
        rev_key = (defect_type2, defect_type1)
        
        if key in interaction_matrix:
            param = interaction_matrix[key]
        elif rev_key in interaction_matrix:
            param = interaction_matrix[rev_key]
        else:
            param = 0.05  # Default
            
        # Simple interaction model (replace with your physics)
        return param / (distance * distance)
    
    def run_kmc_step(self):
        """Perform a single KMC step."""
        if not self.defects:
            return None
            
        # Calculate rates for all possible events
        events = []
        rates = []
        
        # For each defect, consider possible jumps (to nearest neighbor sites)
        for defect_id, defect in self.defects.items():
            current_pos = defect['position']
            
            # Generate possible jump destinations (nearest neighbors on the lattice)
            # This will depend on your lattice structure
            possible_jumps = self._get_possible_jumps(defect_id)
            
            for dest in possible_jumps:
                # Calculate migration energy and rate
                e_mig = self.get_migration_energy(defect_id, dest)
                rate = self._get_rate(e_mig)
                
                if rate > 0:
                    events.append((defect_id, dest))
                    rates.append(rate)
        
        # If no events are possible, return
        if not events:
            return None
            
        # Convert to numpy array for efficient operations
        rates = np.array(rates)
        
        # Calculate cumulative rates
        cumulative_rates = np.cumsum(rates)
        total_rate = cumulative_rates[-1]
        
        # Choose an event with probability proportional to its rate
        rand = np.random.random() * total_rate
        event_index = np.searchsorted(cumulative_rates, rand)
        
        # Execute the chosen event
        chosen_defect, destination = events[event_index]
        self._move_defect(chosen_defect, destination)
        
        # Time increment (exponential distribution)
        dt = -np.log(np.random.random()) / total_rate
        
        return dt
    
    def _get_possible_jumps(self, defect_id):
        """
        Get possible jump destinations for a defect.
        Override this with your lattice structure.
        """
        defect = self.defects[defect_id]
        current_pos = defect['position']
        
        # Example: simple cubic lattice, jumps to 6 nearest neighbors
        jumps = []
        for axis in range(3):
            for direction in [-1, 1]:
                jump = current_pos.copy()
                jump[axis] += direction
                
                # Apply periodic boundary conditions
                jump = jump % self.box
                
                jumps.append(jump)
                
        return jumps
    
    def _get_rate(self, energy):
        """Calculate transition rate using Arrhenius equation."""
        # v0 * exp(-E/kT)
        v0 = 1e13  # Attempt frequency (typical phonon frequency), in Hz
        return v0 * np.exp(-energy / (self.kb * self.temperature))
    
    def _move_defect(self, defect_id, new_position):
        """Move a defect to a new position."""
        if defect_id in self.defects:
            self.defects[defect_id]['position'] = new_position
            
            # Neighbor list is now outdated
            self.neighbor_list = None

# Create simulation
sim = DefectKMC(box_dimensions=(100, 100, 100), cutoff_radius=5.0, temperature=500)

# Add defects
sim.add_defect([25, 25, 25], 'vacancy')
sim.add_defect([30, 30, 30], 'interstitial')

# Run simulation for 1000 steps
time = 0
for _ in range(1000):
    dt = sim.run_kmc_step()
    if dt is not None:
        time += dt
    print(f"Simulation time: {time} seconds")
