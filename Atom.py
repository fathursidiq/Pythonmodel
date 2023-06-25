import numpy as np
import math

class Atom:
    def __init__(self, symbol, mass, x, y, z):
        self.symbol = symbol
        self.mass = mass
        self.x = x
        self.y = y
        self.z = z
        self.temperature = temperature
        self.energy = energy

    def __str__(self):
        return f"{self.symbol} ({self.mass} amu) - ({self.x}, {self.y}, {self.z})"

    def lennard_jones_potential(self, other_atom):
        epsilon = 1.0  # Lennard-Jones parameter
        sigma = 1.0  # Lennard-Jones parameter

        r = np.sqrt((self.x - other_atom.x)**2 + (self.y - other_atom.y)**2 + (self.z - other_atom.z)**2)
        lj_energy = 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6)
        return lj_energy

    def increase_temperature(self, delta_temp):
        self.temperature += delta_temp

    def decrease_temperature(self, delta_temp):
        self.temperature -= delta_temp

    def calculate_entropy(self):
        boltzmann_constant = 1.380649e-23
        entropy = self.energy / (self.temperature * boltzmann_constant)
        return entropy
    def __str__(self):
        return f"{self.symbol}: Energy={self.energy}, Temperature={self.temperature}"
    def update_position(self, time_step):
        self.x += self.velocity[0] * time_step
        self.y += self.velocity[1] * time_step
        self.z += self.velocity[2] * time_step

    def __str__(self):
        return f"{self.symbol} ({self.x:.3f}, {self.y:.3f}, {self.z:.3f})"
    def distance_to(self, other_atom):
        dx = self.x - other_atom.x
        dy = self.y - other_atom.y
        dz = self.z - other_atom.z
        return math.sqrt(dx**2 + dy**2 + dz**2)

    def __str__(self):
        return f"{self.symbol} ({self.x:.3f}, {self.y:.3f}, {self.z:.3f})"
    def create_region(atoms):
        min_x = min(atom.x for atom in atoms)
        max_x = max(atom.x for atom in atoms)
        min_y = min(atom.y for atom in atoms)
        max_y = max(atom.y for atom in atoms)
        min_z = min(atom.z for atom in atoms)
        max_z = max(atom.z for atom in atoms)
        
        region = {
            'min_x': min_x,
            'max_x': max_x,
            'min_y': min_y,
            'max_y': max_y,
            'min_z': min_z,
            'max_z': max_z
            }
        return region
class CarbonDioxide:
    def __init__(self):
        self.atoms = [
            Atom("C", 12.01, 0.0, 0.0, 0.0),
            Atom("O", 16.00, 1.2, 0.0, 0.0),
            Atom("O", 16.00, -1.2, 0.0, 0.0)
        ]

    def __str__(self):
        return "\n".join(str(atom) for atom in self.atoms)


class Boundary:
    def __init__(self, xmin, xmax, ymin, ymax, zmin, zmax):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.zmin = zmin
        self.zmax = zmax

    def is_within_boundary(self, atom):
        return (
            self.xmin <= atom.x <= self.xmax and
            self.ymin <= atom.y <= self.ymax and
            self.zmin <= atom.z <= self.zmax
        )
        import numpy as np
    
class deform:
     def deform(self, displacement):
        self.x += displacement[0]
        self.y += displacement[1]
        self.z += displacement[2]

# Example usage
atoms = [
    Atom("H", 0.0, 0.0, 0.0),
    Atom("O", 1.0, 1.0, 1.0),
    Atom("H", -1.0, -1.0, -1.0)
]

cutoff_distance = 1.5

for atom in atoms:
    neighbors = find_neighbors(atom, atoms, cutoff_distance)
    print("Atom:", atom)
    print("Neighbors:")
    for neighbor in neighbors:
        print("-", neighbor)
    print()
# Example usage
atom = Atom("H", 100, 300)
print(atom)  # Output: H: Energy=100, Temperature=300

atom.increase_temperature(50)
print(atom)  # Output: H: Energy=100, Temperature=350

atom.decrease_temperature(100)
print(atom)  # Output: H: Energy=100, Temperature=250

entropy = atom.calculate_entropy()
print(f"Entropy: {entropy}")


# Example usage
#atom1 = Atom("Ar", 39.95, 0.0, 0.0, 0.0)
#atom2 = Atom("Ar", 39.95, 3.0, 0.0, 0.0)

#energy = atom1.lennard_jones_potential(atom2)
#print(f"Lennard-Jones Energy: {energy}")

# Example usage
#boundary = Boundary(0.0, 10.0, 0.0, 10.0, 0.0, 10.0)

# Example usage
#co2 = CarbonDioxide()
#print(co2)
# Example usage
#atom = Atom("H", 0.0, 0.0, 0.0)
#print("Original Atom:", atom)

#displacement = [0.1, 0.2, 0.3]
#atom.deform(displacement)
#print("Deformed Atom:", atom)
# Example usage
#atom = Atom("H", 0.0, 0.0, 0.0, [0.1, 0.2, 0.3])
#print("Initial Atom:", atom)

#time_step = 0.1

#atom.update_position(time_step)
#print("Atom after time step:", atom)
# Open the file for writing
with open("atoms.xyz", "w") as file:
    # Write the number of atoms
    file.write(str(len(atoms)) + "\n")

    # Write a comment line
    file.write("Atom positions\n")

    # Write atom coordinates to the file
    for atom in atoms:
        file.write(str(atom) + "\n")

print("XYZ file created: atoms.xyz")