from pymatgen.io.vasp.outputs import Chgcar
from pymatgen.io.vasp import Poscar
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# skrypt do rysowania Charge Density z pliku CHGCAR
chgcar_file = "CHGCAR"  
poscar_file = "POSCAR" 
chgcar = Chgcar.from_file(chgcar_file)
poscar = Poscar.from_file(poscar_file)

# chg_cell_size odpowiada za powielenie gęstrości z obliczonej komórki - czasem łatwiejsze to wizualizacji czegoś
chg_cell_size = (2, 2, 1)
charge_density = chgcar.data["total"]
lattice = poscar.structure.lattice.matrix
supercell_charge_density = np.tile(charge_density, chg_cell_size)
charge_density_normalized = supercell_charge_density / np.max(supercell_charge_density)

nx, ny, nz = supercell_charge_density.shape
fx = np.linspace(0, chg_cell_size[0], nx)
fy = np.linspace(0, chg_cell_size[1], ny)
fz = np.linspace(0, chg_cell_size[2], nz)
FX, FY, FZ = np.meshgrid(fx, fy, fz, indexing='ij')


X = FX * lattice[0][0] + FY * lattice[1][0] + FZ * lattice[2][0]
Y = FX * lattice[0][1] + FY * lattice[1][1] + FZ * lattice[2][1]
Z = FX * lattice[0][2] + FY * lattice[1][2] + FZ * lattice[2][2]


fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')


threshold = 0.1  
density_mask = charge_density_normalized > threshold


scatter = ax.scatter(
    X[density_mask], Y[density_mask], Z[density_mask],
    c=charge_density_normalized[density_mask], cmap='viridis', alpha=0.5, s=1
)

fig.colorbar(scatter, ax=ax, label="Normalized Charge Density")
ax.set_axis_off()
plt.show()
