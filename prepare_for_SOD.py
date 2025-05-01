import os
import shutil

# Skrypt do przygotowywania folderów pod uruchamianie SOD
# Potrzebne do tego są pliki początkowe INSOD i SGO z ustawieniami i operatorami symetrii

def prepare_INSOD(path: str, destination: str, cell: tuple, lattice_constants: tuple, n_substitusions: int) -> None: 
    x = n_substitusions/(cell[0]*cell[1]) # concentration of doping, used in some cases
    with open(path, 'r') as file:
        file.seek(0)
        lines = file.readlines()
        lines[1] = f'{cell[0]}x{cell[1]}x{cell[2]} MoS2 substituted with W \n'
        lines[4] = f'{lattice_constants[0]} {lattice_constants[1]} {lattice_constants[2]} {lattice_constants[3]} {lattice_constants[4]} {lattice_constants[5]} \n'
        lines[21] = f'{cell[0]} {cell[1]} {cell[2]} \n'
        lines[27] = f'{n_substitusions} \n'
        with open(destination, 'a+') as destination_file:
            destination_file.seek(0)
            destination_file.writelines(lines)
    return None


lattice_constants = (3.160, 3.160, 19.1699943542, 90.000, 90.000, 120.000) # lattice constants a, b, c, alpha, beta, gamma
cell_sizes = [(4, 4, 1)] # size of the supercell

for cell in cell_sizes:
    folder_path = f'/home/maciej/Documents/2D_Materials/4x4-ClusterExpansion/4x4x1_cell/komórki_3_150'
    os.makedirs(folder_path)
    for i in range(1, cell[0]*cell[1]):
        subfolder_path = folder_path + f'/x-{round(i/(cell[0]*cell[1]), 3)}'
        os.makedirs(subfolder_path)
        SGO_path = '/home/maciej/Documents/2D_Materials/4x4-ClusterExpansion/4x4x1_cell/SGO'
        shutil.copy(SGO_path, subfolder_path + '/SGO')
        INSOD_path = '/home/maciej/Documents/2D_Materials/4x4-ClusterExpansion/4x4x1_cell/INSOD'
        prepare_INSOD(INSOD_path, subfolder_path + '/INSOD', cell, lattice_constants, i)

        
