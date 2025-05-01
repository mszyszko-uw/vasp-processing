import os
import shutil

# skrypt do przygotowywania plików wejściowych do VASP pod serie obliczeń (w tym przypadku na podstawie struktur z SOD)
calc_path = '' 
struct_path = 'komórki_3_150/x-0.5/CALCS'
POTCAR = 'POTCAR'
INCAR = 'INCAR' 
KPOINTS = 'KPOINTS'

i = 0
out_path = os.path.join(calc_path, '3_150_x-0.5')
os.makedirs(out_path)
for file in os.listdir(struct_path):
    if file[-4:] == 'vasp':
        i += 1
        curr_folder = os.path.join(out_path, 'C' + str(i))
        os.makedirs(curr_folder)
        shutil.copy(os.path.join(struct_path, file), curr_folder + '/POSCAR')
        shutil.copy(INCAR, curr_folder + '/INCAR')
        shutil.copy(POTCAR, curr_folder + '/POTCAR')
        shutil.copy(KPOINTS, curr_folder + '/KPOINTS')

