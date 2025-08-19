import shutil
import os
import h5py
import matplotlib.pyplot as plt
import sys
import argparse

# SETTINGS
pseudopotentials_path = 'potpaw_PBE'
#INCAR_templates_path = 'INCAR_templates'
directory_name = 'Calculations_MoTe2'
#system_name = 'bulk MoTe2'
INCAR_settings = {
    'SYSTEM': 'bulk MoTe2',
    'ISMEAR': 0,
    'SIGMA': 0.05,
    'IVDW': 11,
    'LMAXMIX': 4,
    'LASPH': '.TRUE.'
}
dry_settings = {'ALGO': None, 'NELM': 1, 'LWAVE': '.FALSE.', 'LCHARG': '.FALSE.'}

conv_kmesh = [(6, 6, 2), (7, 7, 2)]
conv_encut = [250, 300, 350, 400]



try:
    os.mkdir(directory_name)
except:
    pass


def create_step(step_number, part):
    try:
        os.mkdir(os.path.join(directory_name, step_number))
    except:
        pass

    if step_number == 'step01':
        step01(os.path.join(directory_name, 'step01'), part)
    elif step_number == 'step02':
        step2_geo(os.path.join(directory_name, 'step02'), part)
    elif step_number == 'step03':
        step3_conv_test(os.path.join(directory_name, 'step03'), part, kmesh_list=conv_kmesh, encut_list=conv_encut)

def step01(folder, part):

    if part == 'dry':
    # part 1: dry run
        dry_run = create_folder(folder, '00_t_dry')
        copy_files('', dry_run, ['KPOINTS', 'POSCAR'])
        
        with open('POSCAR', 'r') as f:
            lines = f.readlines()
            elements = lines[5].split()
        
        create_potcar(elements, os.path.join(os.path.join(dry_run, 'POTCAR')))
        with open(os.path.join(dry_run, 'INCAR'), 'w') as f: pass
        modify_incar(os.path.join(dry_run, 'INCAR'), add= INCAR_settings | dry_settings)
        #create_job()
        #run_vasp()
        #check_errors()

    # part 2: scf run
    elif part == 'scf':
        scf_run = create_folder(folder, '01_t_scf')
        dry_run = os.path.join(folder, '00_t_dry')
        copy_files(dry_run, scf_run, ['KPOINTS', 'POSCAR', 'POTCAR', 'INCAR'])
        #os.rename(os.path.join(scf_run, 'CONTCAR'), os.path.join(scf_run, 'POSCAR'))
        NKPTS, NBANDS = get_parallelization_variables(os.path.join(dry_run, 'vaspout.h5'))
        INCAR_settings['NBANDS'] = NBANDS
        INCAR_settings['NPAR'] = 4
        INCAR_settings['KPAR'] = NKPTS
        modify_incar(os.path.join(scf_run, 'INCAR'), remove= dry_settings, add= INCAR_settings)
        #create_job()
        #run_vasp()
        #check_errors()


def step2_geo(folder, part):
    opt_settings = {   
        'IBRION': 2,
        'ISIF': 3,
        'EDIFFG': -1E-3,
        'NSW': 50
    }

    # part 1: dry run
    if part == 'dry':
        dry_run = create_folder(folder, '00_t_dry')
        step01_scf = os.path.join(directory_name, 'step01', '01_t_scf')
        copy_files(step01_scf, dry_run, ['KPOINTS', 'POSCAR', 'POTCAR', 'WAVECAR', 'INCAR'])
        #copy_files(last_step, dry_run, ['KPOINTS', 'POSCAR', 'POTCAR', 'WAVECAR', 'INCAR'])
        modify_incar(os.path.join(dry_run, 'INCAR'), add= opt_settings | dry_settings, remove=['KPAR', 'NPAR'])
        #create_job()
        #run_vasp()
        #check_errors()

    # part 2: conjugate-gradient optimisation
    elif part == 'cg_opt':
        geo1_run = create_folder(folder, '01_t_geo')
        dry_run = os.path.join(folder, '00_t_dry')
        copy_files(dry_run, geo1_run, ['KPOINTS', 'POSCAR', 'POTCAR', 'INCAR'])
        NKPTS, NBANDS = get_parallelization_variables(os.path.join(dry_run, 'vaspout.h5'))
        INCAR_settings['NBANDS'] = NBANDS
        INCAR_settings['NPAR'] = 4
        INCAR_settings['KPAR'] = NKPTS
        modify_incar(os.path.join(geo1_run, 'INCAR'), add= INCAR_settings, remove=dry_settings)
        #create_job()
        #run_vasp()
        #check_errors()

    # part 3: Newton algorithm optimisation
    elif part == 'nw_opt':
        geo2_run = create_folder(folder, '02_t_geo')
        geo1_run = os.path.join(folder, '01_t_geo')
        copy_files(geo1_run, geo2_run, ['KPOINTS', 'CONTCAR', 'POTCAR', 'INCAR', 'WAVECAR'])
        #copy_files(geo1_run, geo2_run, ['KPOINTS', 'CONTCAR', 'POTCAR', 'INCAR', 'WAVECAR'])
        os.rename(os.path.join(geo2_run, 'CONTCAR'), os.path.join(geo2_run, 'POSCAR'))
        opt_settings['IBRION'] = 1
        opt_settings['NSW'] = 100
        modify_incar(os.path.join(geo2_run, 'INCAR'), add= opt_settings)
        #create_job()
        #run_vasp()
        #check_errors()

    # part 4: scf run
    elif part == 'scf':
        scf_run = create_folder(folder, '03_t_scf')
        geo2_run = os.path.join(folder, '02_t_geo')
        copy_files(geo2_run, scf_run, ['KPOINTS', 'CONTCAR', 'POTCAR', 'INCAR', 'WAVECAR'])
        #copy_files(geo2_run, scf_run, ['KPOINTS', 'CONTCAR', 'POTCAR', 'INCAR', 'WAVECAR'])
        os.rename(os.path.join(scf_run, 'CONTCAR'), os.path.join(scf_run, 'POSCAR'))
        modify_incar(os.path.join(scf_run, 'INCAR'), remove= opt_settings)
        #create_job()
        #run_vasp()
        #check_errors()
    elif part == 'report':
        scf_run = os.path.join(folder, '03_t_scf')
        create_report(os.path.join(scf_run, 'vaspout.h5'), os.path.join(scf_run, 'report.pdf'))


def step3_conv_test(folder, part, kmesh_list, encut_list):
    opt_settings = {   
        'IBRION': 2,
        'ISIF': 3,
        'EDIFFG': -1E-3,
        'NSW': 50
    }
    k_folder_paths = []

    # part 1: dry run
    if part == 'dry':
        step02_scf = os.path.join(directory_name, 'step02', '03_t_scf')
        for kmesh in kmesh_list:
            k_folder = create_folder(folder, f'k{kmesh[0]}x{kmesh[1]}x{kmesh[2]}')
            k_folder_paths.append(k_folder)
            dry_run = create_folder(k_folder, '00_t_dry')
            copy_files(step02_scf, dry_run, ['CONTCAR', 'POTCAR', 'INCAR', 'WAVECAR'])
            os.rename(os.path.join(dry_run, 'CONTCAR'), os.path.join(dry_run, 'POSCAR'))
            modify_incar(os.path.join(dry_run, 'INCAR'), remove=['KPAR', 'NPAR'], add= opt_settings | dry_settings)
            with open(os.path.join(dry_run, 'KPOINTS'), 'w') as f:
                f.write(f'{kmesh[0]}x{kmesh[1]}x{kmesh[2]}\n')
                f.write('0\n')
                f.write('Gamma\n')
                f.write(f'{kmesh[0]} {kmesh[1]} {kmesh[2]}\n')
        #create_array_job()
        #run_array_job()
        #check_errors()

    # part 2: conjucate-gradient optimisations
    if part == 'cg_opt':
        for kmesh in kmesh_list:
            k_folder = os.path.join(folder, f'k{kmesh[0]}x{kmesh[1]}x{kmesh[2]}')
            for encut in encut_list:
                e_geo = create_folder(k_folder, f'e{encut}')
                geo1 = create_folder(e_geo, '01_t_geo')
                copy_files(os.path.join(k_folder, '00_t_dry'), geo1, ['POSCAR', 'POTCAR', 'INCAR', 'KPOINTS', 'WAVECAR'])
                NKPTS, NBANDS = get_parallelization_variables(os.path.join(k_folder, '00_t_dry', 'vaspout.h5'))
                INCAR_settings['NBANDS'] = NBANDS
                INCAR_settings['NPAR'] = 4
                INCAR_settings['KPAR'] = NKPTS
                INCAR_settings['ENCUT'] = encut
                modify_incar(os.path.join(geo1, 'INCAR'), remove=dry_settings, add=INCAR_settings)
        #create_array_job()
        #run_array_job()
        #check_errors()

    # part 3: Newton algorithm optimisations
    if part == 'nw_opt':
        for kmesh in kmesh_list:
            k_folder = os.path.join(folder, f'k{kmesh[0]}x{kmesh[1]}x{kmesh[2]}')
            for encut in encut_list:
                geo2 = create_folder(os.path.join(k_folder, f'e{encut}'), '02_t_geo')
                copy_files(os.path.join(k_folder, f'e{encut}', '01_t_geo'), geo2, ['CONTCAR', 'POTCAR', 'INCAR', 'KPOINTS', 'WAVECAR'])
                os.rename(os.path.join(geo2, 'CONTCAR'), os.path.join(geo2, 'POSCAR'))
                opt_settings['IBRION'] = 1
                opt_settings['NSW'] = 100
                modify_incar(os.path.join(geo2, 'INCAR'), add=opt_settings)
        #create_array_job()
        #run_array_job()
        #check_errors()

    # part 4: scf
    if part == 'scf':
        for kmesh in kmesh_list:
            k_folder = os.path.join(folder, f'k{kmesh[0]}x{kmesh[1]}x{kmesh[2]}')
            for encut in encut_list:
                scf = create_folder(os.path.join(k_folder, f'e{encut}'), '03_t_scf')
                copy_files(os.path.join(k_folder, f'e{encut}', '02_t_geo'), scf, ['CONTCAR', 'POTCAR', 'INCAR', 'KPOINTS', 'WAVECAR'])
                os.rename(os.path.join(scf, 'CONTCAR'), os.path.join(scf, 'POSCAR'))
                modify_incar(os.path.join(scf, 'INCAR'), remove=opt_settings)
        #create_array_job()
        #run_array_job()
        #check_errors()
    
    if part == 'report':
        for kmesh in kmesh_list:
            k_folder = os.path.join(folder, f'k{kmesh[0]}x{kmesh[1]}x{kmesh[2]}')
            for encut in encut_list:
                create_report(os.path.join(k_folder, f'e{encut}', '03_t_scf', 'vaspout.h5'), os.path.join(k_folder, f'e{encut}', '03_t_scf', 'report.pdf'))
                #TO DO: zebranie wyników convergence test do .csv/.xlsx


def create_potcar(elements, potcar_path):
    with open(os.path.join(pseudopotentials_path, 'recomended_pseudopotentials.txt'), 'r') as f, open(potcar_path, 'w') as POTCAR:
        lines = f.readlines()
        potcar_paths = {}
        for line in lines:
            if line.split(',')[0].split('_')[0] in elements:
                if 'd' in line.split(',')[2]:
                    INCAR_settings['LMAXMIX'] = 4
                    INCAR_settings['LASPH'] = ".TRUE."
                potcar_paths[line.split(',')[0].split('_')[0]] = os.path.join(pseudopotentials_path, line.split(',')[0], 'POTCAR')
        for element in elements:
            with open(potcar_paths[element], 'r') as file:
                POTCAR.write(file.read())


def copy_files(from_path, to_path, files):
    for file in files:
        shutil.copyfile(os.path.join(from_path, file), os.path.join(to_path, file))


def create_folder(path, folder_name):
    name = os.path.join(path, folder_name)
    try:
        os.mkdir(name)
    except:
        pass
    return name
      

def modify_incar(file, remove=None, add=None):
    remove = remove or []
    add = add or {}

    with open(file, 'r') as f:
        lines = f.readlines()

    incar_dict = {}
    for line in lines:
        if '=' in line:
            key, value = line.split('=', 1)
            incar_dict[key.strip().upper()] = value.strip()

    if isinstance(remove, dict):
        remove_keys = [key.upper() for key in remove.keys()]
    else:
        remove_keys = [key.upper() for key in remove]

    for key in remove_keys:
        incar_dict.pop(key, None)

    for key, value in add.items():
        key_upper = key.upper()
        if isinstance(value, bool):
            incar_dict[key_upper] = f".{str(value).upper()}."
        else:
            incar_dict[key_upper] = str(value)

    with open(file, 'w') as f:
        for key in incar_dict: 
            f.write(f"{key} = {incar_dict[key]}\n")


def get_parallelization_variables(vaspout_file):
    with h5py.File(vaspout_file, "r") as f:
        eigenvalues = f['results/electron_eigenvalues/eigenvalues'][:]
        _ , nkpts, nbands = eigenvalues.shape
        nbands = nbands if nbands % 4 == 0 else ((nbands // 4) + 1) * 4
        nkpts = max(i for i in range(1, 11) if nkpts % i == 0)
        return nkpts, nbands
    

def create_report(vaspout, save_path):
    with h5py.File(vaspout, "r") as f:
        forces = f['intermediate/ion_dynamics/forces'][:]
        stress = f['intermediate/ion_dynamics/stress'][:]
        oszicar = f['intermediate/ion_dynamics/oszicar'][:]

    results = {
        "ENCUT": "250 eV",
        "Stress in kB": f"[{stress[0,0,0]:.4f}, {stress[0, 1, 1]:.4f}, {stress[0, 2, 2]:.4f}]"
    }
    try:
        results['ENCUT'] = f['input/incar/ENCUT'][:]
    except:
        results['ENCUT'] = 'from POTCAR'

    fig = plt.figure(figsize=(8.3, 11.7))  # A4 portrait in inches

    # --- Top: Scalar results and Forces table ---
    ax1 = fig.add_axes([0.05, 0.55, 0.9, 0.4])
    ax1.axis('off')

    # Title
    ax1.text(0.5, 1.05, "Simulation Report", ha="center", va="top", fontsize=16, weight="bold")

    # Scalar results
    y_pos = 0.95
    for key, value in results.items():
        ax1.text(0.05, y_pos, f"{key}:", fontsize=12, weight="bold")
        ax1.text(0.35, y_pos, f"{value}", fontsize=12)
        y_pos -= 0.07

    # Forces table header
    ax1.text(0.05, y_pos - 0.05, "Atomic Forces [eV/Å]:", fontsize=12, weight="bold")

    # Forces table
    col_labels = ["Fx", "Fy", "Fz"]
    table_data = [[f"{f:.4f}" for f in row] for row in forces[0]]

    table = ax1.table(
        cellText=table_data,
        colLabels=col_labels,
        colLoc='center',
        loc='lower center',
        cellLoc='center',
        bbox=[0.05, 0.0, 0.9, 0.65]
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)

    # --- Bottom: Chart ---
    ax2 = fig.add_axes([0.15, 0.08, 0.7, 0.35])
    ax2.plot(oszicar[:, 0], oszicar[:, 2], marker='o', label='dE')
    ax2.plot(oszicar[:, 0], oszicar[:, 3], marker='o', label='d eps')
    ax2.set_yscale('symlog', linthresh=1e-4)
    ax2.set_xlabel("Step")
    ax2.set_ylabel("log(dE/d eps)")
    ax2.legend()
    ax2.grid(True)

    # Save to PDF
    plt.savefig(save_path, format="pdf")
    plt.close()



def create_job():
    pass


def run_vasp():
    pass


def check_errors():
    pass


parser = argparse.ArgumentParser(description="Input agruments for the script")
parser.add_argument("--step", type=str, help="Step number")
parser.add_argument("--part", type=str, help="Part of the step")
args = parser.parse_args()

create_step(step_number=args.step, part=args.part)