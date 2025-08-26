import shutil
import os
import h5py
import matplotlib.pyplot as plt
import sys
import argparse
import numpy as np
import math

# SETTINGS loaded from input (globals: paths, INCAR_settings, etc.)
from input import *

# Create main calculation directory if it does not exist
try:
    os.mkdir(calculation_path)
except:
    pass


# Dispatcher for creating specific calculation steps
# Each step corresponds to a type of VASP step (scf, geo, bandstructure, etc.)

def create_step(step, part):
    try:
        os.mkdir(os.path.join(calculation_path, step))
    except:
        pass

    if step == 't_scf':
        print(t_scf(os.path.join(calculation_path, 't_scf'), part))
    elif step == 't_geo':
        print(t_geo(os.path.join(calculation_path, 't_geo'), part))
    elif step == 't_conv_test':
        paths = t_conv_test(os.path.join(calculation_path, 't_conv_test'), part, kmesh_list=conv_kmesh, encut_list=conv_encut)
        for p in paths:
            print(p)
    elif step == 't_bs':
        print(t_bs(os.path.join(calculation_path, 't_bs'), part))
    elif step == 't_dos':
        print(t_dos(os.path.join(calculation_path, 't_dos'), part))
    elif step == 't_scf_so':
        print(t_scf_so(os.path.join(calculation_path, 't_scf_so'), part))
    elif step == 't_geo_so':
        print(t_geo_so(os.path.join(calculation_path, 't_geo_so'), part))
    elif step == 't_bs_so':
        print(t_bs_so(os.path.join(calculation_path, 't_bs_so'), part))
    elif step == 't_dos_so':
        print(t_dos_so(os.path.join(calculation_path, 't_dos_so'), part))


# ===== Standard SCF step (non-SOC) =====
# Handles dry run and SCF run

def t_scf(folder, part, previous=''):
    if part == 'dry':
        # Dry run: create folder, copy base files, adjust INCAR with dry settings
        dry_run = create_folder(folder, '00_t_dry')
        copy_files(previous, dry_run, ['KPOINTS', 'POSCAR', 'POTCAR', 'INCAR'])
        modify_incar(os.path.join(dry_run, 'INCAR'), add=  dry_settings | INCAR_settings )
        return os.path.abspath(dry_run)

    elif part == 'scf':
        # SCF run: copy files from dry run, set parallelization, update INCAR
        scf_run = create_folder(folder, '01_t_scf')
        dry_run = os.path.join(folder, '00_t_dry')
        copy_files(dry_run, scf_run, ['KPOINTS', 'POSCAR', 'POTCAR', 'INCAR'])
        KPAR, NBANDS = get_parallelization_variables(os.path.join(dry_run, 'vaspout.h5'))
        INCAR_settings['NBANDS'] = NBANDS
        INCAR_settings['NPAR'] = 4
        INCAR_settings['KPAR'] = KPAR
        modify_incar(os.path.join(scf_run, 'INCAR'), remove= dry_settings, add= INCAR_settings)
        return os.path.abspath(scf_run)


# ===== Geometry optimization step (non-SOC) =====
# Handles dry run, Conjugate-Gradient opt, Newton opt, SCF, report

def t_geo(folder, part, previous=os.path.join(calculation_path, 't_scf/01_t_scf')):
    opt_settings = {   
        'IBRION': 2,   # conjugate gradient
        'EDIFFG': -1E-3,
        'NSW': 50,
    }
    if "ISIF" in globals():
        opt_settings["ISIF"] = ISIF
    else:
        opt_settings["ISIF"] = 3

    if part == 'dry':
        # Dry run: copy input, add opt + dry settings
        dry_run = create_folder(folder, '00_t_dry')
        copy_files(previous, dry_run, ['KPOINTS', 'POSCAR', 'POTCAR', 'INCAR'])
        modify_incar(os.path.join(dry_run, 'INCAR'), add= opt_settings | dry_settings | INCAR_settings, remove=['KPAR', 'NPAR'])
        return os.path.abspath(dry_run)

    elif part == 'cg_opt':
        # Conjugate gradient optimization run
        geo1_run = create_folder(folder, '01_t_geo')
        dry_run = os.path.join(folder, '00_t_dry')
        copy_files(dry_run, geo1_run, ['KPOINTS', 'POSCAR', 'POTCAR', 'INCAR'])
        KPAR, NBANDS = get_parallelization_variables(os.path.join(dry_run, 'vaspout.h5'))
        INCAR_settings['NBANDS'] = NBANDS
        INCAR_settings['NPAR'] = 4
        INCAR_settings['KPAR'] = KPAR
        modify_incar(os.path.join(geo1_run, 'INCAR'), add= INCAR_settings, remove=dry_settings)
        return os.path.abspath(geo1_run)

    elif part == 'nw_opt':
        # Newton optimization using CONTCAR → POSCAR
        geo2_run = create_folder(folder, '02_t_geo')
        geo1_run = os.path.join(folder, '01_t_geo')
        copy_files(geo1_run, geo2_run, ['KPOINTS', 'CONTCAR', 'POTCAR', 'INCAR', 'WAVECAR'])
        os.rename(os.path.join(geo2_run, 'CONTCAR'), os.path.join(geo2_run, 'POSCAR'))
        opt_settings['IBRION'] = 1  # switch to quasi-Newton
        opt_settings['NSW'] = 100
        modify_incar(os.path.join(geo2_run, 'INCAR'), add= opt_settings)
        return os.path.abspath(geo2_run)

    elif part == 'scf':
        # Final SCF run after geometry optimization
        scf_run = create_folder(folder, '03_t_scf')
        geo2_run = os.path.join(folder, '02_t_geo')
        copy_files(geo2_run, scf_run, ['KPOINTS', 'CONTCAR', 'POTCAR', 'INCAR', 'WAVECAR'])
        os.rename(os.path.join(scf_run, 'CONTCAR'), os.path.join(scf_run, 'POSCAR'))
        modify_incar(os.path.join(scf_run, 'INCAR'), remove= opt_settings)
        return os.path.abspath(scf_run)
        
    elif part == 'report':
        # Generate PDF report from HDF5 output - checks stress and forces
        scf_run = os.path.join(folder, '03_t_scf')
        create_report(os.path.join(scf_run, 'vaspout.h5'), os.path.join(scf_run, 'report.pdf'))


# ===== Convergence test step =====
# Runs over k-meshes and ENCUT values, performs dry, Conjugate-Gradient, Newton, SCF, report

def t_conv_test(folder, part, kmesh_list, encut_list, previous=os.path.join(calculation_path, 't_geo/03_t_scf')):
    opt_settings = {   
        'IBRION': 2,
        'EDIFFG': -1E-3,
        'NSW': 50
    }
    if "ISIF" in globals():
        opt_settings["ISIF"] = ISIF
    else:
        opt_settings["ISIF"] = 3

    k_folder_paths = []
    calc_paths = []

    if part == 'dry':
        # For each kmesh, create folder, copy files, set corresponding KPOINTS
        for kmesh in kmesh_list:
            k_folder = create_folder(folder, f'k{kmesh[0]}x{kmesh[1]}x{kmesh[2]}')
            k_folder_paths.append(k_folder)
            dry_run = create_folder(k_folder, '00_t_dry')
            calc_paths.append(dry_run)
            copy_files(previous, dry_run, ['CONTCAR', 'POTCAR', 'INCAR', 'WAVECAR'])
            os.rename(os.path.join(dry_run, 'CONTCAR'), os.path.join(dry_run, 'POSCAR'))
            modify_incar(os.path.join(dry_run, 'INCAR'), remove=['KPAR', 'NPAR'], add= opt_settings | dry_settings | INCAR_settings)
            with open(os.path.join(dry_run, 'KPOINTS'), 'w') as f:
                f.write(f'{kmesh[0]}x{kmesh[1]}x{kmesh[2]}\n')
                f.write('0\n')
                f.write('Gamma\n')
                f.write(f'{kmesh[0]} {kmesh[1]} {kmesh[2]}\n')
        return calc_paths

    if part == 'cg_opt':
        # For each encut, perform geometry optimization using CG
        for kmesh in kmesh_list:
            k_folder = os.path.join(folder, f'k{kmesh[0]}x{kmesh[1]}x{kmesh[2]}')
            for encut in encut_list:
                e_geo = create_folder(k_folder, f'e{encut}')
                geo1 = create_folder(e_geo, '01_t_geo')
                calc_paths.append(geo1)
                copy_files(os.path.join(k_folder, '00_t_dry'), geo1, ['POSCAR', 'POTCAR', 'INCAR', 'KPOINTS', 'WAVECAR'])
                KPAR, NBANDS = get_parallelization_variables(os.path.join(k_folder, '00_t_dry', 'vaspout.h5'))
                INCAR_settings['NBANDS'] = NBANDS
                INCAR_settings['NPAR'] = 4
                INCAR_settings['KPAR'] = KPAR
                INCAR_settings['ENCUT'] = encut
                modify_incar(os.path.join(geo1, 'INCAR'), remove=dry_settings, add=INCAR_settings)
        return calc_paths

    if part == 'nw_opt':
        # For each encut, perform geometry optimization using Newton
        for kmesh in kmesh_list:
            k_folder = os.path.join(folder, f'k{kmesh[0]}x{kmesh[1]}x{kmesh[2]}')
            for encut in encut_list:
                geo2 = create_folder(os.path.join(k_folder, f'e{encut}'), '02_t_geo')
                calc_paths.append(geo2)
                copy_files(os.path.join(k_folder, f'e{encut}', '01_t_geo'), geo2, ['CONTCAR', 'POTCAR', 'INCAR', 'KPOINTS', 'WAVECAR'])
                os.rename(os.path.join(geo2, 'CONTCAR'), os.path.join(geo2, 'POSCAR'))
                opt_settings['IBRION'] = 1
                opt_settings['NSW'] = 100
                modify_incar(os.path.join(geo2, 'INCAR'), add=opt_settings)
        return calc_paths

    if part == 'scf':
        # Final SCF run for each encut
        for kmesh in kmesh_list:
            k_folder = os.path.join(folder, f'k{kmesh[0]}x{kmesh[1]}x{kmesh[2]}')
            for encut in encut_list:
                scf = create_folder(os.path.join(k_folder, f'e{encut}'), '03_t_scf')
                calc_paths.append(scf)
                copy_files(os.path.join(k_folder, f'e{encut}', '02_t_geo'), scf, ['CONTCAR', 'POTCAR', 'INCAR', 'KPOINTS', 'WAVECAR'])
                os.rename(os.path.join(scf, 'CONTCAR'), os.path.join(scf, 'POSCAR'))
                modify_incar(os.path.join(scf, 'INCAR'), remove=opt_settings)
        return calc_paths
    
    if part == 'report':
        # Generate convergence report over all kmesh/encut combinations
        for kmesh in kmesh_list:
            k_folder = os.path.join(folder, f'k{kmesh[0]}x{kmesh[1]}x{kmesh[2]}')
            for encut in encut_list:
                create_report(os.path.join(k_folder, f'e{encut}', '03_t_scf', 'vaspout.h5'), os.path.join(k_folder, f'e{encut}', '03_t_scf', 'report.pdf'))
        convergence_report(folder, kmesh_list, encut_list)
        return calc_paths


# ===== Band structure step (non-SOC) =====
# Dry run and band structure calculation

def t_bs(folder, part, previous=os.path.join(calculation_path, f't_conv_test/k{conv_kmesh[-1][0]}x{conv_kmesh[-1][1]}x{conv_kmesh[-1][2]}/e{max(conv_encut)}/03_t_scf')):
    bands_settings = {   
        'ICHARG': 11,
        'LWAVE': '.FALSE.',
        'LCHARG': '.FALSE.',
        'LORBIT': 11
    }

    if part == 'dry':
        # Dry run: set k-path, adjust INCAR
        dry_run = create_folder(folder, '00_t_dry')
        copy_files(previous, dry_run, ['POSCAR', 'POTCAR', 'CHGCAR', 'INCAR'])
        create_kpath_standard(os.path.join(dry_run, 'KPOINTS'), k_points=k_path, points=min_points)
        modify_incar(os.path.join(dry_run, 'INCAR'), add= bands_settings | dry_settings | INCAR_settings, remove=['KPAR', 'NPAR', 'NBANDS'])
        return os.path.abspath(dry_run)

    elif part == 'bs':
        # Band structure run
        bs_run = create_folder(folder, '01_t_bs')
        dry_run = os.path.join(folder, '00_t_dry')
        copy_files(dry_run, bs_run, ['KPOINTS', 'POSCAR', 'POTCAR', 'INCAR', 'CHGCAR'])
        KPAR, NBANDS = get_parallelization_variables(os.path.join(dry_run, 'vaspout.h5'))
        NBANDS = 4 * math.ceil((NBANDS*1.25)/4)  # increase nbands by 25%, round to multiple of 4
        INCAR_settings['NBANDS'] = NBANDS
        INCAR_settings['NPAR'] = 4
        INCAR_settings['KPAR'] = KPAR
        modify_incar(os.path.join(bs_run, 'INCAR'), add= bands_settings | INCAR_settings, remove=dry_settings)
        return os.path.abspath(bs_run)


# ===== Density of states workflow (non-SOC) =====
# Dry run sets finer kmesh, DOS range; second step does calculation

def t_dos(folder, part, previous=os.path.join(calculation_path, f't_conv_test/k{conv_kmesh[-1][0]}x{conv_kmesh[-1][1]}x{conv_kmesh[-1][2]}/e{max(conv_encut)}/03_t_scf')):
    dos_settings = {   
        'ICHARG': 11,
        'LWAVE': '.FALSE.',
        'LCHARG': '.FALSE.',
        'LORBIT': 11
    }

    if part == 'dry':
        # Dry run: generate KPOINTS mesh and set energy range from E-fermi
        dry_run = create_folder(folder, '00_t_dry')
        previous = os.path.join(calculation_path, 'step03', f'k{chosen_kmesh[0]}x{chosen_kmesh[1]}x{chosen_kmesh[2]}', f'e{chosen_encut}', '03_t_scf')
        copy_files(previous, dry_run, ['POSCAR', 'POTCAR', 'CHGCAR', 'INCAR'])
        with open(os.path.join(dry_run, 'KPOINTS'), 'w') as f:
                f.write(f'{dos_mesh_mult[0]*chosen_kmesh[0]}x{dos_mesh_mult[1]*chosen_kmesh[1]}x{dos_mesh_mult[2]*chosen_kmesh[2]}\n')
                f.write('0\n')
                f.write('Gamma\n')
                f.write(f'{dos_mesh_mult[0]*chosen_kmesh[0]} {dos_mesh_mult[1]*chosen_kmesh[1]} {dos_mesh_mult[2]*chosen_kmesh[2]}\n')
        # Determines the energy range for DOS calculation, currently e_fermi +- 6 eV, 1 meV resolution
        with h5py.File(os.path.join(previous, 'vaspout.h5'), "r") as f:
            e_fermi = round(f['results/electron_dos/efermi'][()], 3)
            dos_settings['EMIN'] = e_fermi - 6
            dos_settings['EMAX'] = e_fermi + 6
            dos_settings['NEDOS'] = 12000

        modify_incar(os.path.join(dry_run, 'INCAR'), add= dos_settings | dry_settings | INCAR_settings, remove=['KPAR', 'NPAR', 'NBANDS'])
        return os.path.abspath(dry_run)



    elif part == 'dos':
        # Density of states run 
        bs_run = create_folder(folder, '01_t_dos')
        dry_run = os.path.join(folder, '00_t_dry')
        copy_files(dry_run, bs_run, ['KPOINTS', 'POSCAR', 'POTCAR', 'INCAR', 'CHGCAR'])
        KPAR, NBANDS = get_parallelization_variables(os.path.join(dry_run, 'vaspout.h5'))
        NBANDS = 4 * math.ceil((NBANDS*1.25)/4)
        INCAR_settings['NBANDS'] = NBANDS
        INCAR_settings['NPAR'] = 4
        INCAR_settings['KPAR'] = KPAR
        modify_incar(os.path.join(bs_run, 'INCAR'), add= dos_settings | INCAR_settings, remove=dry_settings)
        return os.path.abspath(bs_run)


# ===== SCF with SOC =====
# Handles dry run and SCF run with LSORBIT enabled
def t_scf_so(folder, part, previous=os.path.join(calculation_path, f't_conv_test/k{conv_kmesh[-1][0]}x{conv_kmesh[-1][1]}x{conv_kmesh[-1][2]}/e{max(conv_encut)}/03_t_scf')):
    soc_settings = {
        'LSORBIT': '.TRUE.',   # enable spin-orbit coupling
        'ENCUT': chosen_encut  # use chosen energy cutoff
    }

    if part == 'dry':
        # Part 1: dry run — prepare folder and adjust INCAR for SOC
        dry_run = create_folder(folder, '00_t_dry_so')
        copy_files(previous, dry_run, ['POSCAR', 'POTCAR', 'CHGCAR', 'KPOINTS', 'INCAR'])
        modify_incar(os.path.join(dry_run, 'INCAR'), add= dry_settings | soc_settings | INCAR_settings)
        return os.path.abspath(dry_run)

    elif part == 'scf':
        # Part 2: SCF run with SOC, using parallelization from dry run
        scf_run = create_folder(folder, '01_t_scf_so')
        dry_run = os.path.join(folder, '00_t_dry_so')
        copy_files(dry_run, scf_run, ['KPOINTS', 'POSCAR', 'POTCAR', 'INCAR', 'CHGCAR'])
        KPAR, NBANDS = get_parallelization_variables(os.path.join(dry_run, 'vaspout.h5'))
        INCAR_settings['NBANDS'] = NBANDS
        INCAR_settings['NPAR'] = 4
        INCAR_settings['KPAR'] = KPAR
        modify_incar(os.path.join(scf_run, 'INCAR'), remove= dry_settings, add= INCAR_settings)
        return os.path.abspath(scf_run)


# ===== Geometry optimization with SOC =====
# Follows dry run → CG opt → Newton opt → SCF → report
def t_geo_so(folder, part, previous=os.path.join(calculation_path, 't_geo_so/01_t_scf_so')):
    opt_settings = {   
        'IBRION': 2,    # CG optimization
        'EDIFFG': -1E-3,
        'NSW': 50
    }
    if "ISIF" in globals:
        opt_settings["ISIF"] = ISIF
    else:
        opt_settings["ISIF"] = 3

    soc_settings = {
        'LSORBIT': '.TRUE.',
        'ENCUT': chosen_encut
    }

    if part == 'dry':
        # Dry run: copy files, adjust INCAR for SOC + optimization
        dry_run = create_folder(folder, '00_t_dry_so')
        copy_files(previous, dry_run, ['KPOINTS', 'POSCAR', 'POTCAR', 'WAVECAR', 'INCAR'])
        modify_incar(os.path.join(dry_run, 'INCAR'), add= opt_settings | dry_settings | soc_settings | INCAR_settings, remove=['KPAR', 'NPAR'])
        return os.path.abspath(dry_run)

    elif part == 'cg_opt':
        # Conjugate-gradient optimization with SOC
        geo1_run = create_folder(folder, '01_t_geo_so')
        dry_run = os.path.join(folder, '00_t_dry_so')
        copy_files(dry_run, geo1_run, ['KPOINTS', 'POSCAR', 'POTCAR', 'WAVECAR', 'INCAR'])
        KPAR, NBANDS = get_parallelization_variables(os.path.join(dry_run, 'vaspout.h5'))
        INCAR_settings['NBANDS'] = NBANDS
        INCAR_settings['NPAR'] = 4
        INCAR_settings['KPAR'] = KPAR
        modify_incar(os.path.join(geo1_run, 'INCAR'), add= soc_settings | INCAR_settings, remove=dry_settings)
        return os.path.abspath(geo1_run)

    elif part == 'nw_opt':
        # Newton optimization step with SOC
        geo2_run = create_folder(folder, '02_t_geo_so')
        geo1_run = os.path.join(folder, '01_t_geo_so')
        copy_files(geo1_run, geo2_run, ['KPOINTS', 'CONTCAR', 'POTCAR', 'INCAR', 'WAVECAR'])
        os.rename(os.path.join(geo2_run, 'CONTCAR'), os.path.join(geo2_run, 'POSCAR'))
        opt_settings['IBRION'] = 1  # Newton method
        opt_settings['NSW'] = 100
        modify_incar(os.path.join(geo2_run, 'INCAR'), add= opt_settings | soc_settings | INCAR_settings)
        return os.path.abspath(geo2_run)

    elif part == 'scf':
        # Final SCF after SOC geometry optimization
        scf_run = create_folder(folder, '03_t_scf_so')
        geo2_run = os.path.join(folder, '02_t_geo_so')
        copy_files(geo2_run, scf_run, ['KPOINTS', 'CONTCAR', 'POTCAR', 'INCAR', 'WAVECAR'])
        os.rename(os.path.join(scf_run, 'CONTCAR'), os.path.join(scf_run, 'POSCAR'))
        modify_incar(os.path.join(scf_run, 'INCAR'), remove= opt_settings, add=soc_settings | INCAR_settings)
        return os.path.abspath(scf_run)
        
    elif part == 'report':
        # Generate PDF report for SOC SCF
        scf_run = os.path.join(folder, '03_t_scf_so')
        create_report(os.path.join(scf_run, 'vaspout.h5'), os.path.join(scf_run, 'report.pdf'))


# ===== Band structure with SOC =====
# Dry run sets k-path, then band structure run
def t_bs_so(folder, part, previous=os.path.join(calculation_path, 't_geo_so/03_t_scf_so')):
    bands_settings = {   
        'ICHARG': 11,
        'LWAVE': '.FALSE.',
        'LCHARG': '.FALSE.',
        'LORBIT': 11,
        'LSORBIT': '.TRUE.',
        'ENCUT': chosen_encut
    }

    if part == 'dry':
        # Dry run: copy files, set k-path, adjust INCAR
        dry_run = create_folder(folder, '00_t_dry_so')
        copy_files(previous, dry_run, ['POSCAR', 'POTCAR', 'CHGCAR', 'INCAR'])
        create_kpath_standard(os.path.join(dry_run, 'KPOINTS'), k_points=k_path, points=min_points)
        modify_incar(os.path.join(dry_run, 'INCAR'), add= bands_settings | dry_settings | INCAR_settings, remove=['KPAR', 'NPAR', 'NBANDS'])
        return os.path.abspath(dry_run)

    elif part == 'bs':
        # Band structure run with SOC
        bs_run = create_folder(folder, '01_t_bs_so')
        dry_run = os.path.join(folder, '00_t_dry_so')
        copy_files(dry_run, bs_run, ['KPOINTS', 'POSCAR', 'POTCAR', 'INCAR', 'CHGCAR'])
        KPAR, NBANDS = get_parallelization_variables(os.path.join(dry_run, 'vaspout.h5'))
        NBANDS = 4 * math.ceil((NBANDS*1.25)/4)
        INCAR_settings['NBANDS'] = NBANDS
        INCAR_settings['NPAR'] = 4
        INCAR_settings['KPAR'] = KPAR
        modify_incar(os.path.join(bs_run, 'INCAR'), add= bands_settings | INCAR_settings, remove=dry_settings)
        return os.path.abspath(bs_run)


# ===== Density of states with SOC =====
# Dry run prepares kmesh + DOS range, second step runs calculation
def t_dos_so(folder, part, previous=os.path.join(calculation_path, 't_geo_so/03_t_scf_so')):
    dos_settings = {   
        'ICHARG': 11,
        'LWAVE': '.FALSE.',
        'LCHARG': '.FALSE.',
        'LORBIT': 11,
        'LSORBIT': '.TRUE.',
        'ENCUT': chosen_encut   
    }

    if part == 'dry':
        # Dry run: generate denser kmesh, set DOS range around E-fermi
        dry_run = create_folder(folder, '00_t_dry_so')
        copy_files(previous, dry_run, ['POSCAR', 'POTCAR', 'CHGCAR', 'INCAR'])
        with open(os.path.join(dry_run, 'KPOINTS'), 'w') as f:
                f.write(f'{dos_mesh_mult[0]*chosen_kmesh[0]}x{dos_mesh_mult[1]*chosen_kmesh[1]}x{dos_mesh_mult[2]*chosen_kmesh[2]}\n')
                f.write('0\n')
                f.write('Gamma\n')
                f.write(f'{2*chosen_kmesh[0]} {2*chosen_kmesh[1]} {2*chosen_kmesh[2]}\n')
        # Determines the energy range for DOS calculation, currently e_fermi +- 6 eV, 1 meV resolution
        with h5py.File(os.path.join(previous, 'vaspout.h5'), "r") as f:
            e_fermi = round(f['results/electron_dos/efermi'][()], 3)
            dos_settings['EMIN'] = e_fermi - 6
            dos_settings['EMAX'] = e_fermi + 6
            dos_settings['NEDOS'] = 12000

        modify_incar(os.path.join(dry_run, 'INCAR'), add= dos_settings | dry_settings | INCAR_settings, remove=['KPAR', 'NPAR', 'NBANDS'])
        return os.path.abspath(dry_run)

    elif part == 'dos':
        # DOS calculation run with SOC
        bs_run = create_folder(folder, '01_t_dos_so')
        dry_run = os.path.join(folder, '00_t_dry_so')
        copy_files(dry_run, bs_run, ['KPOINTS', 'POSCAR', 'POTCAR', 'INCAR', 'CHGCAR'])
        KPAR, NBANDS = get_parallelization_variables(os.path.join(dry_run, 'vaspout.h5'))
        NBANDS = 4 * math.ceil((NBANDS*1.25)/4)
        INCAR_settings['NBANDS'] = NBANDS
        INCAR_settings['NPAR'] = 4
        INCAR_settings['KPAR'] = KPAR
        modify_incar(os.path.join(bs_run, 'INCAR'), add= INCAR_settings | dos_settings, remove=dry_settings)
        return os.path.abspath(bs_run)


# ===== Helper functions =====
# Functions used inside steps for copying files, processing them, extracting information etc.

# Creates a POTCAR file based on the VASP recommended, if POTCAR wasn't supplied by the user
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


# Handles copying files and creates the missing files if possible
def copy_files(from_path, to_path, files):
    for file in files:
        src = os.path.join(from_path, file)
        dst = os.path.join(to_path, file)

        try:
            shutil.copyfile(src, dst)

        except FileNotFoundError:
            if file == 'POSCAR':
                raise FileNotFoundError("POSCAR file is missing. Cannot proceed.")

            elif file == 'POTCAR':
                #print("POTCAR not supplied -> Generating POTCAR from recommended files")
                try:
                    with open(os.path.join(from_path, 'POSCAR'), 'r') as f:
                        lines = f.readlines()
                        elements = lines[5].split()
                    create_potcar(elements, os.path.join(to_path, 'POTCAR'))
                except Exception as e:
                    raise RuntimeError(f"Failed to generate POTCAR: {e}")

            elif file == 'INCAR':
                #print("INCAR not supplied -> Generating from the settings")
                try:
                    incar_path = os.path.join(to_path, 'INCAR')
                    with open(incar_path, 'w') as f:
                        pass
                    modify_incar(incar_path, add=INCAR_settings)
                except Exception as e:
                    raise RuntimeError(f"Failed to generate INCAR: {e}")

            elif file == 'CHGCAR':
                raise FileNotFoundError("CHGCAR file is required for NSCF calculations but was not supplied.")

            elif file == 'WAVECAR':
                #print("WAVECAR file missing, calculations will be started from scratch")
                pass
            else:
                raise FileNotFoundError(f"Required file '{file}' is missing.")

        except PermissionError:
            raise PermissionError(f"Permission denied when copying {file}.")

        except Exception as e:
            raise RuntimeError(f"Unexpected error while handling {file}: {e}")

# Creating folders
def create_folder(path, folder_name):
    name = os.path.join(path, folder_name)
    try:
        os.mkdir(name)
    except:
        pass
    return name
      
# Modifying existing INCAR file based on the settings
# Uses dictionaries as input (can also be the list of tag names to remove)
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

# Extract parallelization variables (KPAR, NBANDS) from vaspout.h5
# Returns NBANDS divisible by 4 to set NPAR = 4 easily
def get_parallelization_variables(vaspout_file):
    with h5py.File(vaspout_file, "r") as f:
        eigenvalues = f['results/electron_eigenvalues/eigenvalues'][:]
        _ , nkpts, nbands = eigenvalues.shape
        nbands = nbands if nbands % 4 == 0 else ((nbands // 4) + 1) * 4
        kpar = max(i for i in range(1, 20) if nkpts % i == 0)
        return kpar, nbands
    
# Creates .pdf reports for SCF steps, that contain basic information about convergence
def create_report(vaspout, save_path):
    with h5py.File(vaspout, "r") as f:
        forces = f['intermediate/ion_dynamics/forces'][:]
        stress = f['intermediate/ion_dynamics/stress'][:]
        oszicar = f['intermediate/ion_dynamics/oszicar'][:]

        results = {
            "Stress in kB": f"[{stress[0,0,0]:.4f}, {stress[0, 1, 1]:.4f}, {stress[0, 2, 2]:.4f}]"
        }
        try:
            results['ENCUT'] = f['input/incar/ENCUT'][()]
        except:
            results['ENCUT'] = 'look in POTCAR'

    fig = plt.figure(figsize=(8.3, 11.7))

    ax1 = fig.add_axes([0.05, 0.55, 0.9, 0.4])
    ax1.axis('off')

    ax1.text(0.5, 1.05, "Simulation Report", ha="center", va="top", fontsize=16, weight="bold")

    y_pos = 0.95
    for key, value in results.items():
        ax1.text(0.05, y_pos, f"{key}:", fontsize=12, weight="bold")
        ax1.text(0.35, y_pos, f"{value}", fontsize=12)
        y_pos -= 0.07

    ax1.text(0.05, y_pos - 0.05, "Atomic Forces [eV/Å]:", fontsize=12, weight="bold")

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

    ax2 = fig.add_axes([0.15, 0.08, 0.7, 0.35])
    ax2.plot(oszicar[:, 0], oszicar[:, 2], marker='o', label='dE')
    ax2.plot(oszicar[:, 0], oszicar[:, 3], marker='o', label='d eps')
    ax2.set_yscale('symlog', linthresh=1e-4)
    ax2.set_xlabel("Step")
    ax2.set_ylabel("log(dE/d eps)")
    ax2.legend()
    ax2.grid(True)

    plt.savefig(save_path, format="pdf")
    plt.close()

# Report for convergence study for different kmeshes and encuts, gives the result in .csv format
def convergence_report(folder, kmesh_list, encut_list, pressure_threshhold=1, forces_threshhold=1e-3):
    data = {'a': {}, 'b': {}, 'c': {}, 'converged': {}}

    for kmesh in kmesh_list:
        kmesh_key = f'k{kmesh[0]}x{kmesh[1]}x{kmesh[2]}'
        data['a'][kmesh_key] = []
        data['b'][kmesh_key] = []
        data['c'][kmesh_key] = []
        data['converged'][kmesh_key] = []

        for encut in encut_list:
            filepath = os.path.join(folder, kmesh_key, f'e{encut}', '03_t_scf', 'vaspout.h5')
            print(filepath)
            with h5py.File(filepath, "r") as h5:
                lattice = h5['results/positions/lattice_vectors'][:]
                forces = h5['intermediate/ion_dynamics/forces'][:]
                stress = h5['intermediate/ion_dynamics/stress'][:]
                a = np.linalg.norm(lattice[0])
                b = np.linalg.norm(lattice[1])
                c = np.linalg.norm(lattice[2])
                converged = (abs(stress) < pressure_threshhold).all() and (abs(forces) < forces_threshhold).all()
                data['a'][kmesh_key].append((encut, a))
                data['b'][kmesh_key].append((encut, b))
                data['c'][kmesh_key].append((encut, c))
                data['converged'][kmesh_key].append((encut, converged))

    with open(os.path.join(folder, 'convergence_report_matrix.csv'), 'w') as f:
        for quantity in ['a', 'b', 'c', 'converged']:
            f.write(f'{quantity}')
            for kmesh in kmesh_list:
                kmesh_key = f'k{kmesh[0]}x{kmesh[1]}x{kmesh[2]}'
                f.write(f', {kmesh_key}')
            f.write('\n')
            for i, encut in enumerate(encut_list):
                f.write(f'e{encut}')
                for kmesh in kmesh_list:
                    kmesh_key = f'k{kmesh[0]}x{kmesh[1]}x{kmesh[2]}'
                    value = data[quantity][kmesh_key][i][1]
                    if quantity == 'converged':
                        f.write(f', {str(value)}')
                    else:
                        f.write(f', {value:.5f}')
                f.write('\n')
            f.write('\n')

# Creates a standard kpath based on the list of high symmetry points with equally spaced points between them
def create_kpath_standard(file, k_points, points):
    with open(file, 'w') as f:
        f.write('bs kpath\n')
        f.write(f' {points}\n')
        f.write('line\n')
        f.write('reciprocal\n')
        for i in range(len(k_points)-1):
            p1 = k_points[i]
            p2 = k_points[i+1]
            f.write(f'  {p1[0]:.5f}  {p1[1]:.5f}  {p1[2]:.5f}    1\n')
            f.write(f'  {p2[0]:.5f}  {p2[1]:.5f}  {p2[2]:.5f}    1\n\n')

# BLANK
def create_job():
    pass

# BLANK
def run_vasp():
    pass

# B:LANK
def check_errors():
    pass

# Exctracts NKPTS and NBANDS from OUTCAR - for use after dry run
def get_NKPTS_NBANDS(file):
    nkpts, nbands = None, None
    with open(file, "r") as f:
        for line in f:
            if "NKPTS" in line and "NBANDS" in line: 
                parts = line.split()
                for i, p in enumerate(parts):
                    if p == "NKPTS":
                        nkpts = int(parts[i+2]) 
                    if "NBANDS" in p:
                        if p == 'NBANDS':
                            nbands = int(parts[i+2])
                        elif p =='NBANDS=':
                            nbands = int(parts[i+1])
                break
    return nkpts, nbands

# Function I used for testing something, maybe usefull later
def find_last_step(path):
    folders = os.listdir(path)
    folder_dict = {f:int(f.split('_')[0]) for f in folders}
    return max(folder_dict, key=folder_dict.get)


# ===== Main working part of the script ======
# Here a corresponding step is launched based on the inputs from console by specifying
# --step (eg. t_scf, t_geo_so ...) and the corresponding --part (eg. dry, bs, cg_opt) 

parser = argparse.ArgumentParser(description="Input agruments for the script")
parser.add_argument("--step", type=str, help="Step")
parser.add_argument("--part", type=str, help="Part of the step")
args = parser.parse_args()

path = create_step(step=args.step, part=args.part)

