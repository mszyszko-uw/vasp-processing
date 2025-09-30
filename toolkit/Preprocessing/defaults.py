# File containing default settings
#if the this setting are not specified in the input.py file, this settings will be used


# ===== Set up variables ======
PSEUDOPOTENTIALS_PATH = 'your_path' # path to the folder containing all VASP pseudopotentials
LOGFILE = 'toolkit.log' # path + name of your log file for toolkit
CALCULATION_PATH = 'Calculations' # where the calculation folders will be created, if toolkit used from terminal
CPU_PER_NODE = 192
CPU_PER_SOCKET = 24



# ===== INCAR variables ======

# Main INCAR variables added to every INCAR file created by toolbox
INCAR_SETTINGS = {
    'EDIFF': 1e-6
}

# INCAR variable responsible for dry runs
DRY_SETTINGS = {
    'ALGO': None, 
    'NELM': 1, 
    'LWAVE': '.FALSE.', 
    'LCHARG': '.FALSE.'}

# INCAR variable responsible for steps that include SOC
SOC_SETTINGS = {
    'LSORBIT': '.TRUE.'
}

# INCAR variable responsible for conjugate-gradient step during structure relaxation
CG_OPT_SETTINGS = {
    'IBRION': 2,   
    'EDIFFG': -1E-3,
    'NSW': 50,
}

# INCAR variable responsible for quasi-Newton algorithm step during structure relaxation
NW_OPT_SETTINGS = {
    'IBRION': 1, 
    'EDIFFG': -1E-3,
    'NSW': 100,
}

# INCAR variable responsible for band structure calculations 
BS_SETTINGS = {   
    'ICHARG': 11,
    'LWAVE': '.FALSE.',
    'LCHARG': '.FALSE.',
    'LORBIT': 11
}

# INCAR variable responsible for density of states calculations
DOS_SETTINGS = {   
    'ICHARG': 11,
    'LWAVE': '.FALSE.',
    'LCHARG': '.FALSE.',
    'LORBIT': 11,
    'NEDOS': 10001
}


# ===== Miscellaneous variables ======
# Variables responsible for different parts of the toolkit functionalities

INCLUDE_DRY_RUN = True # Determines if the dry run is included at every step
ISIF = 3 # Determines the relaxation degrees of freedom - see VASP wiki
DOS_RANGE = 5 # Determines the energy range for DOS calculations -> E_FERMI +- DOS_RANGE eV
KSPACING = 0.5 # Determines KSPACING variable in INCAR in case that KPOINTS file isn't present for k-mesh calculations
DOS_KSPACING = KSPACING / 2 # Determines KSPACING variable in INCAR for DOS calculations in case that KPOINTS file isn't present for k-grid calculations
CONV_KMESH = [None] # Determines the different k-meshes for convergence tests
CONV_ENCUT = [None] # Determines the different ENCUT values for convergence tests
K_PATH = [None] # Determines the high symmetry points for band structure calculations
POINTS_PER_SEGMENT = 10 # Determines the number of points between each high symmetry point on the path
BS_BANDS_MULTIPLIER = 1.25 # The multiplier for the number of bands for band structure calculations
DOS_BANDS_MULTIPLIER = 1.25 # The multiplier for the number of bands for DOS calculations

