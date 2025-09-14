# File containing default settings

INPUT_PATH = ''
LOGFILE = 'toolkit.log'
PSEUDOPOTENTIALS_PATH = '/home/maciej/Documents/Umowa Tomek Woźniak/toolbox/potpaw_PBE'
CALCULATION_PATH = 'Calculations'

INCAR_SETTINGS = {
    'EDIFF': 1e-6
}

DRY_SETTINGS = {
    'ALGO': None, 
    'NELM': 1, 
    'LWAVE': '.FALSE.', 
    'LCHARG': '.FALSE.'}

SOC_SETTINGS = {
    'LSORBIT': '.TRUE.'
}

CG_OPT_SETTINGS = {
    'IBRION': 2,   # conjugate gradient
    'EDIFFG': -1E-3,
    'NSW': 50,
}

NW_OPT_SETTINGS = {
    'IBRION': 1,   # quasi-Newton
    'EDIFFG': -1E-3,
    'NSW': 100,
}

BS_SETTINGS = {   
    'ICHARG': 11,
    'LWAVE': '.FALSE.',
    'LCHARG': '.FALSE.',
    'LORBIT': 11
}

DOS_SETTINGS = {   
    'ICHARG': 11,
    'LWAVE': '.FALSE.',
    'LCHARG': '.FALSE.',
    'LORBIT': 11,
    'NEDOS': 10001
}

ISIF = 3


DOS_RANGE = 5 # E_FERMI +- DOS_RANGE eV
KSPACING = 0.5
DOS_KSPACING = KSPACING / 2

CONV_KMESH = [None]
CONV_ENCUT = [None]
K_PATH = [
    (0.00000,  0.00000,  0.00000),
    (0.50000,  0.00000,  0.00000),  
    (0.33333,  0.33333,  0.00000),
    (0.00000,  0.00000,  0.00000)
]
KPOINTS = 10

CPU_PER_NODE = 192
CPU_PER_SOCKET = 24
BS_BANDS_MULTIPLIER = 1.25
DOS_BANDS_MULTIPLIER = 1.25
