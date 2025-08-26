# INPUT FILE FOR TOOLBOX CONTAINING IMPORTANT SETTINGS FOR VASP CALCULATIONS 
# SUCH AS INCAR SETTINGS or K PATH FOR BAND STRUCTURE CALCULATIONS
#
# the names of the variables should not be changed
# if the variable is not set, the default setting will be used, but this can lead to bad results


pseudopotentials_path = '/net/scratch/hscra/plgrid/plgmszyszko/potpaw_PBE'
calculation_path = 'Calculations'
INCAR_settings = {
    'SYSTEM': 'bulk MoTe2',
    'ISMEAR': 0,
    'SIGMA': 0.05,
    'IVDW': 11,
    'LMAXMIX': 4,
    'LASPH': '.TRUE.',
    'EDIFF': 1e-6
}

ISIF = 2

dry_settings = {
    'ALGO': None, 
    'NELM': 1, 
    'LWAVE': '.FALSE.', 
    'LCHARG': '.FALSE.'}

conv_kmesh = [(6, 6, 2), (7, 7, 2)]
conv_encut = [250, 300, 350]

chosen_kmesh = (7, 7, 2)
chosen_encut = 400
k_path = [
    (0.00000,  0.00000,  0.00000),
    (0.50000,  0.00000,  0.00000),  
    (0.33333,  0.33333,  0.00000),
    (0.00000,  0.00000,  0.00000),
    (0.00000,  0.00000,  0.50000),
    (0.50000,  0.00000,  0.50000),
    (0.33333,  0.33333,  0.50000),
    (0.00000,  0.00000,  0.50000)
]
min_points = 30

dos_mesh_mult = (2, 2, 1)

