from defaults import *
# File containing the settings for varius steps of the calculations
# if the variable is not set, the default setting will be used, but this can lead to bad results

# STEPS parameter is the only on the always has to be specified in this file, the rest is optional
STEPS = {
    't_scf': '',    
    't_geo': 't_scf',
    't_conv_test': 't_geo',
    't_bs': 'converged',
    't_dos': 'converged',
    't_scf_so': 'converged',
    't_geo_so': 't_scf_so',
    't_bs_so': 't_geo_so',
    't_dos_so': 't_geo_so'
}

INCAR_SETTINGS['SYSTEM'] = 'bulk MoTe2'
INCAR_SETTINGS['ISMEAR'] = 0
INCAR_SETTINGS['SIGMA'] = 0.05
INCAR_SETTINGS['IVDW'] = 11
INCAR_SETTINGS['LMAXMIX'] = 4
INCAR_SETTINGS['LASPH'] = '.TRUE.'

CONV_KMESH = [(6, 6, 2), (7, 7, 2), (8, 8, 3)] 
CONV_ENCUT = [250, 300, 350, 400]
K_PATH = [
    (0.00000,  0.00000,  0.00000),
    (0.50000,  0.00000,  0.00000),  
    (0.33333,  0.33333,  0.00000),
    (0.00000,  0.00000,  0.00000),
    (0.00000,  0.00000,  0.50000),
    (0.50000,  0.00000,  0.50000),
    (0.33333,  0.33333,  0.50000),
    (0.00000,  0.00000,  0.50000)
]
