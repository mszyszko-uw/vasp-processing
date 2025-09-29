from defaults import *
# File containing the settings for varius steps of the calculations
# if the variable is not set, the default setting will be used, but this can lead to bad results

# STEPS parameter is the only on the always has to be specified in this file, the rest is optional
STEPS = {
    't_scf': '',    
    't_geo': 't_scf',
    't_bs': 't_geo'
}

INCAR_SETTINGS['SYSTEM'] = 'bulk MoTe2'
INCAR_SETTINGS['ISMEAR'] = 0
INCAR_SETTINGS['SIGMA'] = 0.05
INCAR_SETTINGS['IVDW'] = 11
INCAR_SETTINGS['LMAXMIX'] = 4
INCAR_SETTINGS['LASPH'] = '.TRUE.'


