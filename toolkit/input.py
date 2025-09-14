from defaults import *
# INPUT FILE FOR TOOLBOX CONTAINING IMPORTANT SETTINGS FOR VASP CALCULATIONS 
# SUCH AS INCAR SETTINGS or K PATH FOR BAND STRUCTURE CALCULATIONS

# if the variable is not set, the default setting will be used, but this can lead to bad results


INCAR_SETTINGS['SYSTEM'] = 'bulk MoTe2'
INCAR_SETTINGS['ISMEAR'] = 0
INCAR_SETTINGS['SIGMA'] = 0.05
INCAR_SETTINGS['IVDW'] = 11
INCAR_SETTINGS['LMAXMIX'] = 4
INCAR_SETTINGS['LASPH'] = '.TRUE.'


STEPS = {
    't_scf': '',    
    't_geo': 't_scf',
    't_bs': 't_geo'
}
