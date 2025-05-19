from DFT_Output import DFT_Output
import matplotlib.pyplot as plt

##### Autor: Karol Gałązka
## Cel: przykładowe zastosowania klasy DFT_Output

# Experiment = DFT_Output('VASP')
# Experiment.read_WAVEDERF()
# print( Experiment.Optical_Selection_Rules(1, 52, save = 'VASP_OSR') )
# Experiment.L_Plot(1, 3, 52, plot = 'show', save = 'VASP_L')


Experiment = DFT_Output('VASP')
# Experiment.DOS_Plot()

# EIGENVAL = '/home/karol/Monolayers/VASP/Silicon/EIGENVAL'
# Experiment.set_VASP_files(EIGENVAL=EIGENVAL)

# Experiment.read_EIGENVAL()
# Experiment.band_Plot()

# Experiment.read_WAVEDER()
print( Experiment.Optical_Selection_Rules(1, 26, 2, 2, save = 'VASP_OSR', txt = True) )
print( Experiment.L_Plot(1, 3, 26, 2, 2, plot = 'save', save = False, txt = False)  )

###############

# Experiment = DFT_Output("MoS2.mommat2up", "MoS2")
# # Experiment.read_mommat2up()
# print( Experiment.Optical_Selection_Rules(1, 26, save = 'Wien_OSR', txt = True) )
# Experiment.L_Plot(1, 3, 26, plot = 'show', save = 'Wien_L')

# Experiment.read_P_matrix("MoS2_P_matrix.npy")
# print( Experiment.Optical_Selection_Rules(1, 26, save = 'Wien_OSR') )
# Experiment.L_Plot(1, 3, 26, plot = 'show', save = 'Wien_L')
