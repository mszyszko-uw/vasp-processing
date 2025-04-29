from DFT_Output import DFT_Output

Experiment = DFT_Output('VASP')
Experiment.read_EIGENVAL()
Experiment.band_Plot()