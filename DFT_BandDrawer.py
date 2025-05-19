from DFT_Output import DFT_Output

##### Autor: Karol Gałązka
## Cel: przykładowe zastosowanie klasy DFT_Output

Experiment = DFT_Output('VASP')
Experiment.read_EIGENVAL()
Experiment.band_Plot()
