from vaspout_h5 import vaspout_h5

"""
The way to simply use the postprocessing at the end of calculations
"""

vaspout_h5('T_Wozniak/h5_post/so_bs/vaspout.h5').plot_BS_default()

vaspout_h5('T_Wozniak/h5_post/so_dos/vaspout.h5').plot_DOS_default()