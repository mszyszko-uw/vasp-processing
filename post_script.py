import matplotlib.pyplot as plt
import numpy as np
from vaspout_h5 import vaspout_h5, plot_BS_DOS

BS = vaspout_h5('T_Wozniak/h5_post/so_bs/vaspout.h5')
DOS = vaspout_h5('T_Wozniak/h5_post/so_dos/vaspout.h5')

plot_BS_DOS(BS, DOS, 'Te', bands=np.arange(41,57))
plt.ylim([-2,3])
plt.show()

BS.plot_BS('Te','-1 4 -7', gnuplot=True)
plt.ylim([-6.5,9])
plt.show()

DOS.plot_DOS('2Mo d')
DOS.plot_DOS('3Te px py')
plt.legend()
plt.show()

fig, axes = plt.subplots(ncols=3, sharey=True)
BS = vaspout_h5("Dresden/CrSBr/80BANDS/vaspout.h5")
BS.plot_BS('Cr d down', '1 2 3 4', E0 = 0, ax=axes[0])
BS.plot_BS('S p down', '1 2 3 4', E0 = 0, ax=axes[1])
BS.plot_BS('Br p down', '1 2 3 4', E0 = 0, ax=axes[2])
plt.show()

fig, axes = plt.subplots(ncols=2, sharey=True, gridspec_kw={'width_ratios': [1, 1]})
BS.plot_BS('Cr up d', E0=0, ax=axes[0])
BS.plot_BS('Cr down d', E0=0, ax=axes[1])
plt.show()

BS = vaspout_h5("dCrSBr/50percent/1x1-3/AAFM/band_vaspout.h5")
BS.plot_BS('pure', '1 2 3 4', E0 = 0)
plt.show()