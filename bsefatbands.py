import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os
import stat

##### Autor: Karol Gałązka
## Cel: robienie konkretnych wykresów z plików BSEFATBAND

l_exc  = [31,32]  # List of excitons' indexes to be summed
l_band = [50,51,52,53,54]     # List of bands to draw
n_kz   = 4        # number of k_z points
quart  = 1    # number of quarter (1,2,3 or 4) or 'all'
folder      = '/home/karol/Monolayers/T_Wozniak/CrSBr/'           # PATH to the folder where BSEFATBAND is
plot_folder = '/home/karol/Monolayers/T_Wozniak/CrSBr/BSE_plots/' # Here the plot will be saved


# figsize=(6,6)   # Change the size of the plot, it will be square so it is probably best to keep the dimentions equal
# The function is executed at the end of the script
def ExcitonFatband(l_exc, i_band, i_kz, ax=None, folder = './'):
    for tries in range(2):
        try:
            totmat = []
            print('Looking for small BSE files')
            for i_exc in l_exc:
                with open(folder+f'BSE_{i_exc}_{i_band}_{i_kz}', 'r+') as bsefile:
                    if len(bsefile.readline().split()) > 3:
                        print(f'Cropping BSE_{i_exc}_{i_band}_{i_kz}')
                        bsefile.seek(0)
                        bsemat = bsefile.read()
                        bsemat = bsemat.replace('+i*', ' ', )
                        bsemat = np.fromstring(bsemat, dtype=float, sep=' ').reshape((-1,10))
                        bsemat = bsemat[:,[0,1,5]]
                        bsemat[:len(bsemat)//2,2] += bsemat[len(bsemat)//2:,2]
                        bsemat = bsemat[:len(bsemat)//2,:]
                        kpoint = bsemat[0,0:2]
                        repstep = 0
                        for row in bsemat:
                            if (row[:2] == kpoint).all(): repstep += 1
                            else: break
                        summat = bsemat[0::repstep, :]
                        for rep in range(1,repstep):
                            summat[:,2] += bsemat[rep::repstep,2]
                        bsemat = summat
                        bsefile.seek(0)
                        np.savetxt(bsefile, bsemat, fmt='%.5f \t%.5f \t%.7f')
                        bsefile.truncate()
                    if len(bsefile.readline().split()) == 3:
                        bsefile.seek(0)
                        bsemat = np.fromfile(bsefile, dtype=float, sep=' ').reshape((-1,3))
                if len(totmat) == 0: totmat = bsemat
                else: totmat[:,2] += bsemat[:,2]
                print(f'Sum bse = {np.sum(bsemat[:,2])}')

            print('BSE info retrieved, plotting')
            plotname = 'BSE_'+'-'.join(str(i_exc) for i_exc in l_exc) + f'_{i_band}_{i_kz}.png'
            if  ax!=None:
                ax.set_title(plotname[:-4])
                ax.scatter(totmat[:,0],totmat[:,1],totmat[:,2])
            return totmat

        except FileNotFoundError:
            print('No small BSE files found, opening BSEFATBAND')
            with open(folder+f'BSEFATBAND') as mainfile:
                rank = int(mainfile.readline().split()[0])
                mainfile.readline()
                kzdict = {}
                kznow  = ''
                ikznow = 0
                print('Making k_z values dictionary')
                for line in mainfile:
                    if 'BSE' in line: break
                    kznow = line.split()[2]
                    if kznow not in kzdict.values():
                        ikznow += 1
                        kzdict[ikznow] = kznow
                print(kzdict)
            print('Writing "processing.sh" script')
            with open(folder+'processing.sh', 'w') as processfile:
                processfile.write(f'#!/bin/bash\ncd {folder}\n')
                for i_exc in l_exc:
                    start = 2+(1+rank)*(i_exc-1)+1
                    end   = start + rank - 1
                    processfile.write(f'sed -n "{start},{end}p;{end}q" BSEFATBAND > BSE_{i_exc}\n')
                    processfile.write(f'grep ".\{{18,\}}{kzdict[i_kz]}.\{{46,\}}{i_band}    " ' + \
                                    f'BSE_{i_exc} > BSE_{i_exc}_{i_band}_{i_kz}\n')
                for i_exc in l_exc: processfile.write(f'rm BSE_{i_exc}\n')

            st = os.stat(folder+'processing.sh')
            os.chmod(folder+'processing.sh', st.st_mode | stat.S_IEXEC)
            print('Running "processing.sh" script')
            subprocess.call(f'{folder}processing.sh')

def FatbandGrid(l_exc, l_band, n_kz, quart = 1, folder = './', plot_folder = './', figsize=(10,10)):
    fig = plt.figure(figsize=(2*len(l_band),8))
    gs = fig.add_gridspec(n_kz, len(l_band), 
                          width_ratios = [1]*len(l_band), height_ratios=[1]*n_kz,
                          left=0.05, right=0.98, bottom=0.05, top=0.95,
                          wspace=0.02, hspace=0.2)
    axs = gs.subplots()
    for i_kz, axr in enumerate(axs):
        for i_band, ax in enumerate(axr):
            ax.set_box_aspect(1)
            ExcitonFatband(l_exc, l_band[i_band], i_kz+1, ax, folder=folder)
            if quart != 'all':
                ax.set_xlim((-1)**(quart%3!=1)*np.array([-0.05, 0.55]))
                ax.set_ylim((-1)**(quart>2)*np.array([-0.05, 0.55]))
            if not i_kz+1 == n_kz: 
                ax.set_xticks([])
            if not i_band == 0:
                ax.set_yticks([])
    plt.savefig(plot_folder+f"FatbandGrid{'-'.join(str(x) for x in l_exc)}_{'-'.join(str(x) for x in l_band)}_{quart}")

if __name__ == '__main__':
    # FatbandGrid(l_exc, l_band, n_kz, quart=quart, folder=folder, plot_folder=plot_folder)
    ExcitonFatband([15,16], 54, 1, folder=folder)
