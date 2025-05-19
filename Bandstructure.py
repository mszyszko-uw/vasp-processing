import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from matplotlib import colormaps
from matplotlib.lines import Line2D

##### Autor: Karol Gałązka
## Cel: Wyrysowanie struktury pasmowej o zadanych parametrach
## na podstawie pliku PROCAR

## Powielone czytanie PROCAR, przestarzałe
def read_procar(PROCAR):
    # try: 
        # with open(PROCAR+'.npy', 'rb') as f:
        #     promat = np.load(f, allow_pickle=True)
        #     energs = np.load(f, allow_pickle=True)
        #     E_max  = np.load(f, allow_pickle=True)
        #     E_min  = np.load(f, allow_pickle=True)
        #     E_gap  = np.load(f, allow_pickle=True)
        #     kindexer = np.load(f, allow_pickle=True)
        #     print('Read from PROCAR.npy')
    # except(FileNotFoundError):
        print(f"Checking {PROCAR}")
        nlines = 0
        nprjct = 0
        with open(PROCAR) as procar:
            ion_flag = False
            while True:
                data = procar.readline().split()
                nlines += 1
                if len(data) > 0 :
                    if data[0] == 'ion' and ion_flag == False:
                        norbitals = len(data)-2
                        ion_flag = True
                    elif data[0] == 'tot' and ion_flag == True:
                        nprjct += 1
                    elif data[0] == 'ion' and ion_flag == True:
                        break
            for line in procar:
                if "# of k-points:" in line:
                    nprjct += 1
                nlines += 1

        print(f"Reading {PROCAR}", end=' ')
        with open(PROCAR) as procar:
            procar.readline()
            data = procar.readline().split()
            nkpts, nbands, nions = int(data[3]), int(data[7]), int(data[11])
            kindexer = {}
            proj = 0
            promat = np.zeros((nkpts, nbands, nprjct, nions, norbitals), dtype=float)
            energs = np.zeros((nkpts, nbands), dtype=float)
            vbmax = 0
            for line in tqdm(procar, total=nlines - ((nions+1)*nprjct)*nbands*nkpts - 2):
                if "# of k-points:" in line:
                    proj += 1
                data = line.split()
                if len(data) > 0:
                    if data[0] == 'k-point': 
                        kpoint = int(data[1])-1
                        kindexer[tuple(float(f'{float(num):.4f}') for num in data[3:6])] = kpoint
                    if data[0] == 'band':
                        band = int(data[1])-1
                        energs[kpoint, band] = float(data[4])
                        if vbmax==0 and int(float(data[7]))==0:
                            vbmax = band
                    if data[0] == 'ion':
                        if nprjct != 2:
                            for proj in range(nprjct):
                                for ion in range(nions):
                                    mrow = np.fromstring(procar.readline(),
                                                        dtype=float, sep=' ')[1:-1]
                                    promat[kpoint, band, proj, ion, :] = mrow
                                procar.readline()
                        else:
                            for ion in range(nions):
                                mrow = np.fromstring(procar.readline(),
                                                    dtype=float, sep=' ')[1:-1]
                                promat[kpoint, band, proj, ion, :] = mrow
                            procar.readline()
        print("PROCAR read, promat shape is:", promat.shape, f'vbmax is {vbmax}')
        E_max = np.max(energs[:,vbmax-1])
        E_min = np.min(energs[:,vbmax])
        E_gap = E_min - E_max
        # with open(PROCAR+'.npy', 'wb') as f:
        #     np.save(f, promat, allow_pickle=True)
        #     np.save(f, energs, allow_pickle=True)
        #     np.save(f, E_max, allow_pickle=True)
        #     np.save(f, E_min, allow_pickle=True)
        #     np.save(f, E_gap, allow_pickle=True)
        #     np.save(f, kindexer, allow_pickle=True)
        #     print('Saved to PROCAR.npy')
        return promat, energs, E_max, E_min, E_gap, kindexer


## Rysowanie struktury pasmowej na pojedynczym wykresie
def band_plot(promat, energs, ax,
              kpoints = [], bands = [], projections = [], ions = [], orbits = [], 
              fatbands = True, bandlines = False, sum_orbits = False, sum_ions = True, 
              label = '', bandcolor='black', interpolate = False):
    # x_line, x_ticks = make_xline(klines, kline_nums)
    ebands = np.array([]).reshape(0, energs.shape[1])
    pbands = np.array([]).reshape(0, *promat.shape[1:])
    orbit_names = ['$s$', '$p_y$', '$p_z$', '$p_x$', '$d_{xy}$', 
              '$d_{yz}$', '$d_{z^2}$', '$d_{xz}$', '$d_{x^2-y^2}$']

    if orbits == []: orbits = np.arange(9)
    elif orbits == 's': orbits = [1]
    elif orbits == 'p': orbits = [2,3,4]
    elif orbits == 'd': orbits = [5,6,7,8,9]
    else: orbits = np.array(orbits)-1
    if ions == []: ions = np.arange(promat.shape[3])+1
    if projections == []: projections = np.arange(promat.shape[2])+1
    cmap = colormaps['Set1']
    colors = cmap(np.linspace(0, 1, 9))
    x_line = np.linspace(-0.5,0.5,len(kpoints))
    # if kline_nums == []: kline_nums = np.arange(len(klines))+1
    # if kline_nums[0] > 0:
    #     xlabels = [klines[kline_nums[0]-1].a_name]
    # if kline_nums[0] < 0:
    #     xlabels = [klines[-kline_nums[0]-1].b_name]
    # for num in kline_nums:
    for kpoint in kpoints:
        # si = klines[abs(num)-1].start_index
        # ei = klines[abs(num)-1].end_index
        # if num > 0:
        ebands = np.concatenate([ebands, energs[kpoint:kpoint+1]], axis = 0)
        pbands = np.concatenate([pbands, promat[kpoint:kpoint+1]], axis = 0)
            # xlabels.append(klines[num-1].b_name)
        # if num < 0:
        #     ebands = np.concatenate([ebands, energs[ei:(si-1 if si-1 >= 0 else None):-1]], axis = 0)
        #     pbands = np.concatenate([pbands, promat[ei:(si-1 if si-1 >= 0 else None):-1]], axis = 0)
        #     xlabels.append(klines[-num-1].a_name)
    if len(bands) == 0:
        if bandlines: ax.plot(x_line, ebands, color=bandcolor, linewidth = 0.5)
    elif bandlines: ax.plot(x_line, ebands[:,bands], color=bandcolor, linewidth = 0.5)
    if sum_ions:
        for ion in ions[1:]:
            pbands[:,:,:,ions[0]-1,:] += pbands[:,:,:,ion-1,:]
        ions = [ions[0]]
    if fatbands:
        legend_lines = []
        legend_labels= []
        for ion in ions:
            for proj in projections:
                if sum_orbits:
                    ax.scatter(np.stack([x_line]*energs.shape[1], axis=1), 
                                ebands, np.sum(np.abs([100*pbands[:,:,proj-1,ion-1,j] for j in orbits]), axis=0), 
                                label = label, color = colors[1])
                else:
                    for j in orbits:
                        legend_lines.append(Line2D([0], [0], color=colors[j], lw=4, label = orbit_names[j]))
                        legend_labels.append(orbit_names[j])
                        if interpolate:
                            for band in bands:
                                x_int = np.linspace(x_line[0], x_line[-1], interpolate)
                                e_int = np.interp(x_int, x_line, ebands[:,band])
                                p_int = np.interp(x_int, x_line, 100*pbands[:,band,proj-1,ion-1,j])
                                ax.scatter(x_int, e_int, p_int, 
                                        #    label = orbit_names[j], 
                                           color = colors[j], alpha=0.08)
                        else:
                            if len(bands) == 0:
                                ax.scatter(np.stack([x_line]*energs.shape[1], axis=1), 
                                        ebands, np.abs(100*pbands[:,:,proj-1,ion-1,j]), 
                                        # label = orbit_names[j], 
                                        color = colors[j])
                            else:
                                ax.scatter(np.stack([x_line]*len(bands), axis=1), 
                                        ebands[:,bands], 100*pbands[:,bands,proj-1,ion-1,j], 
                                        # label = orbit_names[j], 
                                        color = colors[j])
        ax.legend(legend_lines, legend_labels)
    # for i in x_ticks[1:-1]:
    #     ax.axvline(i, color="black", linestyle=":")
    ax.axvline(0, color="black", linestyle="--", linewidth = 0.5)
    ax.axhline(0, color="black", linestyle="--", linewidth = 0.5)
    # ax.set_xlim([0, x_line[-1]])
    # ax.set_xticks(x_ticks)
    # ax.set_xticklabels(xlabels, fontsize=18)

# Trzeba podać
folder = '/home/karol/Monolayers/T_Wozniak/BND_AFM/'
title = 'figure'

PROCAR  = folder + 'PROCAR'
promat, energs, E_max, E_min, E_gap, kindexer = read_procar(PROCAR)
print(E_max)
# promat := (nkpts, nbands, nproj, nions, norbitals)
# promat[:,:,0,:,:] += promat[:,:,1,:,:]
energs -= E_max

klist = []
for pos in kindexer.keys():
    if pos[0] == 0 and pos[2] == 0 and pos[1] >= 0 and pos[1] <= 0.5:
        klist.append(kindexer[pos])
klist = klist[-1::-1]
print(klist)
for pos in kindexer.keys():
    if pos[1] == 0 and pos[2] == 0 and pos[0] > 0 and pos[0] <= 0.5:
        klist.append(kindexer[pos])

# Tutaj było wyciąganie odpowiednich k-punktów dla tego pliku ze ścieżką
# klist = []
# klist.append(kindexer[(0.0,0.5,0.0)])
# for pos in list(kindexer.keys())[189:]:
#     if pos[0] == 0 and pos[2] == 0 and pos[1] >= 0 and pos[1] <= 0.5:
#         if kindexer[pos] >= 189:
#             klist.append(kindexer[pos])
# klist.append(kindexer[(0.0,0.0,0.0)])
# for pos in list(kindexer.keys())[189:]:
#     if pos[1] == 0 and pos[2] == 0 and pos[0] > 0 and pos[0] <= 0.5:
#         if kindexer[pos] >= 189:
#             klist.append(kindexer[pos])
# klist.append(kindexer[(0.5,0.0,0.0)])
# for i in [202, 203, 204, 205, 206]:
#     klist.remove(i)

fig = plt.figure(figsize=(5,10))
ax = fig.add_subplot(1,1,1)
band_plot(promat, energs, ax,
          kpoints = klist, 
          bands = np.arange(0,119)+1, # można wybrać mniej, żeby się sprawniej liczyło
          projections = [1], # Faktycznie był wybrany spin dół chyba, góra powinna być z indeksem 0
          ions = [1,2,3,4],  # Wybór jonów, tutaj Chrom 
          orbits = [1,2,3,4,5,6,7,8,9], # Orbitale: S-1, P-2,3,4, D-5,6,7,8,9
          interpolate=1000)
ax.set_ylim([-.7, 3.5])
plt.savefig(folder + title)
