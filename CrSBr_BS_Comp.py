import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from matplotlib import colormaps


class kline():
    def __init__(self, a_pt, b_pt, Npts, start_index, a_name = '', b_name = ''):
        self.a_pt = a_pt
        self.b_pt = b_pt
        self.a_name = a_name
        self.b_name = b_name
        self.Npts = Npts
        self.start_index = start_index
        self.end_index   = start_index + Npts-1
        self.name = self.a_name + ' -> ' + self.b_name

    def __repr__(self):
        return self.name

    def __abs__(self):
        return np.linalg.norm(self.b_pt - self.a_pt)
    

def read_doscar(DOSCAR):
    print(f"Reading {DOSCAR}", end=' ')
    with open(DOSCAR) as dos:
        nions = np.fromstring(dos.readline(), dtype=int, sep=' ')[0]
        for i in range(4):
            dos.readline()
        head = dos.readline()
        data = np.fromstring(head, 
                            dtype=float, 
                            sep = ' ')
        NEDOS, efermi = int(data[2]), data[3]

        first = np.fromstring(dos.readline(), dtype=float, sep=' ')
        totmat = np.zeros((NEDOS, len(first)))
        totmat[0] = first
        totmat[1:,:] = np.fromfile(dos,
                        count = len(first)*(NEDOS-1),
                        sep = ' ', dtype=float).reshape(-1,len(first))

        for ion in range(nions):
            head = dos.readline()
            if ion == 0:
                first = np.fromstring(dos.readline(), dtype=float, sep=' ')
                dosmat = np.zeros((NEDOS, len(first), nions), dtype=float)
                dosmat[0,:,ion] = first
            dosmat[int(ion==0):,:,ion] = np.fromfile(dos,
                                            count = len(first)*(NEDOS-int(ion==0)),
                                            sep = ' ', dtype=float).reshape(-1,len(first))
    print(f"DOSCAR read")
    # if dosmat.shape[1] == 37:
    #     dosmat = np.concatenate([dosmat[:,i:i+1,:] for i in [0]+[1+4*j for j in range(9)]], axis=1)
    print(efermi, totmat.shape, dosmat.shape)
    return efermi, totmat, dosmat
    

def make_klines(KPOINTS, POSCAR):
    print(f"Reading {POSCAR}", end=' ')
    with open(POSCAR) as poscar:
        poscar.readline()
        Const  = np.fromstring(poscar.readline(), sep=' ')
        a1 = np.fromstring(poscar.readline(), sep=' ') * Const
        a2 = np.fromstring(poscar.readline(), sep=' ') * Const
        a3 = np.fromstring(poscar.readline(), sep=' ') * Const
        V = np.dot(a1, np.cross(a2, a3))
        b1 = 2*np.pi*np.cross(a2, a3) / V
        b2 = 2*np.pi*np.cross(a3, a1) / V
        b3 = 2*np.pi*np.cross(a1, a2) / V
    print(f"POSCAR read")
    print(f"Reading {KPOINTS}", end=' ')
    klines = []
    with open(KPOINTS) as kps:
        kps.readline()
        Npts = int(kps.readline().split()[0])
        kps.readline()
        kps.readline()
        start_index = 0
        while True:
            a_line = kps.readline().split()
            b_line = kps.readline().split()
            if a_line == [] or b_line == []: break

            a_pt, a_name = np.array(a_line[0:3], dtype=float), a_line[3]
            b_pt, b_name = np.array(b_line[0:3], dtype=float), b_line[3]
            if a_name == 'GAMMA': a_name = '$\Gamma$'
            if b_name == 'GAMMA': b_name = '$\Gamma$'

            a_pt = a_pt[0]*b1 + a_pt[1]*b2 + a_pt[2]*b3
            b_pt = b_pt[0]*b1 + b_pt[1]*b2 + b_pt[2]*b3

            klines.append(kline(a_pt, b_pt, Npts, start_index, a_name, b_name))
            start_index += Npts

            if kps.readline() == '': break
    print(f"KPOINTS read")
    print(klines)
    return klines
    
    
def read_procar(PROCAR):
    try: 
        with open(PROCAR+'.npy', 'rb') as f:
            promat = np.load(f)
            energs = np.load(f)
            E_max  = np.load(f)
            E_min  = np.load(f)
            E_gap  = np.load(f)
            print('Read from PROCAR.npy')
    except(FileNotFoundError):
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
            proj = 0
            promat = np.zeros((nkpts, nbands, nprjct, nions, 9), dtype=float)
            energs = np.zeros((nkpts, nbands), dtype=float)
            vbmax = 0
            for line in tqdm(procar, total=nlines - ((nions+1)*nprjct)*nbands*nkpts - 2):
                if "# of k-points:" in line:
                    proj += 1
                data = line.split()
                if len(data) > 0:
                    if data[0] == 'k-point': kpoint = int(data[1])-1
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
        with open(PROCAR+'.npy', 'wb') as f:
            np.save(f, promat)
            np.save(f, energs)
            np.save(f, E_max)
            np.save(f, E_min)
            np.save(f, E_gap)
            print('Saved to PROCAR.npy')
    return promat, energs, E_max, E_min, E_gap
    
    
def make_xline(klines, kline_nums = []):
    x_line = np.array([])
    x_ticks = [0]
    x_now   = 0.
    if kline_nums == []: current_klines = klines
    else: current_klines = [klines[abs(i)-1] for i in kline_nums]
    for line in current_klines:
        x_line = np.concatenate([x_line, np.linspace(x_now, x_now+abs(line), line.Npts, endpoint=True)])
        x_now += abs(line)
        x_ticks.append(x_now)
    return x_line, x_ticks


def dos_plot(dosmat, ax, projections = [], ions=[], orbits=[], 
             sum_orbits = False, sum_ions = True, label = '', ylim = [-3, 3]):
    orbit_names = ['$s$', '$p_y$', '$p_z$', '$p_x$', '$d_{xy}$', 
              '$d_{yz}$', '$d_{z^2}$', '$d_{xz}$', '$d_{x^2-y^2}$']
    Nproj = int((dosmat.shape[1]-1)/9)
    cmap = colormaps['Set1']
    colors = cmap(np.linspace(0, 1, 9))
    MAX = 0.
    MIN = 0.
    T_dosmat = dosmat.copy()
    # T_dosmat[:,1:,:] = np.abs(T_dosmat[:,1:,:])
    if orbits == []: orbits = np.arange(9)+1
    if sum_ions:
        for ion in ions[1:]:
            T_dosmat[:,1:,ions[0]-1] += T_dosmat[:,1:,ion-1]
        ions = [ions[0]]
    for ion in ions:
        for proj in projections:
            if sum_orbits:
                data = np.sum([T_dosmat[:,Nproj*(orbit-1)+proj,ion-1] for orbit in orbits], axis=0)
                MAX = np.max([MAX, np.max(data[(T_dosmat[:,0,ion-1]<ylim[1]) & (T_dosmat[:,0,ion-1]>ylim[0])])])
                MIN = np.min([MIN, np.min(data[(T_dosmat[:,0,ion-1]<ylim[1]) & (T_dosmat[:,0,ion-1]>ylim[0])])])
                ax.plot(data, T_dosmat[:,0,ion-1], 
                        label=label, color = colors[1], linewidth=3)
            else:
                for orbit in orbits:
                    data = T_dosmat[:,Nproj*(orbit-1)+proj,ion-1]
                    MAX = np.max([MAX, np.max(data[(T_dosmat[:,0,ion-1]<ylim[1]) & (T_dosmat[:,0,ion-1]>ylim[0])])])
                    MIN = np.min([MIN, np.min(data[(T_dosmat[:,0,ion-1]<ylim[1]) & (T_dosmat[:,0,ion-1]>ylim[0])])])
                    ax.plot(data, T_dosmat[:,0,ion-1], 
                            label=orbit_names[orbit-1], color = colors[orbit-1], linewidth=3)
    ax.set_xlim([MIN,MAX])
    ax.set_ylim(ylim)
    ax.set_xticks([MIN,MAX])


def band_plot(promat, energs, klines, ax,
              kline_nums = [], bands = [], projections = [], ions = [], orbits = [], 
              fatbands = True, bandlines = True, sum_orbits = False, sum_ions = True, 
              label = '', bandcolor='black'):
    x_line, x_ticks = make_xline(klines, kline_nums)
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

    if kline_nums == []: kline_nums = np.arange(len(klines))+1
    if kline_nums[0] > 0:
        xlabels = [klines[kline_nums[0]-1].a_name]
    if kline_nums[0] < 0:
        xlabels = [klines[-kline_nums[0]-1].b_name]
    for num in kline_nums:
        si = klines[abs(num)-1].start_index
        ei = klines[abs(num)-1].end_index
        if num > 0:
            ebands = np.concatenate([ebands, energs[si:ei+1:1]], axis = 0)
            pbands = np.concatenate([pbands, promat[si:ei+1:1]], axis = 0)
            xlabels.append(klines[num-1].b_name)
        if num < 0:
            ebands = np.concatenate([ebands, energs[ei:(si-1 if si-1 >= 0 else None):-1]], axis = 0)
            pbands = np.concatenate([pbands, promat[ei:(si-1 if si-1 >= 0 else None):-1]], axis = 0)
            xlabels.append(klines[-num-1].a_name)
    if len(bands) == 0:
        if bandlines: ax.plot(x_line, ebands, color=bandcolor, linewidth = 0.5)
    elif bandlines: ax.plot(x_line, ebands[:,bands], color=bandcolor, linewidth = 0.5)
    if sum_ions:
        for ion in ions[1:]:
            pbands[:,:,:,ions[0]-1,:] += pbands[:,:,:,ion-1,:]
        ions = [ions[0]]
    if fatbands:
        for ion in ions:
            for proj in projections:
                if sum_orbits:
                    ax.scatter(np.stack([x_line]*energs.shape[1], axis=1), 
                                ebands, np.sum(np.abs([100*pbands[:,:,proj-1,ion-1,j] for j in orbits]), axis=0), 
                                label = label, color = colors[1])
                else:
                    for j in orbits:
                        if len(bands) == 0:
                            ax.scatter(np.stack([x_line]*energs.shape[1], axis=1), 
                                    ebands, np.abs(100*pbands[:,:,proj-1,ion-1,j]), 
                                    label = orbit_names[j], color = colors[j])
                        else:
                            ax.scatter(np.stack([x_line]*len(bands), axis=1), 
                                    ebands[:,bands], 100*pbands[:,bands,proj-1,ion-1,j], 
                                    label = orbit_names[j], color = colors[j])
    for i in x_ticks[1:-1]:
        ax.axvline(i, color="black", linestyle=":")
    ax.axhline(0, color="black", linestyle="--", linewidth = 0.5)
    ax.set_xlim([0, x_line[-1]])
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(xlabels, fontsize=18)


### COMPARISON SPECIFIERS ###
calc   = 'pure-1x1'
mode   = 'SOC'
matr   = 'Cr'
orbt   = 'd'
proj   = 'np'

calc_2 = 'pure-1x2'
mode_2 = 'SOC'
matr_2 = 'Cr'
orbt_2 = 'd'
proj_2 = 'np'

kline_nums   = [3,4,5,6]
kline_nums_2 = [3,4,5,6]
k_colors     = {
                2:(0,1,0,0.07), 4:(0,1,0,0.07), 7:(0,1,0,0.07), 
                1:(1,0,0,0.07), 3:(1,0,0,0.07), 6:(1,0,0,0.07), 8:(1,0,0,0.07), 
                5:(0,0,1,0.07)
                }

one_plot     = False
show         = True
figsize      = (9, 6)
ymin         = -1.1
ymax         =  2.4
T_mod        = '2'
#############################
kline_nums_3   = [3,4,5,6]


### RUN THE COMPARISON ###
proj_dict = {'up':1, 'down':2, 'x':2, 'y':3, 'z':4, 'np':0, 'tot':1}
projection = proj_dict[proj]
projection_2 = proj_dict[proj_2]
projections = [projection]
projections_2 = [projection_2]

fatbands = bool(projection*projection_2)

size = int(calc[-3]) + int(calc[-1])
if size == 2:
    if matr == 'Cr': ions = [1,2,3,4]
    if matr == 'S': ions = [5,6,7,8]
    if matr == 'Br': ions = [9,10,11] + [12]*(calc[:4]=='pure')
    if matr == 'Cl': ions = [12]
elif size == 3:
    if matr == 'Cr': ions = [i for i in range(1,9)]
    if matr == 'S': ions = [i for i in range(9,17)]
    if matr == 'Br': ions = [i for i in range(17,25)]
elif size == 4:
    if matr == 'Cr': ions = [i for i in range(1,17)]
    if matr == 'S': ions = [i for i in range(17,33)]
    if matr == 'Br': ions = [i for i in range(33,49)]

size_2 = int(calc_2[-3]) + int(calc_2[-1])
if size_2 == 2:
    if matr_2 == 'Cr': ions_2 = [1,2,3,4]
    if matr_2 == 'S': ions_2 = [5,6,7,8]
    if matr_2 == 'Br': ions_2 = [9,10,11] + [12]*(calc_2[:4]=='pure')
    if matr_2 == 'Cl': ions_2 = [12]*(calc_2[:4]!='pure')
elif size_2 == 3:
    if matr_2 == 'Cr': ions_2 = [i for i in range(1,9)]
    if matr_2 == 'S': ions_2 = [i for i in range(9,17)]
    if matr_2 == 'Br': ions_2 = [i for i in range(17,25)]
elif size_2 == 4:
    if matr_2 == 'Cr': ions_2 = [i for i in range(1,17)]
    if matr_2 == 'S': ions_2 = [i for i in range(17,33)]
    if matr_2 == 'Br': ions_2 = [i for i in range(33,49)]

orbt_dict = {'s':[1], 'p':[2,3,4], 'd':[5,6,7,8,9]}
orbits = orbt_dict[orbt]
orbits_2 = orbt_dict[orbt_2]


PROCAR  = '/home/karol/Monolayers/VASP/zps/'+calc+'/'+mode+'/PROCAR'
POSCAR  = '/home/karol/Monolayers/VASP/zps/'+calc+'/'+mode+'/POSCAR'
KPOINTS = '/home/karol/Monolayers/VASP/zps/'+calc+'/'+mode+'/KPOINTS'

PROCAR_2  = '/home/karol/Monolayers/VASP/zps/'+calc_2+'/'+mode_2+'/PROCAR'
POSCAR_2  = '/home/karol/Monolayers/VASP/zps/'+calc_2+'/'+mode_2+'/POSCAR'
KPOINTS_2 = '/home/karol/Monolayers/VASP/zps/'+calc_2+'/'+mode_2+'/KPOINTS'

#############
PROCAR_3  = '/home/karol/Monolayers/VASP/zps/pure-2x1/SOC/PROCAR'
POSCAR_3  = '/home/karol/Monolayers/VASP/zps/pure-2x1/SOC/POSCAR'
KPOINTS_3 = '/home/karol/Monolayers/VASP/zps/pure-2x1/SOC/KPOINTS'
#############

title = calc+'-'+mode+('-'+matr+'-'+orbt+'-'+proj)*fatbands+' VS '+\
        calc_2+'-'+mode_2+('-'+matr_2+'-'+orbt_2+'-'+proj_2)*fatbands\
        +' VS '+'pure-2x1-SOC'


klines = make_klines(KPOINTS, POSCAR)
promat, energs, E_max, E_min, E_gap = read_procar(PROCAR)
energs -= E_max

klines_2 = make_klines(KPOINTS_2, POSCAR_2)
promat_2, energs_2, E_max_2, E_min_2, E_gap_2 = read_procar(PROCAR_2)
energs_2 -= E_max_2

############
klines_3 = make_klines(KPOINTS_3, POSCAR_3)
promat_3, energs_3, E_max_3, E_min_3, E_gap_3 = read_procar(PROCAR_3)
energs_3 -= E_max_3
############


width_ratios = [abs(klines[knum-1]) for knum in kline_nums]
width_ratios_2 = [abs(klines_2[knum-1]) for knum in kline_nums_2]

rem = len(kline_nums)
norm = 0
for i in width_ratios: norm += i
for i in width_ratios_2: norm += i
#############
rem2 = len(kline_nums) + len(kline_nums_2) + 1
width_ratios_3 = [abs(klines_3[knum-1]) for knum in kline_nums_3]
for i in width_ratios_3: norm += i
#############
norm /= len(width_ratios) + len(width_ratios_2) + len(width_ratios_3)
kline_nums = kline_nums + [0] + kline_nums_2 + [-1] + kline_nums_2
width_ratios = width_ratios + [norm] + width_ratios_2 + [norm] + width_ratios_3

fig = plt.figure(figsize=figsize)
fig.suptitle(title, size=24)

gs = fig.add_gridspec(1, len(kline_nums_2) if one_plot else len(kline_nums), 
                      width_ratios = width_ratios_2 if one_plot else width_ratios, #height_ratios=(2, 1),
                      left=0.05, right=0.97, bottom=0.06, top=0.92,
                      wspace=0.0, hspace=0.1)

E_gap_now = E_gap
axs = gs.subplots()
for knum, ax in tqdm(np.ndenumerate(axs), total = len(kline_nums)):
    if kline_nums[knum[0]] == 0 : 
        ax.axis('off')
        promat, energs, klines, E_gap_now = promat_2, energs_2, klines_2, E_gap_2
        projections, orbits, ions = projections_2, orbits_2, ions_2
        continue 
    if kline_nums[knum[0]] == -1 : 
        ax.axis('off')
        promat, energs, klines, E_gap_now = promat_3, energs_3, klines_3, E_gap_3
        continue 
    band_plot(promat, energs, klines, 
              kline_nums=[kline_nums[knum[0]]], bands=[], projections=projections, ions=ions, orbits=orbits, 
              ax = ax, sum_orbits=True, fatbands=fatbands, bandcolor='black')
    if one_plot:
        band_plot(promat_2, energs_2, klines_2, 
              kline_nums=[kline_nums[knum[0]]], bands=[], projections=projections, ions=ions, orbits=orbits, 
              ax = ax, sum_orbits=True, fatbands=fatbands, bandcolor='red')
    ax.set(yticks=[])
    ax.set_ylim([ymin, ymax])
    ax.axhline(E_gap_now, color="black", linestyle="--", linewidth = 0.5)
    if kline_nums[knum[0]] in k_colors.keys() :
        ax.set_facecolor(k_colors[kline_nums[knum[0]]]) 
axs[0].set_yticks([ymin,0,E_gap,ymax])
axs[0].set_yticklabels([ymin,0,f'{E_gap:.2f}',ymax])
if not one_plot:
    axs[rem+1].set_yticks([ymin,0,E_gap_2,ymax])
    axs[rem+1].set_yticklabels([ymin,0,f'{E_gap_2:.2f}',ymax])
axs[rem2+1].set_yticks([ymin,0,E_gap_3,ymax])
axs[rem2+1].set_yticklabels([ymin,0,f'{E_gap_3:.2f}',ymax])

if show:
    plt.show()
else:
    name = calc+'-'+mode+('-'+matr+'-'+orbt+'-'+proj)*fatbands+'-VS-'+\
        calc_2+'-'+mode_2+('-'+matr_2+'-'+orbt_2+'-'+proj_2)*fatbands\
           +'pure-2x1-SOC'
    plt.savefig(f"/home/karol/Monolayers/VASP/zps/plots/{name}{T_mod}")



# DOSCAR = '/home/karol/Monolayers/VASP/zps/pure-1x1/SOC/DOSCAR'
# efermi, totmat, dosmat = read_doscar(DOSCAR)
# dosmat[:,0] -= efermi

# ax1 = fig.add_subplot(gs[0, 0:3])
# band_plot(promat, energs, klines, 
#           kline_nums=[4], bands=[], projections=projections, ions=ions, orbits=orbits, 
#           ax = ax1, sum_orbits=True, fatbands=fatbands)
# ax1.set_ylim([-2, 3])
# ax1.axhline(E_gap, color="black", linestyle="--", linewidth = 0.5)
# ax1.set_yticks([-2,0,E_gap,3])
# ax1.set_yticklabels([-2,0,f'{E_gap:.2f}',3])

# ax2 = fig.add_subplot(gs[0, 3:6])
# band_plot(promat, energs, klines, 
#           kline_nums=[5], bands=[], projections=projections, ions=ions, orbits=orbits, 
#           ax = ax2, sum_orbits=True, fatbands=fatbands)
# ax2.set_ylim([-2, 3])
# ax2.set_yticks([])
# ax2.axhline(E_gap, color="black", linestyle="--", linewidth = 0.5)

# klines_2 = make_klines(KPOINTS_2, POSCAR_2)
# promat_2, energs_2, E_max_2, E_min_2, E_gap_2 = read_procar(PROCAR_2)
# energs_2 -= E_max_2

# ax3 = fig.add_subplot(gs[0, 7:10])
# band_plot(promat_2, energs_2, klines_2, 
#           kline_nums=[1,2,3,4], bands=[], projections=projections_2, ions=ions_2, orbits=orbits_2, 
#           ax = ax3, sum_orbits=True, fatbands=fatbands)
# ax3.set_ylim([-2, 3])
# ax3.axhline(E_gap_2, color="black", linestyle="--", linewidth = 0.5)
# ax3.set_yticks([-2,0,E_gap_2,3])
# ax3.set_yticklabels([-2,0,f'{E_gap_2:.2f}',3])

# ax4 = fig.add_subplot(gs[0, 10:13])
# band_plot(promat_2, energs_2, klines_2, 
#           kline_nums=[5,6,7,8], bands=[], projections=projections_2, ions=ions_2, orbits=orbits_2, 
#           ax = ax4, sum_orbits=True, fatbands=fatbands)
# ax4.set_ylim([-2, 3])
# ax4.set_yticks([])
# ax4.axhline(E_gap_2, color="black", linestyle="--", linewidth = 0.5)


# ax3 = fig.add_subplot(gs[0, 6])
# dos_plot(dosmat, ax3, 
#          projections=projections, ions=ions, orbits=orbits,
#          sum_orbits=True, sum_ions=True, ylim = [-3, 3])
# ax3.set_yticks([])
# ax3.set_title('DOS', fontsize=18)

