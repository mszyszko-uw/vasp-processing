import numpy as np
import subprocess
import os
import stat
from tqdm import tqdm

def ReadPROCAR(PROCAR):
    try: 
        with open(PROCAR+'.npy', 'rb') as f:
            promat = np.load(f)
            energs = np.load(f)
            E_max  = np.load(f)
            E_min  = np.load(f)
            E_gap  = np.load(f)
            kindexer = np.load(f)
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
            kindexer = []
            proj = 0
            promat = np.zeros((nkpts, nbands, nprjct, nions, norbitals), dtype=float)
            energs = np.zeros((nkpts, nbands, nprjct), dtype=float)
            vbmax = 0
            for line in tqdm(procar, total=nlines - ((nions+1)*nprjct)*nbands*nkpts - 2):
                if "# of k-points:" in line:
                    proj += 1
                data = line.split()
                if len(data) > 0:
                    if data[0] == 'k-point': 
                        kpoint = int(data[1])-1
                        if proj == 0:
                            kindexer.append([float(num) for num in data[3:6]])
                    if data[0] == 'band': 
                        band = int(data[1])-1
                        energs[kpoint, band, proj] = float(data[4])
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
        kindexer = np.array(kindexer)
        with open(PROCAR+'.npy', 'wb') as f:
            np.save(f, promat)
            np.save(f, energs)
            np.save(f, E_max)
            np.save(f, E_min)
            np.save(f, E_gap)
            np.save(f, kindexer)
            print('Saved to PROCAR.npy')
    return promat, energs, E_max, E_min, E_gap, kindexer


# Ulepszona funkcja do obrobki plikow BSEFATBAND (kiedys ExcitonFatband)
def ProcessBSEFATBAND(l_exc, i_band, 
               i_kpts = [None, None, None],  # Można ograniczyć wybór kpuntow do wybarnych ineksow albo wspolrzednych
               c_kpts = [None, None, None],  # np. dla jednego przekroju kz i_kpts=[None,None,1], c_kpts=[None,None,0.0]
               spins='sum', # spins=('up','down','sum')
               folder = './'):
    
    try:
        with open(folder+'kvecs.npy', 'rb') as f:
            kvecs = [[],[],[]]
            kvecs[0] = np.load(folder+'kvecs.npy')
            kvecs[1] = np.load(folder+'kvecs.npy')
            kvecs[2] = np.load(folder+'kvecs.npy')
            print("kvecs found, updating i_kpts and c_kpts")
            if i_kpts != [None, None, None]:
                c_kpts = [ikpt if ikpt==None else kvecs[dir][ikpt-1] for dir, ikpt in enumerate(i_kpts)]
            if c_kpts[0] != None: 
                for ikxnow, kxnow in enumerate(kvecs[0]): 
                    if f'{c_kpts[0]:9.5f}' == f'{kxnow:9.5f}': i_kpts[0] = ikxnow+1
            if c_kpts[1] != None: 
                for ikynow, kynow in enumerate(kvecs[1]): 
                    if f'{c_kpts[1]:9.5f}' == f'{kynow:9.5f}': i_kpts[1] = ikynow+1
            if c_kpts[2] != None: 
                for ikznow, kznow in enumerate(kvecs[2]): 
                    if f'{c_kpts[2]:9.5f}' == f'{kznow:9.5f}': i_kpts[2] = ikznow+1

    except FileNotFoundError:
        print("kvecs not found, building from BSEFATBAND")
        with open(folder+f'BSEFATBAND') as mainfile:
            mainfile.readline()
            mainfile.readline()
            kvecs = [[], [], []]
            ikxnow, ikynow, ikznow = 0, 0, 0
            print('Making kvecs and updating i_kpts and c_kpts')
            for line in mainfile:
                if 'BSE' in line: break
                kxnow, kynow, kznow = float(line.split()[0]), float(line.split()[1]), float(line.split()[2])
                if kxnow not in kvecs[0]:
                    ikxnow += 1
                    kvecs[0].append(kxnow)
                    if c_kpts[0] != None and f'{c_kpts[0]:9.5f}' == f'{kxnow:9.5f}': i_kpts[0] = ikxnow
                if kynow not in kvecs[1]:
                    ikynow += 1
                    kvecs[1].append(kynow)
                    if c_kpts[1] != None and f'{c_kpts[1]:9.5f}' == f'{kynow:9.5f}': i_kpts[1] = ikynow
                if kznow not in kvecs[2]:
                    ikznow += 1
                    kvecs[2].append(kznow)
                    if c_kpts[2] != None and f'{c_kpts[2]:9.5f}' == f'{kznow:9.5f}': i_kpts[2] = ikznow
            if i_kpts != [None, None, None]:
                c_kpts = [ikpt if ikpt==None else kvecs[dir][ikpt-1] for dir, ikpt in enumerate(i_kpts)]
            np.save(folder+'kvecs.npy', kvecs[0])
            np.save(folder+'kvecs.npy', kvecs[1])
            np.save(folder+'kvecs.npy', kvecs[2])

    print('i_kpts: ', i_kpts, '\tc_kpts: ', c_kpts)
    for tries in range(2):
        try:
            totmat = []
            print('Looking for small BSE files')
            for i_exc in l_exc:
                filename = f'BSE_{spins}_{i_exc}_{i_band}' +\
                          (f'_x{i_kpts[0]}' if i_kpts[0] != None else '') +\
                          (f'_y{i_kpts[1]}' if i_kpts[1] != None else '') +\
                          (f'_z{i_kpts[2]}' if i_kpts[2] != None else '')
                with open(folder+filename, 'r+') as bsefile:
                    if len(bsefile.readline().split()) > 4:
                        print(f'Cropping {filename}')
                        bsefile.seek(0)
                        bsemat = bsefile.read()
                        bsemat = bsemat.replace('+i*', ' ', )
                        bsemat = np.fromstring(bsemat, dtype=float, sep=' ').reshape((-1,10))
                        bsemat = bsemat[:,[0,1,2,5]]
                        if spins=='sum' :
                            bsemat[:len(bsemat)//2,3] += bsemat[len(bsemat)//2:,3]
                            bsemat = bsemat[:len(bsemat)//2,:]
                        kpoint = bsemat[0,0:3]
                        repstep = 0
                        for row in bsemat[:,:]:
                            if (row[:3] == kpoint).all(): repstep += 1
                            else: break
                        summat = bsemat[0::repstep, :]
                        for rep in range(1,repstep):
                            summat[:,3] += bsemat[rep::repstep,3]
                        bsemat = summat
                        bsefile.seek(0)
                        np.savetxt(bsefile, bsemat, fmt='%.5f \t%.5f \t%.5f \t%.7f')
                        bsefile.truncate()
                    if len(bsefile.readline().split()) == 4:
                        bsefile.seek(0)
                        bsemat = np.fromfile(bsefile, dtype=float, sep=' ').reshape((-1,4))
                if len(totmat) == 0: totmat = bsemat
                else: totmat[:,3] += bsemat[:,3]

            print('BSE info retrieved')
            return totmat

        except FileNotFoundError:
            print('No small BSE files found')
            with open(folder+f'BSEFATBAND') as mainfile:
                rank = int(mainfile.readline().split()[0])
            print('Writing "processing.sh" script')
            with open(folder+'processing.sh', 'w') as processfile:
                processfile.write(f'#!/bin/bash\ncd {folder}\n')
                for i_exc in l_exc:
                    start = 2+(1+rank)*(i_exc-1)+1
                    end   = start + rank - 1
                    processfile.write(f'sed -n "{start},{end}p;{end}q" BSEFATBAND > bse_{i_exc}\n')
                    if spins == 'sum':
                        processfile.write(f'cat bse_{i_exc} > BSE_{i_exc}\n')
                    if spins == 'up': 
                        processfile.write(f'sed -n "1,{int(rank/2)}p;{int(rank/2)}q" bse_{i_exc} > BSE_{i_exc}\n')
                    if spins == 'down': 
                        processfile.write(f'sed -n "{int(rank/2)+1},{rank}p;{rank}q" bse_{i_exc} > BSE_{i_exc}\n')
                    kstrs = [f'{ckpt:9.5f}' if ckpt != None else '.\{9\}' for ckpt in c_kpts]
                    outname = f'BSE_{spins}_{i_exc}_{i_band}' +\
                              (f'_x{i_kpts[0]}' if i_kpts[0] != None else '') +\
                              (f'_y{i_kpts[1]}' if i_kpts[1] != None else '') +\
                              (f'_z{i_kpts[2]}' if i_kpts[2] != None else '')
                    processfile.write(f'grep "^{kstrs[0]}{kstrs[1]}{kstrs[2]}.\{{42\}}.*{i_band}.*.\{{29\}}$" ' + \
                                      f'BSE_{i_exc} > {outname}\n')
                for i_exc in l_exc: processfile.write(f'rm bse_{i_exc}\nrm BSE_{i_exc}\n')

            st = os.stat(folder+'processing.sh')
            os.chmod(folder+'processing.sh', st.st_mode | stat.S_IEXEC)
            print('Running "processing.sh" script')
            subprocess.call(f'{folder}processing.sh')


def MakeBSEband(lexc, bands, spins, start_index=1, BSEfolder = './', PROCAR = './PROCAR'):
    promat, energs, E_max, E_min, E_gap, kindexer = ReadPROCAR(PROCAR)
    kindexer = kindexer.round(5)
    with open(BSEfolder + f'bseband_{"-".join([str(exc) for exc in lexc])}_{spins}', 'w') as bseband:
        bseband.write('k_idx\tk_x\tk_y\tk_z\tenergia\t\tAbs(X)/W_k\n')
        for iband in np.arange(bands[0],bands[1]+.1,dtype=int):
            bsemat = ProcessBSEFATBAND(lexc, iband, spins=spins, folder=BSEfolder)
            kptsbse  = bsemat[:,:3].round(5)
            bsevals  = bsemat[:,3]
            for ikpt, kpoint in enumerate(kindexer):
                if ikpt < start_index-1: continue
                try:
                    index = np.where((kptsbse == kpoint).all(1))[0][0]
                    bseval = bsevals[index]
                except IndexError:
                    bseval = 0
                bseband.write(f'{ikpt+1}\t{kpoint[0]}\t{kpoint[1]}\t{kpoint[2]}\t')
                bseband.write(f'{energs[ikpt,iband-1,0]:12.8f}\t{bseval:12.8f}\n')
            bseband.write('\n\n')


if __name__ == '__main__':
    lexc  = [2,3,4]
    bands = (29,52)   # Max i min numer pasma w BSEFATBAND
    spins = 'down'    # Wybór spinu albo 'sum' dla AFM
    start_index = 190 # Punkt z PROCAR od ktorego trzeba zaczac
    BSEfolder = '/home/karol/Monolayers/T_Wozniak/BND_FM/YGX/'
    PROCAR    = '/home/karol/Monolayers/T_Wozniak/BND_FM/YGX/PROCAR'

    MakeBSEband(lexc, bands, spins, start_index, BSEfolder, PROCAR)