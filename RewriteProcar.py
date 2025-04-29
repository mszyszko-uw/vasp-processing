import numpy as np
from tqdm import tqdm

def readprocar(PROCAR):
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


def RewriteProcar(folder, pro_name='PROCAR', re_name = '',
                  spins = 'sum', start_index = 1):
    PROCAR = folder + pro_name
    promat, energs, E_max, E_min, E_gap, kindexer = readprocar(PROCAR)
    promat = promat[:,:,:,:,:9]
    promat = promat[:,:,:,[0,4,8],:]+promat[:,:,:,[1,5,9],:]+\
            promat[:,:,:,[2,6,10],:]+promat[:,:,:,[3,7,11],:]
    promat = np.concatenate((promat[:,:,:,0,:], promat[:,:,:,1,:], promat[:,:,:,2,:]), -1)
    if spins == 'sum':
        promat = promat.sum(2)[:,:,None,:]
        files = ['sum']
    else: files = ['up', 'down']
    if len(re_name) == 0:
        re_name = pro_name + '_from_' + str(start_index)

    kmat = np.hstack((np.arange(0,len(kindexer))[:,None]+1, kindexer))
    for proj, name in enumerate(files):
        with open(folder + re_name + '_' + name, 'w') as procar:
            procar.write('indeks_punktu_k\tk_x\t\tk_y\t\tk_z\t\tenergia\t\t')
            procar.write('Cr s\t\tCr py\t\tCr pz\t\tCr px\t\tCr dxy\t\tCr dyz\t\tCr dz2\t\tCr dxz\t\tCr x2y2\t\t')
            procar.write('S  s\t\tS  py\t\tS  pz\t\tS  px\t\tS  dxy\t\tS  dyz\t\tS  dz2\t\tS  dxz\t\tS  x2y2\t\t')
            procar.write('Br s\t\tBr py\t\tBr pz\t\tBr px\t\tBr dxy\t\tBr dyz\t\tBr dz2\t\tBr dxz\t\tBr x2y2\n')
            for iband in range(promat.shape[1]):
                writemat = np.hstack((kmat, energs[:,iband,proj][:,None], promat[:,iband,proj,:]))
                np.savetxt(procar, writemat[start_index-1:], '%12.8f', '\t')
                procar.write('\n\n')

if __name__ == '__main__':
    folder = '/home/karol/Monolayers/T_Wozniak/BND_AFM/'
    pro_name = 'PROCAR_YGX'
    re_name = ''     # Mozna samemu ustawic nazwe, pusty zachowa stara plus doda indeks startowy
    spins  = 'split' # 'sum' or 'split'
    start_index = 190
    RewriteProcar(folder, pro_name, re_name, spins, start_index)