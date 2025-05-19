import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile

##### Autor: Karol Gałązka
## Cel: Implementacja klasy DFT_Output, czytanie wybranych plików
##      w całości lub tylko w potrzebnym zakresie dla przyspieszenia.
##      Uwspólnienie formatu wczytywania dla różnych plików

class PhysicalConstants:
    pi      = 3.141592653589793     # no unit
    c       = 299792458;            # m/s
    e       = 1.602176565E-19;      # As = C
    me      = 9.10938291E-31;       # kg
    k       = 1.3806488E-23;        # J/K
    h       = 6.62606957E-34;       # Js = Nms = VAs^2
    hbar    = 6.62606957E-34/2/pi;  # Js = Nms = VAs^2
    eps0    = 8.854187817E-12;      # dielectric constant [As/Vm = F/m]
    a0      = 0.52917721092E-10;    # Bohr radius [m]
    Ry      = 2.179872171E-18;      # Rydberg energy [J] 
    u       = 1.660538921E-27;      # atomic mass unit [kg]
    NA      = 6.02214129E23;        # Avogadro constant[1/mol]


class DFT_Output :
    
    def __init__ (self, source, dir = '', name = 'DFT_Output'):
        self.WAVEDER    = None
        self.EIGENVAL   = None
        self.WAVEDERF   = None
        # self.DOSCAR     = None
        self.mommat2up  = None

        if source == 'VASP':
            self.WAVEDER    = dir + "WAVEDER"
            self.EIGENVAL   = dir + "EIGENVAL"
            self.WAVEDERF   = dir + "WAVEDERF"
            # self.DOSCAR     = dir + "DOSCAR"

        else :
            self.mommat2up  = dir + source

        self.source = source
        self.name = name
        self.band_vecs = {}
        self.P_matrix  = np.array([])
        self.eval = np.array([])
        self.occ = np.array([])

    def set_VASP_files(self, WAVEDER = None, EIGENVAL = None, WAVEDERF = None):
        if WAVEDER != None: self.WAVEDER  = WAVEDER
        if EIGENVAL != None: self.EIGENVAL = EIGENVAL
        if WAVEDERF != None: self.WAVEDERF = WAVEDERF

    def set_wien_file(self, file = None):
        self.mommat2up = file

    def read_P_matrix(self, file):
        self.P_matrix = np.load(file, mmap_mode='r')

    def __getitem__ (self, index):
        if self.P_matrix.size > 0 :
            return self.P_matrix[:,index-1,:,:]
        elif index in self.known_bands():
            return self.band_vecs[index]
        else : 
            self.read_bands([index])
            return self.band_vecs[index]
    
    def known_bands (self):
        bands = list(self.band_vecs.keys())
        bands.sort()
        return bands

    def VASP_ThinRead(self, bands):
        vecs = self.band_vecs
        for band in bands:
            vecs[band] = []
        with open(self.WAVEDERF,'r') as fh:
            Q = np.fromstring(fh.readline(), dtype=int, sep = ' ')
            nkpts = Q[1]
            nbands = Q[2]
            ibands = bands.copy()
            for line in fh:
                if len(ibands) == 0: break
                for band in ibands:
                    if f" {band} " in line[:10]:
                        vecs[band].append(np.fromstring(line, dtype=float, sep=' '))
                        if len(vecs[band]) == nkpts*nbands : ibands.remove(band)

        for band in bands:
            vec = np.array(vecs[band])
            vec = np.concatenate([vec[:,6:7]   + vec[:,7:8]*1j,
                                  vec[:,8:9]   + vec[:,9:10]*1j,
                                  vec[:,10:11] + vec[:,11:12]*1j,
                                  -vec[:,4:5]  + vec[:,1:2]     ], axis = 1)
            vec[:,:3] = vec[:,:3] * np.abs(vec[:,3:4]) 
            vecs[band] = vec.reshape((nkpts,nbands,4))

        return vecs
    
    def wien_ThinRead(self, bands):
        nkpnts, nbands = 0, 0
        nlines, skip   = 0, 0     # Optimisation variables
        vecs = self.band_vecs
        for band in bands:
            vecs[band] = []
        with open(self.mommat2up,'r') as fh:
            ibands = bands.copy()
            for line in fh:
                if skip: skip -= 1 ; continue
                if len(ibands) == 0: 
                    skip = nlines
                    ibands = bands.copy()
                if "KP:" in line:
                    nkpnts += 1
                    nbands = int(line.split()[6])
                    nlines = int(0.5*nbands*(nbands+1))
                for band in ibands:
                    if f" {band} " in line:
                        vecs[band].append(np.fromstring(line, dtype=float, sep=' '))
                        if len(vecs[band]) == nkpnts * nbands and band != 1: ibands.remove(band)
                nlines -= 1

        for band in bands:
            vec = np.array(vecs[band])
            vec = np.concatenate([vec[:,2:3]+vec[:,3:4]*1j,
                                vec[:,4:5]+vec[:,5:6]*1j,
                                vec[:,6:7]+vec[:,7:8]*1j,
                                vec[:,8:9]], axis = 1)
            vec[:band-1,0:3] = np.conj(vec[:band-1,0:3])
            vec[:band-1,3] = -vec[:band-1,3]
            vecs[band] = vec.reshape((nkpnts,nbands,4))

        return vecs
    
    def read_bands(self, bands):
        if self.mommat2up != None : return self.wien_ThinRead(bands)
        else : return self.VASP_ThinRead(bands)

    def Optical_Selection_Rules(self, kpoint, vb_max, num_vb, num_cb, save = '', txt = False):
        # calculate the optical matrix elements for circular and linear polarized 
        # light and determine the dipole strengths and degree of circular polarization
        #
        # the k-point is defined by the WAVEDERF file
        # the transitions are given in the vb_list and cb_list

        # find VBM

        vb_list = [vb_max - num_vb + i + 1 for i in range(num_vb)]
        cb_list = [vb_max + i + 1 for i in range(num_cb)]

        bands = vb_list + cb_list

        if self.P_matrix.size == 0 :
            bands_to_read = []
            for band in bands:
                if band not in self.known_bands():
                    bands_to_read.append(band)
            self.read_bands(bands_to_read)

        # pick transitions of interest
        # vb_list = np.array([v-1, v])
        # cb_list = np.array([c, c+1])
        optical = []

        # polarization verctors
        e_plus  = 1/np.sqrt(2)*np.array([[1,  1j, 0]])
        e_minus = 1/np.sqrt(2)*np.array([[1, -1j, 0]])

        vecs = self.band_vecs
        e_flag = 0

        for vb in vb_list:
            for cb in cb_list:
                # momentum matrix element
                # complex valued and not uniquely defined by an arbitray phase

                if self.P_matrix.size > 0 :
                    p = np.array([[self.P_matrix[kpoint-1, cb-1, vb-1, 0], 
                                   self.P_matrix[kpoint-1, cb-1, vb-1, 1], 
                                   self.P_matrix[kpoint-1, cb-1, vb-1, 2]]])

                else :
                    p = np.array([[vecs[cb][kpoint-1, vb-1, 0], 
                                   vecs[cb][kpoint-1, vb-1, 1], 
                                   vecs[cb][kpoint-1, vb-1, 2]]])

                # circular polarization
                p_plus  = e_plus @ p.T  
                p_minus = e_minus @ p.T

                # oscillator strengths
                p_plus_squ  =   p_plus[0,0] * np.conj(p_plus[0,0])
                p_minus_squ =  p_minus[0,0] * np.conj(p_minus[0,0])
                p_x_squ     =        p[0,0] * np.conj(p[0,0])
                p_y_squ     =        p[0,1] * np.conj(p[0,1])
                p_z_squ     =        p[0,2] * np.conj(p[0,2])
                p_squ       = p_x_squ + p_y_squ + p_z_squ
                
                # degree of circular polarization
                polariz = (p_plus_squ - p_minus_squ)/(p_plus_squ + p_minus_squ)
                
                optical += [vb, cb, polariz, p_plus_squ, p_minus_squ, 
                            p_x_squ, p_y_squ, p_z_squ, p_squ]
                
                if self.eval.size > 0 : 
                    optical += [self.eval[cb-1,0] - self.eval[vb-1,0]]
                    e_flag = 1
                
        optical = np.real(np.array(optical)).reshape(-1, 9 + e_flag)
        header = "vb\t\t# cb\t\tcir_pol\t\t|P+|^2\t\t|P-|^2\t\t|Px|^2\t\t|Py|^2\t\t|Pz|^2\t\t|P|^2"
        if e_flag : header += "\t\t  dE"

        if save : 
            if txt : np.savetxt(save+'.txt', optical, fmt = '%.4f', delimiter='\t\t', header=header)
            else :   np.save(save, optical)
                
        return "  vb \tcb \tcir_pol\t|P+|^2\t|P-|^2\t|Px|^2\t|Py|^2\t|Pz|^2\t|P|^2\t  dE\n" + \
                np.array2string(optical, formatter={'float_kind':lambda x: "%.2f" % x}, separator='\t', suppress_small=True)
    
    def L_Plot(self, kpoint, direction, vb_max, num_vb, num_cb, plot = True, save = '', txt = False):
    # Calculates orbital angular momentum (L) and Berry curvature (BC) 
    # from matrix elements 
    # and band energies read from vasp's WAVEDERF
    # for selected k point number
    # direction: 1 - X, 2 - Y, 3 - Z
        if   direction == 1 : beta = 2; gamma = 3
        elif direction == 2 : beta = 1; gamma = 3
        elif direction == 3 : beta = 1; gamma = 2
        else : raise Exception('wrong direction')

        val = [vb_max - num_vb + i + 1 for i in range(num_vb)]
        con = [vb_max + i + 1 for i in range(num_cb)]

        bands_oi = val + con
        nbands_oi = len(bands_oi)

        if self.P_matrix.size == 0 :
            bands_to_read = []
            for band in bands_oi:
                if band not in self.known_bands():
                    bands_to_read.append(band)
            self.read_bands(bands_to_read)
            P_vecs = self.band_vecs
            nbands = P_vecs[vb_max].shape[1]
        else : nbands = self.P_matrix.shape[1]

        L  = np.zeros((nbands, nbands_oi), dtype = complex)

        for count, band in enumerate(bands_oi):
            if self.P_matrix.size > 0 : P_vec = self.P_matrix[kpoint-1, band-1,:,:].copy()
            else :                      P_vec = P_vecs[band][kpoint-1].copy()
            P_vec[abs(P_vec[:,3]) < .00000001] = [0,0,0,1]
            L_iter = 2*np.imag(P_vec[:,beta-1]*np.conj(P_vec[:,gamma-1]))/P_vec[:,3]
            #BC_iter = 2*np.imag(cP[:,beta-1]*np.conj(cP[:,gamma-1]))/cP[:,3]**2
            L[:,count] = np.cumsum(L_iter)

        C = PhysicalConstants
        if self.source == 'VASP': 
            convert =      C.e * 1E-20 * C.me * 4 * C.pi**2 / C.h**2
        else: 
            convert = C.Ry * (C.a0)**2 * C.me * 4 * C.pi**2 / C.h**2

        L = np.real(L*convert)

        if plot:
            plt.figure()
            plt.plot(L, '-')
            dirs = {1: 'X', 2: 'Y', 3: 'Z'}
            plt.title('K-point: ' + str(kpoint+1) + ', ' \
                    + 'Direction: ' + dirs[direction])
            plt.xlabel('number of bands')
            plt.ylabel('L (hbar)')
            plt.legend(['v'+str(i) for i in val]+['c'+str(i) for i in con] )
            if plot == 'show': plt.show()
            if plot == 'save': plt.savefig('L_Plot.png')
            else: plt.savefig(plot)

        if save : 
            if txt : np.savetxt(save+'.txt', np.concatenate([np.array(bands_oi).reshape((1,-1)), L], axis = 0))
            else :   np.save(save, np.concatenate([np.array(bands_oi).reshape((1,-1)), L], axis = 0))

        return L

    def read_EIGENVAL(self):
        with open(self.EIGENVAL) as eig:
            # Get relevant data
            for i in range(5):
                eig.readline()
            data = np.fromstring(eig.readline(), 
                                dtype=int, 
                                sep = ' ')
            eig.readline() #skip next empty line
            vb_max, nkpts, nbands = data

            Q = np.empty((nbands,3,nkpts))
            for kpnt in range(nkpts):
                eig.readline()
                Q[:,:,kpnt] = np.fromfile(eig,
                                        count = 3*nbands,
                                        sep = ' ').reshape(-1,3)
        
        eval, occ = Q[:,1:2,:], Q[:,2:3,:]

        self.eval = eval
        self.occ  = occ

        return eval, occ

    def read_WAVEDER(self):
        with FortranFile(self.WAVEDER, 'r') as file:
            nb_tot, nbands_cder, nkpts, ispin = file.read_record(dtype= np.int32)
            nodesn_i_dielectric_function = file.read_record(dtype= float)
            wplasmon = file.read_record(dtype= float).reshape(3,3)
            cder = file.read_record(dtype= np.complex64).reshape(3, ispin,nkpts,nbands_cder,nb_tot)

        # Format cder matrix to resemble P from WAVEDERF
        P = np.moveaxis(cder, [2, 4, 3, 0], [0, 1, 2, 3])[:,:,:,:,0]
        P = np.concatenate([P, np.zeros_like(P[:,:,:,0:1])], axis = 3)
        eval, occ = self.read_EIGENVAL()
        for kpnt in range(nkpts):
            eval_t = eval[:,:,kpnt]
            z = (np.zeros((nbands_cder, nbands_cder)) - eval_t.T + eval_t).reshape(nbands_cder, nb_tot,1)
            P[kpnt,:,:,:] *= np.abs(z)
            P[kpnt,:,:,3:4] = z

        eval, occ = eval[:,:,-1], occ[:,:,-1]

        np.save(self.name+"_P_matrix.npy", P)
        self.P_matrix = np.load(self.name+"_P_matrix.npy", mmap_mode='r')
        self.eval = eval

    def read_mommat2up(self):
        with open(self.mommat2up, 'r') as fh:
            Q = np.array([])
            fh.readline() 
            fh.readline() # skip first two lines
            nkpnts = 0

            while True:
                try:
                    # Following is a header for a given k-point
                    a = fh.readline().split()
                    nbands = int(a[6])  # number of bands
                    nlines = int(0.5*nbands*(nbands+1))
                    fh.readline() # skip next line
                    Q = np.append(Q, np.fromfile(fh, count = 9*nlines, sep = ' '))
                    nkpnts += 1
                except:
                    break
                
        dat = Q.reshape((-1, 9))
        P = np.zeros((nkpnts, nbands, nbands, 4), dtype=complex)

        for kpnt in range(nkpnts):
            start =  kpnt   *nlines
            stop  = (kpnt+1)*nlines

            mat_vec = np.concatenate(
                    [dat[start:stop, 2:3] + 1j*dat[start:stop, 3:4],
                    dat[start:stop, 4:5] + 1j*dat[start:stop, 5:6],
                    dat[start:stop, 6:7] + 1j*dat[start:stop, 7:8],
                    dat[start:stop, 8:9]],
                    axis = 1 )
            
            indices = np.triu_indices_from(P[kpnt,:,:,0], k =  0)
            P[kpnt,:,:,0][indices] = mat_vec[:,0]
            P[kpnt,:,:,0] += P[kpnt,:,:,0].conj().T - np.conj(np.diag(np.diag(P[kpnt,:,:,0])))

            P[kpnt,:,:,1][indices] = mat_vec[:,1]
            P[kpnt,:,:,1] += P[kpnt,:,:,1].conj().T - np.conj(np.diag(np.diag(P[kpnt,:,:,1])))

            P[kpnt,:,:,2][indices] = mat_vec[:,2]
            P[kpnt,:,:,2] += P[kpnt,:,:,2].conj().T - np.conj(np.diag(np.diag(P[kpnt,:,:,2])))

            P[kpnt,:,:,3][indices] =  mat_vec[:,3]
            P[kpnt,:,:,3] -= P[kpnt,:,:,3].T + np.diag(np.diag(P[kpnt,:,:,3]))

        np.save(self.name+"_P_matrix.npy", P)
        self.P_matrix = np.load(self.name+"_P_matrix.npy", mmap_mode='r')

    def read_WAVEDERF(self):
        with open(self.WAVEDERF, 'r') as fh:

            Q = np.fromstring(fh.readline(), dtype=int, 
                            sep = ' ')

            if Q.shape[0] != 3 : 
                raise Exception('error reading file')

            nkpnts = Q[1]
            nbands = Q[2]

            Q = np.fromfile(fh, sep = ' ')

        if Q.shape[0] != (12*nbands*nbands*nkpnts) : 
            raise Exception('error reading file')

        dat = Q.reshape((-1, 12))

        P = np.zeros((nkpnts, nbands, nbands, 4), dtype = complex)

        for kpnt in range(nkpnts) :
            start =  kpnt   *nbands**2
            stop  = (kpnt+1)*nbands**2

            eval = dat[start:start+nbands, 4:5]
            occ  = dat[start:start+nbands, 5:6]

            mat_vec = np.concatenate(
                    [dat[start:stop,  6: 7] + 1j*dat[start:stop,  7: 8],
                    dat[start:stop,  8: 9] + 1j*dat[start:stop,  9:10],
                    dat[start:stop, 10:11] + 1j*dat[start:stop, 11:12]],
                    axis = 1 )

            # Introduce z matrix for analogous operation to WAVEDER
            z = np.zeros((nbands, nbands)) - eval.T + eval
            
            P[kpnt,:,:,0] = mat_vec[:,0:1].reshape(nbands, nbands) *abs(z)
            P[kpnt,:,:,1] = mat_vec[:,1:2].reshape(nbands, nbands) *abs(z)
            P[kpnt,:,:,2] = mat_vec[:,2:3].reshape(nbands, nbands) *abs(z)
            P[kpnt,:,:,3] = z

        np.save(self.name+"_P_matrix.npy", P)
        self.P_matrix = np.load(self.name+"_P_matrix.npy", mmap_mode='r')
        self.eval = eval

# END
