import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
from Plotting import BandStructurePlot, DensityOfStatesPlot
from InformationCollector import InformationCollector

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

# This script can be used to produce the plots straight from the terminal.
# For the full description go to the bottom.

class vaspout_h5:
    """This class knows how to extract relevant information from vaspout.h5
    and put it into an InformationCollector object
    """
    def __init__(self, vaspout: str = 'vaspout.h5'):
        try: self.vaspout = h5py.File(vaspout, 'r')
        except IsADirectoryError: 
            try: self.vaspout = h5py.File(vaspout+'vaspout.h5', 'r')
            except FileNotFoundError: print('WARNING: No vaspout.h5 file found')
        except FileNotFoundError: print('WARNING: No vaspout.h5 file found')
        self.ic = InformationCollector()

        # Get PATH to directory for conveniance
        try: self.ic.folder = vaspout[:vaspout.rindex('/')+1]
        except ValueError: self.ic.folder = ''

        self.WAVEDER = self.ic.folder+'WAVEDER'
        self.EIGENVAL = self.ic.folder+'EIGENVAL'
        self.source = 'VASP'


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

        self.ic.eval = eval
        # self.occ  = occ
        # return eval, occ


    def read_WAVEDER(self):
        self.read_EIGENVAL()
        try: 
            self.ic.WDR_matrix = np.load(self.ic.folder+"WDR_matrix.npy", mmap_mode='r')
            self.ic.eval = self.ic.eval[:,:,-1]
        except FileNotFoundError:
            with FortranFile(self.WAVEDER, 'r') as file:
                nb_tot, nbands_cder, nkpts, ispin = file.read_record(dtype= np.int32)
                nodesn_i_dielectric_function = file.read_record(dtype= float)
                wplasmon = file.read_record(dtype= float).reshape(3,3)
                cder = file.read_record(dtype= np.complex64).reshape(3, ispin,nkpts,nbands_cder,nb_tot)

            # Format cder matrix to resemble P from WAVEDERF
            P = np.moveaxis(cder, [2, 4, 3, 0], [0, 1, 2, 3])[:,:,:,:,0]
            P = np.concatenate([P, np.zeros_like(P[:,:,:,0:1])], axis = 3)
            # eval, occ = self.read_EIGENVAL()
            eval = self.ic.eval
            for kpnt in range(nkpts):
                eval_t = eval[:,:,kpnt]
                z = (np.zeros((nbands_cder, nbands_cder)) - eval_t.T + eval_t).reshape(nbands_cder, nb_tot,1)
                P[kpnt,:,:,:] *= np.abs(z)
                P[kpnt,:,:,3:4] = z

            # eval, occ = eval[:,:,-1], occ[:,:,-1]

            np.save(self.ic.folder+"WDR_matrix.npy", P)
            self.ic.WDR_matrix = P
            self.ic.eval = eval[:,:,-1]


    def Optical_Selection_Rules(self, kpoint, vb_max, num_vb, num_cb, save_txt = True, save_np = False, **kwargs):
        # calculate the optical matrix elements for circular and linear polarized 
        # light and determine the dipole strengths and degree of circular polarization
        #
        # the k-point is defined by the WAVEDERF file
        # the transitions are given in the vb_list and cb_list

        # find VBM

        vb_list = [vb_max - num_vb + i + 1 for i in range(num_vb)]
        cb_list = [vb_max + i + 1 for i in range(num_cb)]

        bands = vb_list + cb_list

        if len(self.ic.WDR_matrix) == 0 :
            self.read_WAVEDER()
        #     bands_to_read = []
        #     for band in bands:
        #         if band not in self.known_bands():
        #             bands_to_read.append(band)
        #     self.read_bands(bands_to_read)

        optical = []

        # polarization verctors
        e_plus  = 1/np.sqrt(2)*np.array([[1,  1j, 0]])
        e_minus = 1/np.sqrt(2)*np.array([[1, -1j, 0]])

        vecs = {} # self.band_vecs
        e_flag = 0

        for vb in vb_list:
            for cb in cb_list:
                # momentum matrix element
                # complex valued and not uniquely defined by an arbitray phase

                if self.ic.WDR_matrix.size > 0 :
                    p = np.array([[self.ic.WDR_matrix[kpoint-1, cb-1, vb-1, 0], 
                                   self.ic.WDR_matrix[kpoint-1, cb-1, vb-1, 1], 
                                   self.ic.WDR_matrix[kpoint-1, cb-1, vb-1, 2]]])

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
                
                if self.ic.eval.size > 0 : 
                    optical += [self.ic.eval[cb-1,0] - self.ic.eval[vb-1,0]]
                    e_flag = 1
                
        optical = np.real(np.array(optical)).reshape(-1, 9 + e_flag)
        header = "vb\t\t# cb\t\tcir_pol\t\t|P+|^2\t\t|P-|^2\t\t|Px|^2\t\t|Py|^2\t\t|Pz|^2\t\t|P|^2"
        if e_flag : header += "\t\t  dE"

        if save_txt : np.savetxt(self.ic.folder+"OSR.txt", optical, fmt = '%.4f', delimiter='\t\t', header=header)
        if save_np : np.save(self.ic.folder+"OSR.npy", optical)
                
        # return "  vb \tcb \tcir_pol\t|P+|^2\t|P-|^2\t|Px|^2\t|Py|^2\t|Pz|^2\t|P|^2\t  dE\n" + \
        #         np.array2string(optical, formatter={'float_kind':lambda x: "%.2f" % x}, separator='\t', suppress_small=True)
    

    def L_Plot(self, kpoint, direction, vb_max, num_vb, num_cb, save_plot = True, save_txt = False, save_np = False, **kwargs):
    # Calculates orbital angular momentum (L)
    # from matrix elements 
    # and band energies read from vasp's WAVEDERF
    # for selected k point number
    # direction: 1 - X, 2 - Y, 3 - Z
        if isinstance(direction, str): direction = {'x':1,'X':1,'y':2,'Y':2,'z':3,'Z':3}[direction]
        if   direction == 1 : beta = 2; gamma = 3
        elif direction == 2 : beta = 1; gamma = 3
        elif direction == 3 : beta = 1; gamma = 2
        else : raise Exception('wrong direction')

        val = [vb_max - num_vb + i + 1 for i in range(num_vb)]
        con = [vb_max + i + 1 for i in range(num_cb)]

        bands_oi = val + con
        nbands_oi = len(bands_oi)

        if len(self.ic.WDR_matrix) == 0 :
            self.read_WAVEDER()
        nbands = self.ic.WDR_matrix.shape[1]
            # bands_to_read = []
            # for band in bands_oi:
            #     if band not in self.known_bands():
            #         bands_to_read.append(band)
            # self.read_bands(bands_to_read)
            # P_vecs = self.band_vecs
            # nbands = P_vecs[vb_max].shape[1]
        # else : nbands = self.ic.WDR_matrix.shape[1]

        L  = np.zeros((nbands, nbands_oi), dtype = complex)

        for count, band in enumerate(bands_oi):
            P_vec = self.ic.WDR_matrix[kpoint-1, band-1,:,:].copy()
            # if self.ic.WDR_matrix.size > 0 : P_vec = self.ic.WDR_matrix[kpoint-1, band-1,:,:].copy()
            # else :                      P_vec = P_vecs[band][kpoint-1].copy()
            P_vec[abs(P_vec[:,3]) < .00000001] = [0,0,0,1]
            L_iter = 2*np.imag(P_vec[:,beta-1]*np.conj(P_vec[:,gamma-1]))/P_vec[:,3]
            L[:,count] = np.cumsum(L_iter)

        C = PhysicalConstants
        if self.source == 'VASP': 
            convert =      C.e * 1E-20 * C.me * 4 * C.pi**2 / C.h**2
        else: 
            convert = C.Ry * (C.a0)**2 * C.me * 4 * C.pi**2 / C.h**2

        L = np.real(L*convert)

        if save_plot:
            plt.figure()
            plt.plot(L, '-')
            dirs = {1: 'X', 2: 'Y', 3: 'Z'}
            plt.title('K-point: ' + str(kpoint+1) + ', ' \
                    + 'Direction: ' + dirs[direction])
            plt.xlabel('number of bands')
            plt.ylabel('L (hbar)')
            plt.legend(['v'+str(i) for i in val]+['c'+str(i) for i in con] )
            plt.savefig(self.ic.folder+"L_plot.png")

        if save_txt: np.savetxt(self.ic.folder+'L_plot.txt', np.concatenate([np.array(bands_oi).reshape((1,-1)), L], axis = 0))
        if save_np: np.save(self.ic.folder+'L_plot.npy', np.concatenate([np.array(bands_oi).reshape((1,-1)), L], axis = 0))


    def extract_magmom(self, n=3):
        """Get MAGMOM info and put it in POSCAR (as a comment)
        """        
        vaspout = self.vaspout
        ic = self.ic
        ic.ions = vaspout['input/poscar/ion_types'].asstr()[:]
        ic.num_ions = vaspout['input/poscar/number_ion_types'][:]

        magmom_matrix: np.ndarray = vaspout['/intermediate/ion_dynamics/magnetism/spin_moments/values'][:]
        poscar = vaspout['/original/poscar/content'].asstr()[()].split('\n')

        ions = []
        for ion, num_ion in zip(ic.ions, ic.num_ions):
            ions += [ion]*num_ion

        # Get the full MAGMOM matrix and sum over orbitals
        if magmom_matrix.shape[1]==1: # ISPIN = 1
            print("No magmom in this calculation")
        elif magmom_matrix.shape[1]==2: # ISPIN = 2
            magmom = np.sum(magmom_matrix[0,1,...], axis=1)
        else: # lnoncolinear = .true.
            x = np.sum(magmom_matrix[0,1,...], axis=1)[:,None]
            y = np.sum(magmom_matrix[0,2,...], axis=1)[:,None]
            z = np.sum(magmom_matrix[0,3,...], axis=1)[:,None]
            magmom = np.concatenate([x,y,z], axis=1)

        magflag = -1
        with open(f"{self.ic.folder}POSCAR", "w") as f:
            for row in poscar:
                row = row.replace('\r','').replace('\n','')
                f.write(f"{row}")
                if magflag >= 0 and magflag < len(ions):
                    mags = ' '.join(f"{x:9.5f}" for x in magmom[magflag])
                    f.write(f'\t # {ions[magflag]} : \t{mags}')
                    magflag += 1
                f.write('\n')
                if 'Direct' in row or 'Cartesian' in row:
                    magflag = 0


    def read_BS(self):
        """Put all relevant band structure information from vaspout to self.ic: InformationCollector
        """        
        vaspout = self.vaspout
        ic = self.ic
        ic.Pmat = vaspout['results/projectors/par'][:]
        ic.Emat = vaspout['results/electron_eigenvalues/eigenvalues'][:]
        ic.Occp = vaspout['results/electron_eigenvalues/fermiweights'][:]
        try: ic.nSeg = vaspout['input/kpoints/number_line_segments'][()]/2
        except KeyError: ic.nSeg = 1
        ic.Lvec = vaspout['input/poscar/lattice_vectors'][:]
        ic.Scal = vaspout['input/poscar/scale'][()]
        ic.Kpts = vaspout['results/electron_eigenvalues/kpoint_coords'][:]
        ic.Klen = vaspout['input/kpoints/number_kpoints'][()]
        ic.Klab = vaspout['input/kpoints/labels_kpoints'].asstr()[:]
        ic.Efer = vaspout['results/electron_dos/efermi'][()]
        ic.ions = vaspout['input/poscar/ion_types'].asstr()[:]
        ic.num_ions = vaspout['input/poscar/number_ion_types'][:]

        # Rearrange PROCAR and EIGENVAL matrices
        # into the desired format: ( ..., ndirs, nions, norbits )
        # and calculate bandgap
        ic.Pmat = np.moveaxis(ic.Pmat, [-2,-1], [0,1])
        ic.Evalmax = np.max(ic.Emat[ic.Occp>1e-1])
        ic.Econmin = np.min(ic.Emat[ic.Occp<=1e-1])
        ic.Egap = ic.Econmin-ic.Evalmax
        ic.Emat = np.moveaxis(ic.Emat, [-2,-1], [0,1])


    def read_DOS(self):
        """Put all relevant density of states information from vaspout to self.ic: InformationCollector
        """        
        vaspout = self.vaspout
        ic = self.ic

        ic.Vdos = vaspout['results/electron_dos/dos'][:]
        ic.Pdos = vaspout['results/electron_dos/dospar'][:]
        ic.Edos = vaspout['results/electron_dos/energies'][:]
        ic.Efer = vaspout['results/electron_dos/efermi'][()]
        ic.ions = vaspout['input/poscar/ion_types'].asstr()[:]
        ic.num_ions = vaspout['input/poscar/number_ion_types'][:]

        # Rearrange the matrices to desired format
        ic.Vdos = np.moveaxis(ic.Vdos, [-1], [0])
        ic.Pdos = np.moveaxis(ic.Pdos, [-1], [0])

        # Calculate energy range
        ic.Emin = np.min(ic.Edos)
        ic.Emax = np.max(ic.Edos)
        ic.Ediff = ic.Emax - ic.Emin
    

    def plot_BS(self,  
                description: str = '', 
                kpaths: str = '', 
                bands: str = '', 
                **kwargs
                ):
        """Make a band structue plot according to string specifiers 
        and other args that BandStructurePlot accepts

        Args:
            description (str, optional): the sum of which projections should be plotted. Defaults to ''.
            kpaths (str, optional): assuming line mode in band structure calculation 
                                    any ordering of those lines can be selected (including inversion). Defaults to ''.
            bands (str, optional): which bands are to be displayed in the form of an intuitive string. Defaults to ''.

        Kwargs:
            ax (plt.axes): to which axis should the plot be rendered. If None (default) uses plt.gca(). 
            
            E0 (float | str): what energy should be used as 0. str options are ('fermi'). If None (default) will look for valence band maximum. 
            
            folder (str): to which folder should the output(s) be saved. If '' (default) uses vaspout.h5's folder. 
            
            save (bool): should the plot be immediately saved. Defaults to False.
            
            gnuplot (bool): should the plot be written to a file in gnuplot format. Defaults to False, 
            
            bandnums (bool): should band numbers be displayed next to the plot. Defaults to False. 
            
            min_diff (float): if band numbers are to be displayed, how far apart they should be to not group them. Defaults to .1
            
            color (int | float | str | None): same as in all pyplot figures. Defaults to None 
            
            mult (float): the dot size multiplier. Defaults to 1
            
            alpha (float): the dot opacity level in range [0, 1]. Defaults to .1

            skip (int): how many kpoints to skip

            labels (list[(str, str)]): manual labels

            nohlines (bool): disable hlines
        """        
        self.read_BS()
        ic = self.ic
        ic.description, ic.kpaths, ic.bands = description, kpaths, bands

        if 'custom_ions' in kwargs.keys(): ic.ions = kwargs['custom_ions']

        if 'skip' in kwargs.keys():
            skip = kwargs['skip']
            ic.Pmat = ic.Pmat[skip:]
            ic.Emat = ic.Emat[skip:]
            ic.Occp = ic.Occp[skip:]
            ic.Kpts = ic.Kpts[skip:]
            ic.Klen = ic.Klen - skip

        if 'labels' in kwargs.keys():
            labels = kwargs['labels']
            ic.Klab = labels
            ic.nSeg = len(labels)/2
            ic.Klen = int(ic.Klen / ic.nSeg)

        BandStructurePlot(ic, **kwargs)


    def plot_BS_default(self):
        """The default plot, gives overview of the entire band structure
        """        
        plt.figure(figsize=(4,8))
        self.plot_BS('clean', bandnums=True, E0='fermi')
        plt.tight_layout()
        plt.savefig(self.ic.folder+'BSplot')
    

    def plot_DOS(self, 
                 description='', 
                 ax=None, 
                 E0=None,
                 gnuplot = False,
                 **kwargs
                 ):
        """Make a density of states plot according to string specifiers 
        and other args that BandStructurePlot accepts

        Args:
            description (str, optional): the sum of which projections should be plotted. Defaults to ''.
            ax (plt.axes): to which axis should the plot be rendered. If None (default) uses plt.gca(). 
            E0 (float | str): what energy should be used as 0. str options are ('fermi'). If None (default) will look for valence band maximum. 
            gnuplot (bool): should the plot be written to a file in gnuplot format. Defaults to False, 
        """        
        self.read_DOS()
        self.ic.description = description
        DensityOfStatesPlot(self.ic, ax=ax, E0=E0, gnuplot=gnuplot, **kwargs)


    def plot_DOS_default(self):
        """The default plot, shows the entire density of states
        """        
        plt.figure(figsize=(4,8))
        self.plot_DOS('total',E0=0)
        plt.tight_layout()
        plt.savefig(self.ic.folder+'DOSplot')
    

def plot_BS_DOS(BS:vaspout_h5, 
                DOS:vaspout_h5, 
                description='', 
                arrangement='', 
                bands= '',
                ):
    """A simple function to produce an elegant band structure and density of states plot

    Args:
        BS (vaspout_h5): _description_
        DOS (vaspout_h5): _description_
        description (str, optional): _description_. Defaults to ''.
        arrangement (str, optional): _description_. Defaults to ''.
        bands (str, optional): _description_. Defaults to ''.
    """    
    fig, axes = plt.subplots(ncols=2, sharey=True, gridspec_kw={'width_ratios': [3, 1]})
    DOS.read_DOS()
    BS.plot_BS(description, arrangement, bands, 
               min_diff=DOS.ic.Ediff*1e-2, ax=axes[0], bandnums=True, E0 = 'fermi')
    DOS.plot_DOS(description, ax=axes[1])
    plt.tight_layout()


# This script can be used to produce the plots straight from the terminal.
if __name__ == "__main__":
    import argparse

    # Instantiate the parser
    parser = argparse.ArgumentParser(description='Make band stucture or density of states plot')
    
    # Required positional argument
    parser.add_argument('calc_type', type=str,
                        help='Either band or dos')

    parser.add_argument('description', type=str, nargs='?', default='',
                        help='The description of specific projections to be processed, e.g. "1,2 Mo x p dxy"')

    parser.add_argument('--kpaths', type=str, default='',
                        help='Order and direction of kpaths for band, e.g. "1 2 3 -4"')
    
    parser.add_argument('--bands', type=str, default='',
                        help='choice of bands for band')
    
    parser.add_argument("--min_diff", type=float, default=0.1,
                        help="Minimum difference threshold (default: 0.1).")

    parser.add_argument("--ax", default=None,
                        help="Matplotlib Axes object (default: None).")

    parser.add_argument("--E0", type=float, default=None,
                        help="Reference energy level (default: None).")

    parser.add_argument("--gnuplot", action="store_true",
                        help="Enable gnuplot output (default: False).")

    parser.add_argument("--folder", type=str, default="",
                        help="Output folder (default: '').")

    parser.add_argument("--save", action="store_true",
                        help="Save plots to file instead of showing (default: False).")

    parser.add_argument("--bandnums", action="store_true",
                        help="Show band numbers (default: False).")

    parser.add_argument("--color", type=str, default=None,
                        help="Color for plotting (default: None).")

    parser.add_argument("--mult", type=float, default=1,
                        help="Scaling multiplier (default: 1).")

    parser.add_argument("--alpha", type=float, default=0.1,
                        help="Transparency level (default: 0.1).")
    
    parser.add_argument("--kpoint", type=int,
                    help="K-point at which to evaluate the plot.")

    parser.add_argument("--direction", type=str,
                        help="Direction for the L-plot calculation.")

    parser.add_argument("--vb_max", type=int,
                        help="Maximum valence band energy.")

    parser.add_argument("--num_vb", type=int,
                        help="Number of valence bands to include.")

    parser.add_argument("--num_cb", type=int,
                        help="Number of conduction bands to include.")

    parser.add_argument("--save_plot", action="store_true",
                        help="Save the generated plot (default: False).")

    parser.add_argument("--save_txt", action="store_true",
                        help="Save data to a text file (default: False).")

    parser.add_argument("--save_np", action="store_true",
                        help="Save data to a NumPy file (default: False).")

    
    args = parser.parse_args()

    calc = vaspout_h5()

    if args.calc_type == 'band':
        plt.figure(figsize=(4,8))
        calc.plot_BS(**vars(args))
        plt.tight_layout()
        plt.savefig(f'{calc.ic.folder}{args.calc_type}_plot')

    if args.calc_type == 'dos':
        plt.figure(figsize=(4,8))
        calc.plot_DOS(**vars(args))
        plt.tight_layout()
        plt.savefig(f'{calc.ic.folder}{args.calc_type}_plot')

    if args.calc_type == 'OSR':
        calc.Optical_Selection_Rules(**vars(args))

    if args.calc_type == 'L_plot':
        calc.L_Plot(**vars(args))
