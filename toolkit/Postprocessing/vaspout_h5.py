import h5py
import numpy as np
import matplotlib.pyplot as plt
from Plotting import BandStructurePlot, DensityOfStatesPlot
from InformationCollector import InformationCollector


class vaspout_h5:
    """This class knows how to extract relevant information from vaspout.h5
    and put it into InformationCollector object
    """
    def __init__(self, vaspout: str = 'vaspout.h5'):
        self.vaspout = h5py.File(vaspout, 'r')
        self.ic = InformationCollector()
        try: self.ic.folder = vaspout[:vaspout.rindex('/')+1]
        except ValueError: self.ic.folder = ''

    def read_BS(self):
        """_summary_
        """        
        vaspout = self.vaspout
        ic = self.ic
        ic.Pmat = vaspout['results/projectors/par'][:]
        ic.Emat = vaspout['results/electron_eigenvalues/eigenvalues'][:]
        ic.Occp = vaspout['results/electron_eigenvalues/fermiweights'][:]
        ic.nSeg = vaspout['input/kpoints/number_line_segments'][()]/2
        ic.Lvec = vaspout['input/poscar/lattice_vectors'][:]
        ic.Scal = vaspout['input/poscar/scale'][()]
        ic.Kpts = vaspout['results/electron_eigenvalues/kpoint_coords'][:]
        ic.Klen = vaspout['input/kpoints/number_kpoints'][()]
        ic.Klab = vaspout['input/kpoints/labels_kpoints'].asstr()[:]
        ic.ions = vaspout['input/poscar/ion_types'].asstr()[:]
        ic.num_ions = vaspout['input/poscar/number_ion_types'][:]

        ic.Pmat = np.moveaxis(ic.Pmat, [-2,-1], [0,1])
        

        ic.Evalmax = np.max(ic.Emat[ic.Occp>1e-10])
        ic.Econmin = np.min(ic.Emat[ic.Occp<=1e-10])
        ic.Egap = ic.Econmin-ic.Evalmax

        ic.Emat = np.moveaxis(ic.Emat, [-2,-1], [0,1])


    def read_DOS(self):
        """_summary_
        """        
        vaspout = self.vaspout
        ic = self.ic

        ic.Vdos = vaspout['results/electron_dos/dos'][:]
        ic.Vdos = np.moveaxis(ic.Vdos, [-1], [0])

        ic.Pdos = vaspout['results/electron_dos/dospar'][:]
        ic.Pdos = np.moveaxis(ic.Pdos, [-1], [0])

        ic.Edos = vaspout['results/electron_dos/energies'][:]

        ic.Efer = vaspout['results/electron_dos/efermi'][()]

        ic.ions = self.vaspout['input/poscar/ion_types'].asstr()[:]
        ic.num_ions = self.vaspout['input/poscar/number_ion_types'][:]
        
        ic.Emin = np.min(ic.Edos)
        ic.Emax = np.max(ic.Edos)
        ic.Ediff = ic.Emax - ic.Emin
    

    def plot_BS(self,  
                description: str = '', 
                kpaths: str = '', 
                bands: str = '', 
                **kwargs
                ):
        """_summary_

        Args:
            description (str, optional): _description_. Defaults to ''.
            kpaths (str, optional): _description_. Defaults to ''.
            bands (str, optional): _description_. Defaults to ''.
        """        
        self.read_BS()
        ic = self.ic
        ic.description, ic.arrangement, ic.bands = description, kpaths, bands
        BandStructurePlot(ic, **kwargs)


    def plot_BS_default(self):
        """The default plot parameters
        """        
        plt.figure(figsize=(4,8))
        self.plot_BS('clean', bandnums=True)
        plt.tight_layout()
        plt.savefig(self.ic.folder+'BSplot')
    

    def plot_DOS(self, 
                 description='', 
                 ax=None, 
                 E0=None
                 ):
        """_summary_

        Args:
            description (str, optional): _description_. Defaults to ''.
            ax (_type_, optional): _description_. Defaults to None.
            E0 (_type_, optional): _description_. Defaults to None.
        """        
        self.read_DOS()
        self.ic.description = description
        DensityOfStatesPlot(self.ic, ax=ax, E0=E0)


    def plot_DOS_default(self):
        """_summary_
        """        
        plt.figure(figsize=(4,8))
        self.plot_DOS('total',E0=0)
        plt.tight_layout()
        plt.savefig(self.ic.folder+'DOSplot')
    

def plot_BS_DOS(BS:vaspout_h5, 
                DOS:vaspout_h5, 
                description='', 
                arrangement='', 
                bands= ''
                ):
    """_summary_

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
               min_diff=DOS.ic.Ediff*1e-2, ax=axes[0])
    DOS.plot_DOS(description, ax=axes[1], E0=BS.ic.Evalmax)
    plt.tight_layout()