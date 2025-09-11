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


    def extract_magmom(self):
        vaspout = self.vaspout
        Mmat: np.ndarray = vaspout['/intermediate/ion_dynamics/magnetism/spin_moments/values'][:]
        if Mmat.shape[1]==1: print("No magmom in this calculation")
        elif Mmat.shape[1]==2:
            magmom = np.sum(Mmat[0,1,...], axis=1)
        else:
            x = np.sum(Mmat[0,1,...], axis=1)[:,None]
            y = np.sum(Mmat[0,2,...], axis=1)[:,None]
            z = np.sum(Mmat[0,3,...], axis=1)[:,None]
            magmom = np.concatenate([x,y,z], axis=1)

        with open(f"{self.ic.folder}MAGMOM", "w") as f:
            f.write("MAGMOM = ")
            for i, row in enumerate(magmom):
                # Format row with aligned spacing
                row_str = " ".join(f"{x:12.8f}" for x in row)
                if i == 0:
                    f.write(row_str + " \\\n")
                elif i == len(magmom)-1:
                    f.write("         " + row_str + "\n")
                else:
                    f.write("         " + row_str + " \\\n")


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
        
        ic.Efer = vaspout['results/electron_dos/efermi'][()]

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
        """_summary_

        Args:
            description (str, optional): _description_. Defaults to ''.
            ax (_type_, optional): _description_. Defaults to None.
            E0 (_type_, optional): _description_. Defaults to None.
        """        
        self.read_DOS()
        self.ic.description = description
        DensityOfStatesPlot(self.ic, ax=ax, E0=E0, gnuplot=gnuplot, **kwargs)


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


if __name__ == "__main__":
    import argparse

    # Instantiate the parser
    parser = argparse.ArgumentParser(description='Optional app description')
    
    # Required positional argument
    parser.add_argument('calc_type', type=str,
                        help='Either band or dos')

    # Optional positional argument
    parser.add_argument('description', type=str, nargs='?', default='',
                        help='The description of specific projections to be processed')

    # Optional argument
    parser.add_argument('--kpaths', type=str, default='',
                        help='order of kpaths for band')
    
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
    
    args = parser.parse_args()

    calc = vaspout_h5()
    plt.figure(figsize=(4,8))

    if args.calc_type == 'band':
        calc.plot_BS(**vars(args))

    if args.calc_type == 'dos':
        calc.plot_DOS(**vars(args))

    plt.tight_layout()
    plt.savefig(calc.ic.folder+f'{args.calc_type}_plot')